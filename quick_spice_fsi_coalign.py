#!/usr/bin/env python

import argparse
import datetime
import os
import re

from astropy import units as u
from astropy import wcs
from astropy.io import fits
from dateutil.parser import parse as parse_date
import pandas as pd
import scipy.interpolate as si
import yaml

from eui.euiprep import euiprep
from papy.sol.data.solo_eui import EUISelektorClient
import align_images
import spice_stew


class SpiceUtils:
    re_spice_L123_filename = re.compile('''
        solo
        _(?P<level>L[123])
        _spice
            (?P<concat>-concat)?
            -(?P<slit>[wn])
            -(?P<type>(?:ras|sit|exp))
            (?P<db>-db)?
            (?P<int>-int)?
        _(?P<time>\d{8}T\d{6})
        _(?P<version>V\d{2})
        _(?P<SPIOBSID>\d+)-(?P<RASTERNO>\d+)
        \.fits
        ''',
        re.VERBOSE)

    def read_spice_uio_catalog():
        """
        Read UiO text table SPICE FITS files catalog
        http://astro-sdc-db.uio.no/vol/spice/fits/spice_catalog.txt

        Return
        ------
        pandas.DataFrame
            Table

        Example queries that can be done on the result:

        * `df[(df.LEVEL == "L2") & (df["DATE-BEG"] >= "2020-11-17") & (df["DATE-BEG"] < "2020-11-18") & (df.XPOSURE > 60.)]`
        * `df[(df.LEVEL == "L2") & (df.STUDYDES == "Standard dark for cruise phase")]`

        Source: https://spice-wiki.ias.u-psud.fr/doku.php/data:data_analysis_manual:read_catalog_python
        """
        cat_file = os.path.join(
            os.getenv('SOLO_ARCHIVE', '/archive/SOLAR-ORBITER/'),
            'SPICE/fits/spice_catalog.txt')
        columns = list(pd.read_csv(cat_file, nrows=0).keys())
        date_columns = ['DATE-BEG','DATE', 'TIMAQUTC']
        df = pd.read_table(cat_file, skiprows=1, names=columns, na_values="MISSING",
                        parse_dates=date_columns, low_memory=False)
        df.LEVEL = df.LEVEL.apply(lambda string: string.strip())
        df.STUDYTYP = df.STUDYTYP.apply(lambda string: string.strip())
        return df

    def parse_filename(filename):
        m = SpiceUtils.re_spice_L123_filename.match(filename)
        if m is None:
            raise ValueError(f'could not parse SPICE filename: {filename}')
        return m.groupdict()

    def filename_to_date(filename):
        d = SpiceUtils.parse_filename(filename)
        return parse_date(d['time'])

    def ias_fullpath(filename):
        d = SpiceUtils.parse_filename(filename)
        date = parse_date(d['time'])

        fullpath = os.path.join(
            os.getenv('SOLO_ARCHIVE', '/archive/SOLAR-ORBITER/'),
            'SPICE/fits/',
            'level' + d['level'].lstrip('L'),
            f'{date.year:04d}/{date.month:02d}/{date.day:02d}',
            filename)

        return fullpath


class EuiUtils:
    def ias_fullpath(rob_fullpath):
        p = rob_fullpath.lstrip('/data/solo-eui/internal/')
        p = '/archive/SOLAR-ORBITER/EUI/data_internal/' + p
        return p

    def local_L2_path(output_dir, fsi_file_L1):
        base = os.path.basename(fsi_file_L1)
        base = base.replace('L1', 'L2')
        return os.path.join(output_dir, base)


def list_spice_files(start_date, end_date, study_name=None):
    ''' Get list of SPICE files

    Parameters
    ==========
    start_date : str (YYYY-MM-DD)
        Query start date
    end_date : str (YYYY-MM-DD)
        Query end date
    study_name : str or None (default: None)
        Study name. If None, return all studies.

    Returns
    =======
    filenames : list of str
        List of FITS
    '''

    cat = SpiceUtils.read_spice_uio_catalog()
    filters = (
        (cat['DATE-BEG'] > start_date)
        & (cat['DATE-BEG'] <= end_date)
        & (cat['LEVEL'] == 'L2')
        )
    if study_name is not None:
        filters &= (cat['STUDY'] == study_name)
    results = cat[filters]
    return list(results['FILENAME'])


def get_closest_fsi_L1_file_from_selektor(search_date, band, max_t_dist):
    ''' Get FSI L1 file closest to a given date

    Parameters
    ==========
    search_date : datetime.datetime
        Query date
    band : str ('174' or '304')
        Instrument band
    max_t_dist : datetime.timedelta
        Query window size: FSI files are searched within
        [date - max_t_dist, date + max_t_dist].

    Returns
    =======
    filename : str or None
        Closest FSI FITS, if there is one.
    '''
    search_date_min = search_date - max_t_dist
    search_date_max = search_date + max_t_dist

    # Search FSI data
    eui_client = EUISelektorClient()
    base_params = {
        'level[]': 'L1',
        'detector[]': 'FSI',
        'wavelnth[]': band,
        'imgtype[]': 'solar image',
        'limit[]': 100,
        }
    # Because of the results limit, perform two queries:
    # - 'before', in descending from search_date to search_date_min,
    #   yielding the Nth file before search_date
    # - 'after', in ascending from search_date to search_date_max,
    #   yielding the Nth file after search_date
    before_params = {
        'order[]': 'DESC',
        'date_begin_start': search_date_min.isoformat().split('T')[0],
        'date_begin_start_hour': search_date_min.hour,
        'date_begin_start_minute': search_date_min.minute,
        'date_begin_end': search_date.isoformat().split('T')[0],
        'date_begin_end_hour': search_date.hour,
        'date_begin_end_minute': search_date.minute,
        }
    after_params = {
        'order[]': 'ASC',
        'date_begin_start': search_date.isoformat().split('T')[0],
        'date_begin_start_hour': search_date.hour,
        'date_begin_start_minute': search_date.minute,
        'date_begin_end': search_date_max.isoformat().split('T')[0],
        'date_begin_end_hour': search_date_max.hour,
        'date_begin_end_minute': search_date_max.minute,
        }
    before_params.update(base_params)
    after_params.update(base_params)
    res_before = eui_client.search(before_params)
    res_after = eui_client.search(after_params)

    # keep files just before and just after
    if res_before is not None:
        res_before = res_before.iloc[0]
    if res_after is not None:
        res_after = res_after.iloc[0]

    # return closest file, or None
    if res_before is None:
        return res_after  # may be None
    elif res_after is None:
        return res_before  # may be None
    else:
        date_before = parse_date(res_before['date-beg'])
        date_after = parse_date(res_after['date-beg'])
        dt_before = abs((search_date - date_before).total_seconds())
        dt_after = abs((search_date - date_after).total_seconds())
        if dt_before < dt_after:
            return res_before
        else:
            return res_after


def get_fsi_L1(spice_file, band, output_dir, max_t_dist=6):
    ''' Get FSI L1 file to coalign with a SPICE file

    Parameters
    ==========
    spice_file : str
        Path to a SPICE FITS file
    band : str ('174' or '304')
        Instrument band
    output_dir : str
        Output directory
    max_t_dist : float (default: 6)
        Maximum time distance in hour. If no file exists within [date -
        max_t_dist, date + max_t_dist], return None.

    Returns
    =======
    filename : str or None
        Closest FSI FITS, if there is one.
    '''
    os.makedirs(output_dir, exist_ok=True)

    spice_file_base = os.path.splitext(os.path.basename(spice_file))[0]
    cache_file = f'{output_dir}/{spice_file_base}_fsi_info.yml'

    if os.path.isfile(cache_file):
        with open(cache_file, 'r') as f:
            return yaml.safe_load(f)
    else:
        spice_hdu = fits.open(spice_file)
        search_date = parse_date(spice_hdu[0].header['DATE-AVG'])
        max_t_dist = datetime.timedelta(hours=max_t_dist)
        fsi_file = get_closest_fsi_L1_file_from_selektor(
            search_date, band,
            max_t_dist,
            )
        if fsi_file is not None:
            fsi_file = fsi_file.to_dict()
        with open(cache_file, 'w') as f:
            yaml.safe_dump(fsi_file, f, sort_keys=False)
        return fsi_file


def gen_fsi_L2(fsi_file_L1, output_dir):
    fsi_file_L2 = EuiUtils.local_L2_path(output_dir, fsi_file_L1)
    if os.path.isfile(fsi_file_L2):
        print(f'FSI L2 file exists: {fsi_file_L2}, exiting')
    else:
        euimap_L2 = euiprep(
            fsi_file_L1,
            auto=True,
            save_L2=False,
            )
        euimap_L2.save_fits(fsi_file_L2)
    return fsi_file_L2



def get_spice_image_data(filename, window):
    ''' Return SPICE image data

    Parameters
    ==========
    filename : str
        Path to a FITS file containing SPICE intensity maps coaligned with
        spice_stew.
    spice_window : str
        SPICE window name

    Returns
    =======
    img : 2D array
        Image data
    header : astropy.io.fits.Header
        FITS header
    '''
    hdulist = fits.open(filename)
    hdu = hdulist[window]

    img = np.squeeze(hdu.data)  # remove 1-length dimensions, ie t and wvl
    if img.ndim != 2:
        msg = f'invalid image shape {img.shape} for {filename}'
        raise ValueError(msg)

    return img, hdu.header


def get_fsi_image_data(filename):
    ''' Return FSI imge data

    Parameters
    ==========
    filename : dict
        Path to FSI file FITS.

    Returns
    =======
    img : 2D array
        Image data
    header : astropy.io.fits.Header
        FITS header
    '''
    hdulist = fits.open(filename)
    hdu = hdulist[0]
    return hdu.data, hdu.header


def gen_images_to_coalign(spice_file, spice_window, fsi_file, output_dir):
    ''' Generate SPICE and FSI images that can be coaligned together

    Parameters
    ==========
    spice_file : str
        Path to a FITS file containing SPICE intensity maps corrected with
        spice_stew.
    spice_window : str
        SPICE window name
    filename : dict
        Path to FSI file FITS.

    Returns
    =======
    spice_img : str
    fsi_img : str
        Path to FITS images to coalign
    '''
    os.makedirs(output_dir, exist_ok=True)

    basename = os.path.basename(spice_file_aligned).rstrip('_remapped_img.fits')
    new_spice_filename = f'{basename}_spice_img.fits'
    new_spice_filename = os.path.join(output_dir, new_spice_filename)
    new_fsi_filename = f'{basename}_fsi_img.fits'
    new_fsi_filename = os.path.join(output_dir, new_fsi_filename)
    plot_filename = f'{basename}.pdf'
    plot_filename = os.path.join(output_dir, plot_filename)

    if os.path.isfile(new_fsi_filename) and os.path.isfile(new_spice_filename):
        print(f'Coalign images exist for: {basename}, skipping')
        return new_spice_filename, new_fsi_filename

    spice_img, spice_header = get_spice_image_data(
        spice_file,
        spice_window,
        )
    fsi_img, fsi_header = get_fsi_image_data(fsi_file)

    # Cut spice image
    # cut-off sides  FIXME: correct values?
    iymin = 130
    iymax = 670
    spice_img = spice_img[iymin:iymax]

    # SPICE WCS
    wcs_spice = wcs.WCS(spice_header)
    # px coordinates of SPICE cut zone
    ny_spice, nx_spice = spice_img.shape
    x_spice = np.arange(0, nx_spice)
    y_spice = np.arange(iymin, iymax)
    w_spice = np.array([0])
    t_spice = np.array([1])
    assert y_spice.size <= ny_spice
    px_spice = np.meshgrid(x_spice, y_spice, w_spice, t_spice)
    px_spice = np.moveaxis(px_spice, 0, -1)
    # convert to world coordinates
    world_spice = wcs_spice.wcs_pix2world(px_spice.reshape(-1, 4), 0)
    world_spice = world_spice.reshape((y_spice.size, x_spice.size, 1, 1, 4))
    assert px_spice.shape == world_spice.shape
    Tx_spice = u.Quantity(world_spice[:, :, 0, 0, 0], wcs_spice.world_axis_units[0])
    Ty_spice = u.Quantity(world_spice[:, :, 0, 0, 1], wcs_spice.world_axis_units[1])
    # put angles between ]-180, +180] deg
    pi = u.Quantity(180, 'deg')
    Tx_spice = - ((- Tx_spice + pi) % (2*pi) - pi)
    Ty_spice = - ((- Ty_spice + pi) % (2*pi) - pi)
    # convert to arcsec
    Tx_spice = Tx_spice.to('arcsec').value
    Ty_spice = Ty_spice.to('arcsec').value

    # SPICE WCS
    wcs_fsi = wcs.WCS(fsi_header)
    # px coordinates of SPICE cut zone
    ny_fsi, nx_fsi = fsi_img.shape
    x_fsi = np.arange(0, nx_fsi)
    y_fsi = np.arange(0, ny_fsi)
    px_fsi = np.meshgrid(x_fsi, y_fsi)
    px_fsi = np.moveaxis(px_fsi, 0, -1)
    # convert to world coordinates
    world_fsi = wcs_fsi.wcs_pix2world(px_fsi.reshape(-1, 2), 0)
    world_fsi = world_fsi.reshape((y_fsi.size, x_fsi.size, 2))
    assert px_fsi.shape == world_fsi.shape
    Tx_fsi = u.Quantity(world_fsi[:, :, 0], wcs_fsi.world_axis_units[0])
    Ty_fsi = u.Quantity(world_fsi[:, :, 1], wcs_fsi.world_axis_units[1])
    # put angles between ]-180, +180] deg
    pi = u.Quantity(180, 'deg')
    Tx_fsi = - ((- Tx_fsi + pi) % (2*pi) - pi)
    Ty_fsi = - ((- Ty_fsi + pi) % (2*pi) - pi)
    # convert to arcsec
    Tx_fsi = Tx_fsi.to('arcsec').value
    Ty_fsi = Ty_fsi.to('arcsec').value

    # Generate common coordinates
    common_Txy_size = 4  # arcsec
    common_Txy_pad = 200  # arcsec
    Tx_common_1d = np.arange(
        Tx_spice.min() - common_Txy_pad,
        Tx_spice.max() + common_Txy_size + common_Txy_pad,
        common_Txy_size,
        )
    Ty_common_1d = np.arange(
        Ty_spice.min() - common_Txy_pad,
        Ty_spice.max() + common_Txy_size + common_Txy_pad,
        common_Txy_size,
        )
    Tx_common, Ty_common = np.meshgrid(Tx_common_1d, Ty_common_1d)

    # Common WCS
    w_common = wcs.WCS(naxis=2)
    w_common.wcs.cdelt = [common_Txy_size, common_Txy_size]
    w_common.wcs.crpix = [Tx_common_1d.size//2, Ty_common_1d.size//2]
    w_common.wcs.crval = [Tx_common_1d[int(w_common.wcs.crpix[0])],
                          Ty_common_1d[int(w_common.wcs.crpix[1])]]
    w_common.wcs.ctype = ['HPLN-TAN', 'HPLT-TAN']
    w_common.wcs.cunit = ['arcsec', 'arcsec']

    # Pre-cut FSI
    ixmin_fsi = np.where(Tx_fsi > Tx_common.min())[1].min() - 10
    ixmax_fsi = np.where(Tx_fsi < Tx_common.max())[1].max() + 10
    iymin_fsi = np.where(Ty_fsi > Ty_common.min())[0].min() - 10
    iymax_fsi = np.where(Ty_fsi < Ty_common.max())[0].max() + 10
    fsi_img = fsi_img[iymin_fsi:iymax_fsi+1, ixmin_fsi:ixmax_fsi+1]
    Tx_fsi = Tx_fsi[iymin_fsi:iymax_fsi+1, ixmin_fsi:ixmax_fsi+1]
    Ty_fsi = Ty_fsi[iymin_fsi:iymax_fsi+1, ixmin_fsi:ixmax_fsi+1]

    # Remap
    def remap(Tx, Ty, img, new_Tx, new_Ty):
        points = np.array([Tx, Ty]).reshape((2, -1)).T
        values = img.flatten()
        new_points = np.array([new_Tx, new_Ty]).reshape((2, -1)).T
        new_values = si.griddata(points, values, new_points)
        new_values = new_values.reshape(new_Tx.shape)
        return new_values
    new_spice_img = remap(Tx_spice, Ty_spice, spice_img, Tx_common, Ty_common)
    new_fsi_img = remap(Tx_fsi, Ty_fsi, fsi_img, Tx_common, Ty_common)

    # Save FITS
    spice_hdu = fits.PrimaryHDU(new_spice_img, header=w_common.to_header())
    spice_hdu.writeto(new_spice_filename)
    fsi_hdu = fits.PrimaryHDU(new_fsi_img, header=w_common.to_header())
    fsi_hdu.writeto(new_fsi_filename)

    plot_images(new_spice_img, new_fsi_img, w_common, plot_filename)

    return new_spice_filename, new_fsi_filename


def plot_images(spice_img, fsi_img, wcs, filename):
    import sunpy.map

    spice_map = sunpy.map.Map(spice_img, wcs)
    fsi_map = sunpy.map.Map(fsi_img, wcs)
    spice_fov_map = sunpy.map.Map(np.isnan(spice_img).astype(float), wcs)
    fsi_norm = plt.matplotlib.colors.LogNorm(
            vmin=np.max([1, np.nanpercentile(fsi_map.data, 1)]),
            vmax=np.max([10, np.nanpercentile(fsi_map.data, 99.9)]),
            )
    spice_norm = plt.matplotlib.colors.LogNorm(
            vmin=np.max([1, np.nanpercentile(spice_map.data, 0)]),
            vmax=np.max([10, np.nanpercentile(spice_map.data, 99)]),
            )

    plt.clf()
    ax1 = plt.subplot(121, projection=fsi_map)
    fsi_map.plot(
        norm=fsi_norm,
        cmap='viridis',
        )
    ax2 = plt.subplot(122, projection=spice_map)
    spice_map.plot(
        norm=spice_norm,
        cmap='viridis',
        )
    contours = spice_map.contour(
        spice_norm.vmax - (.1*(spice_norm.vmax - spice_norm.vmin)),
        )
    contours_fov = spice_fov_map.contour(.5)
    for ax in [ax1, ax2]:
        for c in contours:
            ax.plot_coord(c, color='w', linewidth=.5)
    if contours_fov:
        ax1.plot_coord(contours_fov[0], color='w', linewidth=.5)
    ax1.set_xlabel('$\\theta_x$')
    ax1.set_ylabel('$\\theta_y$')
    ax2.set_xlabel('$\\theta_x$')
    ax2.set_ylabel(' ')
    ax1.set_title(' ')
    ax2.set_title(' ')
    plt.tight_layout()
    plt.savefig(filename)


def coalign_spice_fsi_images(spice_img, fsi_img, output_dir):
    ''' Coalign SPICE and FSI images

    Parameters
    ==========
    spice_img : str
    fsi_img : str
        Path to FITS images to coalign, written by `gen_images_to_coalign()`
    output_dir : str
        Output directory
    '''
    os.makedirs(output_dir, exist_ok=True)

    basename = os.path.basename(spice_img)
    basename = basename.rstrip('_coalign_spice_img.fits')
    basename = os.path.join(output_dir, basename)
    plot_filename = f'{basename}_coaligned.pdf'
    yml_filename = f'{basename}_coaligned.yml'
    if os.path.isfile(yml_filename):
        return

    spice_hdu = fits.open(spice_img)[0]
    fsi_hdu = fits.open(fsi_img)[0]
    w = wcs.WCS(spice_hdu.header)

    cube = np.stack([spice_hdu.data, fsi_hdu.data])
    shifts, max_cc = align_images.align.track(
        fsi_hdu.data,
        spice_hdu.data,
        missing=np.nan,
        )
    shifts = shifts[::-1]  # y, x to x, y
    aligned_cube = np.squeeze(align_images.align.align_cube(
        cube,
        np.stack([shifts, [0, 0]]),
        ))
    plot_images(aligned_cube[0], aligned_cube[1], w, plot_filename)

    res = dict(
        dx=float(shifts[0]),
        dy=float(shifts[1]),
        max_cc=float(max_cc),
        )
    with open(yml_filename, 'w') as f:
        yaml.safe_dump(res, f, sort_keys=False)



if __name__ == '__main__':

    p = argparse.ArgumentParser()
    p.add_argument('--start-date',
                   help='processing start date (YYYY-MM-DD)')
    p.add_argument('--end-date',
                   help='processing end date (YYYY-MM-DD)')
    p.add_argument('--output-dir', default='./output',
                   help='output directory')
    args = p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print('Listing SPICE files')
    spice_filenames = list_spice_files(
        args.start_date,
        args.end_date,
        study_name='SCI_SYNOPTIC_SC_SL04_60.0S_FF',
        )

    ssp = spice_stew.SpiceSpicePointing()
    for spice_file in spice_filenames:
        spice_file = SpiceUtils.ias_fullpath(spice_file)
        print('\nProcessing', spice_file)

        print('Correcting pointing with SPICE kernels')
        spice_file_aligned = spice_stew.correct_spice_pointing(
            ssp,
            spice_file,
            f'{args.output_dir}/spice_stew',
            overwrite=False,
            plot_results=True,
            sum_wvl=True,
            )

        print('Getting closest FSI image')
        fsi_file_L1 = get_fsi_L1(
            spice_file,
            '304',
            f'{args.output_dir}/fsi_data',
            max_t_dist=3,
            )
        if fsi_file_L1 is None:
            print('No FSI image found for {spice_file}, skipping')
            continue
        fsi_file_L1 = EuiUtils.ias_fullpath(fsi_file_L1['filepath'])
        print(fsi_file_L1)

        print('Generating L2 FSI image')
        fsi_file_L2 = gen_fsi_L2(fsi_file_L1, f'{args.output_dir}/fsi_data')

        print('Generating images to coalign')
        spice_img, fsi_img = gen_images_to_coalign(
            spice_file_aligned,
            'Ly-gamma-CIII group bin (1/4)',
            fsi_file_L2,
            f'{args.output_dir}/coalign_input',
            )

        print('Coaligning images')
        coalign_spice_fsi_images(
            spice_img,
            fsi_img,
            f'{args.output_dir}/coalign_output',
            )
