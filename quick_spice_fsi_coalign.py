#!/usr/bin/env python

import argparse
import datetime
import os
import traceback

from astropy import units as u
from astropy import wcs
from astropy.io import fits
from dateutil.parser import parse as parse_date
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as si
import spiceypy.utils.exceptions
import yaml

from eui.euiprep import euiprep
from eui_selektor_client import EUISelektorClient
from papy.sol.coord import diff_rot
import align_images
import papy.plot
import spice_stew

import common


def list_spice_files(start_date, end_date, study_id):
    """ Get list of SPICE files

    Parameters
    ==========
    start_date : str (YYYY-MM-DD)
        Query start date
    end_date : str (YYYY-MM-DD)
        Query end date
    study_id: str or None (default: None)
        Study id in Miso.

    Returns
    =======
    filenames : list of str
        List of FITS
    """
    if type(study_id) is not str:
        raise ValueError(f'study_id must be str (got {type(study_id)})')

    cat = common.SpiceUtils.read_spice_uio_catalog()
    filters = (
        (cat['DATE-BEG'] > start_date)
        & (cat['DATE-BEG'] <= end_date)
        & (cat['LEVEL'] == 'L2')
        & (cat['MISOSTUD'] == study_id)
        )
    results = cat[filters]
    return list(results['FILENAME'])


def get_closest_fsi_L1_file_from_selektor(search_date, band, max_t_dist):
    """ Get FSI L1 file closest to a given date

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
    """
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
    """ Get FSI L1 file to coalign with a SPICE file

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
    fsi_file : dict or None
        Dictionary containing info about the closest FSI file, if there is one.
    """
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
    fsi_file_L2 = common.EuiUtils.local_L2_path(output_dir, fsi_file_L1)
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


def dummy_stew(
        filename, output_dir, sum_wvl=False,
        windows=None, overwrite=False
        ):
    os.makedirs(output_dir, exist_ok=True)
    # filename operations
    basename = os.path.splitext(os.path.basename(filename))[0]
    if sum_wvl:
        output_fits = f'{output_dir}/{basename}_remapped_img.fits'
    else:
        output_fits = f'{output_dir}/{basename}_remapped.fits'

    if os.path.isfile(output_fits) and not overwrite:
        print(f'Aligned file exists: {output_fits}, exiting')
        return output_fits

    hdulist = fits.open(filename)
    if windows is None:
        windows = [hdu.name for hdu in hdulist
                   if hdu.is_image and (hdu.name != 'WCSDVARR')]
    for win in windows:
        hdu = hdulist[win]
        if sum_wvl:
            cube = np.squeeze(hdu.data)  # remove t axis
            spectral_window = np.any(np.isnan(cube), axis=2)
            iymin, iymax = common.SpiceUtils.vertical_edges_limits(hdu.header)
            spectral_window = spectral_window[:, iymin:iymax + 1]
            valid_columns = ~np.any(spectral_window, axis=1)
            cube = cube[valid_columns]  # columns with no NaNs
            img = np.nansum(cube, axis=0)  # Sum over wavelengths
            hdu.data = img
        hdu.update_header()
        hdu.header.add_history('dummy_stew')
        hdu.add_datasum()
        hdu.add_checksum()

    # save data
    hdulist.writeto(output_fits, overwrite=overwrite)

    return output_fits


def get_spice_image_data(filename, window):
    """ Return SPICE image data

    Parameters
    ==========
    filename : str
        Path to a FITS file containing SPICE intensity maps coaligned with
        spice_stew.
    window : str
        SPICE window name

    Returns
    =======
    img : 2D array
        Image data
    header : astropy.io.fits.Header
        FITS header
    """
    hdulist = fits.open(filename)
    try:
        hdu = hdulist[window]
    except KeyError:
        print(f"Window '{window}' not found in '{filename}'")
        return None

    img = np.squeeze(hdu.data)  # remove 1-length dimensions, ie t and wvl
    if img.ndim != 2:
        msg = f'invalid image shape {img.shape} for {filename}'
        raise ValueError(msg)

    return img, hdu.header


def get_fsi_image_data(filename):
    """ Return FSI imge data

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
    """
    hdulist = fits.open(filename)
    hdu = hdulist[-1]
    return hdu.data, hdu.header


def gen_images_to_coalign(spice_file, spice_window, fsi_file, output_dir):
    """ Generate SPICE and FSI images that can be coaligned together

    Parameters
    ==========
    spice_file : str
        Path to a FITS file containing SPICE intensity maps corrected with
        spice_stew.
    spice_window : str
        SPICE window name
    fsi_file : dict
        FSI FITS info.
    output_dir : str
        Output directory

    Returns
    =======
    spice_img : str,
    fsi_img : str
        Path to FITS images to coalign
    """
    os.makedirs(output_dir, exist_ok=True)

    basename = os.path.basename(spice_file).rstrip('_remapped_img.fits')
    new_spice_filename = f'{basename}_spice_img.fits'
    new_spice_filename = os.path.join(output_dir, new_spice_filename)
    new_fsi_filename = f'{basename}_fsi_img.fits'
    new_fsi_filename = os.path.join(output_dir, new_fsi_filename)
    plot_filename = f'{basename}.pdf'
    plot_filename = os.path.join(output_dir, plot_filename)

    if os.path.isfile(new_fsi_filename) and os.path.isfile(new_spice_filename):
        print(f'Coalign images exist for: {basename}, skipping')
        return new_spice_filename, new_fsi_filename

    spice_img_data = get_spice_image_data(
        spice_file,
        spice_window,
        )
    if spice_img_data is None:
        return None
    spice_img, spice_header = spice_img_data
    fsi_img, fsi_header = get_fsi_image_data(fsi_file)

    iymin, iymax = common.SpiceUtils.vertical_edges_limits(spice_header)
    spice_img = spice_img[iymin:iymax + 1]

    # Correct solar rotation in SPICE header
    # rotation rate (on solar sphere)
    B0 = np.deg2rad(spice_header['SOLAR_B0'])
    band = fsi_header['WAVELNTH']
    omega_car = np.deg2rad(360 / 25.38 / 86400)  # rad s-1
    omega = omega_car + diff_rot(B0, f'EIT {band}')  # rad s-1
    # helioprojective rotation rate for s/c
    Rsun = spice_header['RSUN_REF']  # m
    Dsun = spice_header['DSUN_OBS']  # m
    phi = omega * Rsun / (Dsun - Rsun)  # rad s-1
    phi = np.rad2deg(phi) * 3600  # arcsec s-1
    # time between slit positions
    t = fits.open(spice_file)['VARIABLE_KEYWORDS'].data['TIMAQUTC']
    t = np.array([parse_date(t) for t in np.squeeze(t)])
    dt = [dt.total_seconds() for dt in t[1:] - t[:-1]]  # s
    dt = - np.mean(dt)
    # CDELT correction
    DTx_old = spice_header['CDELT1']
    DTx_new = DTx_old + dt * phi
    spice_header['CDELT1'] = DTx_new
    print(f'changed SPICE CDELT1 from {DTx_old} to {DTx_new}')

    # SPICE WCS
    wcs_spice = wcs.WCS(spice_header)
    # px coordinates of SPICE cut zone
    ny_spice, nx_spice = spice_img.shape
    x_spice = np.arange(0, nx_spice)
    y_spice = np.arange(iymin, iymax + 1)
    w_spice = np.array([0])
    t_spice = np.array([1])
    assert y_spice.size <= ny_spice
    px_spice = np.meshgrid(x_spice, y_spice, w_spice, t_spice)
    px_spice = np.moveaxis(px_spice, 0, -1)
    # convert to world coordinates
    world_spice = wcs_spice.wcs_pix2world(px_spice.reshape(-1, 4), 0)
    world_spice = world_spice.reshape((y_spice.size, x_spice.size, 1, 1, 4))
    assert px_spice.shape == world_spice.shape
    Tx_spice = u.Quantity(
        world_spice[:, :, 0, 0, 0],
        wcs_spice.world_axis_units[0]
        )
    Ty_spice = u.Quantity(
        world_spice[:, :, 0, 0, 1],
        wcs_spice.world_axis_units[1]
        )
    Tx_spice = common.ang2pipi(Tx_spice).to('arcsec').value
    Ty_spice = common.ang2pipi(Ty_spice).to('arcsec').value

    # FSI WCS
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
    Tx_fsi = common.ang2pipi(Tx_fsi).to('arcsec').value
    Ty_fsi = common.ang2pipi(Ty_fsi).to('arcsec').value

    # Generate common coordinates
    common_Txy_size = 4  # arcsec
    common_Txy_pad = 0  # arcsec
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
    w_common.wcs.crpix = [Tx_common_1d.size // 2, Ty_common_1d.size // 2]
    w_common.wcs.crval = [Tx_common_1d[int(w_common.wcs.crpix[0])],
                          Ty_common_1d[int(w_common.wcs.crpix[1])]]
    w_common.wcs.ctype = ['HPLN-TAN', 'HPLT-TAN']
    w_common.wcs.cunit = ['arcsec', 'arcsec']

    # Pre-cut FSI
    corners = [
                  [Tx_common.min(), Ty_common.min()],
                  [Tx_common.min(), Ty_common.max()],
                  [Tx_common.max(), Ty_common.max()],
                  [Tx_common.max(), Ty_common.min()],
                  ] * u.arcsec
    corners_px = np.array(wcs_fsi.world_to_pixel(*corners.T)).T
    ixmin_fsi, iymin_fsi = np.floor(corners_px.min(axis=0)).astype(int)
    ixmax_fsi, iymax_fsi = np.ceil(corners_px.max(axis=0)).astype(int)
    ixmin_fsi = np.clip(ixmin_fsi - 10, 0, None)
    ixmax_fsi = np.clip(ixmax_fsi + 10, None, nx_fsi - 1)
    iymin_fsi = np.clip(iymin_fsi - 10, 0, None)
    iymin_fsi = np.clip(iymin_fsi + 10, None, ny_fsi - 1)
    fsi_img = fsi_img[iymin_fsi:iymax_fsi + 1, ixmin_fsi:ixmax_fsi + 1]
    Tx_fsi = Tx_fsi[iymin_fsi:iymax_fsi + 1, ixmin_fsi:ixmax_fsi + 1]
    Ty_fsi = Ty_fsi[iymin_fsi:iymax_fsi + 1, ixmin_fsi:ixmax_fsi + 1]

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


def plot_images(spice_img, fsi_img, wcs_common, filename):
    assert spice_img.shape == fsi_img.shape
    # assumes CROTA = 0
    ny, nx = spice_img.shape
    Tx, _ = wcs_common.pixel_to_world(np.arange(nx), [0])
    _, Ty = wcs_common.pixel_to_world([0], np.arange(ny))
    pi = u.Quantity(np.pi, 'rad')
    Tx = (Tx + pi) % (2 * pi) - pi
    Ty = (Ty + pi) % (2 * pi) - pi
    Tx = Tx.to('arcsec').value
    Ty = Ty.to('arcsec').value

    plt.clf()

    def nanptp(a):
        return np.nanmax(a) - np.nanmin(a)

    fsi_img = fsi_img - np.nanmin(fsi_img) + 0.01 * nanptp(fsi_img)
    spice_img = spice_img - np.nanmin(spice_img) + 0.01 * nanptp(spice_img)
    fsi_norm = plt.matplotlib.colors.LogNorm(
        vmin=np.nanpercentile(fsi_img, 0),
        vmax=np.nanpercentile(fsi_img, 99.9),
        )
    spice_norm = plt.matplotlib.colors.LogNorm(
        vmin=np.nanpercentile(spice_img, 0),
        vmax=np.nanpercentile(spice_img, 99.9),
        )
    ax1 = plt.subplot(121)
    papy.plot.plot_map(
        plt.gca(),
        fsi_img,
        coordinates=[Tx, Ty],
        regularity_threshold=.2,
        norm=fsi_norm,
        )
    plt.contour(
        Tx, Ty,
        spice_img,
        levels=np.nanpercentile(spice_img, [95, 99]),
        colors='w',
        linewidths=.5,
        )
    plt.contour(
        Tx, Ty,
        np.isnan(spice_img).astype(float),
        levels=[.5],
        colors='w',
        linewidths=.5,
        )
    ax2 = plt.subplot(122)
    papy.plot.plot_map(
        plt.gca(),
        spice_img,
        coordinates=[Tx, Ty],
        regularity_threshold=.2,
        norm=spice_norm,
        )
    ax1.set_xlabel('$\\theta_x$')
    ax1.set_ylabel('$\\theta_y$')
    ax2.set_xlabel('$\\theta_x$')
    ax2.set_ylabel(' ')
    ax1.set_title(' ')
    ax2.set_title(' ')
    plt.tight_layout()
    plt.savefig(filename)


def coalign_spice_fsi_images(spice_img, fsi_img, output_dir, roll=None):
    """ Coalign SPICE and FSI images

    Parameters
    ==========
    spice_img : str
    fsi_img : str
        Path to FITS images to coalign, written by `gen_images_to_coalign()`
    output_dir : str
        Output directory
    roll : float or None (default: None)
        roll in degrees, to save with yml file
    """
    os.makedirs(output_dir, exist_ok=True)

    basename = os.path.basename(spice_img)
    basename = basename.rstrip('_coalign_spice_img.fits')
    basename = os.path.join(output_dir, basename)
    plot_filename = f'{basename}_coaligned.pdf'
    yml_filename = f'{basename}_coaligned.yml'
    if os.path.isfile(yml_filename):
        print(f'Coalign results exist for: {basename}, skipping')
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
    aligned_cube = np.squeeze(
        align_images.align.align_cube(
            cube,
            np.stack([shifts, [0, 0]]),
            )
        )
    plot_images(aligned_cube[0], aligned_cube[1], w, plot_filename)

    res = dict(
        dx=float(shifts[0]),
        dy=float(shifts[1]),
        max_cc=float(max_cc),
        roll=roll,
        wcs=dict(w.to_header())
        )
    with open(yml_filename, 'w') as f:
        yaml.safe_dump(res, f, sort_keys=False)


def process_time_span(
        start_date, end_date, study_id, spec_win,
        output_dir='./output', no_stew=False
        ):

    os.makedirs(output_dir, exist_ok=True)

    print('\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
    print('STUDY:', study_id)
    print('Spectral window:', spec_win)
    print('Listing files...', end=' ')
    spice_filenames = list_spice_files(
        start_date,
        end_date,
        study_id,
        )
    print(len(spice_filenames), 'found')
    for fn in spice_filenames:
        print(fn)

    if not no_stew:
        ssp = spice_stew.SpiceSpicePointing()
    n_tot = len(spice_filenames)
    for i, spice_file in enumerate(spice_filenames):

        # Check if file exists
        if os.path.isfile(
                os.path.join(
                    output_dir, 'coalign_output',
                    f'{os.path.splitext(spice_file)[0]}_coaligned.yml'
                    )
                ):
            print(f'\nSkipping {spice_file} (existing results found)')
            continue

        spice_file = common.SpiceUtils.ias_fullpath(spice_file)
        print('\nProcessing', spice_file, f'{i}/{n_tot}')

        print('Getting closest FSI image')
        fsi_file_L1 = get_fsi_L1(
            spice_file,
            '304',
            f'{output_dir}/fsi_data',
            max_t_dist=3,
            )
        if fsi_file_L1 is None:
            print(f'No FSI image found for {spice_file}, skipping')
            continue
        fsi_file_L1 = common.EuiUtils.ias_fullpath(fsi_file_L1['filepath'])
        print(fsi_file_L1)

        if not no_stew:
            print('Correcting pointing with SPICE kernels')
            try:
                spice_file_aligned = spice_stew.correct_spice_pointing(
                    ssp,
                    spice_file,
                    f'{output_dir}/spice_stew',
                    overwrite=False,
                    plot_results=True,
                    sum_wvl=True,
                    windows=[spec_win],
                    )
            except (spiceypy.utils.exceptions.SpiceNOFRAMECONNECT, ValueError):
                print('An error occurred while running spice_stew:')
                print(traceback.format_exc())
                continue
        else:
            print('Applying dummy SPICE kernel pointing correction')
            spice_file_aligned = dummy_stew(
                spice_file,
                f'{output_dir}/spice_stew_dummy',
                sum_wvl=True,
                overwrite=True,
                windows=[spec_win],
                )
        with fits.open(spice_file) as hdul:
            roll = hdul[0].header['CROTA']

        print('Generating L2 FSI image')
        try:
            fsi_file_L2 = gen_fsi_L2(fsi_file_L1, f'{output_dir}/fsi_data')
        except FileNotFoundError:
            print(f'Could not open {fsi_file_L1}, skipping')
            continue

        print('Generating images to coalign')
        img_to_coalign = gen_images_to_coalign(
            spice_file_aligned,
            spec_win,
            fsi_file_L2,
            f'{output_dir}/coalign_input',
            )
        if img_to_coalign is None:
            print('Could not get images to coalign')
            continue
        spice_img, fsi_img = img_to_coalign

        print('Coaligning images')
        coalign_spice_fsi_images(
            spice_img,
            fsi_img,
            f'{output_dir}/coalign_output',
            roll=roll,
            )


def process_all(conf):
    for time_span in conf['time_spans']:
        process_time_span(
            time_span['start_date'],
            time_span['end_date'],
            time_span['study_id'],
            time_span['spec_win'],
            output_dir=time_span['dir'],
            no_stew=conf['processing']['no_stew'],
            )


if __name__ == '__main__':
    process_all(common.get_conf_from_cli())
