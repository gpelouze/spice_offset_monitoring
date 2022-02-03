#!/usr/bin/env python

import argparse
import datetime
import os
import re

from astropy.io import fits
from dateutil.parser import parse as parse_date
import pandas as pd
import yaml

from papy.sol.data.solo_eui import EUISelektorClient
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
                        parse_dates=date_columns, warn_bad_lines=True)
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


def get_closest_fsi_image(spice_file, band, max_t_dist=6):
    ''' Get FSI image closest to a given date

    Parameters
    ==========
    spice_file : str
        Path to a SPICE FITS file
    band : str ('174' or '304')
        Instrument band
    max_t_dist : float (default: 6)
        Maximum time distance in hour. If no file exists within [date -
        max_t_dist, date + max_t_dist], return None.

    Returns
    =======
    filename : str or None
        Closest FSI FITS, if there is one.
    '''
    # Get SPICE acquisition midpoint date
    spice_hdu = fits.open(spice_file)
    search_date_mid = parse_date(spice_hdu[0].header['DATE-AVG'])

    # Compute date search window
    max_t_dist = datetime.timedelta(hours=max_t_dist)
    search_date_min = search_date_mid - max_t_dist
    search_date_max = search_date_mid + max_t_dist

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
    # - 'before', in descending from search_date_mid to search_date_min,
    #   yielding the Nth images before search_date_mid
    # - 'after', in ascending from search_date_mid to search_date_max,
    #   yielding the Nth images after search_date_mid
    before_params = {
        'order[]': 'DESC',
        'date_begin_start': search_date_min.isoformat().split('T')[0],
        'date_begin_start_hour': search_date_min.hour,
        'date_begin_start_minute': search_date_min.minute,
        'date_begin_end': search_date_mid.isoformat().split('T')[0],
        'date_begin_end_hour': search_date_mid.hour,
        'date_begin_end_minute': search_date_mid.minute,
        }
    after_params = {
        'order[]': 'ASC',
        'date_begin_start': search_date_mid.isoformat().split('T')[0],
        'date_begin_start_hour': search_date_mid.hour,
        'date_begin_start_minute': search_date_mid.minute,
        'date_begin_end': search_date_max.isoformat().split('T')[0],
        'date_begin_end_hour': search_date_max.hour,
        'date_begin_end_minute': search_date_max.minute,
        }
    before_params.update(base_params)
    after_params.update(base_params)
    res_before = eui_client.search(before_params)
    res_after = eui_client.search(after_params)

    # keep images just before and just after
    if res_before is not None:
        res_before = res_before.iloc[0]
    if res_after is not None:
        res_after = res_after.iloc[0]

    # return closest image, or None
    if res_before is None:
        return res_after  # may be None
    elif res_after is None:
        return res_before  # may be None
    else:
        date_before = parse_date(res_before['date-beg'])
        date_after = parse_date(res_after['date-beg'])
        dt_before = abs((search_date_mid - date_before).total_seconds())
        dt_after = abs((search_date_mid - date_after).total_seconds())
        if dt_before < dt_after:
            return res_before
        else:
            return res_after


def coalign_spice_fsi_images(spice_file, fsi_file):
    ''' Coalign SPICE and FSI images

    Parameters
    ==========
    spice_file : str
        Path to a FITS file containing SPICE intensity maps coaligned with
        spice_stew.
    fsi_file : str
        Path to a FSI FITS file.

    Returns
    =======
    coalign : dict
        Coalignment results
    '''
    pass  # TODO
    # return coalign


def save_coalign_results(coalign):
    ''' Save coalignment results

    Parameters
    ==========
    coalign : dict
        Coalignment results
    '''
    pass  # TODO


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

    # List SPICE files to process
    spice_filenames = list_spice_files(
        args.start_date,
        args.end_date,
        study_name='SCI_SYNOPTIC_SC_SL04_60.0S_FF',
        )

    ssp = spice_stew.SpiceSpicePointing()
    for spice_file in spice_filenames:
        spice_file = SpiceUtils.ias_fullpath(spice_file)
        # Correct pointing with SPICE kernels
        spice_file_aligned = spice_stew.correct_spice_pointing(
            ssp,
            spice_file,
            args.output_dir,
            overwrite=False,
            plot_results=True,
            sum_wvl=True,
            )

        # Get closest FSI image
        fsi_file = get_closest_fsi_image(
            spice_file,
            '304',
            max_t_dist=3,
            )

        # Coalign SPICE and FSI image
        coalign = coalign_spice_fsi_images(spice_file_aligned, fsi_file)

        # Write output
        save_coalign_results(coalign)
