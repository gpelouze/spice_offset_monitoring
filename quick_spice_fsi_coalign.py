#!/usr/bin/env python

import argparse
import os
import re

from dateutil.parser import parse as parse_date
import pandas as pd
import yaml

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


def get_closest_fsi_image(date, band, max_t_dist=6):
    ''' Get FSI image closest to a given date

    Parameters
    ==========
    date : str (YYYY-MM-DD)
        Query date
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
    pass  # TODO


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
    return coalign


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
        # Correct pointing with SPICE kernels
        spice_file_aligned = spice_stew.correct_spice_pointing(
            ssp,
            SpiceUtils.ias_fullpath(spice_file),
            args.output_dir,
            overwrite=False,
            plot_results=True,
            sum_wvl=True,
            )

        # Get closest FSI image
        fsi_file = get_closest_fsi_image(spice_file)

        # Coalign SPICE and FSI image
        coalign = coalign_spice_fsi_images(spice_file_aligned, fsi_file)

        # Write output
        save_coalign_results(coalign)
