import argparse
import os
import re

from dateutil.parser import parse as parse_date
import astropy.units as u
import numpy as np
import pandas as pd
import yaml


class SpiceUtils:
    re_spice_L123_filename = re.compile(
        r'''
        solo
        _(?P<level>L[123])
        (?P<level_extension>r?(_quicklook)?)
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
        re.VERBOSE,
        )

    @staticmethod
    def read_spice_uio_catalog():
        """
        Read UiO text table SPICE FITS files catalog
        http://astro-sdc-db.uio.no/vol/spice/fits/spice_catalog.txt

        Return
        ------
        cat : pandas.DataFrame
            Table

        Example queries that can be done on the result:

        * `df[(df.LEVEL == "L2") &
              (df["DATE-BEG"] >= "2020-11-17") &
              (df["DATE-BEG"] < "2020-11-18") &
              (df.XPOSURE > 60.)
              ]`
        * `df[(df.LEVEL == "L2") &
              (df.STUDYDES == "Standard dark for cruise phase")
              ]`

        Source: https://spice-wiki.ias.u-psud.fr/doku.php/\
                data:data_analysis_manual:read_catalog_python
        """
        cat_file = os.path.join(
            os.getenv('SOLO_ARCHIVE', '/archive/SOLAR-ORBITER/'),
            'SPICE/fits/spice_catalog.txt'
            )
        columns = list(pd.read_csv(cat_file, nrows=0).keys())
        date_columns = ['DATE-BEG', 'DATE', 'TIMAQUTC']
        df = pd.read_table(
            cat_file, skiprows=1, names=columns, na_values="MISSING",
            parse_dates=date_columns, low_memory=False
            )
        df.LEVEL = df.LEVEL.apply(lambda string: string.strip())
        df.STUDYTYP = df.STUDYTYP.apply(lambda string: string.strip())
        return df

    @staticmethod
    def parse_filename(filename):
        m = SpiceUtils.re_spice_L123_filename.match(filename)
        if m is None:
            raise ValueError(f'could not parse SPICE filename: {filename}')
        return m.groupdict()

    @staticmethod
    def filename_to_date(filename):
        d = SpiceUtils.parse_filename(filename)
        return parse_date(d['time'])

    @staticmethod
    def ias_fullpath(filename):
        d = SpiceUtils.parse_filename(filename)
        date = parse_date(d['time'])

        fullpath = os.path.join(
            os.getenv('SOLO_ARCHIVE', '/archive/SOLAR-ORBITER/'),
            'SPICE/fits/',
            'level' + d['level'].lstrip('L'),
            f'{date.year:04d}/{date.month:02d}/{date.day:02d}',
            filename
            )

        return fullpath

    @staticmethod
    def slit_px(header):
        """ Compute the first and last pixel of the slit from a FITS header """
        ybin = header['NBIN2']
        h_detector = 1024 / ybin
        if header['DETECTOR'] == 'SW':
            h_slit = 600 / ybin
        elif header['DETECTOR'] == 'LW':
            h_slit = 626 / ybin
        else:
            raise ValueError(f"unknown detector: {header['DETECTOR']}")
        slit_beg = (h_detector - h_slit) / 2
        slit_end = h_detector - slit_beg
        slit_beg = slit_beg - header['PXBEG2'] / ybin + 1
        slit_end = slit_end - header['PXBEG2'] / ybin + 1
        slit_beg = int(np.ceil(slit_beg))
        slit_end = int(np.floor(slit_end))
        return slit_beg, slit_end

    @staticmethod
    def vertical_edges_limits(header):
        iymin, iymax = SpiceUtils.slit_px(header)
        iymin += int(20 / header['NBIN2'])
        iymax -= int(20 / header['NBIN2'])
        return iymin, iymax


class EuiUtils:
    @staticmethod
    def ias_fullpath(rob_fullpath):
        p = rob_fullpath.lstrip('/data/solo-eui/internal/')
        p = '/archive/SOLAR-ORBITER/EUI/data_internal/' + p
        return p

    @staticmethod
    def local_L2_path(output_dir, fsi_file_L1):
        base = os.path.basename(fsi_file_L1)
        base = base.replace('L1', 'L2')
        return os.path.join(output_dir, base)


class Config(dict):
    def __init__(self, filename):
        with open(filename, 'r') as f:
            conf = yaml.safe_load(f)
        super().__init__(conf)


def ang2pipi(ang):
    """ put angle between ]-180, +180] deg """
    pi = u.Quantity(180, 'deg')
    return - ((- ang + pi) % (2 * pi) - pi)


def get_conf_from_cli():
    p = argparse.ArgumentParser()
    p.add_argument(
        'conf_file', metavar='yaml_file',
        help='configuration file',
        )
    args = p.parse_args()
    return Config(args.conf_file)
