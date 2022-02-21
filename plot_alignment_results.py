#!/usr/bin/env python

import argparse
import os

from astropy.io import fits
from dateutil.parser import parse as parse_date
import matplotlib.pyplot as plt
import numpy as np
import yaml

from quick_spice_fsi_coalign import SpiceUtils, list_spice_files


def list_of_dict_to_dict_of_arr(l):
    keys = l[0].keys()
    d_out = {k: [] for k in keys}
    for d in l:
        if d.keys() != keys:
            raise ValueError('mismatched keys')
        for k, v in d.items():
            d_out[k].append(v)
    d_out = {k: np.array(v) for k, v in d_out.items()}
    return d_out


if __name__ == '__main__':

    p = argparse.ArgumentParser()
    p.add_argument('--start-date',
                   help='processing start date (YYYY-MM-DD)')
    p.add_argument('--end-date',
                   help='processing end date (YYYY-MM-DD)')
    p.add_argument('--output-dir', default='./output/coalign_output',
                   help='output directory')
    args = p.parse_args()

    print('Listing SPICE files')
    spice_fnames = list_spice_files(
        args.start_date,
        args.end_date,
        study_name='SCI_SYNOPTIC_SC_SL04_60.0S_FF',
        )

    dat = []
    for spice_fname in spice_fnames:
        base, _ = os.path.splitext(spice_fname)
        yml_fname = f'{args.output_dir}/{base}_coaligned.yml'
        pdf_fname = f'./coalign_output/{base}_coaligned.pdf'
        if os.path.isfile(yml_fname):
            print('opening', spice_fname)
            with open(yml_fname, 'r') as f:
                res = yaml.safe_load(f)
            res['date'] = SpiceUtils.filename_to_date(spice_fname)
            res['plot'] = pdf_fname
            wcs = res.pop('wcs')
            res.update(wcs)
            dat.append(res)
        else:
            print('skipping', spice_fname)

    dat = list_of_dict_to_dict_of_arr(dat)

    # add columns
    dat['dr'] = np.sqrt(np.array(dat['dx'])**2 + np.array(dat['dy'])**2)

    # filter
    m = (dat['max_cc'] > 0.2)
    m &= (dat['dr'] < 50)
    dat_all = dat.copy()
    dat_discard = {k: v[~m] for k, v in dat.items()}
    dat = {k: v[m] for k, v in dat.items()}

    plt.clf()
    plt.plot(dat['dr'], dat['max_cc'], 'ko', ms=3)
    plt.plot(dat_discard['dr'], dat_discard['max_cc'], 'o', color='gray', ms=3)
    plt.xlabel('SPICE-FSI centers distance [arcsec]')
    plt.ylabel('max(CC)')
    plt.savefig('output/coalign_cc_dr.pdf')

    plt.clf()
    plt.plot(dat['date'], dat['dx'], 'ko', ms=3, label='$\\theta_x$')
    plt.plot(dat['date'], dat['dy'], 'rs', ms=3, label='$\\theta_y$')
    plt.xlabel('Date')
    plt.ylabel('SPICE-FSI offset [arcsec]')
    plt.legend()
    plt.gcf().autofmt_xdate()
    # formatter = plt.matplotlib.dates.DateFormatter("%Y-%m")
    # plt.gca().xaxis.set_major_formatter(formatter)
    plt.savefig('./output/coalign.pdf')
    plt.xlim(parse_date('2021-12-02').toordinal(), None)