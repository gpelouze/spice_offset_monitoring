#!/usr/bin/env python

import argparse
import os

from astropy.io import fits
from dateutil.parser import parse as parse_date
import matplotlib.pyplot as plt
import yaml

from quick_spice_fsi_coalign import SpiceUtils, list_spice_files


def list_of_dict_to_dict_of_list(l):
    keys = l[0].keys()
    d_out = {k: [] for k in keys}
    for d in l:
        if d.keys() != keys:
            raise ValueError('mismatched keys')
        for k, v in d.items():
            d_out[k].append(v)
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

    results = []
    for spice_fname in spice_fnames:
        base, _ = os.path.splitext(spice_fname)
        yml_fname = f'{args.output_dir}/{base}_coaligned.yml'
        if os.path.isfile(yml_fname):
            print('opening', spice_fname)
            with open(yml_fname, 'r') as f:
                res = yaml.safe_load(f)
            res['date'] = SpiceUtils.filename_to_date(spice_fname)
            wcs = res.pop('wcs')
            res.update(wcs)
            results.append(res)
        else:
            print('skipping', spice_fname)

    results = list_of_dict_to_dict_of_list(results)

    plt.clf()
    plt.plot(results['date'], results['dx'], 'ko', ms=3, label='$\\theta_x$')
    plt.plot(results['date'], results['dy'], 'rs', ms=3, label='$\\theta_y$')
    plt.xlabel('Date')
    plt.ylabel('SPICE-FSI offset [arcsec]')
    plt.legend()
    plt.gcf().autofmt_xdate()
    # formatter = plt.matplotlib.dates.DateFormatter("%Y-%m")
    # plt.gca().xaxis.set_major_formatter(formatter)
    plt.tight_layout()
    plt.savefig('./output/coalign.pdf')
    plt.xlim(parse_date('2021-12-02').toordinal(), None)
    plt.ylim(-15, 5)
    plt.tight_layout()
    plt.savefig('./output/coalign_zoom.pdf')
