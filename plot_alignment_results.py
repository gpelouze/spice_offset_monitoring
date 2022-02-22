#!/usr/bin/env python

import argparse
import glob
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
    p.add_argument('--output-dir', default='./output',
                   help='data directory')
    args = p.parse_args()

    print('Listing files')
    yml_fnames = glob.glob(f'{args.output_dir}/coalign_output/*_coaligned.yml')

    dat = []
    for yml_fname in yml_fnames:
        if os.path.isfile(yml_fname):
            spice_fname = os.path.basename(yml_fname).rstrip('_coaligned.yml')
            print('opening', spice_fname)
            with open(yml_fname, 'r') as f:
                res = yaml.safe_load(f)
            res['date'] = SpiceUtils.filename_to_date(f'{spice_fname}.fits')
            res['plot'] = f'coalign_output/{spice_fname}_coaligned.pdf'
            wcs = res.pop('wcs')
            res.update(wcs)
            dat.append(res)
        else:
            print('skipping', spice_fname)

    dat = list_of_dict_to_dict_of_arr(dat)

    # add columns
    dat['dr'] = np.sqrt(np.array(dat['dx'])**2 + np.array(dat['dy'])**2)
    r = np.deg2rad(dat['roll'])
    dat['dx_sc'] = dat['dx'] * np.cos(r) - dat['dy'] * np.sin(r)
    dat['dy_sc'] = dat['dx'] * np.sin(r) + dat['dy'] * np.cos(r)

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
    plt.savefig(f'{args.output_dir}/coalign_cc_dr.pdf')

    plt.clf()
    plt.plot(dat['date'], dat['dx'], 'ko', ms=3, label='$\\theta_x$')
    plt.plot(dat['date'], dat['dy'], 'rs', ms=3, label='$\\theta_y$')
    plt.xlabel('Date')
    plt.ylabel('SPICE-FSI WCS offset [arcsec]')
    plt.legend()
    plt.gcf().autofmt_xdate()
    plt.savefig(f'{args.output_dir}/coalign_TxTy_hproj.pdf')
    plt.xlim(parse_date('2021-12-02').toordinal(), None)

    plt.clf()
    plt.plot(dat['date'], dat['dx_sc'], 'ko', ms=3, label='$\\theta_x$ s/c')
    plt.plot(dat['date'], dat['dy_sc'], 'rs', ms=3, label='$\\theta_y$ s/c')
    plt.xlabel('Date')
    plt.ylabel('SPICE-FSI WCS offset [arcsec]')
    plt.legend()
    plt.gcf().autofmt_xdate()
    plt.savefig(f'{args.output_dir}/coalign_TxTy_sc.pdf')
    plt.xlim(parse_date('2021-12-02').toordinal(), None)
