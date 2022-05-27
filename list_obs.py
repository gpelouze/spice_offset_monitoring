#!/usr/bin/env python

import argparse

from astropy.io import fits

from quick_spice_fsi_coalign import SpiceUtils


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--start-date', required=True,
                   help='processing start date (YYYY-MM-DD)')
    p.add_argument('--end-date', required=True,
                   help='processing end date (YYYY-MM-DD)')
    p.add_argument('--study-id', nargs='+',
                   help='study ID in MISO')
    p.add_argument('--spec-win', required=True,
                   help='spectral window')
    p.add_argument('--no-stew', action='store_true',
                   help='skip jitter correction using spice_stew')
    p.add_argument('--output-dir', default='./output',
                   help='output directory')
    p.add_argument('-v', '--verbose', action='store_true',
                   help='print detailled info')
    args = p.parse_args()

    cat = SpiceUtils.read_spice_uio_catalog()
    base_filter = (
        (cat['DATE-BEG'] > args.start_date)
        & (cat['DATE-BEG'] <= args.end_date)
        & (cat['LEVEL'] == 'L2')
        )

    for study_id in args.study_id:
        study_filter = (cat['MISOSTUD'] == study_id)
        res = cat[base_filter & study_filter]

        if len(res) == 0:
            print(study_id, 'none')
            continue

        study = res['STUDY']
        assert len(set(study)) == 1, f'multiple studies: {study}'
        study = study.iloc[0]

        print(
            study_id, study,
            ' ',
            res['DATE-BEG'].min().strftime('%Y-%m-%d'),
            'to',
            res['DATE-BEG'].max().strftime('%Y-%m-%d'),
            ' ',
            len(res['FILENAME']), 'files'
            )

        if args.verbose:

            print('Spectral windows:')
            sample_fits = SpiceUtils.ias_fullpath(res['FILENAME'].iloc[0])
            with fits.open(sample_fits) as hdulist:
                for hdu in hdulist[:-1]:
                    print('   ', hdu.name)

            print('Files:')
            for fn in res['FILENAME']:
                print('   ', fn)

            print()
