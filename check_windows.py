#!/usr/bin/env python


""" Check whether a spectral is present in files of a given study

Example:

```
check_windows.py \
    --start-date "2020-06-14" \
    --end-date "2022-05-27" \
    --study-name "SCI_SYNOPTIC_SC_SL04_60.0S_FF" \
    --spec-win "Ly-gamma-CIII group bin (1/4)"
```

"""

import argparse

from astropy.io import fits

from quick_spice_fsi_coalign import SpiceUtils, list_spice_files


if __name__ == '__main__':

    p = argparse.ArgumentParser()
    p.add_argument('--start-date', required=True,
                   help='processing start date (YYYY-MM-DD)')
    p.add_argument('--end-date', required=True,
                   help='processing end date (YYYY-MM-DD)')
    p.add_argument('--study-name',
                   help='study name')
    p.add_argument('--spec-win', required=True,
                   help='spectral window')
    args = p.parse_args()

    filenames = list_spice_files(args.start_date, args.end_date,
                                 study_name=args.study_name)

    for fn in filenames:
        hdul = fits.open(SpiceUtils.ias_fullpath(fn))
        print(
            fn,
            hdul[0].header['MISOSTUD'],
            hdul[0].header['MISOWIN'],
            args.spec_win in hdul,
            )