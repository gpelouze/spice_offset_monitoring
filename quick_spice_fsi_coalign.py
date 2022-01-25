#!/usr/bin/env python

import argparse

import spice_stew


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
    pass  # TODO


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
            spice_file,
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
