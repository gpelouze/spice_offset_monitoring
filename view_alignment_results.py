#!/usr/bin/env python

import datetime

import argparse
import os

from astropy.io import fits
import bokeh as bk
import bokeh.plotting
import bokeh.models
import numpy as np
import pandas as pd
import yaml

from quick_spice_fsi_coalign import SpiceUtils, list_spice_files


TOOLTIPS_HTML = """
    <div>
        <div>
            <img
                src="@plot" alt="@plot"
                height="480"
                style="float: left; margin: 0px 15px 15px 0px;"
            ></img>
        </div>
        <div>
            <span><b>dx</b>: @dx</span><br/>
            <span><b>dy</b>: @dy</span><br/>
            <span><b>dr</b>: @dr</span><br/>
            <span><b>max_cc</b>: @max_cc</span><br/>
            <span><b>date</b>: @date</span><br/>
            <span><b>CRPIX1</b>: @CRPIX1</span><br/>
            <span><b>CRPIX2</b>: @CRPIX2</span><br/>
            <span><b>CRVAL1</b>: @CRVAL1</span><br/>
            <span><b>CRVAL2</b>: @CRVAL2</span><br/>
        </div>
    </div>
    """


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
        url = 'coalign_output'
        pdf_fname = f'{url}/{base}_coaligned.pdf'
        img_fname = f'{url}/{base}_coaligned.jpg'
        if os.path.isfile(yml_fname):
            print('opening', spice_fname)
            with open(yml_fname, 'r') as f:
                res = yaml.safe_load(f)
            res['date'] = SpiceUtils.filename_to_date(spice_fname)
            res['plot_pdf'] = pdf_fname
            res['plot'] = img_fname
            wcs = res.pop('wcs')
            res.update(wcs)
            dat.append(res)
        else:
            print('skipping', spice_fname)
    dat = pd.DataFrame(dat)

    # add columns
    dat['dr'] = np.sqrt(np.array(dat['dx'])**2 + np.array(dat['dy'])**2)

    # filter data
    m = (dat['max_cc'] > 0.2)
    m &= (dat['dr'] < 50)
    dat = dat[m]

    source = bk.models.ColumnDataSource(data=dat)

    hover_tool = bk.models.HoverTool(
        tooltips=TOOLTIPS_HTML,
        mode='mouse',
        )
    hover_tool.formatters = {"@{date}": "datetime"}
    tap_tool = bk.models.TapTool(
        callback=bk.models.OpenURL(url='@plot'),
        )

    bk.plotting.output_file('output/view.html')
    p = bk.plotting.figure(
        x_axis_label='',
        y_axis_label='',
        tools='pan,box_zoom,wheel_zoom,save,reset',
        plot_height=800,
        plot_width=1000,
        )
    p.add_tools(hover_tool)
    p.add_tools(tap_tool)

    p.dot('date', 'dx', size=20, color='#000000', source=source)
    p.dot('date', 'dy', size=20, color='#ff0000', source=source)
    bk.plotting.save(p)
