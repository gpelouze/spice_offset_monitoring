#!/usr/bin/env python

import datetime

import argparse
import glob
import os

from astropy.io import fits
from dateutil.parser import parse as parse_date
import bokeh as bk
import bokeh.plotting
import bokeh.models
import numpy as np
import pandas as pd
import yaml

from quick_spice_fsi_coalign import SpiceUtils, list_spice_files


class ThompsonData:
    t0 = parse_date('2020-06-19T23:46:38')
    tx, dx = np.array([  # h, arcsec
        [0.3532608695652171, -174.6875],
        [1.3097826086956519, -177.00892857142858],
        [2.0978260869565215, -177.05357142857144],
        [3.0652173913043472, -178.39285714285714],
        [3.847826086956522, -171.96428571428572],
        [4.815217391304348, -167.58928571428572],
        [5.5978260869565215, -175.04464285714286],
        [6.570652173913043, -172.45535714285714],
        ]).T
    ty, dy = np.array([  # h, arcsec
        [0.3453009503695881, -208.8056460369164],
        [1.3167898627243924, -211.5418023887079],
        [2.0929250263991555, -205.33116178067317],
        [3.0644139387539595, -207.37242128121605],
        [3.851108764519536, -217.53528773072745],
        [4.817317845828935, -214.4516829533116],
        [5.604012671594509, -209.2833876221498],
        [6.564941921858502, -210.62975027144407],
        ]).T

    def __init__(self):
        self.date = [self.t0 + datetime.timedelta(hours=np.mean([tx, ty]))
                  for tx, ty in zip(self.tx, self.ty)]
        self.dx = self.dx + 173
        self.dy = self.dy + 205

    def to_df(self):
        return pd.DataFrame(dict(date=self.date, dx=self.dx, dy=self.dy))


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
    p.add_argument('--output-dir', default='./output',
                   help='data directory')
    args = p.parse_args()

    print('Listing files')
    yml_fnames = glob.glob(f'{args.output_dir}/coalign_output/*_coaligned.yml')
    yml_fnames = list(sorted(yml_fnames))

    dat = []
    for yml_fname in yml_fnames:
        if os.path.isfile(yml_fname):
            spice_fname = os.path.basename(yml_fname).rstrip('_coaligned.yml')
            print('opening', spice_fname)
            with open(yml_fname, 'r') as f:
                res = yaml.safe_load(f)
            res['date'] = SpiceUtils.filename_to_date(f'{spice_fname}.fits')
            res['plot_pdf'] = f'coalign_output/{spice_fname}_coaligned.pdf'
            res['plot'] = f'coalign_output/{spice_fname}_coaligned.jpg'
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

    bk.plotting.output_file(f'{args.output_dir}/view.html')
    p = bk.plotting.figure(
        x_axis_type="datetime",
        x_axis_label='Date',
        y_axis_label='SPICE-FSI WCS offset [arcsec]',
        tools='pan,box_zoom,wheel_zoom,save,reset',
        plot_height=800,
        plot_width=1000,
        )
    p.add_tools(hover_tool)
    p.add_tools(tap_tool)

    p.dot('date', 'dx', size=20, color='#000000', source=source, legend_label='X')
    p.dot('date', 'dy', size=20, color='#ff0000', source=source, legend_label='Y')
    if dat['date'].min() <= datetime.datetime(2020, 6, 21):
        source_th = bk.models.ColumnDataSource(data=ThompsonData().to_df())
        p.cross('date', 'dx', size=20, color='#000000', source=source_th, legend_label='X (W. Thompson)')
        p.cross('date', 'dy', size=20, color='#ff0000', source=source_th, legend_label='Y (W. Thompson)')
    bk.plotting.save(p)
