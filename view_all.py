#!/usr/bin/env python

import bokeh as bk
import bokeh.plotting
import bokeh.models
import pandas as pd

import plot_utils

TOOLTIPS_HTML = """
    <div>
        <div>
            <img
                src="@plot_jpg" alt="@plot_jpg"
                height="480"
                style="float: left; margin: 0px 15px 15px 0px;"
            ></img>
        </div>
        <div>
            <span><b>Δx</b>: @dx_sc</span><br/>
            <span><b>Δy</b>: @dy_sc</span><br/>
            <span><b>CC</b>: @max_cc</span><br/>
        </div>
    </div>
    """

if __name__ == '__main__':

    dat = pd.concat([
        plot_utils.get_data('Lyγ CIII SYN', 'output/SYN_Lyg_CIII_new_syn'),
        plot_utils.get_data('Lyγ CIII RSW', 'output/SYN_Lyg_CIII_new_rsw'),
        ])
    dat = plot_utils.Filters.center(dat)
    source = bk.models.ColumnDataSource(data=dat)

    hover_tool = bk.models.HoverTool(
        tooltips=TOOLTIPS_HTML,
        mode='mouse',
        )
    hover_tool.formatters = {"@{date}": "datetime"}
    tap_tool = bk.models.TapTool(
        callback=bk.models.OpenURL(url='@plot_pdf'),
        )

    bk.plotting.output_file('output/coalign_TxTy_sc_all.html')
    p = bk.plotting.figure(
        x_axis_type="datetime",
        x_axis_label='Date',
        y_axis_label='SPICE-FSI WCS offset [arcsec]',
        tools='pan,box_zoom,wheel_zoom,save,reset',
        height=800,
        width=1000,
        )
    p.add_tools(hover_tool)
    p.add_tools(tap_tool)

    p.dot('date', 'dx_sc', size=20, color='#004488', source=source, legend_label='X')
    p.dot('date', 'dy_sc', size=20, color='#bb5566', source=source, legend_label='Y')

    bk.plotting.save(p)
