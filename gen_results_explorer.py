#!/usr/bin/env python
import os

import bokeh as bk
import bokeh.plotting
import bokeh.models
import pandas as pd

import common
import plot_utils


class PointingResultsExplorer:

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

    def __init__(self, conf):
        dat = pd.concat([plot_utils.get_data(conf, time_span)
                         for time_span in conf['time_spans']])
        dat = plot_utils.Filters.center(dat)
        self.conf = conf
        self.source = bk.models.ColumnDataSource(data=dat)

    def gen_explorer(self):
        os.makedirs(self.conf['plot']['dir'], exist_ok=True)

        hover_tool = bk.models.HoverTool(
            tooltips=self.TOOLTIPS_HTML,
            mode='mouse',
            )
        hover_tool.formatters = {"@{date}": "datetime"}
        tap_tool = bk.models.TapTool(
            callback=bk.models.OpenURL(url='@plot_pdf'),
            )

        bk.plotting.output_file(
            f"{self.conf['plot']['dir']}/coalign_TxTy_sc_all.html")
        p = bk.plotting.figure(
            x_axis_type="datetime",
            x_axis_label='Date',
            y_axis_label='Offset [arcsec]',
            tools='pan,box_zoom,wheel_zoom,save,reset',
            height=800,
            width=1000,
            )
        p.add_tools(hover_tool)
        p.add_tools(tap_tool)

        p.dot(
            'date', 'dx_sc', source=self.source,
            legend_label='X',
            size=20, color='#004488',
            )
        p.dot(
            'date', 'dy_sc', source=self.source,
            legend_label='Y',
            size=20,  color='#bb5566',
            )

        bk.plotting.save(p)


if __name__ == '__main__':
    PointingResultsExplorer(common.get_conf_from_cli()).gen_explorer()
