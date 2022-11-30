#!/usr/bin/env python

import common
import quick_spice_fsi_coalign
import gen_results_explorer
import gen_plots


def pipeline(conf):
    quick_spice_fsi_coalign.process_all(conf)
    gen_results_explorer.PointingResultsExplorer(conf).gen_explorer()
    gen_plots.plot_all(conf)


if __name__ == '__main__':
    pipeline(common.get_conf_from_cli())
