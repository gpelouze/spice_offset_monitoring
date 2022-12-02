#!/usr/bin/env python

import offset_determination
import results_explorer
import figures
import utils


def pipeline(conf):
    offset_determination.process_all(conf)
    results_explorer.gen_explorer(conf)
    figures.plot_all(conf)


def main():
    pipeline(utils.get_conf_from_cli())


if __name__ == '__main__':
    main()
