#!/usr/bin/env python
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import plot_utils


def get_legend_handles(df):
    handles = [
        plt.Line2D(
            [], [],
            label='X',
            color='C0',
            marker=plot_utils.marker, mew=0, ms=12, ls=''
            ),
        plt.Line2D(
            [], [],
            label='Y',
            color='C1',
            marker=plot_utils.marker, mew=0, ms=12, ls=''
            ),
        ]
    # add unique combinations of data_source and marker to the legend
    legend_elements = set([
        tuple(v)
        for v in df[['data_source', 'plot_marker']].values
        ])
    for label, marker in legend_elements:
        handles.append(plt.Line2D(
            [], [],
            label=label,
            color='gray',
            marker=marker, fillstyle='none', ls='',
            ))
    return handles


def add_fit(ax, z, dx, dy, fit_func):
    fit_funcs = []
    for i, d in enumerate([dx, dy]):
        if hasattr(fit_func, '__iter__'):
            f = fit_func[i](z, d)
        else:
            f = fit_func(z, d)
        f.fit()
        d_label = {0: 'X', 1: 'Y'}[i]
        fit_label = f.label(f'\\Delta {d_label}', 'z')
        ax.plot(
            f.xopt(), f.yopt(),
            color=f'C{i}', ls='-', lw=1,
            label=fit_label,
            )
        print(fit_label)
        fit_funcs.append(f)
    ax.add_artist(
        plt.legend(
            loc='best',
            fontsize=10,
            )
        )
    return fit_funcs


def plot_dr_cc(df, df_filtered, filename):
    plt.clf()
    ax = plt.gca()
    ax.plot(df['dr'], df['max_cc'], 'o', color='gray', ms=3)
    ax.plot(df_filtered['dr'], df_filtered['max_cc'], 'ko', ms=3)
    ax.set_title('SPICE offset', loc='left')
    ax.set_xlabel('Absolute offset correction [arcsec]')
    ax.set_ylabel('Maximum cross-correlation')
    ax.set_ylim(0, 1)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.tight_layout()
    plt.savefig(filename)


def plot_residuals(df, x_key, x_label, fit_funcs, filename=None):
    plt.clf()
    ax = plt.gca()
    ax.set_title('Residuals', loc='left')
    ax.set_xlabel(x_label)
    ax.set_ylabel('Offset [arcsec]')
    # Points
    for _, r in df.iterrows():
        kw = dict(
            ms=3,
            fillstyle='none',
            ls='',
            marker=r.plot_marker,
            )
        for i, y in enumerate([r['dx_sc'], r['dy_sc']]):
            x = r[x_key]
            f = fit_funcs[i]
            res = y - f(x, *f.popt)
            ax.plot(x, res, color=f'C{i}', **kw)
    # Intervals
    spice_psf_fwhm = 8
    for i, y in enumerate([df['dx_sc'], df['dy_sc']]):
        x = df[x_key]
        f = fit_funcs[i]
        res = y - f(x, *f.popt)
        xlim = ax.get_xlim()
        m = np.mean(res)
        fwhm = np.sqrt(8 * np.log(2)) * np.std(res)
        xy = {0: 'X', 1: 'Y'}[i]
        label = f'${xy}$: {fwhm:.1f}″'
        print('Residuals', label)
        # Percentage within 1 SPICE PSF
        p = np.mean(np.abs(res < spice_psf_fwhm / 2))
        print(f'{p:.2%} of points within 1 SPICE PSF FWHM')
        # Plot regions
        x_ = np.array(xlim)
        m = np.array([m] * 2)
        fwhm = np.array([fwhm] * 2)
        ax.fill_between(
            x_, m + fwhm / 2, m - fwhm / 2,
            color=f'C{i}', alpha=0.2, ec=None,
            label=label,
            )
        ax.set_xlim(*xlim)
    ax.axhline(
        spice_psf_fwhm / 2,
        color=f'k', ls='--', lw=1,
        label=f'SPICE PSF: {spice_psf_fwhm}″',
        )
    ax.axhline(
        -spice_psf_fwhm / 2,
        color=f'k', ls='--', lw=1,
        )
    ax.legend(
        title='FWHMs:', title_fontsize=10, alignment='left',
        loc='lower right', bbox_to_anchor=(1, 1),
        ncol=3,
        fancybox=False,
        fontsize=10,
        )
    ax.set_ylim(-20, +20)
    # remove axes frame
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.tight_layout()
    plt.savefig(filename)


def plot_offsets(
        df, x_key, x_label, filename, date=False,
        title='SPICE offset',
        fit_func=None, filename_residuals=None
        ):
    plt.clf()
    ax = plt.gca()
    ax.set_title(title, loc='left')
    # All data
    for _, r in df.iterrows():
        kw = dict(
            ms=3,
            fillstyle='none',
            ls='',
            marker=r.plot_marker,
            )
        ax.plot(r[x_key], r['dx_sc'], color='C0', **kw)
        ax.plot(r[x_key], r['dy_sc'], color='C1', **kw)
    # fit
    if fit_func is not None:
        print('\nModel for', x_label)
        fit_funcs = add_fit(ax, df[x_key], df['dx_sc'], df['dy_sc'], fit_func)
    # labels and legend
    ax.set_xlabel(x_label)
    ax.set_ylabel('Offset [arcsec]')
    ax.legend(
        handles=get_legend_handles(df),
        ncol=4,
        loc='lower right',
        bbox_to_anchor=(1, 1),
        fancybox=False,
        fontsize=10,
        )
    # remove axes frame
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    # date
    if date:
        ax.xaxis.set_major_locator(
            mpl.dates.MonthLocator(bymonth=[1, 4, 7, 10])
            )
        ax.xaxis.set_minor_locator(
            mpl.dates.MonthLocator()
            )
        ax.xaxis.set_major_formatter(
            mpl.dates.DateFormatter("%b %y")
            )
    plt.tight_layout()
    plt.savefig(filename)

    # plot residuals
    if (fit_func is not None) and (filename_residuals is not None):
        plot_residuals(
            df, x_key, x_label, fit_funcs,
            filename=filename_residuals
            )


def gen_figures(conf):
    dat = pd.concat([plot_utils.get_data(conf, time_span)
                     for time_span in conf['time_spans']])
    dat_filtered = plot_utils.Filters.center(dat)

    os.makedirs(conf['plot']['dir'], exist_ok=True)

    plot_dr_cc(
        dat, dat_filtered,
        f"{conf['plot']['dir']}/coalign_cc_dr.pdf",
        )

    plot_offsets(
        dat_filtered, 'date', 'Date',
        f"{conf['plot']['dir']}/coalign_TxTy_sc_all.pdf", date=True
        )
    plot_offsets(
        dat_filtered, 'DSUN_AU', 'Solar distance [au]',
        f"{conf['plot']['dir']}/coalign_TxTy_sc_all_dsun.pdf",
        fit_func=[plot_utils.FitFunctions.Inverse,
                  plot_utils.FitFunctions.Constant],
        filename_residuals=(
            f"{conf['plot']['dir']}/"
            f"coalign_TxTy_sc_all_dsun_residuals.pdf"),
        )
    plot_offsets(
        dat_filtered, 'CROTA', 'Roll angle [°]',
        f"{conf['plot']['dir']}/coalign_TxTy_sc_all_CROTA.pdf",
        )

    # Detector temperature
    for T_key in ['T_SW', 'T_LW']:
        plot_offsets(
            dat_filtered, T_key, f'{T_key} [°C]',
            f"{conf['plot']['dir']}/coalign_TxTy_sc_all_{T_key}.pdf",
            )

    # Fit on grating and focus mechanism temperatures
    for T_key in ['T_GRAT', 'T_FOCUS']:
        plot_offsets(
            dat_filtered, T_key, f'{T_key} [°C]',
            f"{conf['plot']['dir']}/coalign_TxTy_sc_all_{T_key}.pdf",
            fit_func=[plot_utils.FitFunctions.Linear,
                      plot_utils.FitFunctions.Constant],
            filename_residuals=(
                f"{conf['plot']['dir']}/"
                f"coalign_TxTy_sc_all_{T_key}_residuals.pdf"),
            )

    # Recommended correction
    plot_offsets(
        dat_filtered, 'T_GRAT', f'Grating temperature [°C]',
        f"{conf['plot']['dir']}/coalign_TxTy_sc_all_RECOMM.pdf",
        fit_func=[plot_utils.FitFunctions.RecommX,
                  plot_utils.FitFunctions.RecommY],
        filename_residuals=(
            f"{conf['plot']['dir']}/"
            f"coalign_TxTy_sc_all_RECOMM_residuals.pdf"),
        )

    # Recommended offset for planning
    plot_offsets(
        dat_filtered, 'DSUN_AU', f'Solar distance [au]',
        f"{conf['plot']['dir']}/coalign_TxTy_sc_all_RECOMM_PLANNING.pdf",
        fit_func=[plot_utils.FitFunctions.RecommPlanningX,
                  plot_utils.FitFunctions.RecommPlanningY],
        filename_residuals=(
            f"{conf['plot']['dir']}/"
            f"coalign_TxTy_sc_all_RECOMM_PLANNING_residuals.pdf"),
        )
