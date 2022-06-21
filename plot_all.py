#!/usr/bin/env python

import glob
import os

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import svgpath2mpl
import yaml
from astropy.io import fits
import tqdm

from quick_spice_fsi_coalign import SpiceUtils
from plot_alignment_results import list_of_dict_to_dict_of_arr


marker = svgpath2mpl.parse_path('''
m -49.85316,289.25194 c -2.56128,2.76813 -5.26119,5.3446 -7.66432,8.23714 -2.27869,1.78295 -1.17538,1.70123 -2.14655,3.50822 -3.09091,-3.34488 -6.40062,-6.50127 -9.34717,-9.97289 -1.94049,-4.5167 -5.96962,-1.22916 -7.46765,1.68801 -2.44665,2.23713 -4.29329,4.9939 -6.37245,7.54802 -0.61277,2.85666 3.66246,4.52026 5.14883,6.80331 5.59261,5.76949 11.11525,11.61369 16.66755,17.41598 2.2086,1.83902 5.06987,5.9849 7.52929,1.97188 1.80251,-2.60959 4.23286,-4.45152 5.72762,-7.19003 2.73518,-2.11296 0.15377,-2.16952 -0.11212,-5.331 -1.42302,-4.92029 4.50886,-0.50444 6.808776,-0.29187 3.097384,1.11864 6.018669,-0.63079 7.922927,-3.0082 2.181883,-2.48726 -3.421284,-4.413 -0.04014,-2.63434 4.793558,0.75034 7.013437,0.54105 11.756401,-0.77513 3.179606,-1.2276 5.589153,-3.13953 8.141494,-5.30735 1.880693,-2.15618 5.985017,-7.2305 0.478606,-7.31889 -5.439368,-0.71882 -11.012657,-0.45346 -16.360826,-1.89604 -3.360413,-0.47425 -7.291001,-1.52814 -10.210792,0.89265 -2.620313,1.79523 2.873571,4.67394 -1.432812,1.87326 -2.548494,-2.52714 -6.441884,-7.27211 -9.026664,-6.21273 z
''')
marker.vertices -= marker.vertices.mean(axis=0)


def get_header_data(fname):
    hdul = fits.open(fname)
    header = hdul[0].header
    keywords = [
        'DSUN_AU',
        'T_GRAT',
        'T_FOCUS',
        'T_SW',
        'T_LW',
        ]
    return {kw: header[kw] for kw in keywords}


def get_data(output_dir):
    yml_fnames = glob.glob(f'{output_dir}/coalign_output/*_coaligned.yml')

    dat = []
    for yml_fname in tqdm.tqdm(yml_fnames, desc='Loading data'):
        spice_fname = os.path.basename(yml_fname).rstrip('_coaligned.yml')
        with open(yml_fname, 'r') as f:
            res = yaml.safe_load(f)
        res['date'] = SpiceUtils.filename_to_date(f'{spice_fname}.fits')
        res['plot'] = f'coalign_output/{spice_fname}_coaligned.pdf'
        wcs = res.pop('wcs')
        res.update(wcs)
        header_data = get_header_data(os.path.join(
            output_dir, 'spice_stew', f'{spice_fname}_remapped_img.fits'))
        res.update(header_data)
        dat.append(res)

    if len(dat) == 0:
        return None

    dat = list_of_dict_to_dict_of_arr(dat)

    # add columns
    dat['dr'] = np.sqrt(np.array(dat['dx'])**2 + np.array(dat['dy'])**2)
    r = np.deg2rad(dat['roll'])
    dat['dx_sc'] = dat['dx'] * np.cos(r) - dat['dy'] * np.sin(r)
    dat['dy_sc'] = dat['dx'] * np.sin(r) + dat['dy'] * np.cos(r)
    dat['dx_sc'] -= 83
    dat['dy_sc'] -= 68
    # R_cen: Distance between raster center and Sun center
    R_sun = 959.23 / dat['DSUN_AU']  # arcsec
    R = np.sqrt(dat['CRVAL1']**2 + dat['CRVAL2']**2) * 3600  # arcsec
    dat['R_cen'] = R / R_sun  # R_sun

    return pd.DataFrame(dat)


def filter_data(df):
    print()
    m_cc = (df['max_cc'] > 0.2)
    print(f'Frames with cc < 0.2: {1 - m_cc.mean():.1%}')
    m_dr = (df['dr'] < 50)
    print(f'Frames with dr > 50: {1 - m_dr.mean():.1%}')
    m_R = (df['R_cen'] < 0.2)
    print(f'Frames with R_sun > 0.2: {1 - m_R.mean():.1%}')
    m = m_cc & m_dr & m_R
    print(f'Discarded frames: {1 - m.mean():.1%}')
    return df[m]


if __name__ == '__main__':
    datasets = [
        ('Cal. Lyβ', get_data('output/CAL_COALIGN_Lyb'), dict(marker='o')),
        ('Syn. Lyβ', get_data('output/SYN_Lyb'), dict(marker='s')),
        ('Syn. Lyγ CIII', get_data('output/SYN_Lyg_CIII'), dict(marker='D')),
        ('Syn. CIII', get_data('output/SYN_CIII'), dict(marker='*')),
        ]
    kw_common = dict(fillstyle='none', ls='')
    for _, _, kw in datasets:
        kw.update(kw_common)
    datasets = [(name, dat, kw)
                for (name, dat, kw) in datasets
                if dat is not None]
    datasets_filtered = [(name, filter_data(dat), kw)
                         for (name, dat, kw) in datasets]

    fig1 = plt.figure(1, clear=True)
    fig2 = plt.figure(2, clear=True)
    fig3 = plt.figure(3, clear=True)
    ax1 = fig1.gca()
    ax2 = fig2.gca(sharex=ax1)
    ax3 = fig3.gca(sharex=ax1)
    ax1.set_title('SPICE offset', loc='left')
    ax2.set_title('SPICE offset\n(monthly)', loc='left')
    ax3.set_title('SPICE offset (monthly average)', loc='left')
    # All data
    for name, dat, kw in datasets_filtered:
        ms = 3
        ax1.plot(dat['date'], dat['dx_sc'], color='C0', ms=ms, **kw)
        ax1.plot(dat['date'], dat['dy_sc'], color='C1', ms=ms, **kw)
    # Monthly avg per study group
    for i, (name, dat, kw) in enumerate(datasets_filtered):
        datg = pd.DataFrame(dat).groupby(pd.PeriodIndex(dat['date'], freq="M"))
        datm = datg.mean()
        dats = datg.std()
        t = datm.index.to_timestamp().to_numpy()
        t = t + np.timedelta64(15, 'D') + (i-1)*np.timedelta64(10, 'D')
        ax2.errorbar(t, datm['dx_sc'].values, yerr=dats['dx_sc'].values, color='C0', capsize=3, elinewidth=1, **kw)
        ax2.errorbar(t, datm['dy_sc'].values, yerr=dats['dy_sc'].values, color='C1', capsize=3, elinewidth=1, **kw)
    # Monthly avg all
    dat = pd.concat([pd.DataFrame(d[1]) for d in datasets_filtered])
    datg = pd.DataFrame(dat).groupby(pd.PeriodIndex(dat['date'], freq="M"))
    datm = datg.mean()
    dats = datg.std()
    t = datm.index.to_timestamp().to_numpy()
    t = t + np.timedelta64(15, 'D')
    ax3.errorbar(t, datm['dx_sc'].values, yerr=dats['dx_sc'].values, label='X', color='C0', marker='x', capsize=3, elinewidth=1, **kw_common)
    ax3.errorbar(t, datm['dy_sc'].values, yerr=dats['dy_sc'].values, label='Y', color='C1', marker='x', capsize=3, elinewidth=1, **kw_common)
    #
    for ax in (ax1, ax2, ax3):
        ax.set_xlabel('Date')
        ax.set_ylabel('Offset [arcsec]')
        handles, _ = ax.get_legend_handles_labels()
        if not handles:
            handles = [
                plt.Line2D([], [], label='X', color='C0', marker=marker, mew=0, ms=12, ls=''),
                plt.Line2D([], [], label='Y', color='C1', marker=marker, mew=0, ms=12, ls=''),
                ]
            for name, _, kw in datasets_filtered:
                handles.append(plt.Line2D([], [], label=name, color='gray', **kw))
        legend_prop = {
            ax1: dict(
                ncol=3,
                loc='lower right',
                bbox_to_anchor=(1, 1),
                fancybox=False,
                fontsize=10,
                ),
            ax2: dict(
                ncol=3,
                loc='lower right',
                bbox_to_anchor=(1, 1),
                fancybox=False,
                fontsize=10,
                ),
            ax3: dict(
                ncol=3,
                loc='lower right',
                bbox_to_anchor=(1, 1),
                fancybox=False,
                fontsize=10,
                ),
            }[ax]
        ax.legend(handles=handles, **legend_prop)
        ax.figure.autofmt_xdate()
        ax.xaxis.set_major_locator(mpl.dates.MonthLocator(bymonth=[1, 4, 7, 10]))
        ax.xaxis.set_minor_locator(mpl.dates.MonthLocator())
        ax.xaxis.set_major_formatter(mpl.dates.DateFormatter("%b %y"))
        # etframes.add_range_frame(ax)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.axhline(-83, color='C0', lw=.5)
        ax.axhline(-68, color='C1', lw=.5)
        x0, x1 = ax.get_xlim()
        x = x1 # + (x1 - x0) / 50
        ax.text(x, -83, '−83″', color='C0', fontsize=10, ha='center', va='bottom')
        ax.text(x, -68, '−68″', color='C1', fontsize=10, ha='center', va='bottom')
    fig1.savefig(f'output/coalign_TxTy_sc_all.pdf')
    fig2.savefig(f'output/coalign_TxTy_sc_all_monthly.pdf')
    fig3.savefig(f'output/coalign_TxTy_sc_all_monthly_grouped.pdf')

    for name, dat, kw in datasets_filtered:
        print('\n', name, sep='')
        print('dx:', np.mean(dat['dx']), np.std(dat['dx']))
        print('dy:', np.mean(dat['dy']), np.std(dat['dy']))
        print('dr:', np.mean(dat['dr']), np.std(dat['dr']))
        print('dx s/c:', np.mean(dat['dx_sc']), np.std(dat['dx_sc']))
        print('dy s/c:', np.mean(dat['dy_sc']), np.std(dat['dy_sc']))

    plt.clf()
    ax = plt.gca()
    ax.set_title('SPICE offset', loc='left')
    # All data
    for name, dat, kw in datasets_filtered:
        ms = 3
        ax.plot(dat['DSUN_AU'], dat['dx_sc'], color='C0', ms=ms, **kw)
        ax.plot(dat['DSUN_AU'], dat['dy_sc'], color='C1', ms=ms, **kw)
    ax.set_xlabel('Solar distance [au]')
    ax.set_ylabel('Offset [arcsec]')
    handles, _ = ax.get_legend_handles_labels()
    if not handles:
        handles = [
            plt.Line2D(
                [], [], label='X', color='C0', marker=marker, mew=0, ms=12,
                ls=''
                ),
            plt.Line2D(
                [], [], label='Y', color='C1', marker=marker, mew=0, ms=12,
                ls=''
                ),
            ]
        for name, _, kw in datasets_filtered:
            handles.append(plt.Line2D([], [], label=name, color='gray', **kw))
    legend_prop = dict(
        ncol=3,
        loc='lower right',
        bbox_to_anchor=(1, 1),
        fancybox=False,
        fontsize=10,
        )
    ax.legend(handles=handles, **legend_prop)
    # etframes.add_range_frame(ax)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.axhline(-83, color='C0', lw=.5)
    ax.axhline(-68, color='C1', lw=.5)
    x0, x1 = ax.get_xlim()
    x = x1  # + (x1 - x0) / 50
    ax.text(x, -83, '−83″', color='C0', fontsize=10, ha='center', va='bottom')
    ax.text(x, -68, '−68″', color='C1', fontsize=10, ha='center', va='bottom')
    plt.savefig(f'output/coalign_TxTy_sc_all_dsun.pdf')

    plt.clf()
    ax = plt.gca()
    ax.set_title('SPICE offset\n($d <$ 0.6 au)', loc='left')
    # All data
    for name, dat, kw in datasets:
        ms = 3
        m = (
            (dat['max_cc'] > 0.2)
            & (dat['dr'] < 40)
            & (dat['DSUN_AU'] < .6)
            )
        dat = dat[m]
        ax.plot(dat['R_cen'], dat['dx_sc'], color='C0', ms=ms, **kw)
        ax.plot(dat['R_cen'], dat['dy_sc'], color='C1', ms=ms, **kw)
    ax.set_xlabel('Raster center [$R_\odot$]')
    ax.set_ylabel('Offset [arcsec]')
    handles, _ = ax.get_legend_handles_labels()
    if not handles:
        handles = [
            plt.Line2D(
                [], [], label='X', color='C0', marker=marker, mew=0, ms=12,
                ls=''
                ),
            plt.Line2D(
                [], [], label='Y', color='C1', marker=marker, mew=0, ms=12,
                ls=''
                ),
            ]
        for name, _, kw in datasets_filtered:
            handles.append(plt.Line2D([], [], label=name, color='gray', **kw))
    legend_prop = dict(
        ncol=3,
        loc='lower right',
        bbox_to_anchor=(1, 1),
        fancybox=False,
        fontsize=10,
        )
    ax.legend(handles=handles, **legend_prop)
    # etframes.add_range_frame(ax)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.axhline(-83, color='C0', lw=.5)
    ax.axhline(-68, color='C1', lw=.5)
    x0, x1 = ax.get_xlim()
    x = x1  # + (x1 - x0) / 50
    ax.text(x, -83, '−83″', color='C0', fontsize=10, ha='center', va='bottom')
    ax.text(x, -68, '−68″', color='C1', fontsize=10, ha='center', va='bottom')
    plt.savefig(f'output/coalign_TxTy_sc_all_Rcen.pdf')

    for T_key in ['T_GRAT', 'T_FOCUS', 'T_SW', 'T_LW']:
        plt.clf()
        ax = plt.gca()
        ax.set_title('SPICE offset', loc='left')
        # All data
        for name, dat, kw in datasets_filtered:
            ms = 3
            ax.plot(dat[T_key], dat['dx_sc'], color='C0', ms=ms, **kw)
            ax.plot(dat[T_key], dat['dy_sc'], color='C1', ms=ms, **kw)
        ax.set_xlabel(f'{T_key} [°C]')
        ax.set_ylabel('Offset [arcsec]')
        handles, _ = ax.get_legend_handles_labels()
        if not handles:
            handles = [
                plt.Line2D(
                    [], [], label='X', color='C0', marker=marker, mew=0, ms=12,
                    ls=''
                    ),
                plt.Line2D(
                    [], [], label='Y', color='C1', marker=marker, mew=0, ms=12,
                    ls=''
                    ),
                ]
            for name, _, kw in datasets_filtered:
                handles.append(plt.Line2D([], [], label=name, color='gray', **kw))
        legend_prop = dict(
            ncol=3,
            loc='lower right',
            bbox_to_anchor=(1, 1),
            fancybox=False,
            fontsize=10,
            )
        ax.legend(handles=handles, **legend_prop)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.axhline(-83, color='C0', lw=.5)
        ax.axhline(-68, color='C1', lw=.5)
        plt.savefig(f'output/coalign_TxTy_sc_all_{T_key}.pdf')

    for T_key in ['T_GRAT', 'T_FOCUS', 'T_SW', 'T_LW']:
        plt.clf()
        ax = plt.gca()
        ax.set_title('SPICE offset', loc='left')
        # All data
        handles = []
        for name, dat, kw in datasets_filtered:
            ms = 3
            ax.plot(dat['DSUN_AU'], dat[T_key], color='k', ms=ms, **kw)
            handles.append(plt.Line2D([], [], label=name, color='gray', **kw))
        ax.set_xlabel('Distance [au]')
        ax.set_ylabel(f'{T_key} [°C]')
        legend_prop = dict(
            ncol=3,
            loc='lower right',
            bbox_to_anchor=(1, 1),
            fancybox=False,
            fontsize=10,
            )
        ax.legend(handles=handles, **legend_prop)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        plt.savefig(f'output/coalign_TxTy_sc_all_{T_key}.pdf')