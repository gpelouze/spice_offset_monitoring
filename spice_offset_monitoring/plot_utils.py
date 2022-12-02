import abc
import glob
import os

from astropy.io import fits
import numpy as np
import pandas as pd
import scipy.optimize as sopt
import sklearn.metrics as skm
import svgpath2mpl
import tqdm
import yaml

import utils


marker = svgpath2mpl.parse_path('''
m -49.85316,289.25194 c -2.56128,2.76813 -5.26119,5.3446 -7.66432,8.23714 -2.27869,1.78295 -1.17538,1.70123 -2.14655,3.50822 -3.09091,-3.34488 -6.40062,-6.50127 -9.34717,-9.97289 -1.94049,-4.5167 -5.96962,-1.22916 -7.46765,1.68801 -2.44665,2.23713 -4.29329,4.9939 -6.37245,7.54802 -0.61277,2.85666 3.66246,4.52026 5.14883,6.80331 5.59261,5.76949 11.11525,11.61369 16.66755,17.41598 2.2086,1.83902 5.06987,5.9849 7.52929,1.97188 1.80251,-2.60959 4.23286,-4.45152 5.72762,-7.19003 2.73518,-2.11296 0.15377,-2.16952 -0.11212,-5.331 -1.42302,-4.92029 4.50886,-0.50444 6.808776,-0.29187 3.097384,1.11864 6.018669,-0.63079 7.922927,-3.0082 2.181883,-2.48726 -3.421284,-4.413 -0.04014,-2.63434 4.793558,0.75034 7.013437,0.54105 11.756401,-0.77513 3.179606,-1.2276 5.589153,-3.13953 8.141494,-5.30735 1.880693,-2.15618 5.985017,-7.2305 0.478606,-7.31889 -5.439368,-0.71882 -11.012657,-0.45346 -16.360826,-1.89604 -3.360413,-0.47425 -7.291001,-1.52814 -10.210792,0.89265 -2.620313,1.79523 2.873571,4.67394 -1.432812,1.87326 -2.548494,-2.52714 -6.441884,-7.27211 -9.026664,-6.21273 z
''')
marker.vertices -= marker.vertices.mean(axis=0)


def list_of_dict_to_dict_of_arr(lst):
    keys = lst[0].keys()
    d_out = {k: [] for k in keys}
    for d in lst:
        if d.keys() != keys:
            raise ValueError('mismatched keys')
        for k, v in d.items():
            d_out[k].append(v)
    d_out = {k: np.array(v) for k, v in d_out.items()}
    return d_out


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


def get_data(conf, time_span):
    output_dir = time_span['dir']
    yml_fnames = glob.glob(f'{output_dir}/coalign_output/*_coaligned.yml')

    dat = []
    for yml_fname in tqdm.tqdm(yml_fnames, desc='Loading data'):
        spice_fname = os.path.basename(yml_fname).rstrip('_coaligned.yml')
        with open(yml_fname, 'r') as f:
            res = yaml.safe_load(f)
        res['date'] = utils.SpiceUtils.filename_to_date(f'{spice_fname}.fits')
        plot_dir = os.path.relpath(
            f'{output_dir}/coalign_output',
            conf['plot']['dir'],
            )
        res['plot_pdf'] = f'{plot_dir}/{spice_fname}_coaligned.pdf'
        res['plot_jpg'] = f'{plot_dir}/{spice_fname}_coaligned.jpg'
        res['plot_marker'] = time_span['plot_marker']
        wcs = res.pop('wcs')
        res.update(wcs)
        data_dir = 'spice_L2r'
        header_data = get_header_data(os.path.join(
            output_dir, data_dir, f'{spice_fname}.fits'))
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
    dat['R_cen_x'] = dat['CRVAL1'] * 3600 / R_sun  # arcsec
    dat['R_cen_y'] = dat['CRVAL2'] * 3600 / R_sun  # arcsec

    dat['data_source'] = time_span['title']

    return pd.DataFrame(dat)


class FitFunctions:
    class _Base:
        def __init__(self, x, y):
            self.x = x
            self.y = y
            self.popt = None
            self.pcov = None

        @abc.abstractmethod
        def __call__(self, x, *args):
            pass

        @abc.abstractmethod
        def label(self, y, x, error=True, R2=True, fmt='.2f'):
            pass

        @staticmethod
        def p0():
            return None

        def fit(self):
            self.popt, self.pcov = sopt.curve_fit(
                self,
                self.x,
                self.y,
                p0=self.p0(),
                )

        def xopt(self, n=100):
            # return np.linspace(np.min(self.x), np.max(self.x), n)
            return np.sort(self.x)

        def yopt(self):
            return self(self.xopt(), *self.popt)

        def perr(self):
            return np.sqrt(np.diag(self.pcov))

        def r2(self):
            y_pred = self(self.x, *self.popt)
            return skm.r2_score(self.y, y_pred)

    class Linear(_Base):
        def __call__(self, x, a, b):
            return a * x + b

        def label(self, y, x, error=True, R2=True, fmt='.2f'):
            fmt = f'{{:{fmt}}}'
            a, b = self.popt
            a = fmt.format(a)
            b = fmt.format(b)
            if error:
                aerr, berr = self.perr()
                aerr = fmt.format(aerr)
                berr = fmt.format(berr)
                a = f'({a} ± {aerr})'
                b = f'({b} ± {berr})'
            if R2:
                r2 = '{:.2f}'.format(self.r2())
                return f'${y} = {a} {x} + {b}$ [$R^2 = {r2}$]'
            else:
                return f'${y} = {a} {x} + {b}$'

        def p0(self):
            return [np.ptp(self.y), np.min(self.y)]

    class Constant(_Base):
        def __call__(self, x, b):
            return np.full_like(x, b)

        def label(self, y, x, error=True, R2=True, fmt='.2f'):
            fmt = f'{{:{fmt}}}'
            b, = self.popt
            b = fmt.format(b)
            if error:
                berr, = self.perr()
                berr = fmt.format(berr)
                b = f'({b} ± {berr})'
            if R2:
                r2 = '{:.2f}'.format(self.r2())
                return f'${y} = {b}$ [$R^2 = {r2}$]'
            else:
                return f'${y} = {b}$'

        def p0(self):
            return [np.mean(self.y)]

    class RecommX(Linear):
        def fit(self):
            self.popt = [0.46, -85.0]
            self.pcov = None

        def label(self, *args, **kwargs):
            return f'$\\Delta X = 0.46 T - 85.0$'

    class RecommY(Constant):
        def fit(self):
            self.popt = [-72.2]
            self.pcov = None

        def label(self, *args, **kwargs):
            return f'$\\Delta Y = - 72.2$'

    class InverseSq(_Base):
        def __call__(self, x, a, b):
            return a / x**2 + b

        def label(self, y, x, error=True, R2=True, fmt='.2f'):
            fmt = f'{{:{fmt}}}'
            a, b = self.popt
            a = fmt.format(a)
            b = fmt.format(b)
            if error:
                aerr, berr = self.perr()
                aerr = fmt.format(aerr)
                berr = fmt.format(berr)
                a = f'({a} ± {aerr})'
                b = f'({b} ± {berr})'
            if R2:
                r2 = '{:.2f}'.format(self.r2())
                return f'${y} = {a} / {x}² + {b}$ [$R^2 = {r2}$]'
            else:
                return f'${y} = {a} / {x}² + {b}$'


class Filters:

    @staticmethod
    def center(df):
        print('\nFiltering data')
        m_cc = (df['max_cc'] > 0.2)
        print(f'Frames with cc < 0.2: {1 - m_cc.mean():.1%}')
        m_dr = (df['dr'] < 30)
        print(f'Frames with dr > 30: {1 - m_dr.mean():.1%}')
        m_R = (df['R_cen'] < 0.2)
        print(f'Frames with R_sun > 0.2: {1 - m_R.mean():.1%}')
        m = m_cc & m_dr & m_R
        print(f'Discarded frames: {1 - m.mean():.1%}')
        return df[m]

    @staticmethod
    def disk(df):
        print('\nFiltering data')
        m_cc = (df['max_cc'] > 0.2)
        print(f'Frames with cc < 0.2: {1 - m_cc.mean():.1%}')
        m_dr = (df['dr'] < 30)
        print(f'Frames with dr > 30: {1 - m_dr.mean():.1%}')
        m_dsun = (df['DSUN_AU'] < 0.6)
        print(f'Frames with d_sun > 0.6: {1 - m_dsun.mean():.1%}')
        m = m_cc & m_dr & m_dsun
        print(f'Discarded frames: {1 - m.mean():.1%}')
        return df[m]

    @staticmethod
    def limb(df):
        print('\nFiltering data')
        m_cc = (df['max_cc'] > 0.2)
        print(f'Frames with cc < 0.2: {1 - m_cc.mean():.1%}')
        m_dr = (df['dr'] < 30)
        print(f'Frames with dr > 30: {1 - m_dr.mean():.1%}')
        m_dsun = True
        # m_dsun = (df['DSUN_AU'] > 0.5)
        # m_dsun &= (df['DSUN_AU'] < 0.6)
        # print(f'Frames with d_sun > 0.6: {1 - m_dsun.mean():.1%}')
        m = m_cc & m_dr & m_dsun
        print(f'Discarded frames: {1 - m.mean():.1%}')
        return df[m]
