# SPICE pointing offset monitoring

Monitor the evolution of the pointing offset of the SPICE spectrometer by
coaligning rasters with images of EUI/FSI.

This tool requires access to the internal SPICE and EUI data archives,
[EUI selektor], [`euiprep`], and the [SPICE kernels for Solar
Orbiter][solo-spice-kernels]. See the [installation](#installation) section
for further information.


## Usage

### Command line and Python interfaces

The module has three main commands:

- `process`: computes the offset between SPICE and FSI
- `gen_explorer`: generates an interactive result explorer (bokeh)
- `gen_figures`: generates the figures (matplotlib)

All commands require the path to a configuration file (see section below).

They can be called from the command line:

```shell
spice_offset_monitoring -c config/example.yml process
spice_offset_monitoring -c config/example.yml gen-explorer
spice_offset_monitoring -c config/example.yml gen-figures
```

Commands can also be chained:

```shell
spice_offset_monitoring -c config/example.yml process gen-figures
```

In Python, use:

```python
import spice_offset_monitoring as som
conf = som.Config('config/example.yml')
som.process(conf)
som.gen_explorer(conf)
som.gen_figures(conf)
```

The processing command tries to not regenerate files that already exist. This
means that it can be run several times in a row with the same config file, and
will not rerun unnecessary computations. It allows e.g. to restart it after an
interruption, or to process data newly added to the archive. To regenerate
files, delete them.


### Configuration

Configuration is set in configuration files with the following format:

```yaml
processing:
  jitter_correction: ['fits', 'kernels', 'dummy']
  # List of methods to correct the jitter
  #   fits: use the pointing in the FITS files (files after April 2022)
  #   kernels: compute the pointing from SPICE kernels (requires spice_stew)
  #   dummy: apply no correction, but still generate the files
  # The methods are tried in order, falling back to next one in case of failure

plot:
  dir: 'output/example/plots'


time_spans:
  # List of one or multiple time spans, with the following options:

  - title: 'Synoptic'
    dir: 'output/example/SYN_Lyg_CIII_new_syn/'
    start_date: '2022-05-01'
    end_date: '2022-06-01'
    study_id: '1963'
    spec_win: 'Ly-gamma-CIII group (Merged)'
    plot_marker: 's'
```

See more examples in `config/`.

## Installation

### Package

```shell
pip install git+https://github.com/gpelouze/spice_offset_monitoring
```


### Dependencies

#### Access to the SPICE and EUI data archives

The archive containing the L2 SPICE and L1 FSI FITS files should be mounted
onto your system, and its path stored in the environment variable
`$SOLO_ARCHIVE`, eg:

```shell
export SOLO_ARCHIVE=/archive/SOLAR-ORBITER
```


#### EUI selektor client

This tool uses [EUI selektor] through the [`eui_selektor_client`]. The client
is automatically installed as a dependency, but it requires credentials to
access selektor. Make sure to read the [client's README][`eui_selektor_client`]
if you are not prompted for selektor credentials when running
`spice_offset_monitoring` for the first time.


#### `eui_prep`

[`euiprep`] required to prepare L2 FSI files with an accurate position of the
Sun center, determined by fitting the limb. To install it, clone the
repository and install its dependencies:

```shell
cd /path/to/eui_soft
git clone https://gitlab-as.oma.be/SIDC/SpaceInstruments/eui/
pip install GitPython opencv-python scikit-image
```

and add its location to your `$PYTHONPATH`:

```shell
export PYTHONPATH=/path/to/eui_soft/:$PYTHONPATH
```


#### SPICE kernels

If the `'kernels'` method is present in the `jitter_correction` list, the
[SPICE kernels for Solar Orbiter][solo-spice-kernels] (~2.9 GiB once cloned)
are also required. These are used by [`spice_stew`] through [`spiceypy`]. While
these Python packages are always installed, the kernels are only needed if you
actually use the `'kernels'` method.

To install them:

```shell
SPICE_KERNELS_SOLO="/path/to/spice_kernels/SOLO/"
git clone --depth 1 https://repos.cosmos.esa.int/socci/scm/spice_kernels/solar-orbiter.git "$SPICE_KERNELS_SOLO"
export SPICE_KERNELS_SOLO="$SPICE_KERNELS_SOLO/kernels"
# replace `PATH_VALUES       = ( '..' )` with absolute path:
sed -i "s+\(^[ ]*PATH_VALUES[ ]*=[ ]* ([ ]*'\)\.\.\('[ ]*)\)$+\1$SPICE_KERNELS_SOLO\2+" "$SPICE_KERNELS_SOLO/mk/"*.tm
```


[`eui_selektor_client`]: https://github.com/gpelouze/eui_selektor_client/
[EUI selektor]: https://www.sidc.be/EUI/data_internal/selektor/
[`euiprep`]: https://gitlab-as.oma.be/SIDC/SpaceInstruments/eui/
[solo-spice-kernels]: https://www.cosmos.esa.int/web/spice/solar-orbiter
[`spice_stew`]: https://github.com/gpelouze/spice_stew
[`spiceypy`]: https://pypi.org/project/spiceypy/