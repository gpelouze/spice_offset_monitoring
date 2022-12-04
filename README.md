# SPICE pointing offset monitoring

Monitor the evolution of the pointing offset of the SPICE spectrometer by
coaligning rasters with images of EUI/FSI.

This tool requires access to the internal SPICE and EUI data archives,
[EUI selektor], [`euiprep`], and the [SPICE kernels for Solar
Orbiter][solo-spice-kernels]. See the [installation](#installation) section
for further information.

[EUI selektor]: https://www.sidc.be/EUI/data_internal/selektor/
[`euiprep`]: https://gitlab-as.oma.be/SIDC/SpaceInstruments/eui/
[solo-spice-kernels]: https://www.cosmos.esa.int/web/spice/solar-orbiter


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
conf = 'config/example.yml'
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
pip install git+https://github.com/gpelouze/quick_spice_fsi_coalign
```

### Dependencies (TODO)

- Access to the data archive `$SOLO_ARCHIVE`
- Access to EUI selektor (client installed as dependency, but see documentation
  for credentials)
- Spice_kernels `$SPICE_KERNELS_SOLO`
- eui_prep
