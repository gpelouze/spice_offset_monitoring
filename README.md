# Quick and dirty SPICE-FSI coalignment

The repository requires clean-up before it can be shared.

If not running on IAS sol-calcul, requires environment variable `$SOLO_ARCHIVE`
to point to the directory containing the `EUI` and `SPICE` archives.

Requires and hence the Solar Orbiter SPICE kernels for now, used by 
[`spice_stew`](https://github.com/gpelouze/spice_stew).

Also requires [`euiprep`](https://gitlab-as.oma.be/SIDC/SpaceInstruments/eui/)
to be installed in the `$PYTHONPATH`.