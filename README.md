If not running on IAS sol-calcul, requires environment variable `$SOLO_ARCHIVE`
to point to the directory containing the `EUI` and `SPICE` archives.

Requires [`spice_stew`](https://github.com/gpelouze/spice_stew)

# Troubleshooting

- `SpiceNOFRAMECONNECT: At epoch ..., there is insufficient information
  available to transform from reference frame -144991 (SOLO_SUN_RTN) to
  reference frame -144000 (SOLO_SRF).` is raised when the requested date is not
  in the as-flown SPICE kernels. This happens either because the kernels have
  not been updated, or because a date in the future is being requested.
