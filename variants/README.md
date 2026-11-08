# Additional and earlier versions

The eight reviewed scripts are in the repository root. These versions are retained for comparison and different output options; their assumptions are not all interchangeable.

| File | Purpose and differences |
| --- | --- |
| [red_brace_xz_forces_x_and_z_scipy.py](red_brace_xz_forces_x_and_z_scipy.py) | Earlier SciPy implementation of x/z forces on the rising-x brace. Fixed-crest limits; normal velocity includes u and w. Origin at the brace SWL crossing. |
| [red_brace_xz_force_x_surface_options.py](red_brace_xz_force_x_surface_options.py) | Horizontal x forces with selectable crest, SWL or instantaneous-surface limits and rising/falling orientation. Origin at the brace SWL crossing. |
| [yellow_brace_yz_force_x_from_u_only.py](yellow_brace_yz_force_x_from_u_only.py) | Separate u-only horizontal x model with a common phase along the y-z brace. Drag uses u times abs(u), unlike the combined-speed version in the root. |
| [yellow_brace_yz_forces_y_and_z_from_w_combined_output.py](yellow_brace_yz_forces_y_and_z_from_w_combined_output.py) | Corrected falling-y front brace, default x=0. Both y and z outputs in one script, with normal forces available as reference values. Equivalent components to the two separate w-only root scripts. |
| [yellow_brace_yz_forces_y_and_z_from_w_rising_y_legacy.py](yellow_brace_yz_forces_y_and_z_from_w_rising_y_legacy.py) | Earlier geometry: rising in +y at x=30, positive normal towards -y and +z. Its y signs and reference phase differ from the front-brace scripts. |

Run the surface-options version from the repository root:

```sh
python variants/red_brace_xz_force_x_surface_options.py --surface crest
```

Use `--surface swl` or `--surface instantaneous` only when intentionally comparing those assumptions.
