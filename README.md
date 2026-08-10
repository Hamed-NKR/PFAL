Post-Flame Agglomeration Algorithm (PFAL)
===

PFAL generates soot-like fractal aggregates with a staged Langevin dynamics
workflow. The core simulation pipeline has three stages:

1. `main_LD1_v3.m` generates a raw first-stage aggregate library.
2. `main_scale_v2.m` or `main_scatter_v8.m` calibrates that library to target
   primary-particle and projected-area distributions.
3. `main_LD2_v3.m` runs second-stage post-flame agglomeration on the calibrated
   aggregate population.

`main_tem_analysis_v1.m` is a separate analysis workflow for TEM aggregate and
primary-particle measurements. It produces summary tables, model-input
statistics, and publication figures from configured datasets.

`main_valid_v3.m` compares selected LD2 populations with the processed TEM and
tandem AAC-SMPS effective-density measurements. It produces manuscript-ready
`d_pp` versus `d_a` and `rho_eff` versus `d_m` figures together with pointwise
predictions, log-space validation metrics, and source provenance.

The model represents soot aggregates that first form as compact
diffusion-limited clusters, acquire realistic size distributions, and continue
agglomerating after dilution or cooling. It is intended for studies of hybrid
soot structures with nonuniform primary-particle sizes.

Examples of hybrid soot imaged under transmission electron microscopy:

![Hybrid soot TEM examples](https://github.com/user-attachments/assets/db858f62-1e22-471c-90ad-dc06329de681)

## Main Workflow

```text
main_LD1_v3
first-stage Langevin dynamics
        |
        v
data/main_ld1/ld1_aggregate_library.mat
        |
        +------------------------------+
        |                              |
        v                              v
main_scale_v2                    main_scatter_v8
correlation-based scaling        bivariate scatter sampling (default)
        |                              |
        v                              v
data/main_scale/                 data/main_scatter/
scaled_aggs_for_LD2_             scaled_aggs_for_LD2_
from_main_scale.mat              from_main_scatter.mat
        |                              |
        +--------------+---------------+
                       |
                       v
                 main_LD2_v3
        second-stage Langevin dynamics
                       |
                       v
     results/main_ld2_from_main_*/
       LD2__*__YYYY-MM-DD/
```

### Stage 1: LD1 Aggregate Library

`main_LD1_v3.m` performs first-stage Langevin dynamics. It initializes primary
particles, assigns positions and velocities, evaluates mobility, marches the
particles, applies periodic boundaries, joins collisions, and updates size and
mobility properties after every growth step.

LD1 saves its default output to:

```text
data/main_ld1/ld1_aggregate_library.mat
```

The file contains:

- `pp0`: flat cell array of aggregate primary-particle matrices, directly
  usable by `main_scale_v2.m` and `main_scatter_v8.m`.
- `pp0_n`: number of primaries in each saved aggregate.
- `pars_ld1`: structured version of the saved aggregate library.
- `ensdata0`, `parsdata0`: real-time ensemble summaries and saved aggregate
  properties.
- `cfg_ld1`, `metadata`: config and run provenance.

```text
LD1 config JSON
      |
      v
main_LD1_v3
      |
      +-- pp0                 flat aggregate cell library
      |     |
      |     +--> main_scale_v2
      |     |
      |     `--> main_scatter_v8
      |
      +-- pp0_n               primary-particle counts
      |
      +-- pars_ld1            sizing fields and aggregate structure
      |
      `-- ensdata0, parsdata0 diagnostics and saved snapshots
```

### Stage 2A: Scale Or Scatter

Both calibration scripts start from the LD1 `pp0` library.

`main_scale_v2.m` rescales and filters LD1 aggregates against the Brasil et al.
(1999) aggregate correlation and the Olfert and Rogak (2019) universal
correlation. It writes:

```text
data/main_scale/scaled_aggs_for_LD2_from_main_scale.mat
```

`main_scatter_v8.m` samples a projected-area/primary-particle size
distribution, assigns LD1 aggregate seeds to those targets, and rescales the
selected aggregates. Bivariate sampling is the default. The `sequential` and
`ideal` modes remain available for sensitivity runs through
`options.opt_scale` in the scatter config.

Every scatter run writes the canonical LD2 input:

```text
data/main_scatter/scaled_aggs_for_LD2_from_main_scatter.mat
```

It also writes a mode-specific copy, such as:

```text
data/main_scatter/scaled_aggs_for_LD2_from_main_scatter_bivariate.mat
```

Both calibration workflows save the aggregate payload as `pars_out`, the input
variable expected by LD2. Scatter outputs also contain `scatter_metadata` and
`scatter_summary`.

### Stage 2B: LD2 Post-Flame Agglomeration

`main_LD2_v3.m` loads one of the calibrated `pars_out` files, reinitializes the
aggregate population in a second domain, and performs post-flame Langevin
dynamics. The LD2 config chooses whether the source is the scale output or the
scatter output.

Set `PFAL_MAIN_LD2_CONFIG` before running LD2, for example:

```matlab
setenv('PFAL_MAIN_LD2_CONFIG', ...
    'config/main_ld2/main_ld2_from_main_scatter_config.local.json')
main_LD2_v3
```

LD2 creates a dated run directory under the configured `results.root` and
writes checkpoint and final MAT files using the configured prefixes. An
interrupted checkpoint can be resumed with:

```matlab
UTILS.RESUME_LD2_V3(checkpoint_folder, checkpoint_file)
```

Calling `UTILS.RESUME_LD2_V3` without arguments prompts for both values.

### TEM Measurement Analysis

`main_tem_analysis_v1.m` processes configured TEM datasets independently of
the LD simulation pipeline. Each config entry identifies an aggregate MAT file
(the example expects an `Aggs` variable) and a set of ImageJ primary-particle
area CSV files. The script calculates aggregate- and ensemble-level statistics,
hybridity and collapse summaries, and the four distribution statistics
(`GM_dpp`, `GSD_dpp`, `GM_da`, and `GSD_da`) that correspond to fields in the
scale and scatter configs.

With the example output settings, the main machine-readable files are:

```text
data/main_tem_analysis/tem_analysis_results.mat
data/main_tem_analysis/model_inputs.csv
data/main_tem_analysis/summary_table.csv
```

Additional CSV diagnostics are written to the same directory. When figure
export is enabled, figures are written to `results/main_tem_analysis/`.

### Experimental Validation

The validation workflow consumes processed data and does not rerun TEM or
ODIAS analysis. Keep the original ODIAS workspace at:

```text
data/experimental_effective_density/Effective_Density_Compiled-17_Mar_2025-05_00_44.mat
```

Create the small, self-describing companion used by validation once:

```matlab
source_file = fullfile('data', 'experimental_effective_density', ...
    'Effective_Density_Compiled-17_Mar_2025-05_00_44.mat');
output_file = fullfile('data', 'experimental_effective_density', ...
    'effective_density_validation_data.mat');
UTILS.NORMALIZE_EFFECTIVE_DENSITY_DATA(source_file, output_file, ...
    'ExpectedSourceSHA256', ...
    '3EEF1BC1981EE654BA5FE0198AF50EA977DAEFA5973D37AE27DD787A5844048C');
```

The normalizer loads only `dist_grp` and `test_condition`, copies all four
condition groups without rounding or recalculation, and proves exact numerical
round-trip equality before replacing the companion. The original workspace is
never modified. Its raw tandem distributions, saved figures, and fit objects
are not required by PFAL.

The TEM input is the existing processed artifact:

```text
data/main_tem_analysis/tem_analysis_results.mat
```

`main_valid_v3.m` uses only `aggregate_table.entry_id`, `da_nm`, `dbarpp_nm`,
and the optional `dbarpp_ci95_*` columns. Set the validation config and run:

```matlab
setenv('PFAL_MAIN_VALID_CONFIG', ...
    'config/main_valid/main_valid_config.local.json')
main_valid_v3
```

Low and high agglomeration are enabled by default. Moderate and extensive
collapse remain named, disabled conditions that can be activated after their
TEM and LD2 mappings are supplied. Bayesian fitting, degree, credible level,
posterior samples, priors, plot styling, Segoe UI typography, reference
relations, and exports are controlled independently in JSON.
Legend location, column count, orientation, border visibility, text interpreter,
and text size are configurable under `figures.legend` and
`figures.font.legend_size`. Legend font family, fallback, and weight are also
independent settings; the manuscript profile uses 13-point Segoe UI Light at
normal weight. The plot-frame and tick-line weight is controlled by
`figures.axis_line_width`.
Measurement edge color, marker-face color, outline width, and marker size are
configured per condition through the `conditions[].style.experimental_*`
fields, independently of the numerical-population color. The default marker
face is `none`, leaving the darker measurement outlines transparent.

Stable PDF and PNG figures are written to `results/main_valid/`. Every run also
creates `results/main_valid/runs/<timestamp>/` with the resolved config, source
hashes, pointwise predictions, summary metrics, and a MAT result bundle.

## Configuration and Running

Each configurable workflow has an example JSON file committed under `config/`.
Copy the required example to the corresponding `.local.json` name and edit the
dataset paths and, where present, output settings before running it. Files
matching `.local*.json` are ignored by Git.

```text
config/
|-- main_ld1/
|   `-- main_ld1_config.example.json
|-- main_scale/
|   `-- main_scale_config.example.json
|-- main_scatter/
|   `-- main_scatter_config.example.json
|-- main_ld2/
|   |-- main_ld2_from_main_scale_config.example.json
|   `-- main_ld2_from_main_scatter_config.example.json
|-- main_tem_analysis/
|   `-- main_tem_analysis_config.example.json
`-- main_valid/
    `-- main_valid_config.example.json
```

Config selection follows the script loaders:

| Script | Config selection |
| --- | --- |
| `main_LD1_v3.m` | `PFAL_MAIN_LD1_CONFIG`, or `config/main_ld1/main_ld1_config.local.json` when the environment variable is unset |
| `main_scale_v2.m` | `config/main_scale/main_scale_config.local.json` |
| `main_scatter_v8.m` | `config/main_scatter/main_scatter_config.local.json` |
| `main_LD2_v3.m` | `PFAL_MAIN_LD2_CONFIG` is required |
| `main_tem_analysis_v1.m` | `PFAL_MAIN_TEM_ANALYSIS_CONFIG`, or `config/main_tem_analysis/main_tem_analysis_config.local.json` when the environment variable is unset |
| `main_valid_v3.m` | `PFAL_MAIN_VALID_CONFIG`, or `config/main_valid/main_valid_config.local.json` when the environment variable is unset |

Typical run order:

```matlab
main_LD1_v3
main_scatter_v8      % default bivariate calibration branch
setenv('PFAL_MAIN_LD2_CONFIG', ...
    'config/main_ld2/main_ld2_from_main_scatter_config.local.json')
main_LD2_v3
```

To use correlation-based scaling instead, run `main_scale_v2` in place of
`main_scatter_v8` and select
`config/main_ld2/main_ld2_from_main_scale_config.local.json` for LD2.

The TEM analysis runs separately with `main_tem_analysis_v1` after its local
config and source datasets have been prepared.

Validation also has `local.file`, `local`, and `local.webtest` profiles. The
web-test profile uses compact fixtures created by:

```matlab
UTILS.CREATE_MAIN_VALID_WEBTEST_FIXTURES
```

Serve `data/main_valid_webtest/server/` on localhost port 8765, select
`main_valid_config.local.webtest.json`, and run `main_valid_v3` to exercise the
URL fallback without transferring the full experimental workspace.

## Repository Map

```text
PFAL/
|-- main_LD1_v3.m          first-stage aggregate-library generation
|-- main_scale_v2.m        correlation-based LD1 library scaling
|-- main_scatter_v8.m      LD1 library sampling/scaling; bivariate by default
|-- main_LD2_v3.m          second-stage post-flame agglomeration
|-- main_tem_analysis_v1.m TEM measurement analysis and figure generation
|-- main_valid_v3.m        LD2/TEM/effective-density validation figures
|-- main_*                 validation, shielding, collapse, and utility scripts
|-- post_*                 post-processing scripts for generated results
|-- config/                JSON templates and ignored local configs
|-- data/                  local inputs and generated workflow handoff files
|-- results/               LD2 workspaces and exported TEM figures
|-- inputs/                legacy tab-delimited parameter files
|-- +PAR/                  particle initialization, sizing, projection, geometry
|-- +TRANSP/               transport properties, mobility, marching, boundaries
|-- +COL/                  collision detection, connection, aggregate growth
|-- +UTILS/                config, checkpoint-resume, plotting, fitting, and IO helpers
|-- +VIS/                  visualization routines
|-- @AGG/                  aggregate class utilities
`-- +DEPOT/                archived or experimental scripts, including legacy PFA
```

## Background

Discrete element Langevin dynamics in PFAL follows these repeated operations:

1. Read physical, domain, transport, and output parameters.
2. Initialize primary-particle sizes and aggregate morphology.
3. Randomly initialize particle locations and velocities.
4. Compute mobility and Brownian transport properties.
5. March positions and velocities under Brownian motion, drag, and inertia.
6. Apply periodic boundary conditions.
7. Detect collisions and connect touching aggregates.
8. Update geometry, mobility, and saved diagnostics.

## References

1. Suresh, V., & Gopalakrishnan, R. (2021). Tutorial: Langevin Dynamics
   methods for aerosol particle trajectory simulations and collision rate
   constant modeling. Journal of Aerosol Science, 155, 105746.
2. Heine, M. C., & Pratsinis, S. E. (2007). Brownian coagulation at high
   concentration. Langmuir, 23(19), 9882-9890.
3. Brasil, A. M., Farias, T. L., & Carvalho, M. G. (1999). A recipe for image
   characterization of fractal-like aggregates. Journal of Aerosol Science,
   30(10), 1379-1389.
4. Olfert, J., & Rogak, S. (2019). Universal relations between soot effective
   density and primary particle size for common combustion sources. Aerosol
   Science and Technology, 53(5), 485-492.
5. Baldelli, A., Trivanovic, U., Corbin, J. C., Lobo, P., Gagne, S., Mille,
   J. W., and Rogak, S. (2020). Typical and atypical morphology of
   non-volatile particles from a diesel and natural gas marine engine. Aerosol
   and Air Quality Research, 20(4), 730-740.
6. Nikookar, H., Sipkens, T. A., & Rogak, S. N. (2025). Simulating the effect
   of post-flame agglomeration on the structure of soot. Aerosol Science and
   Technology, 59(1), 1-15.
