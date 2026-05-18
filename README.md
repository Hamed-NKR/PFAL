Post-Flame Agglomeration Algorithm (PFAL)
===

PFAL generates soot-like fractal aggregates with a staged Langevin dynamics
workflow. The current workflow is centered on three scripts:

1. `main_LD1_v3.m` generates a raw first-stage aggregate library.
2. `main_scale_v2.m` or `main_scatter_v8.m` calibrates that library to target
   primary-particle and projected-area distributions.
3. `main_LD2_v3.m` runs second-stage post-flame agglomeration on the calibrated
   aggregate population.

The physical motivation is to model soot aggregates that first form as compact
diffusion-limited clusters, then acquire realistic size distributions, and then
continue agglomerating after dilution/cooling. This gives a controlled way to
study hybrid soot structures with nonuniform primary-particle sizes.

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
correlation-based scaling        bivariate scatter sampling
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
              results/main_ld2_*
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

`main_scale_v2.m` rescales and filters LD1 aggregates against the combined
Brasil et al. and Olfert-Rogak correlations. It writes:

```text
data/main_scale/scaled_aggs_for_LD2_from_main_scale.mat
```

`main_scatter_v8.m` samples a bivariate projected-area/primary-particle size
distribution, assigns LD1 aggregate seeds to those targets, rescales them, and
writes:

```text
data/main_scatter/scaled_aggs_for_LD2_from_main_scatter.mat
```

Both outputs use the variable `pars_out`, which is the input contract expected
by LD2.

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

## Configuration And Running

Each main workflow stage has an example config committed under `config/`.
Machine-specific `.local*.json` configs are ignored by Git.

```text
config/
|-- main_ld1/
|   |-- main_ld1_config.example.json
|   |-- main_ld1_config.local.json
|   |-- main_ld1_config.local.file.json
|   `-- main_ld1_config.local.webtest.json
|-- main_scale/
|   `-- main_scale_config.example.json
|-- main_scatter/
|   `-- main_scatter_config.example.json
`-- main_ld2/
    |-- main_ld2_from_main_scale_config.example.json
    `-- main_ld2_from_main_scatter_config.example.json
```

Typical run order:

```matlab
main_LD1_v3
main_scale_v2        % or main_scatter_v8
main_scatter_v8      % optional alternative branch
setenv('PFAL_MAIN_LD2_CONFIG', ...
    'config/main_ld2/main_ld2_from_main_scale_config.local.json')
main_LD2_v3
```

Use `PFAL_MAIN_LD1_CONFIG` to run LD1 with a non-default config:

```matlab
setenv('PFAL_MAIN_LD1_CONFIG', ...
    'config/main_ld1/main_ld1_config.local.webtest.json')
main_LD1_v3
```

## Repository Map

```text
PFAL/
|-- main_LD1_v3.m          first-stage aggregate-library generation
|-- main_scale_v2.m        correlation-based LD1 library scaling
|-- main_scatter_v8.m      bivariate LD1 library sampling/scaling
|-- main_LD2_v3.m          second-stage post-flame agglomeration
|-- main_*                 validation, shielding, collapse, and utility scripts
|-- post_*                 post-processing scripts for generated results
|-- config/                example and local JSON configs
|-- data/                  generated MAT libraries for workflow handoff
|-- results/               LD2 checkpoints and final workspaces
|-- inputs/                legacy tab-delimited parameter files
|-- +PAR/                  particle initialization, sizing, projection, geometry
|-- +TRANSP/               transport properties, mobility, marching, boundaries
|-- +COL/                  collision detection, connection, aggregate growth
|-- +UTILS/                config loaders, plotting helpers, fitting, IO helpers
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
3. Olfert, J., & Rogak, S. (2019). Universal relations between soot effective
   density and primary particle size for common combustion sources. Aerosol
   Science and Technology, 53(5), 485-492.
4. Baldelli, A., Trivanovic, U., Corbin, J. C., Lobo, P., Gagne, S., Mille,
   J. W., and Rogak, S. (2020). Typical and atypical morphology of
   non-volatile particles from a diesel and natural gas marine engine. Aerosol
   and Air Quality Research, 20(4), 730-740.
5. Nikookar, H., Sipkens, T. A., & Rogak, S. N. (2025). Simulating the effect
   of post-flame agglomeration on the structure of soot. Aerosol Science and
   Technology, 59(1), 1-15.
