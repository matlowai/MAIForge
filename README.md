# MAIForge

Experimental procedural tools for trees, planets and terrain. The repository
contains the upstream generators as well as the refinement code and required
calibration inputs. It does not require a game engine or third-party art packs.

## What is here

- Deterministic tree structure, annual woody growth, geometry and change tracking.
- Native planetary geology, surface processes, hydrology, stratigraphy and climate.
- Regional drainage, relief and detailed terrain-tile generation.
- Geological material models and source-bound generation/checkpoint evidence.
- Small reproducible examples and targeted tests.

These are research and preview models. Supported processes, assumptions and gaps
are recorded in their outputs. Much of the parameter data consists of declared
project priors; a citation does not make a chosen value a measured calibration.

## Quick start

Install a recent stable Rust toolchain, Python 3 and Bash. This candidate was
validated with rustc 1.98.0 and Python 3.14.4 on Linux; older toolchains have not
been qualified for the complete pipeline. Cargo
fetches the dependencies pinned in the lockfiles; the tree tests use only Python's
standard library. Run commands from the repository root.

```sh
bash scripts/check.sh
bash scripts/demo.sh
```

The demo generates a small planet and full column inputs, builds regional drainage,
runs one bounded relief step, refines one tile and writes an engine-neutral
`terrain.obj`. Output goes into a new `generated/demo-*` directory. The example
uses an included recipe and chooses its request from actual generated dry land.
It does not load a saved private world or an art pack. You can also use the native
`pf-lab simulate`, `pf-region-skeleton`, `pf-region-relief` and `pf-detail` commands,
or call the Rust/Python libraries directly.

**Read `demo.json` before using the result.** The current short example completes
its data path but fails several physical/terrain-quality diagnostics. In particular,
refinement is not a calibrated equilibrium solve and drainage/coast/mean-height
checks are not all satisfied. A generated OBJ is not acceptance of those models.
The original diagnostics and logs are retained. To require all reported pass checks:

```sh
python3 scripts/export_terrain.py generated/<your-run> --strict
```

That strict check currently fails for the documented short example. The toolkit
is useful for inspecting and developing the pipeline; it is not a finished game
or a production terrain guarantee.

## Architecture

The core owns generation, stable identities, physical units and model evidence.
Engine adapters consume its geometry/fields, handle coordinate conversion, upload,
rendering and collision representations. Generation rules should not be duplicated
inside an adapter. A Three.js browser setup as an optional Git submodule is planned
for a later phase; it is not included in this initial source release.

Large generated worlds are outputs/caches, not required repository assets. Models,
textures or packs sourced from third parties are not bundled. See DATA_NOTICES.md
for parameter provenance and THIRD_PARTY_NOTICES.md for software dependencies.

## Tests and resource limits

`scripts/check.sh` runs the portable tree, grid, fines, material, river and producer/checkpoint tests.
Some debug producer/checkpoint tests require a 16 MiB thread stack; the script and
local Cargo configuration set RUST_MIN_STACK accordingly. This changes the harness
resource limit, not assertions. Research tests elsewhere may require additional
generated inputs or substantially more compute; no full research-suite pass is
claimed by the quick check. Historical tests tied to unavailable consumer adapters
or captured fixtures are not represented as passing public tests.

CARGO_BUILD_JOBS defaults to 4 in the scripts. CARGO_TARGET_DIR, MAIFORGE_TMP_DIR and
PE_EVIDENCE_DIR can place builds and test evidence elsewhere. Existing demo output
is never silently overwritten. The native pipeline is a preview-tier implementation;
solved-mantle and full game/runtime integration are not demonstrated by this demo.

## License and alternative terms

MAIForge project code is offered under GNU AGPL version 3 only. See LICENSE.
Commercial use is permitted under that license. If AGPL does not fit your project,
[please say hello to **MATLOWAI**](https://github.com/matlowai/MAIForge/issues).

**A separate license may be free of charge or offered on terms we work out together**,
depending on the project, the people involved and the intended use. Individuals,
hobbyists and small teams are welcome. Please do not assume you need a big budget
to start a conversation.

MATLOWAI handles official alternative-license arrangements, including proprietary
integration. Any exception requires a separate written agreement; third-party
rights and dependency licenses remain applicable. See [LICENSING.md](LICENSING.md).


## Contributions

Before submitting a pull request, read [CONTRIBUTING.md](CONTRIBUTING.md) and the
[MATLOWAI Contributor Agreement v1.0](docs/legal/CONTRIBUTOR_AGREEMENT_v1.0.md).
Submitting an acknowledged PR grants MATLOWAI additional licensing permissions
for your contribution, including separate commercial/proprietary licensing;
you retain your copyright. The PR template records that acknowledgement.
