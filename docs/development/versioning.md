# Versioning

MOSE follows [Semantic Versioning 2.0.0](https://semver.org/), as ORION, OSLO, FLINT, ATLAS and ICE do:
a release is a tag `vX.Y.Z` on `main`, and its changes are listed in
[CHANGELOG.md](https://github.com/open-hydra/MOSE/blob/main/CHANGELOG.md). Before 1.0.0 a breaking
change raises `Y`, anything else `Z`.

## How MOSE requires its libraries

`CMakeLists.txt` names the releases MOSE is tested with:

| Library | Requirement | How |
|---|---|---|
| ORION | `MOSE_ORION_VERSION 2.1`: 2.1.0 and later, up to 3.0.0 excluded | `semver_dependency` of `cmake/SemVerDeps.cmake` |
| OSLO | `MOSE_OSLO_VERSION 1.1`: 1.1.0 and later, up to 2.0.0 excluded | `semver_dependency` |
| FLINT | `MOSE_FLINT_VERSION 2.0`: 2.0.0 and later, up to 3.0.0 excluded | `semver_dependency` |
| FiNeR | `MOSE_FINER_TAG v2.3.0` | the submodule `lib/FiNeR`, checked against its `VERSION` file |

`cmake/SemVerDeps.cmake` is FLINT's file, copied unmodified: update it from FLINT, never here.
A requirement accepts the releases from its own up to the next change of its leftmost nonzero
number (`2.1` accepts 2.1.0 to 3.0.0 excluded, `0.3` accepts 0.3.x only). Only `vX.Y.Z` tags are
releases; a tag such as `v2.2.0-rc1` is never selected.

FLINT requires ORION and OSLO as well, so the ranges of MOSE for those two must overlap FLINT's: the
configure stops, naming both projects, when they do not.

FiNeR's own libraries are fetched as FiNeR's author does (`fobis fetch`: their default
branches, `MOSE_FINER_DEPS=FREE`). If an upstream change breaks MOSE, `MOSE_FINER_DEPS=BACKUP`
checks out the releases of `MOSE_FINER_BACKUP_PINS`, the ones ATLAS and ICE share.

## Which release a build gets

The release of ORION, OSLO and FLINT is, in order:

1. the clone given with `-DFETCHCONTENT_SOURCE_DIR_<NAME>=<folder>` (`install.sh --include-orion`,
   `--include-oslo`, `--include-flint`), as it is; the configure stops if it is not a release that
   MOSE accepts (`-DMOSE_PIN_CHECK=WARNING` builds it on purpose);
2. the release in `deps.lock`, if it satisfies the requirement;
3. the newest release on GitHub that satisfies it.

The configure checks each release out in `lib/<NAME>` and writes them to `deps.lock`; both are git
ignored. `-DSEMVER_UPDATE=ALL` (`install.sh build --update-deps`) selects the newest ones again.
`deps.lock` is not committed: in hydra, the top of the build tree, hydra chooses the releases.

## The requirements in hydra

Hydra builds MOSE with its own libraries. When hydra selects ORION, OSLO or FLINT with SemVerDeps,
MOSE's requirement is checked against that release, so:

- a new **major** of one of those libraries in MOSE's requirement is a major (before 1.0.0, minor)
  release of MOSE;
- a later **minor** (from ORION `2.1` to `2.2`) is a minor (before 1.0.0, patch) release;
- another FiNeR release is a major (before 1.0.0, minor) release of MOSE.

## Making a release

1. Choose the number from the `[Unreleased]` entries of `CHANGELOG.md`.
2. Rename `[Unreleased]` to `[X.Y.Z] - YYYY-MM-DD` and open a new, empty `[Unreleased]` above it.
3. Merge into `main` through a pull request.
4. Tag the merge commit and push the tag:
   ```bash
   git fetch origin
   git tag vX.Y.Z origin/main
   git push origin vX.Y.Z
   ```

!!! warning "A published tag stays"
    A pushed tag is never moved or deleted: other projects' lock files name it.
