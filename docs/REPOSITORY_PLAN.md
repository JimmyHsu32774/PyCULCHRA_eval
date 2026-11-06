# Proposed repository organization

Organize shared code by function; record contributors in provenance and commits.
The following is a proposal for team discussion, not a confirmed assignment.
Implemented so far: io/pdb.py, geometry.py, libraries/backbone.py and
reconstruction/backbone.py. Other modules below remain proposed.

```text
src/pyculchra/
  io/pdb.py           # Y4; later add structure adapters/writers
  geometry.py         # Y5/Y6; later add correspondence/RMSD/dihedral helpers
  records.py          # proposed shared record validation
  libraries/          # proposed versioned backbone/side-chain table readers
  reconstruction/     # proposed CA optimization, backbone, candidates,
                      # placement, selection and correction modules
  evaluation/         # proposed metrics, coverage/oracles, summaries
  pipeline.py         # proposed orchestration of named stages
  reference.py        # proposed external C runner, separate from Python engine
tests/                # small unit tests and source-derived stage fixtures
notebooks/            # demonstrations and individual learning notebooks
scripts/              # proposed benchmark and command-line entry points
docs/                 # contracts, source mapping, provenance and decisions
data/fixtures/        # small reviewed fixtures plus metadata
results/              # generated output, ignored by Git
```

## Proposed integration responsibilities

- Yitian: PDB/geometric helpers and backbone reconstruction.
- Jimmy: side-chain library, candidate retrieval/placement and trace capture.
- Caroline: CA optimization, source selection and structure correction.
- Tiffany: C reference runner, benchmark runs and evaluation/reporting.
- Shared: record contracts, stage adapters, end-to-end integration and reviews.

This follows the local exercise contract as a starting point; workload and final
ownership must be agreed by the four members. Do not enforce unconfirmed reviewers
through CODEOWNERS.

## Shared interfaces

Coordinates are finite (N,3) arrays in Angstroms. Structure records align
atom_keys, residue_names (one per atom) and coords. CA records have one CA per
eligible residue; segments explicitly partition row indices in chain order.
Never infer continuity from adjacent array rows alone.

Keep native reference, original CA input, optimized CA, backbone, candidate pools,
selected structure and corrected structure distinct. Reconstruction must not
receive native non-CA coordinates. Preserve candidate IDs and source/frame metadata.
The restricted Y4 result requires an adapter to the full shared structure schema.

Use one declared backbone/CA fit for native evaluation and apply it consistently.
C/Python fidelity uses original-frame coordinates without fitting. Y5 fragment
placement is an internal reconstruction operation, distinct from evaluation fitting.
Y6 frame orientation must be reconciled with each original C placement convention.

## Collaboration and acceptance

Use short feature branches, e.g. member/backbone, with focused pull requests;
merge after another member reviews interfaces, tests and provenance. Keep reusable
implementation in src, with notebooks importing it. Each stage should have a small
normal fixture, boundary/error cases and a source-derived fixture where relevant.

Before integrating new algorithms, agree on PIPELINE_CONTRACT.md and pin the
upstream commit/table versions. Store upstream license and acquisition instructions
when adding source data; do not assign a new project license on behalf of the team.
Record actual human/AI contributions and pending review truthfully.

Implement and validate in order: input/records -> backbone -> candidate retrieval
and placement -> selection/correction -> full pipeline -> evaluation/experiments.
Run unit tests first; source fidelity and real-protein reconstruction need separate
checks. Synthetic tests and empty-input branches are not scientific validation.
Store small reviewed fixtures with input hashes/options/source versions; keep large
proteins and generated benchmark results outside normal Git history.
