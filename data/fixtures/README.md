# Y7 source-derived fixture

Files copied unchanged from the local pinned PULCHRA reference project.
Protein: 1UBQ chain A, 76 CA atoms.
Upstream commit: 23a3e4afc01ab48fe81e5d4a8d5902352d8f0977.
Command: pulchra -s -v y7_1ubq_ca.pdb (CA preserved, side chains disabled).
The JSON records source/library/executable hashes, compiler, native coordinates,
comparison tolerance and the restoration of the blank C-output chain field from
input identity. Original local paths in captured stdout are historical provenance.
See third_party/pulchra/SOURCE.md for the table source and license.

Tests compare Python and recorded C coordinates in the original frame, without
alignment. This does not rerun the C executable. A second protein and C-derived
gap/terminal fixtures remain pending; synthetic segment tests are separate.
