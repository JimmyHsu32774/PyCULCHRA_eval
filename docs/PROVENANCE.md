# Yitian Y4–Y6 provenance

Source: `Yitian_Exercises.ipynb`, functions `parse_pdb_atoms` (Y4),
`fit_rigid_transform` (Y5), and `make_local_frame` (Y6).
Function bodies were extracted without algorithm changes.

The source notebook labels Y5/Y6 implementation as Yitian, with AI assistance
(opencode/DeepSeek); human review is PENDING. Y4 retains a PENDING implementation-author
field in its source cell; its actual implementation authorship needs confirmation.
The original exercise specifications/scaffolds were AI-authored.
This extraction, packaging, documentation and additional Y4 tests were prepared
with OpenAI Codex. Existing Y5/Y6 test cells are labelled Yitian in the notebook;
the notebook also records AI assistance. No human review is claimed here.

These helpers do not establish full PULCHRA equivalence. Y4 is a restricted
ATOM-only parser; Y5 uses SVD rather than the original iterative fitting routine;
Y6 is a generic right-handed frame whose exact source convention needs review.
No upstream PULCHRA source or tables are bundled in this contribution.
