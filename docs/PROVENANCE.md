# Yitian Y4–Y7 provenance

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
Y7 now bundles the original backbone tables and license under third_party/pulchra.

## Y7 addition

Y7 functions and tests were extracted from notebook cells 34–35. The source
labels implementation/test author Yitian, AI-assisted by opencode/DeepSeek;
human reviewer PENDING. Codex performed module extraction, import/path adaptation,
documentation, and a context-manager fix for reading the table hash. Numerical
reconstruction logic is unchanged. The notebook-specific source discovery helper
is replaced by an explicit library path, avoiding machine-dependent discovery.
The C fixture is required by the repository test rather than silently skipped.

Backbone uses descriptor bins, ordinary/PRO tables, SVD placement, end
extrapolation and independent segments. The SVD fit differs from the original C
iterative fitting routine. Segments shorter than five CA are excluded; this is a
documented deviation. The generic Y6 frame is not used by this backbone stage.
The parser is restricted to the pinned table syntax. Input records must already
satisfy the shared identity, finite-coordinate and segment contract; this stage
does not yet provide comprehensive record validation.
