# APOL1 Pysam-Style Proof

`bioscripts/apol1-pysam-proof.py` is the first BioScript assay-shaped script
using:

```python
from bioscript import pysam
```

The current proof opens a CRAM file, fetches the three APOL1 regions, and emits
per-site depth rows. It does not yet produce the same APOL1 `G0/G1/G2`
classification as `bioscripts/apol1.py`.

## Missing Helpers Before Output Parity

To compare the pysam-style proof against the existing high-level lookup assay,
the shim needs read-level helpers equivalent to common pysam workflows:

- base at a reference coordinate
- deletion support across a reference span
- CIGAR-aware query/reference projection
- optional base quality filtering
- clear representation for no-call vs no-coverage

The existing CRAM backend already has SNP and indel pileup logic for the
high-level `GenotypeStore` path. The next implementation step should move or
wrap that logic so `bioscript-libs::pysam` can expose it through read/pileup
objects without duplicating the genomics rules.

