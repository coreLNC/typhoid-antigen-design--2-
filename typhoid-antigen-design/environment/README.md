# Software environments

The analysis used two conda environments on Linux:

| Environment | Used for |
|---|---|
| `antigen` | genome/protein processing, DIAMOND, MAFFT, ranking and classification scripts, Biopython (Shrake-Rupley accessibility), pandas/numpy, matplotlib |
| `sistr` | SISTR serovar prediction (`sistr_cmd`; needs `setuptools<81`) |

Export the exact environments from the machine that produced the results:

```bash
conda env export --from-history -n antigen > environment/antigen.yml     # portable
conda env export -n antigen > environment/antigen.full.yml               # exact builds
conda env export --from-history -n sistr   > environment/sistr.yml
```

Versions of DIAMOND, MAFFT and the NCBI Datasets CLI were not recorded when the analysis was run.
Record them with `diamond --version`, `mafft --version`, `datasets --version` and add them to the README.

The figure scripts only need Python 3, numpy, matplotlib and (for Figure 5) pillow, playwright and 3Dmol.js.
