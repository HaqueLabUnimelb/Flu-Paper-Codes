# Check of the repository

```bash
FLU_RSCRIPT=/path/to/Rscript tests/validate_repository.sh
```

The script checks, without any data: required files, relative links, shell and R syntax, absence of
credentials, e-mail addresses and data objects, the figure manifest and checksums, and that every file of
`original_scripts/` equals the md5 of its manifest. It prints one line per check and changes nothing. The
result of the last run is in [`../logs/validation_summary.txt`](../logs/validation_summary.txt). Tests of the
analyses on the original data: [`../docs/REPRODUCIBILITY.md`](../docs/REPRODUCIBILITY.md).
