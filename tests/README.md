# Repository validation

Run the lightweight validator from any directory:

```bash
FLU_RSCRIPT=/path/to/Rscript tests/validate_repository.sh
```

`FLU_RSCRIPT` is optional. When it is absent, the validator uses `Rscript`
from `PATH`; if neither is available, R parsing is reported as `NOT TESTED`.
The validator does not read research objects or run an analysis.
