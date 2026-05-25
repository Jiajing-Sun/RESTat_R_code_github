.PHONY: check smoke requirements

requirements:
	Rscript requirements.R --check

check:
	Rscript tests/parse_all_R_files.R
	Rscript tests/check_no_private_artifacts.R

smoke: check
	Rscript tests/smoke_critical_values.R
	Rscript tests/smoke_alt_critical_values.R
