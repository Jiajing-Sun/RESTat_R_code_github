#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"

Rscript "$ROOT/tests/parse_all_R_files.R"
Rscript "$ROOT/tests/check_no_private_artifacts.R"
Rscript "$ROOT/tests/smoke_critical_values.R"
Rscript "$ROOT/tests/smoke_alt_critical_values.R"
