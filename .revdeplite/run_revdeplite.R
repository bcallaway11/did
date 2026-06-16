# =============================================================================
# Title: revdeplite reverse dependency check
# Description: Local reverse dependency check using revdeplite (faster than
#              revdepcheck). Run from the package root before CRAN submission.
#              CI uses revdepcheck (full) via .github/workflows/revdep-check.yml.
# Author: Brant Callaway
# Last update: 2026-06-13
# Date created: 2026-06-13
# =============================================================================

revdeplite::revdeplite(
  github_deps = character(0)
)
