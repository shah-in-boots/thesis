library(targets)
library(tarchetypes)

# Functions
source("R/options.R")

# Set target-specific options such as packages.
tar_option_set(
	packages = c(
		# Personal
		"card", "marksman",
		# Tidy
		"tidyverse", "tidymodels", "gt", "gtsummary", "ggdag",
		# Stats
		"lme4", "Hmisc"
	),
	error = "save"
)

# Define targets
targets <- list(
	# Project overview
	tar_render(overview, "R/overview.Rmd")
)

# End with a call to tar_pipeline() to wrangle the targets together.
# This target script must return a pipeline object.
tar_pipeline(targets)
