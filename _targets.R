library(targets)
library(tarchetypes)

# Functions
source("R/packages.R")

# Set target-specific options such as packages.
tar_option_set(
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
