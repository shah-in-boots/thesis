library(targets)
library(tarchetypes)

# Functions
source("R/packages.R")

# Set target-specific options such as packages.
tar_option_set()

# Define targets
targets <- list(
	# Make project plan
)

# End with a call to tar_pipeline() to wrangle the targets together.
# This target script must return a pipeline object.
tar_pipeline(targets)
