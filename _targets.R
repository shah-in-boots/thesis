library(targets)
library(tarchetypes)

# Functions
source("R/options.R")
source("R/biobank.R")
source("R/mims.R")
source("R/twins.R")
source("R/thesis.R")

# Set target-specific options such as packages.
tar_option_set(

	packages = c(
		# Personal
		"card", "octomod",
		# Tidy statistics
		"tidyverse", "tidymodels", "gt", "gtsummary", "labelled",
		# Publication/presentation
		"knitr", "kableExtra", "RefManageR", "bibtex",
		# Stats
		"lme4", "Hmisc", "survival", "multilevelmod", "broom.mixed",
		# Graphing
		"ggdag", "survminer",
		# Helpers
		"magrittr", "magick", "equatiomatic"
	),

	error = "continue"
)

# Define targets
targets <- list(

	### Research Data and Reference files ---------------------------------------

	# Bibliography
	tar_file(file_bib, "../bibliography/Neurocardiology.bib"),

	# Biobank
	tar_file(file_biobank_clinical, "../biobank/_targets/objects/clinical"),
	tar_target(biobank_clinical, readRDS(file_biobank_clinical)),
	tar_file(file_biobank_labels, "../biobank/_targets/objects/labels"),
	tar_target(biobank_labels, readRDS(file_biobank_labels)),
	tar_file(file_biobank_ecg, "../biobank/_targets/objects/ecg"),
	tar_target(biobank_ecg, readRDS(file_biobank_ecg)),

	# MIMS
	tar_file(file_mims_raw, "../mims/_targets/objects/proc"),
	tar_target(mims_raw, readRDS(file_mims_raw)),
	tar_file(file_mims_clinical, "../mims/_targets/objects/clinical"),
	tar_target(mims_clinical, readRDS(file_mims_clinical)),
	tar_file(file_mims_outcomes, "../mims/_targets/objects/outcomes"),
	tar_target(mims_outcomes, readRDS(file_mims_outcomes)),

	# Twins
	tar_file(file_twins_clinical, "../twins/_targets/objects/clinical"),
	tar_target(twins_clinical, readRDS(file_twins_clinical)),
	tar_file(file_twins_ecg, "../twins/_targets/objects/ecg"),
	tar_target(twins_ecg, readRDS(file_twins_ecg)),
	tar_file(file_twins_outcomes, "../twins/_targets/objects/outcomes"),
	tar_target(twins_outcomes, readRDS(file_twins_outcomes)),
	tar_file(file_twins_cosinors, "../twins/_targets/objects/cosinors"),
	tar_target(twins_cosinors, readRDS(file_twins_cosinors)),

	### Analysis ---------------------------------------------------------------

	# Overview
	tar_target(diagrams, draw_diagrams()),

	# Biobank
	tar_target(biobank_tables, make_biobank_tables(biobank_clinical, biobank_ecg, biobank_labels)),

	# MIMS
	tar_target(mims_tables, make_mims_tables(mims_clinical, mims_raw)),
	tar_target(mims_figures, make_mims_figures(mims_clinical, mims_outcomes)),
	tar_target(mims_models, make_mims_models(mims_clinical)),
	tar_target(mims_survival, make_mims_survival(mims_clinical, mims_outcomes)),
	tar_target(mims_reports, report_mims_models(mims_models, mims_survival)),

	# Twins
	tar_target(twins_tables, make_twins_tables(twins_clinical, twins_ecg)),
	tar_target(twins_models, make_twins_models(twins_clinical, twins_ecg)),
	tar_target(twins_survival, make_twins_survival(twins_clinical, twins_ecg, twins_outcomes, twins_cosinors)),
	tar_target(twins_circadian, make_twins_circadian(twins_clinical, twins_cosinors)),
	tar_target(twins_reports, report_twins_models(twins_models, twins_survival, twins_circadian)),

	### Publication / Presentation ---------------------------------------------

	# Project overview
	#tar_render(overview, "R/overview.Rmd"),           # Latest Version: 01/01/21

	# Thesis Presentation
	#tar_render(defense, "R/defense.Rmd"),             # Latest Version: 02/08/21

	# Dissertation
	tar_file(index, "./index.Rmd"),
	tar_target(paths, list.files(path = "./thesis/", pattern = "*.Rmd", full.names = TRUE)),
	tar_target(chapters, paths, format = "file", pattern = map(paths)),
	tar_target(thesis, write_dissertation(index, chapters, diagrams))

)
