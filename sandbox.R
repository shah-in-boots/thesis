library(tidymodels)
library(tidyverse)
library(multilevelmod)
library(broom.mixed)
library(card)
library(octomod)

# Make simple cosinors
tar_load(twins_clinical)
tar_load(twins_cosinors)
df <- inner_join(twins_clinical, twins_cosinors$single, by = "patid")

om <-
	octomod() %>%
	core(df) %>%
	arm(
		title = "test",
		plan = global_cfr ~ mesor + amp1 + phi1,
		approach = linear_reg() %>% set_engine("lm"),
		pattern = "sequential",
		strata = "outcomes"
	) %>%
	equip()
