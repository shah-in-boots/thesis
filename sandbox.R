library(tidymodels)
library(tidyverse)
library(multilevelmod)
library(broom.mixed)
library(card)
library(octomod)

tar_load(twins_models)
models <- twins_models
models$equipment %>%
	bind_rows(.id = "arm") %>%
	select(outcomes, test_num, level, tidied) %>%
	filter(level == 7) %>%
	unnest(tidied) %>%
	filter(term == "rr")
