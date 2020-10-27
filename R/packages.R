# Load all your packages before calling make().

# Setup / basics
library(tidyverse)
library(tidymodels)
library(gt)

# Statistical tools
library(Hmisc)
library(lme4)

# Tidying 
library(magrittr)
library(janitor)
library(readxl)

# Personal
library(card)

# Conflicts
library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("summarize", "dplyr")
