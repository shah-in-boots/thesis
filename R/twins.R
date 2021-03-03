# Summary tables
make_twins_tables <- function(clinical, ecg) {

	# Table 1
	one <-
		clinical %>%
		select(study, age, bmi, race, smoking, prevchd, chf, hptn, dm, ptsd, sad_bin, pet_bin) %>%
		mutate(study = factor(study, levels = c("THS1", "SAVEIT", "THS2", "ETSF"))) %>%
		tbl_summary(
			by = study,
			missing = "no",
			value = list(c(smoking, prevchd, chf, hptn, dm, ptsd, sad_bin, pet_bin) ~ "0")
		) %>%
		modify_header(update = list(
			label ~ '**Characteristic**',
			stat_1 ~ '**ETSF** <br> N = 279',
			stat_2 ~ '**SAVEIT** <br> N = 206',
			stat_3 ~ '**THS1** <br> N = 360',
			stat_4 ~ '**THS2** <br> N = 164'
		))

	# Table on HRV findings
	hrv <-
		ecg %>%
		select(study, n_nmean, sdnn, rmssd, pnn50, ulf, vlf, lf, hf, lfhf, ttlpwr, ac, dc, samp_en, ap_en, dyx) %>%
		tbl_summary(
			by = study,
			missing = "no"
		) %>%
		modify_header(update = list(
			label ~ '**ECG Metric**',
			stat_1 ~ '**ETSF** <br> N = 279',
			stat_2 ~ '**SAVEIT** <br> N = 206',
			stat_3 ~ '**THS1** <br> N = 360',
			stat_4 ~ '**THS2** <br> N = 164'
		))

	# Findings
	tables <- list(
		one = one,
		hrv = hrv
	)

	# Return
	tables

}

# Models
make_twins_models <- function(clinical, ecg) {

	# Hourly HRV and psych/ischemia ---------------------------------------

	df <-
		full_join(clinical, ecg, by = c("vetrid", "study", "patid", "pair")) %>%
		filter(hour %in% c(6:7))

	octobeast <-
		octomod() %>%
		core(df) %>%
		# Logistic - MPI, Depression, PTSD
		arm(
			title = "log_lf_adjusted",
			plan = ptsd + sad_bin + pet_bin ~ lf + age + bmi + smoking + prevchd + chf + hptn + dm + (1 | vetrid) + (1 | pair),
			exposure = c("(1 | vetrid)", "(1 | pair)"),
			pattern = "sequential",
			approach = logistic_reg() %>% set_engine("glmer"),
			strata = "hour"
		) %>%
		arm(
			title = "log_hf_adjusted",
			plan = ptsd + sad_bin + pet_bin ~ hf + age + bmi + smoking + prevchd + chf + hptn + dm + (1 | vetrid) + (1 | pair),
			exposure = c("(1 | vetrid)", "(1 | pair)"),
			pattern = "sequential",
			approach = logistic_reg() %>% set_engine("glmer"),
			strata = "hour"
		) %>%
		arm(
			title = "log_vlf_adjusted",
			plan = ptsd + sad_bin + pet_bin ~ vlf + age + bmi + smoking + prevchd + chf + hptn + dm + (1 | vetrid) + (1 | pair),
			exposure = c("(1 | vetrid)", "(1 | pair)"),
			pattern = "sequential",
			approach = logistic_reg() %>% set_engine("glmer"),
			strata = "hour"
		) %>%
		arm(
			title = "log_dyx_adjusted",
			plan = ptsd + sad_bin + pet_bin ~ dyx + age + bmi + smoking + prevchd + chf + hptn + dm + (1 | vetrid) + (1 | pair),
			exposure = c("(1 | vetrid)", "(1 | pair)"),
			pattern = "sequential",
			approach = logistic_reg() %>% set_engine("glmer"),
			strata = "hour"
		) %>%
		arm(
			title = "log_ac_adjusted",
			plan = ptsd + sad_bin + pet_bin ~ ac + age + bmi + smoking + prevchd + chf + hptn + dm + (1 | vetrid) + (1 | pair),
			exposure = c("(1 | vetrid)", "(1 | pair)"),
			pattern = "sequential",
			approach = logistic_reg() %>% set_engine("glmer"),
			strata = "hour"
		) %>%
		# Linear models - CFR
		arm(
			title = "lin_lf_adjusted",
			plan = global_cfr ~ lf + age + bmi + smoking + prevchd + chf + hptn + dm + (1 | vetrid) + (1 | pair),
			exposure = c("(1 | vetrid)", "(1 | pair)"),
			pattern = "sequential",
			approach = linear_reg() %>% set_engine("lmer"),
			strata = "hour"
		) %>%
		arm(
			title = "lin_hf_adjusted",
			plan = global_cfr ~ hf + age + bmi + smoking + prevchd + chf + hptn + dm + (1 | vetrid) + (1 | pair),
			exposure = c("(1 | vetrid)", "(1 | pair)"),
			pattern = "sequential",
			approach = linear_reg() %>% set_engine("lmer"),
			strata = "hour"
		) %>%
		arm(
			title = "lin_vlf_adjusted",
			plan = global_cfr ~ vlf + age + bmi + smoking + prevchd + chf + hptn + dm + (1 | vetrid) + (1 | pair),
			exposure = c("(1 | vetrid)", "(1 | pair)"),
			pattern = "sequential",
			approach = linear_reg() %>% set_engine("lmer"),
			strata = "hour"
		) %>%
		arm(
			title = "lin_dyx_adjusted",
			plan = global_cfr ~ dyx + age + bmi + smoking + prevchd + chf + hptn + dm + (1 | vetrid) + (1 | pair),
			exposure = c("(1 | vetrid)", "(1 | pair)"),
			pattern = "sequential",
			approach = linear_reg() %>% set_engine("lmer"),
			strata = "hour"
		) %>%
		arm(
			title = "lin_ac_adjusted",
			plan = global_cfr ~ ac + age + bmi + smoking + prevchd + chf + hptn + dm + (1 | vetrid) + (1 | pair),
			exposure = c("(1 | vetrid)", "(1 | pair)"),
			pattern = "sequential",
			approach = linear_reg() %>% set_engine("lmer"),
			strata = "hour"
		) %>%
		equip()

	# Return  -----------------------------------------------------------
	octobeast

}

make_twins_circadian <- function(clinical, cosinors) {

	octobeast <-
		octomod() %>%
		core(left_join(cosinors$single, clinical, by = "patid")) %>%
		arm(
			title = "single_log",
			plan = ptsd + sad_bin + pet_bin ~ mesor + amp1 + phi1 + age + bmi + race + smoking + prevchd + chf + hptn + dm + (1 | vetrid) + (1 | pair),
			approach = logistic_reg() %>% set_engine("glmer"),
			exposure = c("(1 | pair)", "(1 | vetrid)"),
			pattern = "direct",
			strata = "outcomes"
		) %>%
		arm(
			title = "single_lin",
			plan = beck_total + global_cfr ~ mesor + amp1 + phi1 + age + bmi + race + smoking + prevchd + chf + hptn + dm + (1 | vetrid) + (1 | pair),
			approach = linear_reg() %>% set_engine("lmer"),
			exposure = c("(1 | pair)", "(1 | vetrid)"),
			pattern = "direct",
			strata = "outcomes"
		) %>%
		equip(which_arms = c("single_lin", "single_log")) %>%
		change_core(left_join(cosinors$multiple, clinical, by = "patid")) %>%
		arm(
			title = "multiple_log",
			plan = ptsd + sad_bin + pet_bin ~ mesor + amp1 + phi1 + amp2 + phi2 + age + bmi + race + smoking + prevchd + chf + hptn + dm + (1 | vetrid) + (1 | pair),
			approach = logistic_reg() %>% set_engine("glmer"),
			exposure = c("(1 | pair)", "(1 | vetrid)"),
			pattern = "direct",
			strata = "outcomes"
		) %>%
		arm(
			title = "multiple_lin",
			plan = beck_total + global_cfr ~ mesor + amp1 + phi1 + amp2 + phi2 + age + bmi + race + smoking + prevchd + chf + hptn + dm + (1 | vetrid) + (1 | pair),
			approach = linear_reg() %>% set_engine("lmer"),
			exposure = c("(1 | pair)", "(1 | vetrid)"),
			pattern = "direct",
			strata = "outcomes"
		) %>%
		equip(which_arms = c("multiple_log", "multiple_lin"))

	# Return octobeast
	octobeast

}

# Outcomes
make_twins_survival <- function(clinical, ecg, outcomes) {

	# Clinical data
	df <- inner_join(ecg, clinical, by = "vetrid")

	octobeast <-
		octomod() %>%
		# Overall mortality
		core(inner_join(df, outcomes$death, by = c("vetrid" = "id"))) %>%
		arm(
			title = "hf_death",
			plan = Surv(stop, status) ~ hf + pet_bin + age + bmi + race + smoking + prevchd + chf + hptn + dm + sad_bin + ptsd,
			approach = cox_reg() %>% set_engine("survival", method = "breslow"),
			pattern = "sequential",
			strata = "hour"
		) %>%
		arm(
			title = "lf_death",
			plan = Surv(stop, status) ~ lf + pet_bin + age + bmi + race + smoking + prevchd + chf + hptn + dm + sad_bin + ptsd,
			approach = cox_reg() %>% set_engine("survival", method = "breslow"),
			pattern = "sequential",
			strata = "hour"
		) %>%
		arm(
			title = "vlf_death",
			plan = Surv(stop, status) ~ vlf + pet_bin + age + bmi + race + smoking + prevchd + chf + hptn + dm + sad_bin + ptsd,
			approach = cox_reg() %>% set_engine("survival", method = "breslow"),
			pattern = "sequential",
			strata = "hour"
		) %>%
		arm(
			title = "dyx_death",
			plan = Surv(stop, status) ~ dyx + pet_bin + age + bmi + race + smoking + prevchd + chf + hptn + dm + sad_bin + ptsd,
			approach = cox_reg() %>% set_engine("survival", method = "breslow"),
			pattern = "sequential",
			strata = "hour"
		) %>%
		arm(
			title = "ac_death",
			plan = Surv(stop, status) ~ ac + pet_bin + age + bmi + race + smoking + prevchd + chf + hptn + dm + sad_bin + ptsd,
			approach = cox_reg() %>% set_engine("survival", method = "breslow"),
			pattern = "sequential",
			strata = "hour"
		) %>%
		arm(
			title = "rr_death",
			plan = Surv(stop, status) ~ rr + pet_bin + age + bmi + race + smoking + prevchd + chf + hptn + dm + sad_bin + ptsd,
			approach = cox_reg() %>% set_engine("survival", method = "breslow"),
			pattern = "sequential",
			strata = "hour"
		) %>%
		equip(which_arms = c("hf_death", "lf_death", "vlf_death", "dyx_death", "ac_death", "rr_death")) %>%
		# CVD mortality
		change_core(inner_join(df, outcomes$cvd, by = c("vetrid" = "id"))) %>%
		arm(
			title = "hf_cvd",
			plan = Surv(stop, status) ~ hf + pet_bin + age + bmi + race + smoking + prevchd + chf + hptn + dm + sad_bin + ptsd,
			approach = cox_reg() %>% set_engine("survival", method = "breslow"),
			pattern = "sequential",
			strata = "hour"
		) %>%
		arm(
			title = "lf_cvd",
			plan = Surv(stop, status) ~ lf + pet_bin + age + bmi + race + smoking + prevchd + chf + hptn + dm + sad_bin + ptsd,
			approach = cox_reg() %>% set_engine("survival", method = "breslow"),
			pattern = "sequential",
			strata = "hour"
		) %>%
		arm(
			title = "vlf_cvd",
			plan = Surv(stop, status) ~ vlf + pet_bin + age + bmi + race + smoking + prevchd + chf + hptn + dm + sad_bin + ptsd,
			approach = cox_reg() %>% set_engine("survival", method = "breslow"),
			pattern = "sequential",
			strata = "hour"
		) %>%
		arm(
			title = "dyx_cvd",
			plan = Surv(stop, status) ~ dyx + pet_bin + age + bmi + race + smoking + prevchd + chf + hptn + dm + sad_bin + ptsd,
			approach = cox_reg() %>% set_engine("survival", method = "breslow"),
			pattern = "sequential",
			strata = "hour"
		) %>%
		arm(
			title = "ac_cvd",
			plan = Surv(stop, status) ~ ac + pet_bin + age + bmi + race + smoking + prevchd + chf + hptn + dm + sad_bin + ptsd,
			approach = cox_reg() %>% set_engine("survival", method = "breslow"),
			pattern = "sequential",
			strata = "hour"
		) %>%
		arm(
			title = "rr_cvd",
			plan = Surv(stop, status) ~ rr + pet_bin + age + bmi + race + smoking + prevchd + chf + hptn + dm + sad_bin + ptsd,
			approach = cox_reg() %>% set_engine("survival", method = "breslow"),
			pattern = "sequential",
			strata = "hour"
		) %>%
		equip(which_arms = c("hf_cvd", "lf_cvd", "vlf_cvd", "dyx_cvd", "ac_cvd", "rr_cvd"))

	# Return
	octobeast
}

# Report out twin models
report_twins_models <- function(models, survival, circadian) {

	# Set Up Hourly Data ------------------------------------------------------
	unadj <- 2 # HRV
	adj_demo <- 4 # + age + BMI
	adj_cv <- 9 # + smoking, CAD, CHF, HTN, DM

	hourly <-
		models$equipment %>%
		bind_rows(.id = "arm") %>%
		select(outcomes, test_num, level, tidied) %>%
		unnest(tidied) %>%
		filter(level %in% 6:7) %>%
		filter(term %in% c("lf", "hf", "vlf", "ac", "dyx")) %>%
		filter(test_num %in% c(unadj, adj_demo, adj_cv)) %>%
		group_by(outcomes, test_num, term) %>%
		arrange(-abs(statistic), .by_group = TRUE) %>%
		slice(1) %>%
		ungroup() %>%
		select(-c(group, x, level, effect, std.error, p.value)) %>%
		pivot_wider(names_from = "term", values_from = c("estimate", "conf.low", "conf.high", "statistic"), names_glue = "{term}_{.value}") %>%
		mutate(test_num = factor(test_num, labels = c("Model 1", "Model 2", "Model 3"))) %>%
		mutate(outcomes = factor(outcomes, levels = c("global_cfr", "pet_bin", "ptsd", "sad_bin"), labels = c("Coronary Flow Reserve", "Abnormal MPI", "PTSD", "Depression")))

	# Psych Models -------------------------------
	psych <-
		hourly %>%
		filter(outcomes %in% c("PTSD", "Depression")) %>%
		gt(rowname_col = "test_num", groupname_col = "outcomes") %>%
		cols_merge(columns = starts_with("hf_"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("lf_"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("vlf_"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("dyx_"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("ac"), pattern = "{1} ({2}, {3})") %>%
		fmt_number(columns = everything(), decimals = 2) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
			locations = cells_body(columns = "lf_estimate", rows = abs(lf_statistic) > 2.0)
		) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
			locations = cells_body(columns = "hf_estimate", rows = abs(hf_statistic) > 2.0)
		) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
			locations = cells_body(columns = "vlf_estimate", rows = abs(vlf_statistic) > 2.0)
		) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
			locations = cells_body(columns = "ac_estimate", rows = abs(ac_statistic) > 2.0)
		) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
			locations = cells_body(columns = "dyx_estimate", rows = abs(dyx_statistic) > 2.0)
		) %>%
		cols_label(
			lf_estimate = "LF",
			hf_estimate = "HF",
			vlf_estimate = "VLF",
			ac_estimate = "AC",
			dyx_estimate = md("*Dyx*")
		) %>%
		tab_footnote(
			footnote = "Model 1 = HRV",
			locations = cells_stub(rows = c(1:2))
		) %>%
		tab_footnote(
			footnote = "Model 2 = Model 1 + Age + BMI",
			locations = cells_stub(rows = c(3:4))
		) %>%
		tab_footnote(
			footnote = "Model 3 = Model 2 + Smoking + HTN + Cardiovascular Disease",
			locations = cells_stub(rows = c(5:6))
		)

	# Ischemia Models  ----------------------------------------------------
	mpi <-
		hourly %>%
		filter(outcomes %in% c("Coronary Flow Reserve", "Abnormal MPI")) %>%
		gt(rowname_col = "test_num", groupname_col = "outcomes") %>%
		cols_merge(columns = starts_with("hf_"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("lf_"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("vlf_"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("dyx_"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("ac"), pattern = "{1} ({2}, {3})") %>%
		fmt_number(columns = everything(), decimals = 2) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
			locations = cells_body(columns = "lf_estimate", rows = abs(lf_statistic) > 2.0)
		) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
			locations = cells_body(columns = "hf_estimate", rows = abs(hf_statistic) > 2.0)
		) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
			locations = cells_body(columns = "vlf_estimate", rows = abs(vlf_statistic) > 2.0)
		) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
			locations = cells_body(columns = "ac_estimate", rows = abs(ac_statistic) > 2.0)
		) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
			locations = cells_body(columns = "dyx_estimate", rows = abs(dyx_statistic) > 2.0)
		) %>%
		cols_label(
			lf_estimate = "LF",
			hf_estimate = "HF",
			vlf_estimate = "VLF",
			ac_estimate = "AC",
			dyx_estimate = md("*Dyx*")
		) %>%
		tab_footnote(
			footnote = "Model 1 = HRV",
			locations = cells_stub(rows = c(1:2))
		) %>%
		tab_footnote(
			footnote = "Model 2 = Model 1 + Age + BMI",
			locations = cells_stub(rows = c(3:4))
		) %>%
		tab_footnote(
			footnote = "Model 3 = Model 2 + Smoking + HTN + Cardiovascular Disease",
			locations = cells_stub(rows = c(5:6))
		)

	# Survival -------------------------------------------------------------
	unadj <- 1
	adj_mpi <- 2
	adj_demo <- 6
	adj_cvd <- 10
	adj_sad <- 11

	outcomes <-
		survival$equipment %>%
		bind_rows(.id = "arm") %>%
		filter(str_detect(arm, "death")) %>%
		select(level, test_num, tidied) %>%
		unnest(tidied) %>%
		filter(level %in% c(6:9)) %>%
		filter(term %in% c("ac", "dyx", "hf", "lf", "vlf", "rr")) %>%
		filter(test_num %in% c(unadj, adj_mpi, adj_demo, adj_cvd, adj_sad)) %>%
		mutate(test_num = factor(test_num, labels = paste("Model", 1:5))) %>%
		group_by(test_num, term) %>%
		arrange(-abs(statistic), .by_group = TRUE) %>%
		slice(1) %>%
		ungroup() %>%
		select(test_num, term, estimate, conf.low, conf.high, p.value) %>%
		pivot_wider(
			names_from = term,
			values_from = c(estimate, conf.low, conf.high, p.value),
			names_glue = "{term}_{.value}"
		) %>%
		gt(rowname_col = "test_num") %>%
		cols_merge(columns = starts_with("ac"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("dyx"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("hf"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("lf"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("vlf"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("rr"), pattern = "{1} ({2}, {3})") %>%
		fmt_number(
			columns = contains(c("estimate", "conf")),
			decimals = 2,
			drop_trailing_zeros = TRUE
		) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
			locations = cells_body(columns = starts_with("ac"),
														 rows = ac_p.value < .05)
		) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
			locations = cells_body(columns = starts_with("dyx"),
														 rows = dyx_p.value < .05)
		) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
			locations = cells_body(columns = starts_with("hf"),
														 rows = hf_p.value < .05)
		) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
			locations = cells_body(columns = starts_with("lf"),
														 rows = lf_p.value < .05)
		) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
			locations = cells_body(columns = starts_with("vlf"),
														 rows = vlf_p.value < .05)
		) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
			locations = cells_body(columns = starts_with("rr"),
														 rows = rr_p.value < .05)
		) %>%
		cols_label(
			ac_estimate = "Acceleration Capacity",
			dyx_estimate = "Dyx",
			hf_estimate = "High Frequency HRV",
			lf_estimate = "Low Frequency HRV",
			vlf_estimate = "Very Low Frequency HRV",
			rr_estimate = "RR Intervals"
		)

	# Set up cosinor setup --------------------------------------------------------

	circ <-
		circadian$equipment %>%
		bind_rows(.id = "arm") %>%
		filter(str_detect(arm, "single")) %>%
		filter(level %in% c("hf", "lf", "vlf", "ac", "dyx", "rr")) %>%
		select(outcomes, level, tidied) %>%
		unnest(tidied) %>%
		select(-c(group, effect, p.value, std.error)) %>%
		filter(term %in% c("mesor", "amp1", "phi1")) %>%
		pivot_wider(names_from = "term", values_from = c("estimate", "conf.low", "conf.high", "statistic"), names_glue = "{term}_{.value}") %>%
		mutate(outcomes = factor(outcomes, levels = c("global_cfr", "pet_bin", "ptsd", "sad_bin", "beck_total"), labels = c("Coronary Flow Reserve", "Abnormal MPI", "PTSD", "Depression", "BDI Score"))) %>%
		mutate(level = case_when(
			level == "hf" ~ "High Frequency HRV",
			level == "lf" ~ "Low Frequency HRV",
			level == "vlf" ~ "Very Low Frequency HRV",
			level == "dyx" ~ "Dyx",
			level == "ac" ~ "Acceleration Capacity",
			level == "rr" ~ "RR Intervals"
		))

	# Psych Cosinor --------------------------------------------------------
	psych_cosinor <-
		circ %>%
		group_by(outcomes) %>%
		filter(outcomes %in% c("Depression", "PTSD")) %>%
		gt(rowname_col = "level") %>%
		cols_merge(columns = starts_with("mesor_"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("amp1_"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("phi1_"), pattern = "{1} ({2}, {3})") %>%
		fmt_number(
			columns = starts_with(c("mesor_", "amp1", "phi1")),
			decimals = 2,
			drop_trailing_zeros = TRUE
		) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
			locations = cells_body(columns = starts_with("mesor"),
														 rows = abs(mesor_statistic) > 2.0)
		) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
			locations = cells_body(columns = starts_with("amp1"),
														 rows = abs(amp1_statistic) > 2.0)
		) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
			locations = cells_body(columns = starts_with("phi1"),
														 rows = abs(phi1_statistic) > 2.0)
		) %>%
		cols_label(
			mesor_estimate = "MESOR",
			amp1_estimate = "Amplitude",
			phi1_estimate = "Phi"
		)


	# Ischemia Cosinor --------------------------------------------------------
	ischemia_cosinor <-
		circ %>%
		group_by(outcomes) %>%
		filter(outcomes %in% c("Abnormal MPI", "Coronary Flow Reserve")) %>%
		gt(rowname_col = "level") %>%
		cols_merge(columns = starts_with("mesor_"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("amp1_"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("phi1_"), pattern = "{1} ({2}, {3})") %>%
		fmt_number(
			columns = starts_with(c("mesor_", "amp1", "phi1")),
			decimals = 2,
			drop_trailing_zeros = TRUE
		) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
			locations = cells_body(columns = starts_with("mesor"),
														 rows = abs(mesor_statistic) >= 2.0)
		) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
			locations = cells_body(columns = starts_with("amp1"),
														 rows = abs(amp1_statistic) >= 2.0)
		) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
			locations = cells_body(columns = starts_with("phi1"),
														 rows = abs(phi1_statistic) >= 2.0)
		) %>%
		cols_label(
			mesor_estimate = "MESOR",
			amp1_estimate = "Amplitude",
			phi1_estimate = "Phi"
		)

	# Return -------------------------------------------------------------
	tables <- list(
		psych = psych,
		mpi = mpi,
		outcomes = outcomes,
		psych_cosinor = psych_cosinor,
		ischemia_cosinor = ischemia_cosinor
	)

	# Returning
	tables

}
