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
			value = list(c(smoking, prevchd, chf, hptn, dm, ptsd, sad_bin, pet_bin) ~ "1")
		) %>%
		modify_header(update = list(
			label ~ "**Characteristic**",
			stat_1 ~ "**THS1**, N = 361",
			stat_2 ~ "**SAVEIT**, N = 206",
			stat_3 ~ "**THS2**, N = 165",
			stat_4 ~ "**ETSF**, N = 280"
		))

	# Table on HRV findings
	hrv <-
		ecg %>%
		select(study, n_nmean, sdnn, rmssd, pnn50, ulf, vlf, lf, hf, lfhf, ttlpwr, ac, dc, samp_en, ap_en, dyx) %>%
		mutate(study = factor(study, levels = c("THS1", "SAVEIT", "THS2", "ETSF"))) %>%
		tbl_summary(
			by = study,
			missing = "no"
		) %>%
		modify_header(update = list(
			label ~ "**ECG/HRV Metric**",
			stat_1 ~ "**THS1**, N = 361",
			stat_2 ~ "**SAVEIT**, N = 206",
			stat_3 ~ "**THS2**, N = 165",
			stat_4 ~ "**ETSF**, N = 280"
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
		filter(hour %in% c(6:7)) %>%
		mutate(bpm = (1 / (rr * 1/1000 * 1 / 60))/10)

	octobeast <-
		octomod() %>%
		core(df) %>%
		# Logistic - MPI, Depression, PTSD
		arm(
			title = "log_lf_adjusted",
			plan = ptsd + sad_bin + pet_bin ~ lf + age + bmi + race + smoking + prevchd + chf + hptn + dm + (1 | vetrid) + (1 | pair),
			exposure = c("(1 | vetrid)", "(1 | pair)"),
			pattern = "sequential",
			approach = logistic_reg() %>% set_engine("glmer"),
			strata = "hour"
		) %>%
		arm(
			title = "log_hf_adjusted",
			plan = ptsd + sad_bin + pet_bin ~ hf + age + bmi + race + smoking + prevchd + chf + hptn + dm + (1 | vetrid) + (1 | pair),
			exposure = c("(1 | vetrid)", "(1 | pair)"),
			pattern = "sequential",
			approach = logistic_reg() %>% set_engine("glmer"),
			strata = "hour"
		) %>%
		arm(
			title = "log_vlf_adjusted",
			plan = ptsd + sad_bin + pet_bin ~ vlf + age + bmi + race + smoking + prevchd + chf + hptn + dm + (1 | vetrid) + (1 | pair),
			exposure = c("(1 | vetrid)", "(1 | pair)"),
			pattern = "sequential",
			approach = logistic_reg() %>% set_engine("glmer"),
			strata = "hour"
		) %>%
		arm(
			title = "log_dyx_adjusted",
			plan = ptsd + sad_bin + pet_bin ~ dyx + age + bmi + race + smoking + prevchd + chf + hptn + dm + (1 | vetrid) + (1 | pair),
			exposure = c("(1 | vetrid)", "(1 | pair)"),
			pattern = "sequential",
			approach = logistic_reg() %>% set_engine("glmer"),
			strata = "hour"
		) %>%
		arm(
			title = "log_ac_adjusted",
			plan = ptsd + sad_bin + pet_bin ~ ac + age + bmi + race + smoking + prevchd + chf + hptn + dm + (1 | vetrid) + (1 | pair),
			exposure = c("(1 | vetrid)", "(1 | pair)"),
			pattern = "sequential",
			approach = logistic_reg() %>% set_engine("glmer"),
			strata = "hour"
		) %>%
		arm(
			title = "log_bpm_adjusted",
			plan = ptsd + sad_bin + pet_bin ~ bpm + age + bmi + race + smoking + prevchd + chf + hptn + dm + (1 | vetrid) + (1 | pair),
			exposure = c("(1 | vetrid)", "(1 | pair)"),
			pattern = "sequential",
			approach = logistic_reg() %>% set_engine("glmer"),
			strata = "hour"
		) %>%
		# Linear models - CFR
		arm(
			title = "lin_lf_adjusted",
			plan = global_cfr ~ lf + age + bmi + race + smoking + prevchd + chf + hptn + dm + (1 | vetrid) + (1 | pair),
			exposure = c("(1 | vetrid)", "(1 | pair)"),
			pattern = "sequential",
			approach = linear_reg() %>% set_engine("lmer"),
			strata = "hour"
		) %>%
		arm(
			title = "lin_hf_adjusted",
			plan = global_cfr ~ hf + age + bmi + race + smoking + prevchd + chf + hptn + dm + (1 | vetrid) + (1 | pair),
			exposure = c("(1 | vetrid)", "(1 | pair)"),
			pattern = "sequential",
			approach = linear_reg() %>% set_engine("lmer"),
			strata = "hour"
		) %>%
		arm(
			title = "lin_vlf_adjusted",
			plan = global_cfr ~ vlf + age + bmi + race + smoking + prevchd + chf + hptn + dm + (1 | vetrid) + (1 | pair),
			exposure = c("(1 | vetrid)", "(1 | pair)"),
			pattern = "sequential",
			approach = linear_reg() %>% set_engine("lmer"),
			strata = "hour"
		) %>%
		arm(
			title = "lin_dyx_adjusted",
			plan = global_cfr ~ dyx + age + bmi + race + smoking + prevchd + chf + hptn + dm + (1 | vetrid) + (1 | pair),
			exposure = c("(1 | vetrid)", "(1 | pair)"),
			pattern = "sequential",
			approach = linear_reg() %>% set_engine("lmer"),
			strata = "hour"
		) %>%
		arm(
			title = "lin_ac_adjusted",
			plan = global_cfr ~ ac + age + bmi + race + smoking + prevchd + chf + hptn + dm + (1 | vetrid) + (1 | pair),
			exposure = c("(1 | vetrid)", "(1 | pair)"),
			pattern = "sequential",
			approach = linear_reg() %>% set_engine("lmer"),
			strata = "hour"
		) %>%
		arm(
			title = "lin_bpm_adjusted",
			plan = global_cfr ~ bpm + age + bmi + race + smoking + prevchd + chf + hptn + dm + (1 | vetrid) + (1 | pair),
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
			plan = ptsd + sad_bin + pet_bin ~ mesor + amp1 + phi1 + (1 | pair),
			approach = logistic_reg() %>% set_engine("glmer"),
			exposure = c("(1 | pair)"),
			pattern = "direct",
			strata = "outcomes"
		) %>%
		arm(
			title = "single_lin",
			plan = beck_total + global_cfr ~ mesor + amp1 + phi1 + (1 | pair),
			approach = linear_reg() %>% set_engine("lmer"),
			exposure = c("(1 | pair)"),
			pattern = "direct",
			strata = "outcomes"
		) %>%
		equip(which_arms = c("single_lin", "single_log")) %>%
		change_core(left_join(cosinors$multiple, clinical, by = "patid")) %>%
		arm(
			title = "multiple_log",
			plan = ptsd + sad_bin + pet_bin ~ mesor + amp1 + phi1 + (1 | pair),
			approach = logistic_reg() %>% set_engine("glmer"),
			exposure = c("(1 | pair)"),
			pattern = "direct",
			strata = "outcomes"
		) %>%
		arm(
			title = "multiple_lin",
			plan = beck_total + global_cfr ~ mesor + amp1 + phi1 + amp2 + phi2 + (1 | pair),
			approach = linear_reg() %>% set_engine("lmer"),
			exposure = c("(1 | pair)"),
			pattern = "direct",
			strata = "outcomes"
		) %>%
		equip(which_arms = c("multiple_log", "multiple_lin"))

	# Return octobeast
	octobeast

}

# Outcomes
make_twins_survival <- function(clinical, ecg, outcomes, cosinors) {

	# With cosinor outcomes
	df <-
		inner_join(ecg, clinical, by = "vetrid") %>%
		filter(hour %in% c(6:7)) %>%
		mutate(bpm = (1 / (rr * 1/1000 * 1 / 60))/10)
	cos_df <- inner_join(cosinors$single, clinical, by = "patid")

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
			title = "bpm_death",
			plan = Surv(stop, status) ~ bpm + pet_bin + age + bmi + race + smoking + prevchd + chf + hptn + dm + sad_bin + ptsd,
			approach = cox_reg() %>% set_engine("survival", method = "breslow"),
			pattern = "sequential",
			strata = "hour"
		) %>%
		equip(which_arms = c("hf_death", "lf_death", "vlf_death", "dyx_death", "ac_death", "bpm_death")) %>%
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
			title = "bpm_cvd",
			plan = Surv(stop, status) ~ bpm + pet_bin + age + bmi + race + smoking + prevchd + chf + hptn + dm + sad_bin + ptsd,
			approach = cox_reg() %>% set_engine("survival", method = "breslow"),
			pattern = "sequential",
			strata = "hour"
		) %>%
		equip(which_arms = c("hf_cvd", "lf_cvd", "vlf_cvd", "dyx_cvd", "ac_cvd", "bpm_cvd")) %>%
		# Cosinor analysis
		change_core(inner_join(cos_df, outcomes$death, by = c("vetrid" = "id"))) %>%
		arm(
			title = "cosinor_death",
			plan = Surv(stop, status) ~ mesor + amp1 + phi1 + pet_bin + age + bmi + race + smoking + prevchd + chf + hptn + dm + sad_bin + ptsd,
			approach = cox_reg() %>% set_engine("survival", method = "breslow"),
			pattern = "sequential",
			strata = "outcomes"
		) %>%
		equip(which_arms = "cosinor_death") %>%
		change_core(inner_join(cos_df, outcomes$cvd, by = c("vetrid" = "id"))) %>%
		arm(
			title = "cosinor_cvd",
			plan = Surv(stop, status) ~ mesor + amp1 + phi1 + pet_bin + age + bmi + race + smoking + prevchd + chf + hptn + dm + sad_bin + ptsd,
			approach = cox_reg() %>% set_engine("survival", method = "breslow"),
			pattern = "sequential",
			strata = "outcomes"
		) %>%
		equip(which_arms = "cosinor_cvd")

	# Return
	octobeast
}

# Report out twin models
report_twins_models <- function(models, survival, circadian) {

	teal <- "#487F84"

	# Set Up ------------------------------------------------------
	unadj <- 2 # HRV
	adj_demo <- 5 # + age + BMI
	adj_cv <- 10 # + smoking, CAD, CHF, HTN, DM

	hourly <-
		models$equipment %>%
		bind_rows(.id = "arm") %>%
		select(outcomes, test_num, level, tidied) %>%
		unnest(tidied) %>%
		filter(level %in% 6:7) %>%
		filter(term %in% c("lf", "hf", "vlf", "ac", "dyx", "bpm")) %>%
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
		cols_merge(columns = starts_with("bpm"), pattern = "{1} ({2}, {3})") %>%
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
			locations = cells_body(columns = "bpm_estimate", rows = abs(bpm_statistic) > 2.0)
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
			dyx_estimate = md("*Dyx*"),
			bpm_estimate = "BPM"
		) %>%
		cols_hide("bpm_estimate") %>%
		tab_footnote(
			footnote = "Model 1 = HRV",
			locations = cells_stub(rows = c(1:2))
		) %>%
		tab_footnote(
			footnote = "Model 2 = Model 1 + Age + BMI + Race",
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
		cols_merge(columns = starts_with("bpm"), pattern = "{1} ({2}, {3})") %>%
		fmt_number(columns = everything(), decimals = 2) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
			locations = cells_body(columns = "bpm_estimate", rows = abs(bpm_statistic) > 2.0)
		) %>%
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
			dyx_estimate = md("*Dyx*"),
			bpm_estimate = "BPM"
		) %>%
		cols_hide("bpm_estimate") %>%
		tab_footnote(
			footnote = "Model 1 = HRV",
			locations = cells_stub(rows = c(1:2))
		) %>%
		tab_footnote(
			footnote = "Model 2 = Model 1 + Age + BMI + Race",
			locations = cells_stub(rows = c(3:4))
		) %>%
		tab_footnote(
			footnote = "Model 3 = Model 2 + Smoking + HTN + Cardiovascular Disease",
			locations = cells_stub(rows = c(5:6))
		)

	# Survival -------------------------------------------------------------
	unadj <- 1
	adj_mpi <- 2
	adj_demo <- 5
	adj_cvd <- 10
	adj_sad <- 12

	outcomes <- survival$equipment[1:12] %>%
		bind_rows(.id = "arm") %>%
		separate(arm, into = c("hrv", "outcome"), sep = "_") %>%
		select(outcome, test_num, level, tidied) %>%
		filter(level %in% c(6:7)) %>%
		unnest(tidied) %>%
		filter(term %in% c("ac", "dyx", "hf", "lf", "vlf", "bpm")) %>%
		filter(test_num %in% c(unadj, adj_mpi, adj_demo, adj_cvd, adj_sad)) %>%
		group_by(outcome, term, test_num) %>%
		arrange(-abs(statistic), .by_group = TRUE) %>%
		slice(1) %>%
		ungroup() %>%
		select(-c(level, std.error, statistic)) %>%
		pivot_wider(
			names_from = term,
			values_from = c(estimate, conf.low, conf.high, p.value),
			names_glue = "{term}_{.value}"
		) %>%
		mutate(outcome = case_when(
			outcome == "cvd" ~ "Cardiovascular Death",
			outcome == "death" ~ "All Cause Mortality"
		)) %>%
		mutate(test_num = paste("Model", c(1:5, 1:5))) %>%
		gt(rowname_col = "test_num", groupname_col = "outcome") %>%
		cols_merge(columns = starts_with("ac"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("dyx"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("hf"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("lf"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("vlf"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("bpm"), pattern = "{1} ({2}, {3})") %>%
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
			locations = cells_body(columns = starts_with("bpm"),
														 rows = bpm_p.value < .05)
		) %>%
		cols_label(
			ac_estimate = "Acceleration Capacity",
			dyx_estimate = "Dyx",
			hf_estimate = "High Frequency HRV",
			lf_estimate = "Low Frequency HRV",
			vlf_estimate = "Very Low Frequency HRV",
			bpm_estimate = "Heart Rate"
		) %>%
		cols_hide("bpm_estimate") %>%
		tab_footnote(
			footnote = "Model 1 = HRV",
			locations = cells_stub(rows = c(1,6))
		) %>%
		tab_footnote(
			footnote = "Model 2 = Model 1 + Myocardial Perfusion Imaging",
			locations = cells_stub(rows = c(2,7))
		) %>%
		tab_footnote(
			footnote = "Model 3 = Model 2 + Age + BMI + Race",
			locations = cells_stub(rows = c(3,8))
		) %>%
		tab_footnote(
			footnote = "Model 4 = Model 3 + Cardiovascular Disease + Hypertension + Diabetes + Smoking",
			locations = cells_stub(rows = c(4,9))
		) %>%
		tab_footnote(
			footnote = "Model 5 = Model 4 + Depression + PTSD",
			locations = cells_stub(rows = c(5,10))
		)


	# Set up cosinor setup -----------------------------------------------------

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
			level == "rr" ~ "RR Interval"
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

	# Cosinor survival -----------------------------------------------
	unadj <- 3
	adj_mpi <- 4
	adj_demo <- 7
	adj_cvd <- 12
	adj_sad <- 14

	outcomes_cosinor <-
		survival$equipment[13:14] %>%
		bind_rows(.id = "arm") %>%
		separate(arm, into = c("trash", "outcome")) %>%
		select(outcome, test_num, level, tidied) %>%
		filter(test_num %in% c(adj_sad)) %>%
		unnest(tidied) %>%
		filter(term %in% c("mesor", "amp1", "phi1")) %>%
		select(-std.error, -statistic, -test_num) %>%
		pivot_wider(names_from = term, values_from = c(estimate, conf.low, conf.high, p.value), names_glue = "{term}_{.value}") %>%
		mutate(outcome = case_when(
			outcome == "cvd" ~ "Cardiovascular Death",
			outcome == "death" ~ "All Cause Mortality"
		)) %>%
		mutate(level = case_when(
			level == "hf" ~ "High Frequency HRV",
			level == "lf" ~ "Low Frequency HRV",
			level == "vlf" ~ "Very Low Frequency HRV",
			level == "dyx" ~ "Dyx",
			level == "ac" ~ "Acceleration Capacity",
			level == "rr" ~ "RR"
		)) %>%
		gt(rowname_col = "level", groupname_col = "outcome") %>%
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
														 rows = mesor_p.value < 0.05)
		) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
			locations = cells_body(columns = starts_with("amp1"),
														 rows = amp1_p.value < 0.05)
		) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
			locations = cells_body(columns = starts_with("phi1"),
														 rows = phi1_p.value < 0.05)
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
		ischemia_cosinor = ischemia_cosinor,
		outcomes_cosinor = outcomes_cosinor
	)

	# Returning
	tables

}
