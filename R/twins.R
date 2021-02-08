# Summary tables
make_twins_tables <- function(clinical, ecg) {

	# Table 1
	one <-
		clinical %>%
		select(study, age, bmi, race, smoking, prevchd, chf, hptn, dm, ptsd, sad_bin, pet_bin) %>%
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
		)) %>%
		as_gt() %>%
		tab_header(title = "Twins: Cohort Descriptions") %>%
		tab_source_note("ETSF = Emory Twins Study Follow-up, SAVEIT = Stress and Vascular Evaluation in Twins, THS = Twins Heart Study") %>%
		tab_options(table.font.size = "11px") %>%
		as_raw_html()

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
		)) %>%
		as_gt() %>%
		tab_header(title = "HRV In Twin Cohorts") %>%
		tab_source_note("ETSF = Emory Twins Study Follow-up, SAVEIT = Stress and Vascular Evaluation in Twins, THS = Twins Heart Study") %>%
		tab_options(table.font.size = "11px") %>%
		as_raw_html()

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

# Outcomes
make_twins_survival <- function(clinical, ecg, outcomes) {

	octobeast <-
		octomod() %>%
		core(inner_join(ecg, outcomes$death, by = c("vetrid" = "id"))) %>%
		arm(
			title = "death_outcomes",
			plan = Surv(stop, status) ~ hf + lf + vlf + dyx + ac + rr,
			approach = cox_reg() %>% set_engine("survival", method = "breslow"),
			pattern = "parallel",
			strata = "hour"
		) %>%
		equip(which_arms = "death_outcomes") %>%
		change_core(inner_join(ecg, outcomes$cvd, by = c("vetrid" = "id"))) %>%
		arm(
			title = "cvd_outcomes",
			plan = Surv(stop, status) ~ hf + lf + vlf + dyx + ac + rr,
			approach = cox_reg() %>% set_engine("survival", method = "breslow"),
			pattern = "parallel",
			strata = "hour"
		) %>%
		equip(which_arms = "cvd_outcomes")

	# Return
	octobeast
}

# Report out twin models
report_twins_models <- function(models, survival) {

	# Set Up All Models ------------------------------------------------------
	unadj <- 2 # HRV
	adj_demo <- 4 # + age + BMI
	adj_cv <- 9 # + smoking, CAD, CHF, HTN, DM

	df <-
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
		df %>%
		filter(outcomes %in% c("PTSD", "Depression")) %>%
		gt(rowname_col = "test_num", groupname_col = "outcomes") %>%
		cols_merge(columns = starts_with("hf_"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("lf_"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("vlf_"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("dyx_"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("ac"), pattern = "{1} ({2}, {3})") %>%
		fmt_number(columns = everything(), decimals = 2) %>%
		tab_style(
			style = list(cell_fill(color = "#487f84"), cell_text(color = "white")),
			locations = cells_body(columns = "lf_estimate", rows = abs(lf_statistic) > 2.0)
		) %>%
		tab_style(
			style = list(cell_fill(color = "#487f84"), cell_text(color = "white")),
			locations = cells_body(columns = "hf_estimate", rows = abs(hf_statistic) > 2.0)
		) %>%
		tab_style(
			style = list(cell_fill(color = "#487f84"), cell_text(color = "white")),
			locations = cells_body(columns = "vlf_estimate", rows = abs(vlf_statistic) > 2.0)
		) %>%
		tab_style(
			style = list(cell_fill(color = "#487f84"), cell_text(color = "white")),
			locations = cells_body(columns = "ac_estimate", rows = abs(ac_statistic) > 2.0)
		) %>%
		tab_style(
			style = list(cell_fill(color = "#487f84"), cell_text(color = "white")),
			locations = cells_body(columns = "dyx_estimate", rows = abs(dyx_statistic) > 2.0)
		) %>%
		tab_header(
			title = "Early Morning HRV and Chronic Psychological Stress",
			subtitle = "Emory Twins Study"
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
		) %>%
		tab_source_note("Depression is measured as a binary outcome with Beck Depression Inventory score > 14. Green highlight signifies p-value < 0.05. PTSD = Post-Traumatic Stress Disorder, HRV = Heart Rate Variability, LF = Low Frequency HRV, HF = High Frequency HRV, VLF = Very Low Frequency HRV, AC = Acceleration Capacity" )

	# Ischemia Models  ----------------------------------------------------
	mpi <-
		df %>%
		filter(outcomes %in% c("Coronary Flow Reserve", "Abnormal MPI")) %>%
		gt(rowname_col = "test_num", groupname_col = "outcomes") %>%
		cols_merge(columns = starts_with("hf_"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("lf_"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("vlf_"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("dyx_"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("ac"), pattern = "{1} ({2}, {3})") %>%
		fmt_number(columns = everything(), decimals = 2) %>%
		tab_style(
			style = list(cell_fill(color = "#487f84"), cell_text(color = "white")),
			locations = cells_body(columns = "lf_estimate", rows = abs(lf_statistic) > 2.0)
		) %>%
		tab_style(
			style = list(cell_fill(color = "#487f84"), cell_text(color = "white")),
			locations = cells_body(columns = "hf_estimate", rows = abs(hf_statistic) > 2.0)
		) %>%
		tab_style(
			style = list(cell_fill(color = "#487f84"), cell_text(color = "white")),
			locations = cells_body(columns = "vlf_estimate", rows = abs(vlf_statistic) > 2.0)
		) %>%
		tab_style(
			style = list(cell_fill(color = "#487f84"), cell_text(color = "white")),
			locations = cells_body(columns = "ac_estimate", rows = abs(ac_statistic) > 2.0)
		) %>%
		tab_style(
			style = list(cell_fill(color = "#487f84"), cell_text(color = "white")),
			locations = cells_body(columns = "dyx_estimate", rows = abs(dyx_statistic) > 2.0)
		) %>%
		tab_header(title = "Early Morning HRV and Myocardial Perfusion") %>%
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
		) %>%
		tab_source_note("Green highlight signifies p-value < 0.05. HRV = Heart Rate Variability, LF = Low Frequency HRV, HF = High Frequency HRV, VLF = Very Low Frequency HRV, AC = Acceleration Capacity" )

	# Survival -------------------------------------------------------------
	outcomes <-
		survival$equipment %>%
		bind_rows(.id = "arm") %>%
		select(arm, level, test_num, tidied) %>%
		unnest(tidied) %>%
		filter(level %in% c(6:9)) %>%
		group_by(arm, term) %>%
		arrange(p.value, .by_group = TRUE) %>%
		slice(1) %>%
		ungroup() %>%
		select(arm, term, estimate, conf.low, conf.high, p.value) %>%
		pivot_wider(
			names_from = arm,
			values_from = c(estimate, conf.low, conf.high, p.value),
			names_glue = "{arm}_{.value}"
		) %>%
		gt(rowname_col = "term") %>%
		cols_merge(columns = starts_with("death"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("cvd"), pattern = "{1} ({2}, {3})") %>%
		fmt_number(
			columns = contains("outcomes"),
			decimals = 2,
			drop_trailing_zeros = TRUE
		) %>%
		tab_style(
			style = list(cell_fill(color = "#487f84"), cell_text(color = "white")),
			locations = cells_body(columns = starts_with("death"),
														 rows = death_outcomes_p.value < .05)
		) %>%
		tab_style(
			style = list(cell_fill(color = "#487f84"), cell_text(color = "white")),
			locations = cells_body(columns = starts_with("cvd"),
														 rows = cvd_outcomes_p.value < .05)
		) %>%
		tab_header(
			title = "Outcomes by HRV",
			subtitle = "Emory Twins Study"
		)

	# Return -------------------------------------------------------------
	tables <- list(
		psych = psych,
		mpi = mpi,
		outcomes = outcomes
	)

	# Returning
	tables

}
