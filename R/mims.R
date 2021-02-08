# Basic tables
make_mims_tables <- function(clinical) {

	# Making table 1
	one <-
		clinical %>%
		select(c(studygroup, age_bl, female_bl, race_bl, bmi, smk_now1, cath_cad70binary, starts_with("hx_"), rdr_msi_bl, rdr_psi2_bl, scid_depression_bl, scid_ptsd_bl)) %>%
		select(-c(hx_revasc_bl)) %>%
		tbl_summary(
			by = rdr_msi_bl,
			type = list(c(rdr_psi2_bl, scid_ptsd_bl, scid_depression_bl, female_bl) ~ "dichotomous", c(race_bl) ~ "categorical"),
			value = list(female_bl ~ "Female"),
			missing = "no",
			label = list(rdr_psi2_bl ~ "PSIMI", female_bl ~ "Sex (Female)"),
		) %>%
		add_p() %>%
		add_strata(
			strata = studygroup,
			method = "merge",
		) %>%
		modify_header(update = list(
			label ~ 'Characteristic',
			stat_1_1 ~ '**MSIMI = 0**, <br> N = 256',
			stat_2_1 ~ '**MSIMI = 1**, <br> N = 50',
			stat_1_2 ~ '**MSIMI = 0**, <br> N = 454',
			stat_2_2 ~ '**MSIMI = 1**, <br> N = 193'
		)) %>%
		modify_spanning_header(ends_with("_1") ~ "MIMS") %>%
		modify_spanning_header(ends_with("_2") ~ "MIPS") %>%
		as_gt() %>%
		tab_header(
			title = "MIMS and MIPS Cohorts",
			subtitle = "Clinical Characteristics"
		) %>%
		tab_source_note("MSIMI = Mental Stress Induced Myocardial Ischemia; PSIMI = Physical Stress Induced Myocardial Ischemia, MIMS = Myocardial Infarction and Mental Stress, MIPS = Mental Stress Ischemia Mechanisms and Prognosis Study")

	# Distribution of HRV
	paired <-
		octomod() %>%
		core(clinical) %>%
		arm(
			title = "paired_hr",
			plan = hr_rest ~ hr_stress + hr_recovery,
			pattern = "parallel",
			approach = "t.test",
			paired = TRUE
		) %>%
		arm(
			title = "paired_hf",
			plan = hf_rest ~ hf_stress + hf_recovery,
			pattern = "parallel",
			approach = "t.test",
			paired = TRUE
		) %>%
		arm(
			title = "paired_lf",
			plan = lf_rest ~ lf_stress + lf_recovery,
			pattern = "parallel",
			approach = "t.test",
			paired = TRUE
		) %>%
		arm(
			title = "paired_twa",
			plan = twa_rest ~ twa_stress + twa_recovery,
			pattern = "parallel",
			approach = "t.test",
			paired = TRUE
		) %>%
		equip() %>%
		pluck("equipment") %>%
		bind_rows(.id = "arm") %>%
		select(outcomes, vars, tidied) %>%
		unnest(tidied) %>%
		unnest(vars) %>%
		separate(vars, into = c("hrv", "phase"), sep = "_") %>%
		select(c(hrv, phase, estimate, conf.low, conf.high, statistic, p.value)) %>%
		mutate(hrv = case_when(
			hrv == "hr" ~ "Heart Rate",
			hrv == "hf" ~ "High Frequency HRV",
			hrv == "lf" ~ "Low Frequency HRV",
			hrv == "twa" ~ "T Wave Area"
		)) %>%
		mutate(phase = case_when(
			phase == "stress" ~ "Stress",
			phase == "recovery" ~ "Recovery"
		)) %>%
		gt(rowname_col = "phase", groupname_col = "hrv") %>%
		cols_merge(
			columns = c("estimate", "conf.low", "conf.high"),
			pattern = "{1} ({2}, {3})"
		) %>%
		cols_hide("p.value") %>%
		fmt_number(columns = everything(), decimals = 1) %>%
		tab_style(
			style = list(cell_fill(color = "#487f84"), cell_text(color = "white")),
			locations = cells_body(rows = p.value < 0.05)
		) %>%
		cols_label(estimate = "Mean (95% CI)", statistic = "T-statistic") %>%
		tab_header(
			title = "Distribution of HRV by Stress Phase",
			subtitle = "Within Subject Testing in MIMS/MIPS Cohorts"
		) %>%
		tab_stubhead("ECG/HRV Metric") %>%
		tab_source_note("Each metric is tested against corresponding resting value. Highlighted cells signify p-value < 0.05. HRV = Heart Rate Variability")

	# MSIMI and HRV
	hrv <-
		clinical %>%
		select(starts_with(c("lf", "hf", "twa", "hr", "rdr_msi"))) %>%
		select(-contains(c("rest_", "stress_"))) %>%
		tbl_summary(
			by = rdr_msi_bl,
			type = list(starts_with(c("lf", "hf", "lh", "twa", "hr", "rdr_msi")) ~ "continuous"),
			label = list(contains("rest") ~ "Rest", contains("stress") ~ "Stress", contains("recovery") ~ "Recovery"),
			missing = "no"
		) %>%
		add_p() %>%
		modify_header(update = list(
			label ~ '**Characteristic**',
			stat_1 ~ '**MSIMI = 0**, <br> N = 710',
			stat_2 ~ '**MSIMI = 1**, <br> N = 243',
			p.value ~ '**p-value**'
		)) %>%
		as_gt() %>%
		tab_header(
			title = "HRV and Mental Stress",
			subtitle = "MIMS and MIPS Cohorts"
		) %>%
		tab_stubhead("ECG/HRV Metrics") %>%
		tab_row_group(group = "Low Frequency HRV", rows = 1:3) %>%
		tab_row_group(group = "High Frequency HRV", rows = 4:6) %>%
		tab_row_group(group = "T Wave Area", rows = 7:9) %>%
		tab_row_group(group = "Heart Rate", rows = 10:12) %>%
		tab_source_note("MSIMI = Mental Stress Induced Myocardial Ischemia") %>%
		cols_hide(columns = "p.value") %>%
		tab_style(
			style = list(cell_fill(color = "#487f84"), cell_text(color = "white")),
			locations = cells_body(rows = p.value < 0.05)
		)

	# Tables
	tables <- list(
		one = one,
		paired = paired,
		hrv = hrv
	)

	# Return
	tables

}

# Figures
make_mims_figures <- function(clinical, outcomes) {

}

# Make cross-sectional models
make_mims_models <- function(clinical) {

	octobeast <-
		octomod() %>%
		core(clinical) %>%
		arm(
			title = "psych",
			plan =  scid_depression_bl + scid_ptsd_bl ~ lf_rest + lf_stress + lf_recovery + hf_rest + hf_stress + hf_recovery + twa_rest + twa_stress + twa_recovery + hr_rest + hr_stress + hr_recovery,
			approach = logistic_reg() %>% set_engine("glm"),
			pattern = "parallel"
		) %>%
		arm(
			title = "ischemia",
			plan =  rdr_combined + rdr_msi_bl + rdr_psi2_bl ~ lf_rest + lf_stress + lf_recovery + hf_rest + hf_stress + hf_recovery + twa_rest + twa_stress + twa_recovery + hr_rest + hr_stress + hr_recovery,
			approach = logistic_reg() %>% set_engine("glm"),
			pattern = "parallel"
		) %>%
		arm(
			title = "adjusted_lf_stress",
			plan = rdr_msi_bl ~ lf_stress + age_bl + bmi + female_bl + race_bl + smk_now1 + hx_diabetes_bl + hx_hypertension_bl + hx_hbchol_bl + hx_revasc_bl + hx_ptca_bl + hx_cabg_bl + scid_depression_bl + scid_ptsd_bl,
			approach = logistic_reg() %>% set_engine("glm"),
			pattern = "sequential"
		) %>%
		arm(
			title = "adjusted_lf_rest",
			plan = rdr_msi_bl ~ lf_rest + age_bl + bmi + female_bl + race_bl + smk_now1 + hx_diabetes_bl + hx_hypertension_bl + hx_hbchol_bl + hx_revasc_bl + hx_ptca_bl + hx_cabg_bl + scid_depression_bl + scid_ptsd_bl,
			approach = logistic_reg() %>% set_engine("glm"),
			pattern = "sequential"
		) %>%
		arm(
			title = "adjusted_hf_stress",
			plan = rdr_msi_bl ~ hf_stress + age_bl + bmi + female_bl + race_bl + smk_now1 + hx_diabetes_bl + hx_hypertension_bl + hx_hbchol_bl + hx_revasc_bl + hx_ptca_bl + hx_cabg_bl + scid_depression_bl + scid_ptsd_bl,
			approach = logistic_reg() %>% set_engine("glm"),
			pattern = "sequential"
		) %>%
		arm(
			title = "adjusted_hf_rest",
			plan = rdr_msi_bl ~ hf_rest + age_bl + bmi + female_bl + race_bl + smk_now1 + hx_diabetes_bl + hx_hypertension_bl + hx_hbchol_bl + hx_revasc_bl + hx_ptca_bl + hx_cabg_bl + scid_depression_bl + scid_ptsd_bl,
			approach = logistic_reg() %>% set_engine("glm"),
			pattern = "sequential"
		) %>%
		equip() %>%
		sharpen()

}

# Outcomes
make_mims_survival <- function(clinical, outcomes) {

	octobeast <-
		octomod() %>%
		core(inner_join(clinical, outcomes$death, by = c("patid" = "id"))) %>%
		arm(
			title = "death_outcomes",
			plan = Surv(stop, status) ~ rdr_msi_bl + lf_stress + lf_rest + hf_stress + hf_rest,
			approach = cox_reg() %>% set_engine("survival", method = "breslow"),
			pattern = "parallel",
			exposure = c("rdr_msi_bl")
		) %>%
		equip(which_arms = "death_outcomes") %>%
		change_core(inner_join(clinical, outcomes$death_cv, by = c("patid" = "id"))) %>%
		arm(
			title = "cv_death_outcomes",
			plan = Surv(stop, status) ~ rdr_msi_bl + lf_stress + lf_rest + hf_stress + hf_rest,
			approach = cox_reg() %>% set_engine("survival", method = "breslow"),
			pattern = "parallel",
			exposure = c("rdr_msi_bl")
		) %>%
		equip(which_arms = "cv_death_outcomes") %>%
		change_core(inner_join(clinical, outcomes$mace_marginal, by = c("patid" = "id"))) %>%
		arm(
			title = "marginal_outcomes",
			plan = Surv(start, stop, status) ~ rdr_msi_bl + lf_stress + lf_rest + hf_stress + hf_rest + cluster(patid),
			approach = cox_reg() %>% set_engine("survival", method = "breslow"),
			pattern = "parallel",
			exposure = c("rdr_msi_bl", "cluster(patid)")
		) %>%
		equip(which_arms = "marginal_outcomes") %>%
		change_core(inner_join(clinical, outcomes$mace_pwptt, by = c("patid" = "id"))) %>%
		arm(
			title = "pwptt_outcomes",
			plan = Surv(start, stop, status) ~ rdr_msi_bl + lf_stress + lf_rest + hf_stress + hf_rest + cluster(patid) + strata(strata),
			approach = cox_reg() %>% set_engine("survival", method = "breslow"),
			pattern = "parallel",
			exposure = c("rdr_msi_bl", "cluster(patid)", "strata(strata)")
		) %>%
		equip(which_arms = "pwptt_outcomes") %>%
		change_core(inner_join(clinical, outcomes$mace_pwpgt, by = c("patid" = "id"))) %>%
		arm(
			title = "pwpgt_outcomes",
			plan = Surv(start, stop, status) ~ rdr_msi_bl + lf_stress + lf_rest + hf_stress + hf_rest + cluster(patid) + strata(strata),
			approach = cox_reg() %>% set_engine("survival", method = "breslow"),
			pattern = "parallel",
			exposure = c("rdr_msi_bl", "cluster(patid)", "strata(strata)")
		) %>%
		equip(which_arms = "pwpgt_outcomes") %>%
		change_core(inner_join(clinical, outcomes$mace_ag, by = c("patid" = "id"))) %>%
		arm(
			title = "ag_outcomes",
			plan = Surv(start, stop, status) ~ rdr_msi_bl + lf_stress + lf_rest + hf_stress + hf_rest,
			approach = cox_reg() %>% set_engine("survival", method = "breslow", robust = TRUE),
			pattern = "parallel",
			exposure = c("rdr_msi_bl", "cluster(patid)", "strata(strata)")
		) %>%
		equip(which_arms = "ag_outcomes")

	# Return
	octobeast

}

# Report
report_mims_models <- function(models, survival) {

	# Psych ---------------------------------------------------------------

	psych <-
		models$equipment$psych %>%
		bind_rows(.id = "arm") %>%
		unnest(tidied) %>%
		filter(term != "(Intercept)") %>%
		select(outcomes, term, estimate, conf.low, conf.high, p.value, metric) %>%
		pivot_wider(names_from = outcomes, values_from = c(estimate, conf.low, conf.high, p.value, metric), names_glue = "{outcomes}_{.value}") %>%
		separate(term, into = c("hrv", "phase"), sep = "_") %>%
		mutate(phase = case_when(
			phase == "rest" ~ "Rest",
			phase == "stress" ~ "Stress",
			phase == "recovery" ~ "Recovery"
		)) %>%
		gt(rowname_col = "phase") %>%
		tab_stubhead("ECG/HRV Metric") %>%
		tab_header(
			title = "HRV and Chronic Psychological Stress",
			subtitle = "MIMS/MIPS Cohorts"
		) %>%
		cols_merge(
			columns = starts_with("scid_depression_bl"),
			pattern = "{1} ({2}, {3}) <br> AUC {5}"
		) %>%
		cols_merge(
			columns = starts_with("scid_ptsd_bl"),
			pattern = "{1} ({2}, {3}) <br> AUC {5}"
		) %>%
		fmt_number(
			columns = starts_with("scid"),
			drop_trailing_zeros = TRUE,
			decimals = 2
		) %>%
		tab_style(
			style = list(cell_fill(color = "#487f84"), cell_text(color = "white")),
			locations = cells_body(
				columns = vars(scid_ptsd_bl_estimate),
				rows = scid_ptsd_bl_p.value < 0.05
			)
		) %>%
		tab_style(
			style = list(cell_fill(color = "#487f84"), cell_text(color = "white")),
			locations = cells_body(
				columns = vars(scid_depression_bl_estimate),
				rows = scid_depression_bl_p.value < 0.05
			)
		) %>%
		tab_row_group(group = "Low Frequency HRV", rows = 1:3) %>%
		tab_row_group(group = "High Frequency HRV", rows = 4:6) %>%
		tab_row_group(group = "T Wave Area", rows = 7:9) %>%
		tab_row_group(group = "Heart Rate", rows = 10:12) %>%
		cols_label(
			scid_depression_bl_estimate = "SCID Depression",
			scid_ptsd_bl_estimate = "SCID PTSD"
		) %>%
		tab_source_note(
			source_note = md("Highlighted boxes represent significant findings. SCID = Structured Clinical Interview for the DSM-IV, PTSD = Post-Traumatic Stress Disorder")
		) %>%
		tab_footnote(
			footnote = "Logistic regression model, OR with 95% CI and concordance statistic.",
			locations = cells_column_labels(columns = starts_with("scid"))
		) %>%
		cols_hide("hrv")

	# Unadjusted Ischemia Models ---------------------------------------------

	stress <-
		models$equipment %>%
		bind_rows(.id = "arm") %>%
		filter(arm == "ischemia") %>%
		unnest(tidied) %>%
		filter(term != "(Intercept)") %>%
		select(outcomes, term, estimate, conf.low, conf.high, p.value, metric) %>%
		pivot_wider(names_from = outcomes, values_from = c(estimate, conf.low, conf.high, p.value, metric), names_glue = "{outcomes}_{.value}") %>%
		separate(term, into = c("hrv", "phase"), sep = "_") %>%
		mutate(phase = case_when(
			phase == "rest" ~ "Rest",
			phase == "stress" ~ "Stress",
			phase == "recovery" ~ "Recovery"
		)) %>%
		gt(rowname_col = "phase") %>%
		tab_stubhead("ECG/HRV Metric") %>%
		tab_header(
			title = "HRV and Myocardial Ischemia",
			subtitle = "MIMS and MIPS Cohorts"
		) %>%
		cols_merge(
			columns = starts_with("rdr_msi_bl"),
			pattern = "{1} ({2}, {3}) <br> AUC {5}"
		) %>%
		cols_merge(
			columns = starts_with("rdr_psi2_bl"),
			pattern = "{1} ({2}, {3}) <br> AUC {5}"
		) %>%
		cols_merge(
			columns = starts_with("rdr_combined"),
			pattern = "{1} ({2}, {3}) <br> AUC {5}"
		) %>%
		fmt_number(
			columns = starts_with("rdr"),
			drop_trailing_zeros = TRUE,
			decimals = 2
		) %>%
		tab_style(
			style = list(cell_fill(color = "#487f84"), cell_text(color = "white")),
			locations = cells_body(
				columns = vars(rdr_combined_estimate),
				rows = rdr_combined_p.value < 0.05
			)
		) %>%
		tab_style(
			style = list(cell_fill(color = "#487f84"), cell_text(color = "white")),
			locations = cells_body(
				columns = vars(rdr_msi_bl_estimate),
				rows = rdr_msi_bl_p.value < 0.05
			)
		) %>%
		tab_style(
			style = list(cell_fill(color = "#487f84"), cell_text(color = "white")),
			locations = cells_body(
				columns = vars(rdr_psi2_bl_estimate),
				rows = rdr_psi2_bl_p.value < 0.05
			)
		) %>%
		tab_row_group(group = "Low Frequency HRV", rows = 1:3) %>%
		tab_row_group(group = "High Frequency HRV", rows = 4:6) %>%
		tab_row_group(group = "T Wave Area", rows = 7:9) %>%
		tab_row_group(group = "Heart Rate", rows = 10:12) %>%
		cols_label(
			rdr_combined_estimate = "Combined MSIMI/PSIMI",
			rdr_msi_bl_estimate = "MSIMI",
			rdr_psi2_bl_estimate = "PSIMI"
		) %>%
		tab_source_note(
			source_note = md("Highlighted boxes represent significant findings. MSIMI = Mental Stress-Induced Myocardial Ischemia, PSIMI = Physical Stress-Induced Myocardial Ischemia, HRV = Heart Rate Variability")
		) %>%
		tab_footnote(
			footnote = "Logistic regression model, OR with 95% CI and concordance statistic.",
			locations = cells_column_labels(columns = starts_with("rdr"))
		) %>%
		cols_hide("hrv")


	# Adjusted Ischemia Models ---------------------------------------------

	unadj <- 1
	adj_demo <- 5
	adj_risk <- 9
	adj_cad <- 12
	adj_psych <- 14

	# Get and Shape the models
	mpi <-
		models$equipment %>%
		bind_rows(.id = "arm") %>%
		filter(str_detect(arm, "adjusted")) %>%
		unnest(tidied) %>%
		filter(term %in% c("lf_stress", "lf_rest", "hf_stress", "hf_rest")) %>%
		select(test_num, term, estimate, conf.low, conf.high, p.value, metric) %>%
		filter(test_num %in% c(unadj, adj_demo, adj_risk, adj_cad, adj_psych)) %>%
		pivot_wider(names_from = term, values_from = c(estimate, conf.low, conf.high, p.value, metric), names_glue = "{term}_{.value}") %>%
		mutate(test_num = c("Model 1", "Model 2", "Model 3", "Model 4" , "Model 5")) %>%
		gt(rowname_col = "test_num") %>%
		cols_merge(
			columns = starts_with("lf_stress"),
			pattern = "{1} ({2}, {3}) <br> AUC {5}"
		) %>%
		cols_merge(
			columns = starts_with("lf_rest"),
			pattern = "{1} ({2}, {3}) <br> AUC {5}"
		) %>%
		cols_merge(
			columns = starts_with("hf_stress"),
			pattern = "{1} ({2}, {3}) <br> AUC {5}"
		) %>%
		cols_merge(
			columns = starts_with("hf_rest"),
			pattern = "{1} ({2}, {3}) <br> AUC {5}"
		) %>%
		fmt_number(
			columns = everything(),
			drop_trailing_zeros = TRUE,
			decimals = 2
		) %>%
		tab_style(
			style = list(cell_fill(color = "#487f84"), cell_text(color = "white")),
			locations = cells_body(
				columns = "lf_stress_estimate",
				rows = lf_stress_p.value < .01
			)
		) %>%
		tab_style(
			style = list(cell_fill(color = "#487f84"), cell_text(color = "white")),
			locations = cells_body(
				columns = "lf_rest_estimate",
				rows = lf_rest_p.value < .05
			)
		) %>%
		tab_style(
			style = list(cell_fill(color = "#487f84"), cell_text(color = "white")),
			locations = cells_body(
				columns = "hf_stress_estimate",
				rows = hf_stress_p.value < .05
			)
		) %>%
		tab_style(
			style = list(cell_fill(color = "#487f84"), cell_text(color = "white")),
			locations = cells_body(
				columns = "hf_rest_estimate",
				rows = hf_rest_p.value < .05
			)
		) %>%
		cols_label(
			lf_stress_estimate = "Stress LF",
			lf_rest_estimate = "Rest LF",
			hf_stress_estimate = "Stress HF",
			hf_rest_estimate = "Rest HF"
		) %>%
		tab_stubhead(label = "Adjusted Models") %>%
		tab_header(
			title = "Mental Stress-Induced Myocardial Ischemia and HRV",
			subtitle = "Sequential Models in MIMS/MIPS"
		) %>%
		tab_source_note("Highlighted boxes signify p.value < 0.05. MSIMI = Mental Stress-Induced Myocardial Ischemia, LF = Low Frequency HRV, HF = High Frequency HRV") %>%
		tab_footnote(
			footnote = "Model 1 = MSIMI ~ HRV",
			locations = cells_stub(rows = 1)
		) %>%
		tab_footnote(
			footnote = "Model 2 = Model 1 + Age + BMI + Sex + Race",
			locations = cells_stub(rows = 2)
		) %>%
		tab_footnote(
			footnote = "Model 3 = Model 2 + Smoking + Diabetes + Hypertension + Hyperlipidemia",
			locations = cells_stub(rows = 3)
		) %>%
		tab_footnote(
			footnote = "Model 4 = Model 3 + Known Coronary/Peripheral Artery Disease",
			locations = cells_stub(rows = 4)
		) %>%
		tab_footnote(
			footnote = "Model 5 = Model 4 + Depression + Post-Traumatic Stress Disorder",
			locations = cells_stub(rows = 5)
		)

	# Outcomes -----------------------------------------------------------

	outcomes <-
		survival$equipment %>%
		bind_rows(.id = "arm") %>%
		filter(str_detect(arm, "_outcomes")) %>%
		unnest(tidied) %>%
		select(c(arm, test_num, term, estimate, conf.low, conf.high, p.value)) %>%
		pivot_wider(
			names_from = arm,
			values_from = c(estimate, conf.low, conf.high, p.value),
			names_glue = "{arm}_{.value}"
		) %>%
		mutate(term = case_when(
			term == "rdr_msi_bl1" ~ "MSIMI",
			term == "lf_stress" ~ "Stress LF HRV",
			term == "lf_rest" ~ "Rest LF HRV",
			term == "hf_stress" ~ "Stress HF HRV",
			term == "hf_rest" ~ "Rest HF HRV"
		)) %>%
		gt() %>%
		cols_merge(columns = starts_with("death"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("cv_"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("marg"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("pwptt"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("pwpgt"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("ag"), pattern = "{1} ({2}, {3})") %>%
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
			locations = cells_body(columns = starts_with("cv_death"),
														 rows = cv_death_outcomes_p.value < .05)
		) %>%
		tab_style(
			style = list(cell_fill(color = "#487f84"), cell_text(color = "white")),
			locations = cells_body(columns = starts_with("marg"),
														 rows = marginal_outcomes_p.value < .05)
		) %>%
		tab_style(
			style = list(cell_fill(color = "#487f84"), cell_text(color = "white")),
			locations = cells_body(columns = starts_with("pwptt"),
														 rows = pwptt_outcomes_p.value < .05)
		) %>%
		tab_style(
			style = list(cell_fill(color = "#487f84"), cell_text(color = "white")),
			locations = cells_body(columns = starts_with("pwpgt"),
														 rows = pwpgt_outcomes_p.value < .05)
		) %>%
		tab_style(
			style = list(cell_fill(color = "#487f84"), cell_text(color = "white")),
			locations = cells_body(columns = starts_with("ag"),
														 rows = ag_outcomes_p.value < .05)
		) %>%
		cols_label(
			death_outcomes_estimate = "Death",
			cv_death_outcomes_estimate = "Cardiovascular Death",
			marginal_outcomes_estimate = "Marginal",
			pwptt_outcomes_estimate = "PWP Total Time",
			pwpgt_outcomes_estimate = "PWP Gap Time",
			ag_outcomes_estimate = "Anderson Gill",
			term = "Parameters"
		) %>%
		tab_row_group(rows = 1, group = "Model 1") %>%
		tab_row_group(rows = 2:3, group = "Model 2") %>%
		tab_row_group(rows = 4:5, group = "Model 3") %>%
		tab_row_group(rows = 6:7, group = "Model 4") %>%
		tab_row_group(rows = 8:9, group = "Model 5") %>%
		row_group_order(groups = c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5")) %>%
		tab_style(
			style = cell_borders(sides = c("top", "bottom")),
			locations = cells_row_groups(groups  = TRUE)
		) %>%
		tab_style(
			style = cell_borders(sides = "right"),
			locations = cells_body(columns = "term")
		) %>%
		cols_hide(columns = "test_num") %>%
		cols_width(vars(term) ~ px(150)) %>%
		cols_align(align = "right", columns = vars(term)) %>%
		tab_header(
			title = "Outcomes Analysis for Mental Stress and HRV",
			subtitle = "Traditional and Recurrent Event Models in MIMS/MIPS"
		) %>%
		tab_source_note("Estimates = HR (95% CI). Highlighted boxes signify p-value < 0.05. PWP = Prentice, Williams, and Peterson models, MSIMI = Mental Stress-Induced Myocardial Ischemia, LF = Low Frequency, HF = High Frequency, HRV = Heart Rate Variability")

	# Return -------------------------------------------------------

	# List
	tables <- list(
		psych = psych,
		stress = stress,
		mpi = mpi,
		outcomes = outcomes
	)

	# Return
	tables
}
