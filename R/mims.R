# Basic tables
make_mims_tables <- function(clinical, raw) {

	# Making table 1
	one <-
		clinical %>%
		select(c(studygroup, age_bl, female_bl, race_bl, bmi, smk_now1, cath_cad70binary, starts_with("hx_"), rdr_msi_bl, rdr_psi2_bl, scid_depression_bl, scid_ptsd_bl)) %>%
		select(-c(hx_revasc_bl)) %>%
		filter(!is.na(studygroup)) %>%
		tbl_strata(
			strata = studygroup,
			.tbl_fun =
				~.x %>%
				tbl_summary(
					by = rdr_msi_bl,
					type = list(c(rdr_psi2_bl, scid_ptsd_bl, scid_depression_bl, female_bl) ~ "dichotomous", c(race_bl) ~ "categorical"),
					value = list(female_bl ~ "Female"),
					missing = "no",
					label = list(rdr_psi2_bl ~ "PSIMI", female_bl ~ "Sex (Female)"),
				),
			.combine_with = "tbl_merge"
		) %>%
		modify_header(update = list(
			label ~ 'Characteristic',
			stat_1_1 ~ '**MSIMI = 0**, N = 256',
			stat_2_1 ~ '**MSIMI = 1**, N = 50',
			stat_1_2 ~ '**MSIMI = 0**, N = 440',
			stat_2_2 ~ '**MSIMI = 1**, N = 188'
		)) %>%
		modify_spanning_header(ends_with("_1") ~ "MIMS") %>%
		modify_spanning_header(ends_with("_2") ~ "MIPS")

	# Distribution of HRV
	paired <-
		octomod() %>%
		core(raw) %>%
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
			hrv == "hr" ~ "Heart Rate (beats/minute)",
			hrv == "hf" ~ "High Frequency HRV (ln ms^2)",
			hrv == "lf" ~ "Low Frequency HRV (ln ms^2)",
			hrv == "twa" ~ "T Wave Area (µVs)"
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
		fmt_number(columns = everything(), decimals = 2) %>%
		cols_label(estimate = "Mean (95% CI)", statistic = "T-statistic")

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
			stat_1 ~ '**MSIMI = 0**, N = 710',
			stat_2 ~ '**MSIMI = 1**, N = 243',
			p.value ~ '**p-value**'
		)) %>%
		as_gt() %>%
		tab_row_group(group = "Low Frequency HRV", rows = 1:3) %>%
		tab_row_group(group = "High Frequency HRV", rows = 4:6) %>%
		tab_row_group(group = "T Wave Area", rows = 7:9) %>%
		tab_row_group(group = "Heart Rate", rows = 10:12)

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

	violin <-
		clinical %>%
		select(c(starts_with(c("lf_", "hf_", "twa_", "hr_")), rdr_msi_bl)) %>%
		select(-contains(c("rest_stress", "stress_rec"))) %>%
		pivot_longer(cols = starts_with(c("lf_", "hf_", "twa_", "hr_")), names_to = "measures") %>%
		separate(col = "measures", into = c("hrv", "phase"), sep = "_") %>%
		filter(value < 150) %>%
		mutate(
			phase = factor(phase, levels = c("rest", "stress", "recovery"), labels = c("Rest", "Stress", "Recovery")),
			hrv = factor(hrv, levels = c("hf", "lf", "twa", "hr"), labels = c("High Frequency", "Low Frequency", "T Wave Area", "Heart Rate"))
		) %>%
		ggplot(aes(x = phase, y = value, fill = phase)) +
		facet_wrap(~hrv, scales = "free") +
		geom_violin() +
		geom_jitter(size = 0.2, width = 0.1, alpha = 0.2, color = "white")

	# List
	figures <- list(
		violin = violin
	)

	# Return
	figures

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
			title = "lf_death",
			plan = Surv(stop, status) ~ lf_stress + rdr_msi_bl + age_bl + bmi + female_bl + race_bl + smk_now1 + hx_diabetes_bl + hx_hypertension_bl + hx_hbchol_bl + hx_revasc_bl + hx_ptca_bl + hx_cabg_bl + scid_depression_bl + scid_ptsd_bl,
			approach = cox_reg() %>% set_engine("survival", method = "breslow"),
			pattern = "sequential",
			exposure = c("lf_stress")
		) %>%
		equip(which_arms = "lf_death") %>%
		change_core(inner_join(clinical, outcomes$death_cv, by = c("patid" = "id"))) %>%
		arm(
			title = "lf_cv_death",
			plan = Surv(stop, status) ~ lf_stress + rdr_msi_bl + age_bl + bmi + female_bl + race_bl + smk_now1 + hx_diabetes_bl + hx_hypertension_bl + hx_hbchol_bl + hx_revasc_bl + hx_ptca_bl + hx_cabg_bl + scid_depression_bl + scid_ptsd_bl,
			approach = cox_reg() %>% set_engine("survival", method = "breslow"),
			pattern = "sequential",
			exposure = c("lf_stress")
		) %>%
		equip(which_arms = "lf_cv_death") %>%
		change_core(inner_join(clinical, outcomes$mace_marginal, by = c("patid" = "id"))) %>%
		arm(
			title = "lf_marg",
			plan = Surv(stop, status) ~ lf_stress + rdr_msi_bl + age_bl + bmi + female_bl + race_bl + smk_now1 + hx_diabetes_bl + hx_hypertension_bl + hx_hbchol_bl + hx_revasc_bl + hx_ptca_bl + hx_cabg_bl + scid_depression_bl + scid_ptsd_bl + cluster(patid),
			approach = cox_reg() %>% set_engine("survival", method = "breslow"),
			pattern = "sequential",
			exposure = c("lf_stress", "cluster(patid)")
		) %>%
		equip(which_arms = "lf_marg") %>%
		change_core(inner_join(clinical, outcomes$mace_pwptt, by = c("patid" = "id"))) %>%
		arm(
			title = "lf_pwptt",
			plan = Surv(stop, status) ~ lf_stress + rdr_msi_bl + age_bl + bmi + female_bl + race_bl + smk_now1 + hx_diabetes_bl + hx_hypertension_bl + hx_hbchol_bl + hx_revasc_bl + hx_ptca_bl + hx_cabg_bl + scid_depression_bl + scid_ptsd_bl + cluster(patid) + strata(strata),
			approach = cox_reg() %>% set_engine("survival", method = "breslow"),
			pattern = "sequential",
			exposure = c("lf_stress", "cluster(patid)", "strata(strata)")
		) %>%
		equip(which_arms = "lf_pwptt") %>%
		change_core(inner_join(clinical, outcomes$mace_pwpgt, by = c("patid" = "id"))) %>%
		arm(
			title = "lf_pwpgt",
			plan = Surv(stop, status) ~ lf_stress + rdr_msi_bl + age_bl + bmi + female_bl + race_bl + smk_now1 + hx_diabetes_bl + hx_hypertension_bl + hx_hbchol_bl + hx_revasc_bl + hx_ptca_bl + hx_cabg_bl + scid_depression_bl + scid_ptsd_bl + cluster(patid) + strata(strata),
			approach = cox_reg() %>% set_engine("survival", method = "breslow"),
			pattern = "sequential",
			exposure = c("lf_stress", "cluster(patid)", "strata(strata)")
		) %>%
		equip(which_arms = "lf_pwpgt") %>%
		change_core(inner_join(clinical, outcomes$mace_ag, by = c("patid" = "id"))) %>%
		arm(
			title = "lf_ag",
			plan = Surv(stop, status) ~ lf_stress + rdr_msi_bl + age_bl + bmi + female_bl + race_bl + smk_now1 + hx_diabetes_bl + hx_hypertension_bl + hx_hbchol_bl + hx_revasc_bl + hx_ptca_bl + hx_cabg_bl + scid_depression_bl + scid_ptsd_bl + cluster(patid),
			approach = cox_reg() %>% set_engine("survival", method = "breslow"),
			pattern = "sequential",
			exposure = c("lf_stress", "cluster(patid)")
		) %>%
		equip(which_arms = "lf_ag") %>%
		change_core(inner_join(clinical, outcomes$death, by = c("patid" = "id"))) %>%
		arm(
			title = "hf_death",
			plan = Surv(stop, status) ~ hf_stress + rdr_msi_bl + age_bl + bmi + female_bl + race_bl + smk_now1 + hx_diabetes_bl + hx_hypertension_bl + hx_hbchol_bl + hx_revasc_bl + hx_ptca_bl + hx_cabg_bl + scid_depression_bl + scid_ptsd_bl,
			approach = cox_reg() %>% set_engine("survival", method = "breslow"),
			pattern = "sequential",
			exposure = c("hf_stress")
		) %>%
		equip(which_arms = "hf_death") %>%
		change_core(inner_join(clinical, outcomes$death_cv, by = c("patid" = "id"))) %>%
		arm(
			title = "hf_cv_death",
			plan = Surv(stop, status) ~ hf_stress + rdr_msi_bl + age_bl + bmi + female_bl + race_bl + smk_now1 + hx_diabetes_bl + hx_hypertension_bl + hx_hbchol_bl + hx_revasc_bl + hx_ptca_bl + hx_cabg_bl + scid_depression_bl + scid_ptsd_bl,
			approach = cox_reg() %>% set_engine("survival", method = "breslow"),
			pattern = "sequential",
			exposure = c("hf_stress")
		) %>%
		equip(which_arms = "hf_cv_death") %>%
		change_core(inner_join(clinical, outcomes$mace_marginal, by = c("patid" = "id"))) %>%
		arm(
			title = "hf_marg",
			plan = Surv(stop, status) ~ hf_stress + rdr_msi_bl + age_bl + bmi + female_bl + race_bl + smk_now1 + hx_diabetes_bl + hx_hypertension_bl + hx_hbchol_bl + hx_revasc_bl + hx_ptca_bl + hx_cabg_bl + scid_depression_bl + scid_ptsd_bl + cluster(patid),
			approach = cox_reg() %>% set_engine("survival", method = "breslow"),
			pattern = "sequential",
			exposure = c("hf_stress", "cluster(patid)")
		) %>%
		equip(which_arms = "hf_marg") %>%
		change_core(inner_join(clinical, outcomes$mace_pwptt, by = c("patid" = "id"))) %>%
		arm(
			title = "hf_pwptt",
			plan = Surv(stop, status) ~ hf_stress + rdr_msi_bl + age_bl + bmi + female_bl + race_bl + smk_now1 + hx_diabetes_bl + hx_hypertension_bl + hx_hbchol_bl + hx_revasc_bl + hx_ptca_bl + hx_cabg_bl + scid_depression_bl + scid_ptsd_bl + cluster(patid) + strata(strata),
			approach = cox_reg() %>% set_engine("survival", method = "breslow"),
			pattern = "sequential",
			exposure = c("hf_stress", "cluster(patid)", "strata(strata)")
		) %>%
		equip(which_arms = "hf_pwptt") %>%
		change_core(inner_join(clinical, outcomes$mace_pwpgt, by = c("patid" = "id"))) %>%
		arm(
			title = "hf_pwpgt",
			plan = Surv(stop, status) ~ hf_stress + rdr_msi_bl + age_bl + bmi + female_bl + race_bl + smk_now1 + hx_diabetes_bl + hx_hypertension_bl + hx_hbchol_bl + hx_revasc_bl + hx_ptca_bl + hx_cabg_bl + scid_depression_bl + scid_ptsd_bl + cluster(patid) + strata(strata),
			approach = cox_reg() %>% set_engine("survival", method = "breslow"),
			pattern = "sequential",
			exposure = c("hf_stress", "cluster(patid)", "strata(strata)")
		) %>%
		equip(which_arms = "hf_pwpgt") %>%
		change_core(inner_join(clinical, outcomes$mace_ag, by = c("patid" = "id"))) %>%
		arm(
			title = "hf_ag",
			plan = Surv(stop, status) ~ hf_stress + rdr_msi_bl + age_bl + bmi + female_bl + race_bl + smk_now1 + hx_diabetes_bl + hx_hypertension_bl + hx_hbchol_bl + hx_revasc_bl + hx_ptca_bl + hx_cabg_bl + scid_depression_bl + scid_ptsd_bl + cluster(patid),
			approach = cox_reg() %>% set_engine("survival", method = "breslow"),
			pattern = "sequential",
			exposure = c("hf_stress", "cluster(patid)")
		) %>%
		equip(which_arms = "hf_ag")

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
		cols_merge(
			columns = starts_with("scid_depression_bl"),
			pattern = "{1} ({2}, {3})"
		) %>%
		cols_merge(
			columns = starts_with("scid_ptsd_bl"),
			pattern = "{1} ({2}, {3})"
		) %>%
		fmt_number(
			columns = starts_with("scid"),
			drop_trailing_zeros = TRUE,
			decimals = 2
		) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
			locations = cells_body(
				columns = vars(scid_ptsd_bl_estimate),
				rows = scid_ptsd_bl_p.value < 0.05
			)
		) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
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
			scid_depression_bl_estimate = "Depression",
			scid_ptsd_bl_estimate = "PTSD"
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
			pattern = "{1} ({2}, {3})"
		) %>%
		cols_merge(
			columns = starts_with("rdr_psi2_bl"),
			pattern = "{1} ({2}, {3})"
		) %>%
		cols_merge(
			columns = starts_with("rdr_combined"),
			pattern = "{1} ({2}, {3})"
		) %>%
		fmt_number(
			columns = starts_with("rdr"),
			drop_trailing_zeros = TRUE,
			decimals = 2
		) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
			locations = cells_body(
				columns = vars(rdr_combined_estimate),
				rows = rdr_combined_p.value < 0.05
			)
		) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
			locations = cells_body(
				columns = vars(rdr_msi_bl_estimate),
				rows = rdr_msi_bl_p.value < 0.05
			)
		) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
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
			pattern = "{1} ({2}, {3})"
		) %>%
		cols_merge(
			columns = starts_with("lf_rest"),
			pattern = "{1} ({2}, {3})"
		) %>%
		cols_merge(
			columns = starts_with("hf_stress"),
			pattern = "{1} ({2}, {3})"
		) %>%
		cols_merge(
			columns = starts_with("hf_rest"),
			pattern = "{1} ({2}, {3})"
		) %>%
		fmt_number(
			columns = everything(),
			drop_trailing_zeros = TRUE,
			decimals = 2
		) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
			locations = cells_body(
				columns = "lf_stress_estimate",
				rows = lf_stress_p.value < .05
			)
		) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
			locations = cells_body(
				columns = "lf_rest_estimate",
				rows = lf_rest_p.value < .05
			)
		) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
			locations = cells_body(
				columns = "hf_stress_estimate",
				rows = hf_stress_p.value < .05
			)
		) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
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
		tab_stubhead(label = "Sequential Models") %>%
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

	unadj <- 1
	adj_ms <- 2
	adj_demo <- 6
	adj_risk <- 10
	adj_cv <- 13
	adj_psych <- 15

	outcomes <-
		survival$equipment %>%
		bind_rows(.id = "arm") %>%
		separate(arm, into = c("hrv", "model"), sep = "_", extra = "merge") %>%
		unnest(tidied) %>%
		select(c(hrv, model, test_num, term, estimate, conf.low, conf.high, p.value)) %>%
		filter(test_num %in% c(unadj, adj_ms, adj_demo, adj_risk, adj_cv, adj_psych)) %>%
		filter(str_detect(term, pattern = "stress")) %>%
		select(-term) %>%
		pivot_wider(
			names_from = model,
			values_from = c(estimate, conf.low, conf.high, p.value),
			names_glue = "{model}_{.value}"
		) %>%
		mutate(test_num = paste("Model", c(1:6, 1:6))) %>%
		mutate(hrv = case_when(
			hrv == "lf" ~ "Stress Low Frequency HRV (ln ms^2)",
			hrv == "hf" ~ "Stress High Frequency HRV (ln ms^2)",
		)) %>%
		gt(rowname_col = "test_num", groupname_col = "hrv") %>%
		cols_merge(columns = starts_with("death"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("cv_"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("marg"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("pwptt"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("pwpgt"), pattern = "{1} ({2}, {3})") %>%
		cols_merge(columns = starts_with("ag"), pattern = "{1} ({2}, {3})") %>%
		fmt_number(
			columns = contains("_"),
			decimals = 2,
			drop_trailing_zeros = TRUE
		) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
			locations = cells_body(columns = starts_with("death"),
														 rows = death_p.value < .01)
		) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
			locations = cells_body(columns = starts_with("cv_death"),
														 rows = cv_death_p.value < .01)
		) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
			locations = cells_body(columns = starts_with("marg"),
														 rows = marg_p.value < .01)
		) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
			locations = cells_body(columns = starts_with("pwptt"),
														 rows = pwptt_p.value < .01)
		) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
			locations = cells_body(columns = starts_with("pwpgt"),
														 rows = pwpgt_p.value < .01)
		) %>%
		tab_style(
			style = list(cell_text(weight = "bold")),
			locations = cells_body(columns = starts_with("ag"),
														 rows = ag_p.value < .01)
		) %>%
		cols_label(
			death_estimate = "Death",
			cv_death_estimate = "Cardiovascular Death",
			marg_estimate = "Marginal",
			pwptt_estimate = "PWP Total Time",
			pwpgt_estimate = "PWP Gap Time",
			ag_estimate = "Anderson Gill",
		) %>%
		tab_footnote(
			footnote = "Model 1 = MSIMI ~ HRV",
			locations = cells_stub(rows = c(1, 7))
		) %>%
		tab_footnote(
			footnote = "Model 2 = Model 1 + MSIMI",
			locations = cells_stub(rows = c(2, 8))
		) %>%
		tab_footnote(
			footnote = "Model 3 = Model 2 + Age + BMI + Sex + Race",
			locations = cells_stub(rows = c(3, 9))
		) %>%
		tab_footnote(
			footnote = "Model 4 = Model 3 + Smoking + Diabetes + Hypertension + Hyperlipidemia",
			locations = cells_stub(rows = c(4, 10))
		) %>%
		tab_footnote(
			footnote = "Model 5 = Model 4 + Known Coronary/Peripheral Artery Disease",
			locations = cells_stub(rows = c(5, 11))
		) %>%
		tab_footnote(
			footnote = "Model 6 = Model 5 + Depression + Post-Traumatic Stress Disorder",
			locations = cells_stub(rows = c(6, 12))
		)

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
