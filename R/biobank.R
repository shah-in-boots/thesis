# Biobank Tables
make_biobank_tables <- function(clinical, ecg, labels) {

	# Clean up ECG data
	labs <- list(
		n_nmean ~ "RR Interval",
		sdnn ~ "SDNN",
		rmssd ~ "RMSSD",
		pnn50 ~ "PNN50",
		ulf ~ "Ultra Low Frequency",
		vlf ~ "Very Low Frequency",
		lf ~ "Low Frequency",
		hf ~ "High Frequency",
		lfhf ~ "Low/High Frequency Ratio",
		ttlpwr ~ "Total Power",
		ac ~ "Acceleration Capacity",
		dc ~ "Deceleration Capacity",
		samp_en ~ "Sample Entropy",
		ap_en ~ "Approximate Entropy",
		dyx ~ "Dyx"
	)


	# Table 1 ---------------------------------------------------------------
	one <-
		clinical %>%
		mutate(
			race = labels$Race,
			gend = labels$Gender,
			adm_reason = labels$reason,
			setting = labels$Setting
		) %>%
		labelled::set_variable_labels(
			race = "Race",
			gend = "Sex",
			setting = "Setting",
			adm_reason = "Reason for Catheterization",
		) %>%
		select(c(age, race, blbmi, gend, setting, adm_reason, phq, sad, gensini, stenosis, cass70)) %>%
		tbl_summary(
			missing = "no",
			value = list(c(stenosis, sad) ~ "1")
		)

	# Depression -----------------------------------------------------------
	sad <-
		ecg$merged %>%
		group_by(patid) %>%
		summarise(across(c(n_nmean:ap_en, dyx), ~ mean(.x, na.rm = TRUE)), .groups = "keep") %>%
		left_join(clinical, ., by = "patid") %>%
		select(c(sad, n_nmean:dyx)) %>%
		tbl_summary(
			by = sad,
			missing = "no",
			value = list(c(sad) ~ "1"),
			label = labs
		) %>%
		add_p() %>%
		modify_header(update = list(
			label ~ '**HRV Metric**',
			stat_1 ~ '**No Depression**, N = 38',
			stat_2 ~ '**Depression**, N = 10',
			p.value ~ '**p-value**'
		))


	# Revascularization Table -----------------------------------------
	revasc <-
		ecg$merged %>%
		group_by(patid) %>%
		summarise(across(c(n_nmean:ap_en, dyx), ~ mean(.x, na.rm = TRUE)), .groups = "keep") %>%
		left_join(clinical, ., by = "patid") %>%
		select(c(pci_adm, n_nmean:dyx)) %>%
		tbl_summary(
			by = pci_adm,
			missing = "no",
			label = labs
		) %>%
		add_p() %>%
		modify_header(update = list(
			label ~ '**Characteristic**',
			stat_1 ~ '**No Revascularization**, N = 34',
			stat_2 ~ '**Revascularization**, N = 22',
			p.value ~ '**p-value**'
		))

	# HRV by Obstructive CAD --------------------------------------------------
	cad <-
		ecg$merged %>%
		group_by(patid) %>%
		summarise(across(c(n_nmean:ap_en, dyx), ~ mean(.x, na.rm = TRUE)), .groups = "keep") %>%
		left_join(clinical, ., by = "patid") %>%
		mutate(stenosis = replace_na(stenosis, 0)) %>%
		select(c(stenosis, n_nmean:dyx)) %>%
		tbl_summary(
			by = stenosis,
			missing = "no",
			label = labs
		) %>%
		add_p() %>%
		modify_header(update = list(
			label ~ '**Characteristic**',
			stat_1 ~ '**Nonobstructive CAD**, N = 22',
			stat_2 ~ '**Obstructive CAD**, N = 34',
			p.value ~ '**p-value**'
		))


	# HRV by Timing Context --------------------------------------------------
	timing <-
		ecg$timed %>%
		select(patid, context, n_nmean:ap_en, -pnn50) %>%
		filter(context %in% c("start", "balloon")) %>%
		left_join(clinical[c("patid", "pci_adm")], ., by = "patid") %>%
		filter(!is.na(pci_adm) & !is.na(context)) %>%
		mutate(pci_adm = factor(pci_adm, levels = c(0, 1), labels = c("No Revascularization", "Revascularization"))) %>%
		select(-patid) %>%
		tbl_strata(
			strata = pci_adm,
			.tbl_fun =
				~.x %>%
				tbl_summary(
					by = context,
					digits = all_continuous() ~ 1,
				) %>%
				add_p(),
			.combine_with = "tbl_merge"
		) %>%
		modify_header(update = list(
			label ~ '**Characteristic**',
			stat_1_1 ~ '**Angiography**, N = 13',
			stat_2_1 ~ '**Before**, N = 14',
			p.value_1 ~ '**p-value**',
			stat_1_2 ~ '**Balloon**, N = 10',
			stat_2_2 ~ '**Before**, N = 13',
			p.value_2 ~ '**p-value**'
		))

	# Make list
	tables <- list(
		one = one,
		sad = sad,
		revasc = revasc,
		cad = cad,
		timing = timing
	)

	# Return
	tables
}
