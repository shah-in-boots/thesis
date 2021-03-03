# Biobank Tables
make_biobank_tables <- function(clinical, ecg, labels) {

	# Table 1
	one <-
		clinical %>%
		mutate(
			race = labels$Race,
			gend = labels$Gender
		) %>%
		labelled::set_variable_labels(
			race = "Race",
			gend = "Sex"
		) %>%
		select(c(age, race, blbmi, gend, phq, sad, gensini, stenosis, cass70)) %>%
		tbl_summary(
			missing = "no",
			value = list(c(stenosis, sad) ~ "1")
		)

	# HRV Table
	hrv <-
		ecg$merged %>%
		group_by(patid) %>%
		summarise(across(c(n_nmean:ap_en, dyx), ~ mean(.x, na.rm = TRUE)), .groups = "keep") %>%
		left_join(clinical, ., by = "patid") %>%
		select(c(stenosis, n_nmean:dyx)) %>%
		tbl_summary(
			by = stenosis,
			missing = "no",
			value = list(c(stenosis) ~ "1")
		) %>%
		add_p() %>%
		modify_header(update = list(
			label ~ '**HRV Metric**',
			stat_1 ~ '**No Revascularization** <br> N = 14',
			stat_2 ~ '**Revascularization** <br> N = 34',
			p.value ~ '**p-value**'
		))

	# Comparing HRV by context
	x <- ecg$timed %>%
		select(patid, context, n_nmean:ap_en) %>%
		filter(context %in% c("pre", "sedation", "balloon", "end")) %>%
		pivot_longer(cols = n_nmean:ap_en, names_to = "hrv") %>%
		pivot_wider(names_from = context, values_from = "value") %>%
		left_join(clinical[c("patid", "stenosis")], ., by = "patid") %>%
		select(-patid) %>%
		filter(!is.na(stenosis) & !is.na(hrv)) %>%
		group_by(stenosis, hrv) %>%
		summarise(across(balloon:sedation, list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ sd(.x, na.rm = TRUE)), .names = "{.col}_{.fn}"), .groups = "keep") %>%
		ungroup() %>%
		group_by(stenosis) %>%
		gt(rowname_col = "hrv")

	# HRV by Context
	timing <-
		ecg$timed %>%
		select(patid, context, n_nmean:ap_en) %>%
		filter(context %in% c("start", "balloon")) %>%
		left_join(clinical[c("patid", "stenosis")], ., by = "patid") %>%
		filter(!is.na(stenosis) & !is.na(context)) %>%
		mutate(stenosis = factor(stenosis, levels = c(0, 1), labels = c("No Revascularization", "Revascularization"))) %>%
		select(-patid) %>%
		tbl_strata(
			strata = stenosis,
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
			label ~ '**ECG Metrics**',
			stat_1_1 ~ '**Balloon** <br> N = 6',
			stat_2_1 ~ '**Start** <br> N = 5',
			p.value_1 ~ '**p-value**',
			stat_1_2 ~ '**Balloon** <br> N = 15',
			stat_2_2 ~ '**Start** <br> N = 20',
			p.value_2 ~ '**p-value**'
		))

	# Make list
	tables <- list(
		one = one,
		hrv = hrv,
		timing = timing
	)

	# Return
	tables
}
