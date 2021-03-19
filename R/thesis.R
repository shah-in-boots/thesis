# Write dissertation
write_dissertation <- function(index, ...) {

	# Load all the chapters by targets to see all chapter sections
	# Set up book rendering
	bookdown::render_book(
		input = index,
		output_format = "bookdown::pdf_book",
		config_file = "_bookdown.yml"
	)

}

# Supporting schematics and diagrams
draw_diagrams <- function() {

	# Overview DAG -------------------------------------------------------
	dag <- dagify(
		psych ~ ans,
		ans ~ stress,
		ihd ~ ans,
		mace ~ ans + ihd + psych,
		labels = c(
			"ihd" = "Ischemic Heart Disease",
			"ans" = "Autonomic Nervous System",
			"stress" = "Stress",
			"psych" = "Psychological Distress",
			"mace" = "Major Adverse Cardiovascular Events"
		),
		outcome = c("ihd", "psych", "mace"),
		exposure = c("ans")
	) %>%
		tidy_dagitty(layout = "gem", seed = 16) %>%
		node_status() %>%
		ggplot(aes(x = x, y = y, xend = xend, yend = yend, color = status)) +
		geom_dag_point() +
		geom_dag_edges_link() +
		geom_dag_label_repel(
			aes(label = label), size = 3,
			nudge_y = -5, nudge_x = -5
		) +
		theme_dag()


	# Example of Cosinor ----------------------------------------------------
	data(triplets)
	single_lf <- cosinor(LF ~ hour, tau = 24, data = triplets)
	multi_lf <- cosinor(LF ~ hour, tau = c(24, 12), data = triplets)
	gg <- ggcosinor(single_lf)

	# Returning list --------------------------------------------------------
	diagrams <- list(
		dag = dag,
		cosinor = gg
	)

	# Return
	diagrams
}
