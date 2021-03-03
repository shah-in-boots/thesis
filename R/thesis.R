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

	# Overview DAG
	dag <- dagify(
		ans ~ cad,
		out ~ cad,
		ans ~~ ms + psych + circ,
		stress ~ psych + ms + circ + ps,
		cad ~ ps + ms,
		out ~ stress + cad + ans,
		labels = c(
			"cad" = "Coronary Artery Disease",
			"ans" = "Autonomic Nervous System",
			"circ" = "Circadian Disruption",
			"ms" = "Acute Mental Stress",
			"ps" = "Acute Physical Stress",
			"psych" = "Psychological Distress",
			"stress" = "Stress Reactivity",
			"out" = "Outcomes"
		),
		outcome = "out",
		exposure = c("ans", "stress")
	) %>%
		tidy_dagitty(layout = "kk") %>%
		node_status() %>%
		ggplot(aes(x = x, y = y, xend = xend, yend = yend, color = status)) +
		geom_dag_point() +
		geom_dag_edges_link() +
		geom_dag_label_repel(
			aes(label = label), size = 3,
			nudge_y = -5, nudge_x = -5
		) +
		theme_dag()

	# List
	diagrams <- list(
		dag = dag
	)

	# Return
	diagrams
}
