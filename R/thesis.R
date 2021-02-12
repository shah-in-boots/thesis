# Write dissertation
write_dissertation <- function(index, chapters) {

	# Load all the chapters by targets to see all chapter sections
	# Set up book rendering
	bookdown::render_book(
		input = index,
		output_format = "bookdown::pdf_book",
		config_file = "_bookdown.yml"
	)

}
