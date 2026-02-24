library(rmarkdown)

dir.create("Docs", showWarnings = FALSE)

rmds <- list.files("Code/", pattern = "\\.Rmd$", full.names = TRUE)

for (f in rmds) {
  output_file <- file.path("Docs", tools::file_path_sans_ext(basename(f)))
  
  render(
    f,
    output_format = github_document(html_preview = FALSE), # extract relative filename path without .Rmd/.md extension
    output_file = output_file,
    output_dir = dirname(output_file),  # Just "Docs"
    quiet = TRUE
  )
}
