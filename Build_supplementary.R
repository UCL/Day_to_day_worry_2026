library(rmarkdown)

dir.create("Docs", showWarnings = FALSE)

rmds <- list.files("Code/", pattern = "\\.Rmd$", full.names = TRUE)

for (f in rmds) {
  render(f, output_dir = "Docs/")
}
