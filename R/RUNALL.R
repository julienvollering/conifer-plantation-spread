Sys.setenv(RSTUDIO_PANDOC="/usr/lib/rstudio/bin/pandoc")

rmarkdown::render("./reports/ClimateData.Rmd", knit_root_dir = "../")
rmarkdown::render("./reports/locality-WALD.Rmd", knit_root_dir = "../")
rmarkdown::render("./reports/wildlings-NiN.Rmd", knit_root_dir = "../")
rmarkdown::render("./reports/locality-maps.Rmd", knit_root_dir = "../")
rmarkdown::render("./reports/exploringData.Rmd", knit_root_dir = "../")
rmarkdown::render("./reports/modelingData.Rmd", knit_root_dir = "../")
rmarkdown::render("./reports/hypothesisComparisons.Rmd", knit_root_dir = "../")