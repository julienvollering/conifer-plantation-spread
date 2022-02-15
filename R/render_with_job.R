# CC-BY https://github.com/sellorm
# From https://gist.github.com/sellorm/20525aff669acafb975b7495b8f2066e#gistcomment-3707855

render_with_job <- function() {
  rstudioapi::verifyAvailable()
  
  job_file <- tempfile(fileext = ".R")
  
  active_doc_ctx <- rstudioapi::getSourceEditorContext()
  rmd_path <- active_doc_ctx$path
  
  if (identical(rmd_path, "")) {
    rstudioapi::showDialog(
      "Cannot Render Unsaved R Markdown Document",
      "Please save the current document before rendering."
    )
    return(invisible())
  }
  
  rstudioapi::documentSave(active_doc_ctx$id)
  
  rmd_path <- normalizePath(rmd_path, mustWork = TRUE)
  
  cat(
    'res <- rmarkdown::render("', basename(rmd_path), '")\n',
    'unlink("', job_file, '")\n',
    'rstudioapi::viewer(res)',
    sep = "",
    file = job_file
  )
  
  rstudioapi::jobRunScript(
    path = job_file,
    name = basename(rmd_path),
    workingDir = dirname(rmd_path),
    importEnv = FALSE
  )
} 