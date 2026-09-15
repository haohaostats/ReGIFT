.comparator_root <- local({
  paths <- lapply(sys.frames(), function(frame) frame$ofile)
  paths <- Filter(function(x) is.character(x) && length(x) == 1L, paths)
  if (!length(paths)) stop("Load this file with source().")
  dirname(normalizePath(tail(paths, 1L)[[1L]], winslash = "/", mustWork = TRUE))
})
.comparator_script <- local({
  root <- .comparator_root
  function(name) file.path(root, "python", name)
})
for (.comparator_file in list.files(file.path(.comparator_root, "R"),
                                   pattern = "\\.R$", full.names = TRUE)) {
  sys.source(.comparator_file, envir = environment())
}
rm(.comparator_file)
