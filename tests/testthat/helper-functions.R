skip_if_not_travis_or_bioc <- function() {
    if (identical(Sys.getenv("TRAVIS"), "true") ||
        !identical(Sys.getenv("BBS_HOME"), "")) {
        return(invisible(TRUE))
    }
    
    skip("Not on Travis or Bioconductor")
}

skip_if_not_in_path <- function(program) {
    tryCatch(
        {
            .checkDep(program)
            return(invisible(TRUE))
        }, error=function(e) {
            skip(paste(program, 'not in path'))
        })
}