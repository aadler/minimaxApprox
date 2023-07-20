.onAttach <- function(libname, pkgname) {
  # Warn Mac M1 users
  # See https://stackoverflow.com/questions/69913603/skip-test-on-m1mac-in-r
  problemSystems <- c("aarch64", "aarch32", "ppc", "ppc64")
  if ((tolower(Sys.info()[["sysname"]]) == "darwin" &&
       R.version[["arch"]] %in% problemSystems)) {
    packageStartupMessage("As per R Installation and Administration, the ",
                          "native build for Apple Silicon or Power PC may ",
                          "give different numerical results from other ",
                          "systems, or even fail where they succeed, ",
                          "as that hardware either lacks extended-precision ",
                          "floating-point operations or handles them ",
                          "differently.")
  }
}
