# Function to create and load the virtual environment to use the .py gotcha script
.onAttach <- function(libname, pkgname){
  packageStartupMessage("Checking if r-reticulate-gotcha virtual environment is available...")
  if(!reticulate::virtualenv_exists(envname = "r-reticulate-gotcha")){
    packageStartupMessage("r-reticulate-gotcha virtual environment not found, setting up environment...")

    # Create isolated virtual environment using reticulate to install python modules
    reticulate::virtualenv_create(envname = "r-reticulate-gotcha",
                                  python = Sys.which("python3"),
                                  version = "python3",
                                  packages = c("pip","wheel"),
                                  pip_version = "pip-22.1.1")

    # Set up reticulate to use the virtual environment
    reticulate::use_virtualenv("r-reticulate-gotcha")

    # Install required python modules
    reticulate::py_install(packages = c("pandas","matplotlib","seaborn","numpy","sklearn"),
                           python_version = "python3",
                           pip = T,
                           envname = "r-reticulate-gotcha")


  }else{

    # Use the r-reticulate-gotcha virtual environment
    reticulate::use_virtualenv("r-reticulate-gotcha")
    packageStartupMessage("Virutal environment r-reticulate-gotcha is available | use_virtualenv(r-reticulate-gotcha)")

  }

  # Source labeling functions from GitHub
  packageStartupMessage("Sourcing python script from GitHub...")
  reticulate::source_python(file = "https://github.com/landau-lab/Gotcha/blob/main/Gotcha/Python/gotcha_labeling.py?raw=TRUE", envir = globalenv())
}