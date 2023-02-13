
library(GenABEL.data)
library(GenABEL)
library(melonnpan)
library(ggplot2)

###Installation of MelonnPan####

#first, need Rtools
devtools::install_version("Rtools",repos = "http://cran.us.r-project.org")

# then need to link Rtools filepath 
Sys.getenv("PATH")
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "C:\\rtools40\\usr\\bin", "C:\\rtools40\\mingw64\\bin", sep = ";")) #need to replace these with filepaths for usr and mingw64 which should be where you installed Rtools
Sys.setenv(BINPREF = "C:/rtools40/mingw64/bin/") #need to replace this with filepath for mingw64
Sys.which("make")#should now show you filepath for make.exe
Sys.which("g++")#should now show you filepath for g++.exe

#then need to install Melonnpan Package dependencies before installing melonnpan
devtools::install_version("GenABEL.data", version = "1.0.0", repos = "http://cran.us.r-project.org")
devtools::install_version("GenABEL", version = "1.8-0", repos = "http://cran.us.r-project.org")
devtools::install_github("biobakery/melonnpan")

###Running of MelonnPan
predict_metabolites()
