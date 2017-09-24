This package constitutes an interactive R problem set based on the RTutor package (https://github.com/skranz/RTutor). 

The interactive R problem set "Pollution Reductions from Renewable Energies: Social Benefits and Policy Implications" guides users through the paper "Valuing the Wind: Renewable Energy Policies and Air Pollution Avoided" by Kevin Novan (2015). Novan estimates the impact of renewable energies on pollution from conventional power plants and discusses policy implications from his findings. The problem set reproduces the main analysis while familiarising the user with necessary R functions and the econometric background. The paper can be found at http://dx.doi.org/10.1257/pol.20130268.

## 1. Installation

RTutor and this package is hosted on Github. To install everything, run the following code in your R console.
```s
if (!require(devtools))
  install.packages("devtools")
source_gist("gist.github.com/skranz/fad6062e5462c9d0efe4")
install.rtutor(update.github=TRUE)

devtools::install_github("asbara/RTutorPollutionReductions", upgrade_dependencies=FALSE)
```

## 2. Show and work on the problem set
To start the problem set first create a working directory in which files like the data sets and your solution will be stored. Then adapt and run the following code.
```s
library(RTutorPollutionReductions)

# Adapt your working directory to an existing folder
setwd("C:/problemsets/RTutorPollutionReductions")
# Adapt your user name
run.ps(user.name="Jon Doe", package="RTutorPollutionReductions",
       load.sav=TRUE, sample.solution=FALSE)
```
If everything works fine, a browser window should open, in which you can start exploring the problem set.
