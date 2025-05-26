# HPV-COMPARE
R-code to evaluate the cost-effectiveness profile of 9v vs 2v HPV vaccination in the Netherlands. \
Used in _Projected health and economic effects of nonavalent versus bivalent HPV vaccination in the Netherlands: a data-driven analysis_ (link to paper).

# .R files
_HPVcompare_functions.R_ \
_HPVcompare_RiskReduction_functions.R_ \
_HPVcompare_ParametersData.R_ \
_HPVcompare_RiskReductions.R_ \
_HPVcompare_Computation.R_ \
_AnogenitalWarts_Incidence.R_ 

# Data
The _Inc__ files contain cancer incidence data from the Netherlands Cancer Registry for the various cancer types \
The _Surv__ files contain cancer survival data from the Netherlands Cancer Registry for the various cancer types \
The _TablePop__ files contain overall survival data of the general population in 2020 from the Statistics Netherlands database \
The _Pop_ files contain data on the population size stratified by age group from the Statistics Netherlands database \
The _BirthsPer1000__ file contains data on the average number of births in the population (per age of the mother) from the Statistics Netherlands database

# Transmission model
The folder _TransmissionModel_ contains all the output from the transmission model needed to compute the HPV risk reductions in the code. 

# Example
_HPVcompare_Example.R_ contains code with which an example scenario can be run:
- It first loads all the required functions by loading _HPVcompare_functions.R_ and _HPVcompare_RiskReduction_functions.R_
- It inputs all the data and simulates the parameters by loading _HPVcompare_ParametersData.R_
- It defines the example scenario
- It computes the corresponding risk reductions by loading _HPVcompare_RiskReductions.R_
- It runs the analysis by loading _HPVcompare_Computation.R_
- Various results can be found in the list _Results_


To run the example; save the folders _Rcode_, _Data_ and _TransmissionModel_ (including the subfolders and content) in a main folder.
Change your R working directory to this main folder and run the _HPVcompare_Example.R_ from this main folder.
