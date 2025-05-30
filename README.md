# SurveyTogether #
 
This repository serves as the the compliment to ``Integrating representitive and non-representitive surveys for efficient inferenece" to reproduce the results and figures.

Scripts
----------------
Includes the R scripts to reproduce the results.

**proofofconcept.R**: includes the functions to generate the parameters, data , and the ability to fit the model using parallel cores (Section 2.2).

**par_survey_together_<t>binom2.R**: are three independent scripts for the simulation study. This includes generating the data for T = 5, 10, 15 time-points and functions to fit the data over the 2,000 repetitions. These scripts then save the data including with the corresponding MCSE. (These simulations were run a Digital Alliance of Canada computing cluster) (Simulation study Section 3)


**extendeddata.R**: This script includes cleaning some of the data for the vaccination application, formatting it so that it may be used in the modeling, and fitting the model according to different parameritizations of $\phi_{kt}$. The JAGS models for the models are found in JagsMod.R (Section 4).

**NowCastperSurvey.R**: Is a repetition of the extended data, but the modeling is performed 48 times to get the now-casting results (Section 4.3). Includes calculating the $n_{iid}$ for the synthesis method.

**prior_quantiles.R**: a script that reviews the interpretations and quantiles of the priors (Section 2).

**other scripts**: old scripts exploring the identifiability of $\phi_{kt}$.

Some of the data in the list format can he found in the main folder. 

Includes data presented in <br> BRADLEY, V. C., KURIWAKI, S., ISAKOV, M., SEJDINOVIC, D., MENG, X.-L. and FLAXMAN, S.<br>  (2021).
 \"**Unrepresentative big surveys significantly overestimated US vaccine uptake**"\ _Nature_ 600 695–700. Dec 2021, doi:10.1038/s41586-021-04198-4 .

and the extended data.

