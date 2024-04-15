$$ SurveyTogether$$
 
This repository serves as the the compliment to ``Integrating representitive and non-representitive surveys for efficient inferenece" to reproduce the results and figures.

$$Scripts$$
Includes the R scripts to reproduce the results.

proofofconcept.R: includes the functions to generate the parameters, data , and the ability to fit the model using parallel cores.

par_survey_together_<t>binom2.R: are three independent scripts for the simulation study. This includes generating the data for that T = 5, 10, 15 time-points and functions to fit the data over the 2,000 repetitions. Saves the data including MCSE. (These simulations were run a Digital Alliance of Canada computing cluster)


extendeddata.R: This script includes cleaning some of the data for the vaccination application, formatting it so that it may be used in the modeling and fitting the model according to different parameritizations of $\phi_{kt}$. The JAGS models for the models are found in JagsMod.R .

NowCastperSurvey.R: Is a repetition of the extended data but repeating the modeling 48 times to get the now-casting results. Includes calculating the $n_{iid}$ for the synthesis method.

test_sigmasq.R: a script that reviews the interpretations and quantiles of the priors (methods section).
