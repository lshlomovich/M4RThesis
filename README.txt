README

MultitaperEst.m
   File containing code to generate MA and AR processes. The multitapering function is used to estimate and compare to the known sdf. This code tests the multitapering estimating function for p=2.
Spectrum can be plotted and for just a poisson process, should be flat at lambda*T/m

TestingMultitaper.m
   File containing code to generate QQ plots for coherence vs Goodman Distribution (code for cdf in Goodman_QQ_Plots.m) 

TestingMultitaperPingPartial.m
   File containing code to generate QQ plots for partial coherence vs Beta Distribution, plus a test against the F distribution. Beta performed better. 

regularpulsesims.m
   This file contains code needed to simulate figures for this section of the thesis. That is, looking at point processes formed of a Poisson process and a regular pulse/`impulse' process.

myfejer.m 
   Code to generate a plot of Fejer's kernel. 

ShiftedPingPoisson.m
   Code used for the simulation where we investigate shifting the pulse process. 

BernoulliTrialledCts.m
   Contains code to generate simulations when Bernoulli trialling. User can change variables as desired in order to test different probabilities, periodicities and rates.





