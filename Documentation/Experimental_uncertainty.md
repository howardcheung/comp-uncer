#Uncertainty calculation of steady state observation

This document describes how the uncertainty of a steady state data point is obtained.

A steady state data point is obtained by averaging data collected for a period of time when the experimental setup is under steady state. There are two contributions of its uncertainty: zero-order uncertainty and first-order uncertainty.

##Zero-order uncertainty

Zero-order uncertainty is the uncertainty of a sensor reading in a standardized environment. Being imperfect, a sensor in a standardized environment may give a reading different from the true reading of the standard environment. For instance, when a thermocouple placed in an percetly 0C ice water bath, it may still give a temperature 24.8C because of its deficiency. Zero-order uncertainty is the possible range of readings that a measurement apparatus can give under this situation.

This uncertainty is usually given by the manufacturer specification or professional standard. For instance, the zero-order uncertainty of a T-type thermocouple is given in the 2006 ASHRAE standard 41.2 as the instrument accuracy at 0.5K. The other uncertainty due to thermal expansion or change of humidity in the environment can also be counted towards this type of uncertainty. In this case, each measurement in a time-series data carries its own zero-order uncertainty.

##First-order uncertainty

First-order uncertainty is the uncertainty contributed by the noise in the environment or the surroundings. It is estimated by taking the a reading in an environment through a period of time. The data should fluctutate around a mean reading due to the noise in the environment, and the first-order uncertainty is estimated by interpreting the fluctuation quantitatively with statistics.

The fluctuation of data with time is usually examined by the sample standard deviation of the time-series data because the data are collected in a finite period of time. Population standard deviation is only used when the data are taken in a inifinite period of time which is impossible for experiments. The sample standard deviation is given by

![equation](http://www.sciweavers.org/tex2img.php?eq=%5Csigma%3D%5Csqrt%7B%5CSigma_i%20%5Cfrac%7B%28x_i%20-%20%5Coverline%7Bx%7D%29%5E2%7D%7Bn-1%7D%7D&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)

where x<sub>i</sub> is the i<sup>th</sup> observation in the time-series data, x(bar) is the mean value of the time-series data and n is the number of observations in the time-series data.

As we are looking at the first-order uncertainty of a steady-state measurement observed by averaging data in a finite period of time, the sample standard deviation of the mean, instead of sample standard deviation, should be used. The sample standard deviation of the mean is given by

![equation](http://www.sciweavers.org/tex2img.php?eq=s%20%3D%20%5Cfrac%7B%5Csigma%7D%7B%5Csqrt%7Bn%7D%7D&bc=White&fc=Black&im=png&fs=12&ff=arev&edit=0)

Assuming that the noise follows a Student T-distribution and its confidence level is 95%, the first-order uncertainty of the steady-state observation is given by

![equation](http://www.sciweavers.org/tex2img.php?eq=u_%7B1st%7D%3Dt_%7B1-0.95%2Cn-1%7Ds&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)

##Final uncertainty of steady state observation

The uncertainty of a steady state observation is given by

![equation](http://www.sciweavers.org/tex2img.php?eq=u_%7B%5Coverline%7Bx%7D%7D%3D%5Csqrt%7B%5Cfrac%7B1%7D%7Bn%7D%5CSigma_i%28u_%7B0th%2Ci%7D%5E2%29%20%2B%20u_%7B1st%7D%5E2%20%7D&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)