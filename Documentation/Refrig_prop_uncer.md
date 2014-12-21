#Uncertainty calculation of dewpoint temperature

This document describes how the uncertainty of the dewpoint temperature of the equation of state of refrigerant and the uncertainty of the dewpoint pressure measurement propagate to the dewpoint temperature.

##Background

Compressor maps require dewpoint temperature as inputs, and the temperature is calculated from pressure readings and the equation of state. But the equation of state, just like the pressure readings, has its own accuracy. This document describes how the two uncertainties combine to form the uncertainty of the dewpoint temperature using the pseudo-fluid equation of state of R410A as an example.

##Uncertainty from pressure reading

Consider the calculation of dewpoint temperature from refrigerant pressure as 

![equation](http://www.sciweavers.org/tex2img.php?eq=T_%7Bdew%7D%3Df%28P%29&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)

The uncertainty propagated from the pressure reading can be estimated by 

![equation](http://www.sciweavers.org/tex2img.php?eq=u_%7BT_%7Bdew%7D%2CP%7D%3D%5Cfrac%7B%5Cpartial%20f%7D%7B%5Cpartial%20P%7D%20u_%7BP%7D&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)

where u<sub>P</sub> is the uncertainty of the pressure reading.

The partial derivatives can be estimated by finite difference method as

![equation](http://www.sciweavers.org/tex2img.php?eq=%5Cfrac%7B%5Cpartial%20f%7D%7B%5Cpartial%20P%7D%20%5Capprox%20%5Cfrac%7Bf%28P%2Bh%29-f%28P-h%29%7D%7B2h%7D&bc=White&fc=Black&im=png&fs=12&ff=arev&edit=0)

where h is a small number.

##Uncertainty from the equation of state

The pseudo-fluid equation of state of R410A only states the uncertainty as the accuracy of the dewpoint pressure as 0.5% but not the accuracy of the dewpoint temperature. To examine how it is propagated to the result of the dewpoint temperature, an equation with dewpoint pressure as the independent variable is used as

![equation](http://www.sciweavers.org/tex2img.php?eq=P%3Dg%28T_%7Bdew%7D%29&bc=White&fc=Black&im=png&fs=12&ff=arev&edit=0)

After calculation, the uncertainty of the dewpoint temperature can be expressed as

![equation](http://www.sciweavers.org/tex2img.php?eq=u_%7BT_%7Bdew%7D%2CEOS%7D%3D%28%5Cfrac%7B%5Cpartial%20g%7D%7B%5Cpartial%20T_%7Bdew%7D%7D%29%5E%7B-1%7D%20u_%7BP%2C%20EOS%7D&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)

where u<sub>P,EOS</sub> is the absolute uncertainty of the equation of state to estimate dewpoint pressure.

##Final calculation

The uncertainty components are combined together to estimate the uncertainty of the dewpoint temperature as

![equation](http://www.sciweavers.org/tex2img.php?eq=u_%7BT_%7Bdew%7D%7D%3D%5Csqrt%7Bu_%7BT_%7Bdew%7D%2C%20P%7D%5E2%2Bu_%7BT_%7Bdew%7D%2C%20EOS%7D%5E2%7D&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)