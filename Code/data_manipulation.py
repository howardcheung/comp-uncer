#!/usr/bin/python

"""
    This file contains functions that are used to predict measurement
    data and its uncertainty based on expected readings and its
    experimental methods
"""

from math import floor, sqrt

from CoolProp.CoolProp import PropsSI as Props
import numpy as np
from random import normalvariate
from scipy.stats import norm

import exp_uncer
import misc_func


class EXP_METHOD:
    """
        This class contains information about how a steady state
        data point was obtained in an experiment

        freq: float
            frequency of data acquisition in Hz
        period: float
            length of data acquisition period in seconds

    """

    def __init__(self, freq=0.0, period=0.0):
        self._freq = freq
        self._period = period

    def set_freq(self, freq):
        self._freq = freq

    def get_freq(self):
        return self._freq

    def set_period(self, period):
        self._period = period

    def get_period(self):
        return self._period

    def get_num(self):
        """
            Return the number of time-series data points as an integer
        """
        #get the number of data points in the time-series data
        return int(floor(self._period*self._freq))

    def exp_data_generator(
        self, exp_val, apparatus_uncer, first_abs=0.0,
        first_rel=0.0, conf=0.95
    ):
        """
            Predict experimental data according to the experimental
            method and apparatus uncertainty information with the
            assumption that the apparatus uncertainty is described
            by a 95% confidence interval under a normal distribution

            Parameters:
            ===========
            exp_val: float
                expected reading of the sensor
            apparatus_uncer: class exp_uncer.APPARATUS_UNCER
                uncertainty information of the measurement apparatus
            first_abs: float
                extra absolute uncertainty to the measurement
                due to noise, in the same engineering unit as the
                measurement
            first_rel: float
                extra relative uncertainty to the measurement
                due to noise
            conf: float, optional
                the confidence level of the apparatus uncertainty. Default 95%

            Returns:
            ===========
            data: numpy array
                readings in time-series

        """

        num = self.get_num()  # number of data points
        zero_order = apparatus_uncer.zero_order_uncer(exp_val)
        first_order = sqrt(
            first_abs**2+(first_rel*exp_val)**2
        )

        norm_std = sqrt(
            zero_order**2+first_order**2
        )/norm.interval(0.95)[1]
        data = []
        for ii in xrange(num):
            data.append(normalvariate(exp_val, norm_std))

        return np.array(data)


def sat_temp_uncer_cal(pres, uncer_pres, refri, full_output=0):
    """
        Calculate the uncertainty of the saturation temperature
        as a result of the refrigerant property calculation
        from the mean pressure value averaged from a time-series data

        Parameters:
        ============
        pres:   float, list or numpy array
            pressure readings averaged from time-series data in kPa
        uncer_pres: float, list or numpy array
            uncertainty due to sensor and time-fluctuation
            of pressure readings averaged from time-series data in kPa
        refri:  string
            name of refrigerant
        full_output: int, optional
            return an array of showing uncertainty propagated from
            zero- and first-order uncertainty of pressure readings
            and the uncertainty due to the refrigerant property equation

        Returns:
        ============
        Tsat:   float or numpy array
            dewpoint temperature in K
        uncer_Tsat:   float or numpy array
            uncertainty of Tsat in K
        uncer_read:   float or numpy array, optional
            uncertainty from time-fluctuation and sensor of
            pressure readings in K
        uncer_theo:   float or numpy array, optional
            uncertainty from refrigerant property calculation in K
    """

    #make all pressure and uncer_pressure readings to be in numpy arrays
    pres_ary = misc_func.check_size(pres)
    uncer_ary = misc_func.check_size(uncer_pres)
    thres = 1.e-6

    #calculate the uncertainty propagated to the saturation temperature result
    #from the sensors and time fluctuation

    def _PropsTP(p):
        return Props('T', 'P', p*1000.0, 'Q', 1, refri)

    dTdP_ary = [
        misc_func.finite_difference(
            p, _PropsTP, p*thres, thres
        ) for p in pres_ary
    ]
    uncer_read = [dTdP*uncer for dTdP, uncer in zip(dTdP_ary, uncer_ary)]

    # NEED MORE DIFFERENT REFRIGERANTS LATER
    # use the R410A default that the dewpoint pressure
    # uncertainty is 0.5% in Lemmon(2003) to calculate
    # the uncertainty due to the equation

    if refri is 'R22' or refri is 'R134a':
        # for R22, use the relative uncertainty from Kamei et al. (1995)
        # for R134a, use the relative unceratinty from Tillner-Roth et al. (1994)
        p_uncer_rel = 0.002
    else:  # use R410A ones as default
        p_uncer_rel = 0.005

    Tsat_ary = [Props('T', 'P', p*1000.0, 'Q', 1, refri) for p in pres_ary]

    uncer_theo = [
        dTdP*p*p_uncer_rel for dTdP, uncer in zip(dTdP_ary, pres_ary)
    ]

    # calculate total uncertainy of the saturation temperature
    uncer_total = np.sqrt(
        (np.array(uncer_theo))**2+(np.array(uncer_read))**2
    )
    if len(uncer_total) == 1:
        if full_output:
            return Tsat_ary[0], uncer_total[0], uncer_read[0], uncer_theo[0]
        else:
            return Tsat_ary[0], uncer_total[0]
    else:
        if full_output:
            return Tsat_ary, uncer_total, uncer_read, uncer_theo
        else:
            return Tsat_ary, uncer_total


if __name__ == "__main__":

    """
        Testing scripts
    """

    #import libraries
    from math import sqrt
    import random

    import exp_uncer

    # define experimental method and apparatus
    test = EXP_METHOD(0.1, 600)
    power_meter = exp_uncer.APPARATUS_UNCER(0.0, 0.005)
    p_trans = exp_uncer.APPARATUS_UNCER(0.0, 0.008)

    # expected readings
    power = 2000.0
    p_suc = 500.0
    p_dischg = 1000.0
    refri = 'R410A'
    random.seed(10)  # define seed

    # generate random time-series data from expected readings with exp_method
    # and apparatus information
    power_data = test.exp_data_generator(power, power_meter)
    p_suc_data = test.exp_data_generator(p_suc, p_trans, first_abs=0.9)
    p_dischg_data = test.exp_data_generator(p_dischg, p_trans, first_abs=0.4)

    # return the mean values of the time-series data
    print('Power meter reading:')
    print(
        'Mean of time series [W], Uncertainty [W], ' +
        'Uncertainty from sensor [W], Uncertainty from time fluctuation [W]:'
    )
    print(power_meter.measure_av_result(power_data, full_output=1))

    p_suc_info = p_trans.measure_av_result(p_suc_data, full_output=1)
    print('Compressor suction pressure reading:')
    print(
        'Mean of time series [kPa], Uncertainty [kPa],' +
        ' Uncertainty from sensor [kPa], ' +
        'Uncertainty from time fluctuation [kPa]:'
    )
    print(p_suc_info)

    p_dischg_info = p_trans.measure_av_result(p_dischg_data, full_output=1)
    print('Compressor discharge pressure reading:')
    print(
        'Mean of time series [kPa], Uncertainty [kPa],' +
        ' Uncertainty from sensor [kPa], ' +
        'Uncertainty from time fluctuation [kPa]:'
    )
    print(p_dischg_info)

    print('Compressor suction saturation temperature reading:')
    print(
        'Dewpoint temperature [K], Uncertainty [K],' +
        ' Uncertainty from pressure reading [K], ' +
        'Uncertainty from equation [K]:'
    )
    print(sat_temp_uncer_cal(
        p_suc_info[0], p_suc_info[1], refri, full_output=1
    ))

    print('Compressor discharge saturation temperature reading:')
    print(
        'Dewpoint temperature [K], Uncertainty [K],' +
        ' Uncertainty from pressure reading [K], ' +
        'Uncertainty from equation [K]:'
    )
    print(sat_temp_uncer_cal(
        p_dischg_info[0], p_dischg_info[1], refri, full_output=1
    ))
