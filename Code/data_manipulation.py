#/usr/bin/python

"""
	This file contains functions that are used to predict measurement
	data and its uncertainty based on expected readings and its
	experimental methods
"""

from math import floor

import numpy as np
from random import normalvariate
from scipy.stats import norm

import exp_uncer

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
            Return the number of time-series data points
        """
        #get the number of data points in the time-series data
        return int(floor(self._period*self._freq))


def exp_data_generator(exp_val, exp_method, apparatus_uncer, conf=0.95):
    """
        Predict experimental data according to the experimental
        method and apparatus uncertainty information with the
        assumption that the apparatus uncertainty is described
        by a 95% confidence interval under a normal distribution
        
        Parameters:
        ===========
        exp_val: float
            expected reading of the sensor
        exp_method: class EXP_METHOD
            experimental method information
        apparatus_uncer: class exp_uncer.APPARATUS_UNCER
            uncertainty information of the measurement apparatus
        conf: float, optional
            the confidence level of the apparatus uncertainty. Default 95%
            
        Returns:
        ===========
        data: numpy array
            readings in time-series
        
    """
    
    num = exp_method.get_num()  #number of data points
    zero_order = apparatus_uncer.zero_order_uncer(exp_val)
    
    norm_std = zero_order/norm.interval(0.95)[1]
    data = []
    for ii in xrange(num):
        data.append(normalvariate(exp_val, norm_std))
        
    return np.array(data)


if __name__ == "__main__":

    """
        Testing scripts
    """
    
    #import libraries
    from math import sqrt
    
    from CoolProp.CoolProp import PropsSI
    
    import exp_uncer
    
    # define experimental method and apparatus
    test = EXP_METHOD(0.1, 600)
    power_meter = exp_uncer.APPARATUS_UNCER(0.0, 0.005)
    p_trans = exp_uncer.APPARATUS_UNCER(0.0, 0.05)
    
    # expected readings
    power = 500.0
    p_suc = 500.0
    p_dischg = 1000.0
    refri = 'R410A'
    
    # generate random time-series data from expected readings with exp_method
    # and apparatus information
    power_data = exp_data_generator(power, test, power_meter)
    p_suc_data = exp_data_generator(p_suc, test, p_trans)
    p_dischg_data = exp_data_generator(p_dischg, test, p_trans)
    
    # return the mean values of the time-series data
    print('Power meter reading:')
    print(power_meter.measure_av_result(power_data, full_output=1))
    
    print('Compressor suction pressure reading:')
    print(p_trans.measure_av_result(p_suc_data, full_output=1))
    
    print('Compressor discharge pressure reading:')
    print(p_trans.measure_av_result(p_dischg_data, full_output=1))
    
    