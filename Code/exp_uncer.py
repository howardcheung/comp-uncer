#/usr/bin/python

"""
    This file contains functions and classes that help the calculation of
    measurement uncertainty of a variable obtained by averaging
    data collected in a time series
"""

from math import sqrt

import numpy as np
from scipy.stats import t


class APPARATUS_UNCER:
    """
        This class stores the information of measurement uncertainty
        of different apparatus

        _abs_uncer: float
                absolute uncertainty
        _rel_uncer: float
                relative uncertainty
    """

    def __init__(self, abs_uncer=0.0, rel_uncer=0.0):
        self._abs_uncer = abs_uncer
        self._rel_uncer = rel_uncer

    def set_abs_uncer(self, abs_uncer):
        self._abs_uncer = abs_uncer

    def get_abs_uncer(self):
        return self._abs_uncer

    def set_rel_uncer(self, rel_uncer):
        self._rel_uncer = rel_uncer

    def get_rel_uncer(self):
        return self._rel_uncer

    def zero_order_uncer(self, reading):
        """
            This function returns the zero-order uncertainty of the mean value
            of a time-series data based on the information in the class

            Parameters:
            ===========
            reading: array
                readings of the measurement in time-series

            Returns:
            ===========
            uncertainty of the measurement
        """

        try:
            num = len(reading)
            if self.get_abs_uncer() > 0.0:
                return sqrt(self.get_abs_uncer()**2/num)
            else:
                return np.sqrt(np.sum(
                    (self.get_rel_uncer()*np.array(reading))**2
                )/num)
        except TypeError:  #not list nor numpy array
            if self.get_abs_uncer() > 0.0:
                return self.get_abs_uncer()
            else:
                return self.get_rel_uncer()*np.array(reading)
            

    def measure_av_result(self, reading, conf=0.95, full_output=0):
        """
            This function calculates the mean reading and its uncertainty
            based on its reading in a time-series and the apparatus used
            to collect the readings

            Parameters:
            ===========
            reading: array
                readings of the measurement
            conf:   float
                Confidence level. Default is 95%
            full_output: int
                if the number is bigger than one, return zero and first
                order uncertainty for reference

            Returns:
            ===========
            mean_val: float
                mean of the time-series data
            mean_uncer: float
                uncertainty of the measurement
            zero_order: float, optional
                magnitude of zero order uncertainty
            first_order: float, optional
                magnitude of first order uncertainty
        """

        #calculate the mean value
        mean_val = np.mean(reading)

        #calculate the uncertainty of the mean value
        zero_uncer = self.zero_order_uncer(reading)
        first_uncer = first_order_uncer(reading, conf)
        mean_uncer = sqrt(zero_uncer**2+first_uncer**2)

        if full_output <= 0:
            return mean_val, mean_uncer
        else:
            return mean_val, mean_uncer, zero_uncer, first_uncer


def first_order_uncer(reading, conf=0.95):
    """
        This function calculates the first order uncertainty
        of the mean of data in a time series

        Parameters:
        ===========
        reading: array
            readings of the measurement in time-series
        conf:   float
            Confidence level. Default is 95%

        Returns:
        ===========
        first order uncertainty of the mean of the time series data

    """

    #sample standard deviation
    num = len(reading)
    sample_sigma = np.std(reading, ddof=1)

    #standard deviation of mean
    mean_sigma = sample_sigma/sqrt(num-1)

    #t-statistics
    k = t.interval(conf, num-1)[1]

    return k*mean_sigma


if __name__ == "__main__":
    """
        Testing scripts
    """

    #example time-series data
    data = [10, 20, 30, 20, 10]
    print("Data:")
    print(data)

    #0.5K absolute uncertainty and 0% relative uncertainty
    thermocouple = APPARATUS_UNCER(0.5, 0.0)

    print("Mean, Uncertainty, Zero-order, First-order:")
    print(thermocouple.measure_av_result(data, full_output=1))
