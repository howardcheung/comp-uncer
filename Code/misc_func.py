#!/usr/bin/python

"""
    This file consists of supportive functions that help
    calculation in all python scripts
"""

import numpy as np


class OperatingPoint:
    """
        This class is a class to store a compressor operating point
        defined by the condensing temperature and evaporating
        temperature
    """

    def __init__(self, CondTempInF=float('inf'), EvapTempInF=float('-inf')):
        self._condtempinF = CondTempInF
        self._evaptempinF = EvapTempInF

    def get_CondTempInF(self):
        return self._condtempinF

    def set_CondTempInF(self, CondTempInF):
        self._condtempinF = CondTempInF

    def get_EvapTempInF(self):
        return self._evaptempinF

    def set_EvapTempInF(self, EvapTempInF):
        self._evaptempinF = EvapTempInF

    def __eq__(self, SecondOP):
        if self._evaptempinF == SecondOP._evaptempinF and \
                self._condtempinF == SecondOP._condtempinF:
            return True
        else:
            return False

    def is_(self, SecondOP):
        return self.__eq__(SecondOP)

    def __ne__(self, SecondOP):
        if self == SecondOP:
            return False
        else:
            return True

    def is_not(self, SecondOP):
        return self.__ne__(SecondOP)

    def __str__(self):
        # function for printing
        return "(CondTempInF: "+str(self._condtempinF) +\
            ", EvapTempInF: "+str(self._evaptempinF)+")"


def check_size(entry):
    """
        This function returns a list or numpy array of an entry even when
        it is an integer or float

        Parameters:
        ===========
        x   : entry
            variable to be checked if it is a float, a list or a numpy array

        Returns:
        ==========
        ary : list of numpy array
            float that is converted to list or the original variable

    """
    if type(entry) is list or type(entry) is np.ndarray:
        ary = entry
    else:
        ary = []
        ary.append(entry)
    return ary


def central_difference(x, func, thres):
    """
        This function returns the finite difference estimate of
        the partial derivative of function func at x by a threshold thres

        Parameters:
        ===========
        x   : float
            variables that the derivative is taken with
        func: function
            function that the derivative is taken with
        thres: float
            increment that used to evaluate the finite difference

        Returns:
        ==========
            Derivative of func at x by central difference

    """

    return (func(x+thres)-func(x-thres))/(thres*2.)


def finite_difference(x, func, thres_initial, thres_check):
    """
        This function keeps decreasing the size of threshold for the
        central difference estimate of the partial derivative of
        func at x until the difference between the iterations is
        smaller than threshold thres_check.

        Parameters:
        ===========
        x   : float
            variables that the derivative is taken with
        func: function
            function that the derivative is taken with
        thres_initial: float
            intial increment that used to evaluate the finite difference
        thres_check: float
            threshold to be checked with the difference between
            derivatives

        Returns:
        ==========
            Derivative of func at x that is converged with
            decreasing size of increment

    """

    deriv = central_difference(x, func, thres_initial)
    if abs(deriv) < thres_check:
        deriv_base = thres_check
    else:
        deriv_base = deriv

    iter = 1
    iter_max = 20
    deriv_diff = float('inf')
    while abs(deriv_diff) > thres_check and iter < iter_max:
        deriv_next = central_difference(x, func, thres_initial/2.)
        deriv_diff = (deriv_next-deriv)/deriv_base
        deriv = deriv_next
        thres_initial = thres_initial/2.
        iter = iter+1

    return deriv


def F2K(T_F):
    """
        Convert temperature in Fahrenheit to Kelvin, code from ACHP

        Parameters:
        ===========
        T_F: float
            temperature in degrees Fahrenheit

        Returns:
        ===========
        T_K: float
            temperature in Kelvin

    """
    return R2K(T_F+459.67)


def R2K(T_R):
    """
        Convert temperature in Rankine to Kelvin, code from ACHP

        Parameters:
        ===========
        T_R: float
            temperature in degrees Rankine

        Returns:
        ===========
        T_K: float
            temperature in Kelvin

    """
    return 5./9.*(T_R)


def K2F(T_K):
    """
        Convert temperature in Kelvin to Fahrenheit , code from ACHP

        Returns:
        ===========
        T_K: float
            temperature in Kelvin

        Parameters:
        ===========
        T_F: float
            temperature in degrees Fahrenheit

    """
    return K2R(T_K)-459.67


def K2R(T_K):
    """
        Convert temperature in Kelvin to Rankine , code from ACHP

        Returns:
        ===========
        T_K: float
            temperature in Kelvin

        Parameters:
        ===========
        T_R: float
            temperature in degrees Rankine

    """
    return T_K*9./5.
