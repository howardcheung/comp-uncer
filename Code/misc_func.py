#!/usr/bin/python

"""
    This file consists of supportive functions that help
    calculation in all python scripts
"""

import numpy as np
import pylab as plt
import csv

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


def F2C(T_F):
    """
        Convert temperature in Fahrenheit to Celcius, code from ACHP

        Parameters:
        ===========
        T_F: float
            temperature in degrees Fahrenheit

        Returns:
        ===========
        T_C: float
            temperature in Celcius

    """
    return R2K(T_F-32.0)


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

        Parameters:
        ===========
        T_K: float
            temperature in Kelvin

        Returns:
        ===========
        T_R: float
            temperature in degrees Rankine

    """
    return T_K*9./5.

def parity_plot(
        x_data, uncer_x_data, x_label,
        y_data, uncer_y_data, y_label,
        filepath, showit = False
        ):
    """
        This function prints a parity plot with deviation
        lines and save it as pdf

        Parameters:
        ===========
        x_data: list or numpy array
            data to be plotted on x-axis

        uncer_x_data: list or numpy array
            error bar on x-axis

        x_label: string
            label text on x-axis

        y_data: list or numpy array
            data to be plotted on y-axis

        uncer_y_data: list or numpy array
            error bar on y-axis

        y_label: string
            label text on y-axis

        path: string
            path and filename where the figure is saved

    """

    fig = plt.figure(1)
    plt.scatter(
        x_data, y_data, s=30, c='b', marker='o'
    )
    plt.legend(loc=2)
    plt.errorbar(
        x_data, y_data, xerr=uncer_x_data,
        yerr=uncer_y_data, fmt='bo'
    )
    plt.xlabel(
        x_label, fontsize='x-large'
    )
    plt.ylabel(
        y_label,
        fontsize='x-large', multialignment='center'
    )
    plt.rc('xtick', labelsize='x-large')
    plt.rc('ytick', labelsize='x-large')
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.gcf().subplots_adjust(left=0.15)
    datamax = max([
        max(x_data), max(y_data),
    ])*1.1
    datamin = min([
        min(x_data), min(y_data),
    ])*0.9
    degreeline = [datamin, datamax]
    plt.plot(degreeline, degreeline, 'r--')
    a = max(abs(np.array(y_data)-np.array(x_data)))
    tenpercentline = [datamin-a, datamax-a]
    plt.plot(degreeline, tenpercentline, 'b:')
    s = '- %2.2f' % a
    plt.figtext(
        0.55, 0.45, s+' W',
        fontdict=None, size='x-large', color='k',
        horizontalalignment='left'
    )
    tenpercentline = [datamin+a, datamax+a]
    plt.plot(degreeline, tenpercentline, 'b:')
    s = '+ %2.2f' % a
    plt.figtext(
        0.5, 0.65, s+' W',
        fontdict=None, size='x-large', color='k',
        horizontalalignment='right')
    plt.axis([datamin, datamax, datamin, datamax])
    graph_filename = filepath
    plt.savefig(graph_filename, dpi=300)
    if not showit:
        plt.close(fig)
    else:
        plt.show()

def xy_plot(
        x_data, uncer_x_data, x_label,
        y_data, uncer_y_data, y_label,
        filepath, showit = False
        ):
    """
        This function prints a xy plot with and save it

        Parameters:
        ===========
        x_data: list or numpy array
            data to be plotted on x-axis

        uncer_x_data: list or numpy array
            error bar on x-axis

        x_label: string
            label text on x-axis

        y_data: list or numpy array
            data to be plotted on y-axis

        uncer_y_data: list or numpy array
            error bar on y-axis

        y_label: string
            label text on y-axis

        path: string
            path and filename where the figure is saved

    """

    fig = plt.figure(1)
    plt.scatter(
        x_data, y_data, s=30, c='b', marker='o'
    )
    plt.legend(loc=2)
    plt.errorbar(
        x_data, y_data, xerr=uncer_x_data,
        yerr=uncer_y_data, fmt='bo'
    )
    plt.xlabel(
        x_label, fontsize='x-large'
    )
    plt.ylabel(
        y_label,
        fontsize='x-large', multialignment='center'
    )
    plt.rc('xtick', labelsize='x-large')
    plt.rc('ytick', labelsize='x-large')
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.gcf().subplots_adjust(left=0.15)
    plt.axis([min(min(x_data),-25), max(max(x_data),60),
              min(min(y_data),75), max(max(y_data),155)])
    graph_filename = filepath
    plt.savefig(graph_filename, dpi=300)
    #save the data
    ofile = open(filepath[:-3]+"csv", 'wb')
    writersummary = csv.writer(ofile)
    writersummary.writerow(list(x_data))
    writersummary.writerow(list(y_data))
    ofile.close()
    if not showit:
        plt.close(fig)
    else:
        plt.show()

def plot_map(map_input,filepath,showit=True):
    #plots the data from a map as defined in main.py
    x_data = map_input.MeaEvapTempInF
    uncer_x_data = map_input.UncerMeaEvapTempInF
    x_label = "Evaporation Temperature [F]"
    y_data = map_input.MeaCondTempInF
    uncer_y_data = map_input.UncerMeaCondTempInF
    y_label = "Condensing Temperature [F]"
    xy_plot(x_data, uncer_x_data, x_label,
                y_data, uncer_y_data, y_label,
                filepath,showit=showit)
