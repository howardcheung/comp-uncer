#!/usr/bin/python
"""
    This file contains scripts that create coefficients
    of the 10-coefficient compressor map and the related
    paramters for the calculation of model uncertainty
"""

import copy
from math import sqrt

import numpy as np
from scipy.stats import t


class MAP_PARA:
    """
        This structure contains information related to
        a compressor map, including the coefficients
        and parameters to caculate the coefficients
    """

    def __init__(self, num_data=10):
        """
            Parameters:
            ===========
            num_data: int
                number of data points to train the map
        """
        # a numpy array of coefficients
        self._coeff = np.matrix(np.zeros([1, 10]))
        # matrix X in training data
        self._X = np.matrix(np.zeros([num_data, 10]))
        # uncertainty of matrix X
        self._uncer_X = np.matrix(np.zeros([num_data, 10]))
        # matrix y in training data
        self._y = np.matrix(np.zeros([num_data, 1]))
        # uncertainty of vector y
        self._uncer_y = np.matrix(np.zeros([num_data, 1]))
        self._sigma = 0.0  # standard deviation of y
        # quantity (X^T*X)^-1
        self._X_inverse_prod = np.matrix(np.zeros([10, 10]))
        # uncertainty of coefficients because of training data
        # a vector with number of entries equal to number of coefficients
        self._dBdydeltay = np.matrix(np.zeros([num_data, 10]))
        self._dBdXdeltaX = np.matrix(np.zeros([num_data*2, 10]))

    def get_coeff(self):
        return self._coeff

    def get_X(self):
        return self._X

    def get_uncer_X(self):
        return self._uncer_X

    def get_y(self):
        return self._y

    def get_uncer_y(self):
        return self._uncer_y

    def get_sigma(self):
        return self._sigma

    def get_dBdydeltay(self):
        return self._dBdydeltay

    def get_dBdXdeltaX(self):
        return self._dBdXdeltaX

    def get_X_inverse_prod(self):
        return self._X_inverse_prod

    def set_coeff(self, coeff):
        self._coeff = coeff

    def set_X(self, X):
        self._X = np.matrix(X)

    def set_uncer_X(self, uncer_X):
        self._uncer_X = np.matrix(uncer_X)

    def set_y(self, y):
        self._y = np.matrix(y)
        if self._y.shape[0] == 1:
            self._y = self._y.transpose()

    def set_uncer_y(self, uncer_y):
        self._uncer_y = np.matrix(uncer_y)
        if self._uncer_y.shape[0] == 1:
            self._uncer_y = self._uncer_y.transpose()

    def set_sigma(self, sigma):
        self._sigma = sigma

    def set_dBdydeltay(self, y):
        self._dBdydeltay = y

    def set_dBdXdeltaX(self, X):
        self._dBdXdeltaX = X

    def set_X_inverse_prod(self, X_inverse_prod):
        self._X_inverse_prod = X_inverse_prod


def set_regression_ind(
    CondTemp, EvapTemp, UncerCondTemp, UncerEvapTemp, para=MAP_PARA()
):
    """
        This function transforms the array of condensing temperature
        and evaporating temperature into the regression matrix X that
        trains an ARI 10-coefficient map. It also appends the uncertainty
        array into it

        Parameters:
        ===========
        CondTemp: list or numpy array
            Condensing temperature in F

        EvapTemp: list or numpy array
            Evaporating temperature in F

        UncerCondTemp: list or numpy array
            Uncertainty of condensing temperature in F

        UncerEvapTemp: list or numpy array
            Uncertainty of evaporating temperature in F

        para: MAP_PARA()
            structure of map parameters. defaulted to be an
            empty variable

        Returns:
        ===========
        para: MAP_PARA()
            structure of map parameters with X defined according
            to the input condensing temperature and evaporating
            temperature

    """

    num_data = len(CondTemp)
    X = []
    X.append(np.ones(num_data))
    X.append(EvapTemp)
    X.append(CondTemp)
    X.append([et**2 for et in EvapTemp])
    X.append([et*ct for et, ct in zip(EvapTemp, CondTemp)])
    X.append([ct**2 for ct in CondTemp])
    X.append([et**3 for et in EvapTemp])
    X.append([et**2*ct for et, ct in zip(EvapTemp, CondTemp)])
    X.append([et*ct**2 for et, ct in zip(EvapTemp, CondTemp)])
    X.append([ct**3 for ct in CondTemp])
    para.set_X(np.transpose(np.matrix(X)))

    # include uncertainty
    X = []
    X.append(np.zeros(num_data))
    X.append(UncerEvapTemp)
    X.append(UncerCondTemp)
    X.append([2.*et*uet for et, uet in zip(EvapTemp, UncerEvapTemp)])
    X.append([sqrt(
        (uet*ct)**2+(uct*et)**2
    ) for et, ct, uet, uct in zip(
        EvapTemp, CondTemp, UncerEvapTemp, UncerCondTemp
    )])
    X.append([2.*ct*uct for ct, uct in zip(CondTemp, UncerCondTemp)])
    X.append([3.*et**2*uet for et, uet in zip(EvapTemp, UncerEvapTemp)])
    X.append([sqrt(
        (2.*et*ct*uet)**2+(et**2*uct)**2
    ) for et, ct, uet, uct in zip(
        EvapTemp, CondTemp, UncerEvapTemp, UncerCondTemp
    )])
    X.append([sqrt(
        (2.*ct*et*uct)**2+(ct**2*uet)**2
    ) for et, ct, uet, uct in zip(
        EvapTemp, CondTemp, UncerEvapTemp, UncerCondTemp
    )])
    X.append([3.*ct**2*uct for ct, uct in zip(CondTemp, UncerCondTemp)])
    para.set_uncer_X(np.transpose(np.matrix(X)))

    return para


def set_regression_power(power, uncer_power, para=MAP_PARA()):
    """
        This function transforms the array of compressor power consumption
        into the regression matrix y that trains an ARI 10-coefficient map

        Parameters:
        ===========
        power: list or numpy array
            Compressor power consumption in W

        uncer_power: list or numpy array
            Uncertainty of compressor power consumption in W

        para: MAP_PARA()
            structure of map parameters. defaulted to be an
            empty variable

        Returns:
        ===========
        para: MAP_PARA()
            structure of map parameters with y defined according
            to the input power consumption

    """

    num_data = len(power)
    para.set_y(np.transpose(np.matrix(np.array(power))))
    para.set_uncer_y(np.transpose(
        np.matrix(np.array(uncer_power))
    ))

    return para


def inverse_X_prod(X):
    """
        Calculates and return (X^T*X)^-1 matrix

        Parameters:
        ===========
        X: numpy matrix
            matrix of training data

    """

    return (X.transpose()*X).getI()


def set_regression_coeff(para):
    """
        This function gives a 10-coefficient map based on the data
        stored in the input MAP_PARA()

        Parameters:
        ===========
        para: MAP_PARA()
            structure of map parameters. Should have X, y and
            the uncertainty defined beforehand

        Returns:
        ===========
        para: MAP_PARA()
            structure of map parameters containing the
            coefficients and sigma

    """

    # calculate coefficients
    para.set_X_inverse_prod(inverse_X_prod(para.get_X()))
    para.set_coeff(
        para.get_X_inverse_prod()*np.transpose(
            para.get_X()
        )*para.get_y()
    )

    # calculate standard deviation
    y_est = para.get_X()*para.get_coeff()
    y_diff = y_est-para.get_y()
    para.set_sigma(
        sqrt((
            np.multiply(y_diff, y_diff)
        ).sum()*1.0/(len(y_est)-10-1))
    )

    # calculate uncertainty progated to
    # coefficients from training data
    m = 10  # number of coefficients
    n = len(para.get_y())  # number of data pts
    dBdy = (para.get_X_inverse_prod()*np.transpose(
        para.get_X()
    )).transpose()  # a nxm matrix
    deltay = para.get_uncer_y()
    for ii in range(m):
        for jj in range(n):
            dBdy[jj, ii] = dBdy[jj, ii]*deltay[jj]
    para.set_dBdydeltay(dBdy)

    X_ori = para.get_X()
    coeff_ori = np.array(para.get_coeff().transpose().tolist()[0])
    dbdtdeltaet = []
    uncer_t = para.get_uncer_X()
    for ii in range(n):
        X_plus = copy.deepcopy(X_ori)
        X_plus[ii, 1] = X_plus[ii, 1]*1.001
        X_plus[ii, 3] = X_plus[ii, 3]*1.001**2
        X_plus[ii, 4] = X_plus[ii, 4]*1.001
        X_plus[ii, 6] = X_plus[ii, 6]*1.001**3
        X_plus[ii, 7] = X_plus[ii, 7]*1.001**2
        X_plus[ii, 8] = X_plus[ii, 8]*1.001
        beta_plus = np.array((
            inverse_X_prod(X_plus)*np.transpose(
                X_plus
            )*para.get_y()
        ).transpose().tolist()[0])
        dbdtdeltaet.append(
            (
                beta_plus-coeff_ori
            )/X_ori[ii, 1]/0.001*uncer_t[ii, 1]
        )
        del X_plus

    for ii in range(n):
        X_plus = copy.deepcopy(X_ori)
        X_plus[ii, 2] = X_plus[ii, 2]*1.001
        X_plus[ii, 4] = X_plus[ii, 4]*1.001
        X_plus[ii, 5] = X_plus[ii, 5]*1.001**2
        X_plus[ii, 7] = X_plus[ii, 7]*1.001
        X_plus[ii, 8] = X_plus[ii, 8]*1.001**2
        X_plus[ii, 9] = X_plus[ii, 9]*1.001**3
        beta_plus = np.array((
            inverse_X_prod(X_plus)*np.transpose(
                X_plus
            )*para.get_y()
        ).transpose().tolist()[0])
        dbdtdeltaet.append(
            (
                beta_plus-coeff_ori
            )/X_ori[ii, 2]/0.001*uncer_t[ii, 2]
        )
        del X_plus
    para.set_dBdXdeltaX(np.matrix(dbdtdeltaet))

    return para


def cal_regression_power(
    t_evap, t_cond, uncer_t_evap, uncer_t_cond,
    rel_uncer_power, abs_uncer_power,
    para, full_output=False
):
    """
        Estimate the compressor power and
        its uncertainty based on evaporating
        and condensing temperature

        Parameters:
        ===========
        t_evap: float
            Evaporating temperature in F

        t_cond: float
            Condensing temperature in F

        uncer_t_evap: float
            Uncertainty of evaporating temperature in F

        uncer_t_cond: float
            Uncertainty of condensing temperature in F

        rel_uncer_power: float
            Relative uncertainty of measured
            power consumption in %

        abs_uncer_power: float
            Absolute uncertainty of measured
            power consumption in W

        para: MAP_PARA() object
            Object containing the coefficients

        full_output: boolean
            Whether to output all other uncertainties.
            Default false.

        Returns:
        ===========
        power: float
            Estimated power in W

        uncer: float
            Uncertainty of the estimation in W

        uncer_input: float
            Uncertainty from inputs in W. Only output
            when full_output=True

        uncer_output: float
            Uncertainty from output in W. Only output
            when full_output=True

        uncer_train: float
            Uncertainty from training data in W. Only
            output when full_output=True

        uncer_dev: float
            Uncertainty from deviation in W. Only
            output when full_output=True

        uncer_cov: float
            Uncertainty from covariance in W. Only
            output when full_output=True

    """

    # form x vector
    coeff = para.get_coeff()
    x = (np.matrix([
        1.0, t_evap, t_cond, t_evap**2,
        t_evap*t_cond, t_cond**2, t_evap**3,
        t_evap**2*t_cond, t_evap*t_cond**2,
        t_cond**3
    ])).transpose()
    dyestdet = (np.matrix([
        0.0, 1.0, 0.0, 2.*t_evap,
        t_cond, 0.0, 3.*t_evap**2,
        2.*t_evap*t_cond, t_cond**2,
        0.0
    ]))*coeff
    dyestdct = (np.matrix([
        0.0, 0.0, 1.0, 0.0,
        t_evap, 2.*t_cond, 0.0,
        t_evap**2, 2.*t_cond*t_evap,
        3.*t_cond**2
    ]))*coeff

    # estimate power
    power = (x.transpose()*coeff)[0, 0]

    # estimate uncer_input
    uncer_input = sqrt((
        (dyestdet*uncer_t_evap).sum()
    )**2+(
        (dyestdct*uncer_t_cond).sum()
    )**2)

    # estimate uncer_output
    uncer_output = sqrt(
        abs_uncer_power**2 +
        (rel_uncer_power*power)**2
    )

    # estimate uncer_train
    train_x_entry = para.get_dBdXdeltaX()*x
    train_y_entry = para.get_dBdydeltay()*x
    uncer_train = sqrt((np.multiply(
        train_x_entry, train_x_entry
    )).sum()+(np.multiply(
        train_y_entry, train_y_entry
    )).sum())

    # estimate uncer_dev
    m = len(para.get_y())
    t_stat = t.interval(0.95, m-10-1)[1]
    uncer_dev = t_stat*para.get_sigma()

    # estimate uncer_cov
    uncer_cov = t_stat*sqrt(
        x.transpose() *
        para.get_X_inverse_prod() *
        x
    )*para.get_sigma()

    # estimate uncer
    uncer = sqrt(
        uncer_input**2+uncer_output**2 +
        uncer_train**2+uncer_dev**2 +
        uncer_cov**2
    )

    if full_output:
        return power, uncer, uncer_input, uncer_output, \
            uncer_train, uncer_dev, uncer_cov
    else:
        return power, uncer


if __name__ == '__main__':
    """
        Testing
    """

    #import libraries
    from math import sqrt
    import random

    from CoolProp.CoolProp import PropsSI as Props

    import data_manipulation
    import exp_uncer
    import misc_func
    import read_perform_data

    # define experimental method and apparatus
    test = data_manipulation.EXP_METHOD(0.1, 600)
    power_meter = exp_uncer.APPARATUS_UNCER(0.0, 0.005)
    p_trans = exp_uncer.APPARATUS_UNCER(0.0, 0.008)

    # read data
    filename = "..//Data//mfg_data_sheets//compressor//H23A463DBL_data.csv"
    print("Reading "+filename)
    df = read_perform_data.parse_one_sh_data(filename)

    # expected readings
    refri = 'R22'
    df['CondPInkPa'] = [Props(
        'P', 'T', misc_func.F2K(T), 'Q', 1, refri
    )/1000. for T in df.CondTempInF]
    df['EvapPInkPa'] = [Props(
        'P', 'T', misc_func.F2K(T), 'Q', 1, refri
    )/1000. for T in df.EvapTempInF]

    # defind space for uncertainty and new values
    random.seed(10)  # fix the random number generator seed for consistency
    df['MeaPowerInW'] = 0.0
    df['MeaCondPInkPa'] = 0.0
    df['MeaEvapPInkPa'] = 0.0
    df['MeaCondTempInF'] = 0.0
    df['MeaEvapTempInF'] = 0.0
    df['UncerMeaPowerInW'] = 0.0
    df['UncerMeaCondPInkPa'] = 0.0
    df['UncerMeaEvapPInkPa'] = 0.
    df['UncerMeaCondTempInF'] = 0.0
    df['UncerMeaEvapTempInF'] = 0.0

    # generate random time-series data from expected readings with exp_method
    # and apparatus information
    for ind in df.index:
        power_data = test.exp_data_generator(df.PowerInW[ind], power_meter)
        p_suc_data = test.exp_data_generator(
            df.EvapPInkPa[ind], p_trans, first_abs=0.9
        )
        p_dischg_data = test.exp_data_generator(
            df.CondPInkPa[ind], p_trans, first_abs=0.4
        )

        df.MeaPowerInW[ind], df.UncerMeaPowerInW[ind] = \
            power_meter.measure_av_result(power_data)
        df.MeaEvapPInkPa[ind], df.UncerMeaEvapPInkPa[ind] = \
            p_trans.measure_av_result(p_suc_data)
        df.MeaCondPInkPa[ind], df.UncerMeaCondPInkPa[ind] = \
            p_trans.measure_av_result(p_dischg_data)
        df.MeaEvapTempInF[ind], df.UncerMeaEvapTempInF[ind] = \
            data_manipulation.sat_temp_uncer_cal(
                df.MeaEvapPInkPa[ind], df.UncerMeaEvapPInkPa[ind], refri
            )
        df.MeaEvapTempInF[ind] = misc_func.K2F(df.MeaEvapTempInF[ind])
        df.UncerMeaEvapTempInF[ind] = misc_func.K2R(
            df.UncerMeaEvapTempInF[ind]
        )
        df.MeaCondTempInF[ind], df.UncerMeaCondTempInF[ind] = \
            data_manipulation.sat_temp_uncer_cal(
                df.MeaCondPInkPa[ind], df.UncerMeaCondPInkPa[ind], refri
            )
        df.MeaCondTempInF[ind] = misc_func.K2F(df.MeaCondTempInF[ind])
        df.UncerMeaCondTempInF[ind] = misc_func.K2R(
            df.UncerMeaCondTempInF[ind]
        )

    # create dataset for regression
    para = set_regression_ind(
        df.MeaCondTempInF, df.MeaEvapTempInF,
        df.UncerMeaCondTempInF, df.UncerMeaEvapTempInF
    )

    para = set_regression_power(
        df.MeaPowerInW, df.UncerMeaPowerInW, para
    )

    para = set_regression_coeff(para)

    # check power and uncertainty
    print(
        'Estimated power with 30F evaporating temperature' +
        ' and 100F condensing temperature:'
    )
    print(
        'Power [W], Uncertainty [W],' +
        ' Uncertainty from inputs [W], ' +
        ' Uncertainty from outputs [W], ' +
        ' Uncertainty from training data [W], ' +
        ' Uncertainty from deviation [W], ' +
        ' Uncertainty from covariance [W]:'
    )
    print(cal_regression_power(
        t_evap=30., t_cond=100.,
        uncer_t_evap=.3, uncer_t_cond=.3,
        rel_uncer_power=0.005, abs_uncer_power=0.0,
        para=para, full_output=True
    ))
    print(
        'Estimated power with 40F evaporating temperature' +
        ' and 120F condensing temperature:'
    )
    print(
        'Power [W], Uncertainty [W],' +
        ' Uncertainty from inputs [W], ' +
        ' Uncertainty from outputs [W], ' +
        ' Uncertainty from training data [W], ' +
        ' Uncertainty from deviation [W], ' +
        ' Uncertainty from covariance [W]:'
    )
    print(cal_regression_power(
        t_evap=40., t_cond=120.,
        uncer_t_evap=.3, uncer_t_cond=.3,
        rel_uncer_power=0.005, abs_uncer_power=0.0,
        para=para, full_output=True
    ))
    print(
        'Estimated power with -25F evaporating temperature' +
        ' and 170F condensing temperature:'
    )
    print(
        'Power [W], Uncertainty [W],' +
        ' Uncertainty from inputs [W], ' +
        ' Uncertainty from outputs [W], ' +
        ' Uncertainty from training data [W], ' +
        ' Uncertainty from deviation [W], ' +
        ' Uncertainty from covariance [W]:'
    )
    print(cal_regression_power(
        t_evap=-25., t_cond=170.,
        uncer_t_evap=.3, uncer_t_cond=.3,
        rel_uncer_power=0.005, abs_uncer_power=0.0,
        para=para, full_output=True
    ))
