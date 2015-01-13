#!/usr/bin/python

"""
    This file contains functions that read performance
    data from the csv file, return a raw data file of
    the performance data, and filter the data according
    to the user-defined test matrix
"""

import copy

import pandas as pd

from misc_func import OperatingPoint


def parse_one_sh_data(filename):
    """
        Read the data from a file and return a pandas
        dataframe of performance data that has the same
        superheat. If there are multiple superheat
        in the data, the number of data points per superheat
        will be counted and only the data with the superheat
        corresponding to the largest number of data points
        will be returned

        Parameters:
        ===========
        filename: string
            name of file, or path of the CSV file that contains
            the data with the column EvapTempInF, CondTempInF,
            PowerInF, MassFlowRateInlbhr and SuperheatInF

        Returns:
        ===========
        df: pandas dataframe
            a pandas dataframe containing the performance data

    """

    # pass the data in the csv file into a pandas dataframe
    df = pd.read_csv(filename)

    # count the number of data points per superheat
    sh_counts = df.SuperheatInF.value_counts()  # dictionary
    num_val = 0
    sh_val = sh_counts.keys()[0]
    for sh, num_sh in zip(sh_counts.keys(), sh_counts.values):
        if num_sh > num_val:
            sh_val = sh

    # form a column of OperatingPoint
    df['OperatingPoint'] = OperatingPoint()  # define type
    for ind in df.index:
        df.OperatingPoint[ind] = OperatingPoint(
            df.CondTempInF[ind], df.EvapTempInF[ind]
        )  # relink to a new pointer

    return df[df.SuperheatInF == sh_val]


def data_filter(
        df, CondTempRange=[float('-inf'), float('inf')],
        EvapTempRange=[float('-inf'), float('inf')],
        RemovalPoint=[OperatingPoint()],
        AddPoint=[OperatingPoint()]
        ):
    """
        This file reads the raw file and return a pandas
        dataframe that only contains the data within the
        defined condensing temperature and evaporating
        temperature defined by the user. If the user defines
        extra operating points to be removed, the pandas
        dataframe will not contain the defined data point.

        Parameters:
        ============
        df: pandas dataframe
            performance data obtained from the function parse_one_sh_data

        CondTempRange: list
            range of condensing temperature that should be included
            in the new dataframe. So the resultant dataframe has
            condensing temperature >= CondTempRange[0] and condensing
            temperature <= CondTempRand[1].

        EvapTempRange: list
            range of evaporating temperature that should be included
            in the new dataframe So the resultant dataframe has
            evaporating temperature >= EvapTempRange[0] and evaporating
            temperature <= EvapTempRand[1].

        RemovalPoint: list
            list of OperatingPoint() that should be excluded in the new
            dataframe

        AddPoint: list
            list of OperatingPoint() that should be included in the new
            dataframe

        Returns:
        ============
        df_new: pandas dataframe
            performance data after filtering
    """

    # copy new dataframe
    df_new = copy.deepcopy(df)

    # condition list
    cond = []
    cond.append(df.CondTempInF >= CondTempRange[0])
    cond.append(df.CondTempInF <= CondTempRange[1])
    cond.append(df.EvapTempInF >= EvapTempRange[0])
    cond.append(df.EvapTempInF <= EvapTempRange[1])
    for point in RemovalPoint:
        cond.append(df.OperatingPoint != point)
    addcond = []
    for point in AddPoint:
        addcond.append(df.OperatingPoint == point)

    # Apply AND to all conditions
    final_condition = cond[0]
    for ii in xrange(1, len(cond)):
        final_condition = final_condition*cond[ii]

    # Apply OR to all conditions
    for ii in xrange(0, len(addcond)):
        final_condition = final_condition+addcond[ii]

    # Return the data that satisfy all conditions
    return df_new[final_condition]


if __name__ == '__main__':
    """
        for Testing
    """

    # read data
    filename = "..//Data//mfg_data_sheets//compressor//H23A463DBL_data.csv"
    print("Reading "+filename)
    df = parse_one_sh_data(filename)
    print(df)

    # filter data such that the new dataframe contains data with condensing
    # temperature within 90F and 130F, with evaporating temperature within -5F
    # and 10F and without the operating points (Condensing temp.,
    # Evaporating temp.) at (90F, -5F) and (100F, 0F) and with the operating
    # point at (80, 0)
    df_new = data_filter(
        df, CondTempRange=[90, 130],
        EvapTempRange=[-5, 10],
        RemovalPoint=[OperatingPoint(90, -5), OperatingPoint(100, 0)],
        AddPoint=[OperatingPoint(80, 0)]
    )
    print(df_new)
