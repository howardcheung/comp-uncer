#!/usr/bin/python

"""
    This file contains functions that read performance
    data from the csv file, return a raw data file of
    the performance data, and filter the data according
    to the user-defined test matrix
"""

import pandas as pd


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

    return df[df.SuperheatInF == sh_val]


if __name__ == '__main__':
    """
        for testing
    """

    filename = "..//Data//mfg_data_sheets//compressor//H23A463DBL_data.csv"
    print("Reading "+filename)
    df = parse_one_sh_data(filename)
    print(df)
