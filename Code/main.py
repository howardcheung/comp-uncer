#!/usr/bin/python

"""
    This file is the main function that carries out
    the regression of 6 different maps
"""

#import libraries
import copy
import csv
from math import sqrt
import pdb
import random

from CoolProp.CoolProp import PropsSI as Props
import numpy as np
import pylab as plt

import calib_proc
import data_manipulation
import exp_uncer
from info_class import MAP_INFO
import misc_func
from misc_func import OperatingPoint
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
random.seed(80)  # fix the random number generator seed for consistency
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
    power_data = test.exp_data_generator(df.PowerInW[ind], 
        power_meter, first_rel=0.03
    )
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

# generate new dataframes for all tested maps
# need to recode to read the conditions from files
# including testing cases and ARI cases
df_maps = {}
maps = ["Map 1", "Map 2", "Map 3", "Map 4", "Map 5", "Map 6", "Map 7", "Map 8"]
# map 1
df_maps[maps[0]] = read_perform_data.data_filter(
    df, CondTempRange=[80, 150], EvapTempRange=[5, 55],
    RemovalPoint=[
        OperatingPoint(110, 20), OperatingPoint(120, 25),
        OperatingPoint(100, 15), OperatingPoint(130, 15),
        OperatingPoint(140, 30), OperatingPoint(120, 10),
        OperatingPoint(90, 25), OperatingPoint(110, 25),
        OperatingPoint(120, 20), OperatingPoint(120, 30),
        OperatingPoint(140, 50), OperatingPoint(120, 50),
        OperatingPoint(100, 50), OperatingPoint(130, 40),
        OperatingPoint(110, 40), OperatingPoint(90, 40),
        OperatingPoint(100, 35)
    ],
)
# map 2
df_maps[maps[1]] = read_perform_data.data_filter(
    df, CondTempRange=[80, 150], EvapTempRange=[0, 50],
    RemovalPoint=[
        OperatingPoint(110, 20), OperatingPoint(120, 25),
        OperatingPoint(100, 15), OperatingPoint(130, 15),
        OperatingPoint(140, 30), OperatingPoint(120, 10),
        OperatingPoint(90, 25), OperatingPoint(110, 25),
        OperatingPoint(120, 20), OperatingPoint(120, 30),
        OperatingPoint(130, 5), OperatingPoint(100, 5),
        OperatingPoint(130, 40), OperatingPoint(110, 40),
        OperatingPoint(90, 40), OperatingPoint(100, 35)
    ],
)
# map 3
df_maps[maps[2]] = read_perform_data.data_filter(
    df, CondTempRange=[80, 140], EvapTempRange=[-10, 45],
    RemovalPoint=[
        OperatingPoint(110, 20), OperatingPoint(120, 25),
        OperatingPoint(110, 25), OperatingPoint(120, 20),
        OperatingPoint(120, 30), OperatingPoint(100, 10),
        OperatingPoint(100, 0), OperatingPoint(120, 0),
        OperatingPoint(120, 10), OperatingPoint(100, 20),
        OperatingPoint(100, 30)
    ],
)
# map 4
df_maps[maps[3]] = read_perform_data.data_filter(
    df, CondTempRange=[80, 130], EvapTempRange=[-15, 45],
    RemovalPoint=[
        OperatingPoint(110, 20), OperatingPoint(120, 25),
        OperatingPoint(110, 25), OperatingPoint(120, 20),
        OperatingPoint(120, 30)
    ],
)
# map 5
df_maps[maps[4]] = read_perform_data.data_filter(
    df, CondTempRange=[80, 130], EvapTempRange=[-20, 45],
    RemovalPoint=[
        OperatingPoint(110, 20), OperatingPoint(120, 25),
        OperatingPoint(100, 15), OperatingPoint(100, 10),
        OperatingPoint(110, 15), OperatingPoint(110, 25),
        OperatingPoint(120, 20), OperatingPoint(120, 30)
    ],
)
# map 6
df_maps[maps[5]] = read_perform_data.data_filter(
    df, CondTempRange=[130, 130], EvapTempRange=[45, 45],
    AddPoint=[
        OperatingPoint(130, 45), OperatingPoint(110, 45),
        OperatingPoint(100, 45), OperatingPoint(110, 30),
        OperatingPoint(90, 5), OperatingPoint(80, 45),
        OperatingPoint(90, 35), OperatingPoint(120, 45),
        OperatingPoint(80, 55), OperatingPoint(150, 55),
        OperatingPoint(150, 10), OperatingPoint(100, -20),
        OperatingPoint(80, -20), OperatingPoint(130, -5),
        OperatingPoint(130, 15), OperatingPoint(150, 25),
        OperatingPoint(80, 0), OperatingPoint(130, 25),
        OperatingPoint(110, -5), OperatingPoint(80, -10),
        OperatingPoint(110, 5), OperatingPoint(80, 10),
        OperatingPoint(90, 15), OperatingPoint(150, 40),
        OperatingPoint(140, 10), OperatingPoint(140, 20),
        OperatingPoint(140, 30), OperatingPoint(140, 50),
        OperatingPoint(120, -10), OperatingPoint(120, 0),
        OperatingPoint(120, 10), OperatingPoint(120, 35),
        OperatingPoint(120, 55), OperatingPoint(100, -10),
        OperatingPoint(100, 0), OperatingPoint(100, 20),
        OperatingPoint(100, 30), OperatingPoint(100, 55),
        OperatingPoint(90, -15), OperatingPoint(90, -5),
        OperatingPoint(90, 45), OperatingPoint(140, 0),
        OperatingPoint(140, 40), OperatingPoint(130, 5),
        OperatingPoint(130, 35), OperatingPoint(130, 55),
        OperatingPoint(120, 15), OperatingPoint(110, -15),
        OperatingPoint(110, 15), OperatingPoint(110, 50),
        OperatingPoint(100, 10), OperatingPoint(100, 40),
        OperatingPoint(90, 25), OperatingPoint(90, 55),
        OperatingPoint(80, 20), OperatingPoint(80, 30),
        OperatingPoint(80, 40), OperatingPoint(150, 20),
        OperatingPoint(150, 35), OperatingPoint(150, 50),
        OperatingPoint(110, 35), OperatingPoint(100, 35),
        OperatingPoint(100, 25), OperatingPoint(100, 15),
        OperatingPoint(130, 20), OperatingPoint(130, 30),
        OperatingPoint(110, 10), OperatingPoint(90, 50),
        OperatingPoint(120, 40), OperatingPoint(130, 50)
    ],
)
# testing points
df_maps["testing"] = read_perform_data.data_filter(
    df, CondTempRange=[80, 80], EvapTempRange=[-20, -20],
    AddPoint=[
        OperatingPoint(80, -20), OperatingPoint(100, -20),
        OperatingPoint(120, -10), OperatingPoint(150, 10),
        OperatingPoint(150, 55), OperatingPoint(120, 55),
        OperatingPoint(80, 55), OperatingPoint(110, 20),
        OperatingPoint(110, 25), OperatingPoint(120, 20),
        OperatingPoint(120, 30), OperatingPoint(80, 15),
        OperatingPoint(150, 30), OperatingPoint(120, 25)
    ],
)

# ARI points
df_maps["ARI"] = read_perform_data.data_filter(
    df, CondTempRange=[80, 80], EvapTempRange=[-20, -20],
    AddPoint=[
        OperatingPoint(130, 45), OperatingPoint(110, 45),
        OperatingPoint(100, 45), OperatingPoint(110, 30),
        OperatingPoint(90, 5), OperatingPoint(80, 45),
        OperatingPoint(90, 35), OperatingPoint(120, 45)
    ],
)

# All points
df_maps["All"] = read_perform_data.data_filter(
    df, CondTempRange=[80, 150], EvapTempRange=[-20, 55]
)

# use the same range as map 1 but much fewer data points to
# create map 7
df_maps[maps[6]] = read_perform_data.data_filter(
    df, CondTempRange=[80, 80], EvapTempRange=[5, 5],
    AddPoint=[
        OperatingPoint(100, 5), OperatingPoint(120, 5),
        OperatingPoint(140, 5), OperatingPoint(150, 20),
        OperatingPoint(130, 20), OperatingPoint(110, 20),
        OperatingPoint(90, 20), OperatingPoint(80, 35),
        OperatingPoint(100, 35), OperatingPoint(120, 35),
        OperatingPoint(140, 35), OperatingPoint(150, 55),
        OperatingPoint(130, 55), OperatingPoint(110, 55),
        OperatingPoint(90, 55),
        OperatingPoint(130, 45), OperatingPoint(110, 45),
        OperatingPoint(100, 45), OperatingPoint(110, 30),
        OperatingPoint(90, 5), OperatingPoint(80, 45),
        OperatingPoint(90, 35), OperatingPoint(120, 45)
    ],
)

# ARI map and with 1 dof
df_maps[maps[7]] = read_perform_data.data_filter(
    df, CondTempRange=[80, 80], EvapTempRange=[5, 5],
    AddPoint=[
        OperatingPoint(130, 55), OperatingPoint(110, 55),
        OperatingPoint(90, 55),
        OperatingPoint(130, 45), OperatingPoint(110, 45),
        OperatingPoint(100, 45), OperatingPoint(110, 30),
        OperatingPoint(90, 5), OperatingPoint(80, 45),
        OperatingPoint(90, 35), OperatingPoint(120, 45)
    ],
)

# regression for each map
map_infos = {}
for map in df_maps:
    if map != "ARI" and map != "testing" and map != "All":
        para = calib_proc.set_regression_ind(
            df_maps[map].MeaCondTempInF, df_maps[map].MeaEvapTempInF,
            df_maps[map].UncerMeaCondTempInF, df_maps[map].UncerMeaEvapTempInF
        )
        para = calib_proc.set_regression_power(
            df_maps[map].MeaPowerInW, df_maps[map].UncerMeaPowerInW, para
        )
        para = calib_proc.set_regression_coeff(para)
        map_infos[map] = MAP_INFO(
            name=map, map_data=df_maps[map], map_para=copy.deepcopy(para)
        )  # deepcopy to avoid assigning to the same pointer
        del para

# calculate and plot the parity plot for each map
result_dict = {}
for map in map_infos:
    power_pre = []
    uncer_power_pre = []
    df_temp = map_infos[map].map_data
    rel_uncer_power = np.mean(
        df_temp.UncerMeaPowerInW/df_temp.MeaPowerInW
    )

    # at training data points
    ss_tot_train = 0.0
    ss_res_train = 0.0
    for ind in map_infos[map].map_data.index:
        power, uncer_power = \
            calib_proc.cal_regression_power(
                t_evap=df_temp.MeaEvapTempInF[ind],
                t_cond=df_temp.MeaCondTempInF[ind],
                uncer_t_evap=df_temp.UncerMeaEvapTempInF[ind],
                uncer_t_cond=df_temp.UncerMeaCondTempInF[ind],
                rel_uncer_power=rel_uncer_power,
                abs_uncer_power=0.0,
                para=map_infos[map].map_para, full_output=False
            )
        power_pre.append(power)
        uncer_power_pre.append(uncer_power)
        ss_tot_train = ss_tot_train+(power-np.mean(
            df_temp.MeaPowerInW
        ))**2
        ss_res_train = ss_res_train+(power-df_temp.MeaPowerInW[ind])**2

    misc_func.parity_plot(
        df_temp.MeaPowerInW, df_temp.UncerMeaPowerInW,
        'Measured \nPower Consumption [W]',
        power_pre, uncer_power_pre,
        'Predicted \nPower Consumption [W]',
        '..//Results//'+map+'_scatter.pdf'
    )

    # at ARI and testing conditions
    cases = ["ARI", "testing", "All"]
    result = []
    result.append([
        "Map", "MeaPower", "UncerMeaPower", "CondTempInF", "EvapTempInF", "EstPower",
        "UncerOverall", "UncerInput", "UncerOutput",
        "UncerTrain", "UncerDev", "UncerCov", "PowerDiff"
    ])
    ss_tot_all = 0.0
    ss_res_all = 0.0
    for case in cases:
        power_pre = []
        uncer_power_pre = []
        for ind in df_maps[case].index:
            power, uncer_power, uncer_input, uncer_output, \
                uncer_train, uncer_dev, uncer_cov = \
                    calib_proc.cal_regression_power(
                        t_evap=df_maps[case].MeaEvapTempInF[ind],
                        t_cond=df_maps[case].MeaCondTempInF[ind],
                        uncer_t_evap=df_maps[case].UncerMeaEvapTempInF[ind],
                        uncer_t_cond=df_maps[case].UncerMeaCondTempInF[ind],
                        rel_uncer_power=rel_uncer_power,
                        abs_uncer_power=0.0,
                        para=map_infos[map].map_para, full_output=True
                    )
            power_pre.append(power)
            uncer_power_pre.append(uncer_power)

            # write results to file
            result.append([
                case, "%.5e"%df_maps[case].MeaPowerInW[ind],
                "%.5e"%df_maps[case].UncerMeaPowerInW[ind],
                "%.5e"%df_maps[case].MeaCondTempInF[ind],
                "%.5e"%df_maps[case].MeaEvapTempInF[ind],
                "%.5e"%power, "%.5e"%uncer_power,
                "%.5e"%uncer_input, "%.5e"%uncer_output,
                "%.5e"%uncer_train, "%.5e"%uncer_dev,
                "%.5e"%uncer_cov,
                "%.5e"%(power-df_maps[case].MeaPowerInW[ind])
            ])
            
            if case is "All":
                ss_tot_all = ss_tot_all+(power-np.mean(
                    df_maps[case].MeaPowerInW
                ))**2
                ss_res_all = ss_res_all+(
                    power-df_maps[case].MeaPowerInW[ind]
                )**2
                

        misc_func.parity_plot(
            df_maps[case].MeaPowerInW, df_maps[case].UncerMeaPowerInW,
            'Measured \nPower Consumption [W]',
            power_pre, uncer_power_pre,
            'Predicted \nPower Consumption [W]',
            '..//Results//'+map+'_'+case+'_scatter.pdf'
        )

    result_dict[map] = {
        'r2_train': 1-ss_res_train/ss_tot_train,
        'r2_all': 1-ss_res_all/ss_tot_all,
        'cov_train': map_infos[map].map_para.get_sigma()/np.mean(
            df_temp.MeaPowerInW
        ),
        'cov_all': sqrt(
            ss_res_all/(len(df_maps["All"].index)-10.)
        )/np.mean(
            df_maps["All"].MeaPowerInW
        )
    }

    ofile =  open('..//Results//'+map+'_all_result.csv', 'wb')
    writersummary = csv.writer(ofile)
    for row in result:
        writersummary.writerow(row)
    ofile.close()
    
    # detailed analysis of different uncertainty components
    case_points = {}
    case_points["ARI_points"] = [
        OperatingPoint(130, 45), OperatingPoint(110, 45),
        OperatingPoint(100, 45), OperatingPoint(110, 30),
        OperatingPoint(90, 5), OperatingPoint(80, 45),
        OperatingPoint(90, 35), OperatingPoint(120, 45)
    ]
    case_points["testing_points"] = [
        OperatingPoint(80, -20), OperatingPoint(100, -20),
        OperatingPoint(120, -10), OperatingPoint(150, 10),
        OperatingPoint(150, 55), OperatingPoint(120, 55),
        OperatingPoint(80, 55), OperatingPoint(110, 20),
        OperatingPoint(110, 25), OperatingPoint(120, 20),
        OperatingPoint(120, 30), OperatingPoint(80, 15),
        OperatingPoint(150, 30), OperatingPoint(120, 25)
    ]
    result = []
    result.append([
        "Map", "Case", "CondTempInF", "EvapTempInF", "EstPower",
        "UncerOverall", "UncerInput", "UncerOutput",
        "UncerTrain", "UncerDev", "UncerCov"
    ])
    for cases in case_points:
        for point in case_points[cases]:
            power, uncer_power, uncer_input, uncer_output, \
                uncer_train, uncer_dev, uncer_cov = \
                    calib_proc.cal_regression_power(
                        t_evap=point.get_EvapTempInF(),
                        t_cond=point.get_CondTempInF(),
                        uncer_t_evap=0.9,
                        uncer_t_cond=0.9,
                        rel_uncer_power=rel_uncer_power,
                        abs_uncer_power=0.0,
                        para=map_infos[map].map_para,
                        full_output=True
                    )
            result.append([
                map, cases, "%.5e"%(point.get_CondTempInF()),
                "%.5e"%(point.get_EvapTempInF()),
                "%.5e"%power, "%.5e"%uncer_power,
                "%.5e"%uncer_input, "%.5e"%uncer_output,
                "%.5e"%uncer_train, "%.5e"%uncer_dev,
                "%.5e"%uncer_cov
            ])
    ofile =  open('..//Results//'+map+'_result.csv', 'wb')
    writersummary = csv.writer(ofile)
    for row in result:
        writersummary.writerow(row)
    ofile.close()

# compare the accuracies of different maps
result = [[
    'Map', 'Mean square error in data', 'r2_train',
    'r2_all', 'cov_train', 'cov_all'
]]
coeff = []
for map in map_infos:
    result.append([
        map, "%.5e"%(map_infos[map].map_para.get_sigma()),
        "%.5e"%(result_dict[map]['r2_train']),
        "%.5e"%(result_dict[map]['r2_all']),
        "%.5e"%(result_dict[map]['cov_train']),
        "%.5e"%(result_dict[map]['cov_all']),
    ])
    coeff_entries = [map]
    for coeff_ind in map_infos[map].map_para.get_coeff():
        coeff_entries.append("%.5e"%coeff_ind)
    coeff.append(coeff_entries)

ofile =  open('..//Results//overall_result.csv', 'wb')
writersummary = csv.writer(ofile)
for row in result:
    writersummary.writerow(row)
ofile.close()

ofile =  open('..//Results//coeff.csv', 'wb')
writersummary = csv.writer(ofile)
for row in coeff:
    writersummary.writerow(row)
ofile.close()
    

