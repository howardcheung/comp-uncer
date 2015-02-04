#!/usr/bin/python

"""
    This file contains a function to calculate the distance of
    a point to the nearest point in the training data set.
"""

import numpy as np
import pylab as plt


class point_distance:
    """
        This class stores condensing and evaporating temperatures and
        calculates the distance of a point to that data.
    """

    def __init__(self, CondTempInF=float('inf'), EvapTempInF=float('-inf'),
                 test_pt_CTF=float('-inf'), test_pt_ETF=float('inf')):

        self.condtempinF = CondTempInF
        self.evaptempinF = EvapTempInF
        self.test_pt_CTF = test_pt_CTF
        self.test_pt_ETF = test_pt_ETF

        if len(self.condtempinF) != len(self.evaptempinF):  #check if inputs OK
            raise()

    def get_min_DistanceF(self):
        #returns the distance of a point to the closest point in the
        #training data set
        distance_ctf = (self.test_pt_CTF - self.condtempinF)
        distance_etf = (self.test_pt_ETF - self.evaptempinF)
        min_distance = np.min(distance_ctf * distance_ctf +
                              distance_etf * distance_etf)
        min_distance = min_distance**0.5
        self.min_distance = min_distance

        return self.min_distance

    def get_max_DistanceF(self):
        #returns the distance of the point that is furthest away
        distance_ctf = (self.test_pt_CTF - self.condtempinF)
        distance_etf = (self.test_pt_ETF - self.evaptempinF)
        max_distance = np.max(distance_ctf * distance_ctf +
                              distance_etf * distance_etf)
        max_distance = max_distance**0.5
        self.max_distance = max_distance

        return self.max_distance

    def get_min_norm_DistanceF(self):
        #returns the minimum normalized distance
        distance_ctf = (self.test_pt_CTF - self.condtempinF)
        distance_etf = (self.test_pt_ETF - self.evaptempinF)
        max_c_distance = np.max(distance_ctf)
        max_e_distance = np.max(distance_etf)
        min_c_distance = np.min(distance_ctf)
        min_e_distance = np.min(distance_etf)
        #normalizing constants
        c_norm = max_c_distance - min_c_distance
        e_norm = max_e_distance - min_e_distance
        self.min_norm_distance = np.min((distance_ctf**2) / (c_norm**2) +
                                        (distance_etf**2) / (e_norm**2))
        self.min_norm_distance = np.sqrt(self.min_norm_distance)

        return self.min_norm_distance

if __name__ == "__main__":

    #generate some randome data
    np.random.seed(seed=20)  #set generator seed
    cond_temp = np.random.rand(5)*100.
    evap_temp = np.random.rand(5)*100.

    #some specific points for testing max/min distance
    cond_temp[0] = 100.
    evap_temp[0] = 10.
    cond_temp[1] = 0.
    evap_temp[1] = 100.

    #point that we want to have the distance to
    test_pt_cnd = 90.
    test_pt_evap = 20.

    #use our class to get the distance
    kwargs = {
        'CondTempInF': cond_temp,
        'EvapTempInF': evap_temp,
        'test_pt_CTF': test_pt_cnd,
        'test_pt_ETF': test_pt_evap
        }
    pt = point_distance(**kwargs)
    print pt.get_min_DistanceF(), pt.get_max_DistanceF(),\
          pt.get_min_norm_DistanceF()

    #plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.plot(cond_temp, evap_temp, 'x', markersize=10)
    plt.plot([test_pt_cnd], [test_pt_evap], 'd')
    plt.xlim([-1, 101.])
    plt.ylim([-1, 101.])
    ax.set_xlabel('Condensing Temperature')
    ax.set_ylabel('Evaporating Temperature')
    plt.show()
