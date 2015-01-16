#!/usr/bin/python

"""
    This file contains the script for the
	class structure containing the map
	training data and map parameters
"""

import pandas as pd

from calib_proc import MAP_PARA

class MAP_INFO:
    """
        This structure contains information related to
        parameters in a compressor map and pandas
        dataframe for training data
    """

    def __init__(
            self, name="", map_data=pd.DataFrame(),
            map_para = MAP_PARA()
        ):
        """
            Parameters:
            ===========
            num_data: int
                number of data points to train the map
        """
        # name of map
        self.name = name
        # pandas dataframe for training data
        self.map_data = map_data
        # MAP_PARA object for compressor map parameters
        self.map_para = map_para
