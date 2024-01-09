#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 10:00:00 2024

@author: ghart

This script is designed to aggregate the outputs from the branching process when
simulations have been run on COMPS the experiment ID (and name for clearity)
should be changed for each COMPS experiment that is run with run_onCOMPS.py
"""

from idmtools.core.platform_factory import Platform
from idmtools.analysis.platform_anaylsis import PlatformAnalysis
import os
from io import BytesIO
from typing import Dict, Any, Union
from idmtools.entities.ianalyzer import IAnalyzer as BaseAnalyzer

import pandas as pd
from idmtools.entities.iworkflow_item import IWorkflowItem
from idmtools.entities.simulation import Simulation


# Experiment ID for the aggregation and name of process
experiment_ID = ['a8aefc66-95ae-ee11-aa0f-9440c9be2c51']
analysis_title = 'Aggregate big sweep - HIV branching process'


file_name = os.path.join('output','parameter-sweep-results.csv')
class NNDistAnalyzer(BaseAnalyzer):

    def __init__(self, title='idm'):
        super().__init__(filenames=[file_name], parse=False)
        print(title)

    def initialize(self):
        """
        Initialize our Analyzer. At the moment, this just creates our output folder
        Returns:

        """
        if not os.path.exists(os.path.join(self.working_dir, "output")):
            os.mkdir(os.path.join(self.working_dir, "output"))

    def map(self, data: Dict[str, Any], item: Union[IWorkflowItem, Simulation]) -> Any:
        """
        Extracts the Statistical Population, Data channel from InsetChart.
        Called for Each WorkItem/Simulation.

        Args:
            data: Data mapping str to content of file
            item: Item to Extract Data from (Usually a Simulation)

        Returns:

        """
    
        # Get clustering by attribute information
        dataframe = pd.read_csv(BytesIO(data[file_name]))
        dataframe['sim_num'] = item.tags['simNum']

        return dataframe

    def reduce(self, all_data: Dict[Union[IWorkflowItem, Simulation], Any]) -> Any:
        """
        Create the Final Population JSON and Plot

        Args:
            all_data: Populate data from all the Simulations

        Returns:
            None
        """

        output_dir = os.path.join(self.working_dir, "output")

        print('Creating output dataframe')
        out_data = pd.DataFrame()
        
        for key, data in all_data.items():
            out_data = pd.concat([out_data, data])
        # out_data = pd.concat(out_data, ignore_index=True)
        out_data.to_csv(os.path.join(output_dir, 'aggregateData.csv'), index=False)
        
        return


if __name__ == "__main__":
    docker_image = "docker-staging.packages.idmod.org/greg/idmtools_greg:1.0.0"
    platform = Platform('Calculon', docker_image=docker_image)
    analysis = PlatformAnalysis(
        platform=platform, experiment_ids=experiment_ID,
        analyzers=[NNDistAnalyzer], analyzers_args=[{'title': 'idm'}],
        analysis_name=analysis_title,
        # You can pass any additional arguments needed to AnalyzerManager through the extra_args parameter
        extra_args=dict(max_workers=8)
    )

    analysis.analyze(check_status=True)
    wi = analysis.get_work_item()
    print(wi)
