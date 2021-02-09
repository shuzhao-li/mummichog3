# Licensed under the BSD 3-Clause License.
#
# mummichog - pathway and network analysis for metabolomics
# Online documentation: http://mummichog.org
#
#


## dev v3, overhaul 

VERSION = '3.0.4'
RELEASE = False

import json

import logging
import random
import itertools

#import pandas as pd
#from scipy import stats

from .io.models import *
from .io.get_user_data import *

from .algorithms.pathwayAnalysis import PathwayAnalysis
from .algorithms.modularAnalysis import ModularAnalysis
from .algorithms.activityNetwork import ActivityNetwork

from .report.reporting import json_export_all


fishlogo = '''     
    --------------------------------------------
    
             oO                      ooooooooo
           oOO   OOOOO  ooooo       ooo oooo
     oOO   O       ooooo  oooooo ooooo
    oooO           oooooo         oooo ooooo
        Oooo   o      OOOOOO   oooo   oooooooo
            ooooo  oooo      
                 o
    
    --------------------------------------------
    '''


def main():
    
    print (fishlogo)
    print ( "mummichog version %s \n" %VERSION )
    optdict = dispatcher()

    print("Started @ %s\n" %time.asctime())
    userData = InputUserData(optdict)
    
    #specify which metabolic model 
    if userData.paradict['network'] in ['human', 'hsa', 'Human', 'human_mfn', 'hsa_mfn', '']:
        theoreticalModel = metabolicNetwork(metabolicModels[ 'human_model_mfn' ])
    elif userData.paradict['network'] in ['worm', 'C. elegans', 'icel1273', 'Caenorhabditis elegans']:
        theoreticalModel = metabolicNetwork(metabolicModels[ 'worm_model_icel1273' ])
        
    else:
        raise KeyError( "Unsupported species/model. Pls contact author." )
    
    mixedNetwork = DataMeetModel(theoreticalModel, userData)

    # getting a list of Pathway instances, with p-values, in PA.resultListOfPathways
    PA = PathwayAnalysis(mixedNetwork.model.metabolic_pathways, mixedNetwork)
    PA.cpd_enrich_test()
    
    # Module analysis, getting a list of Mmodule instances
    MA = ModularAnalysis(mixedNetwork)
    MA.dispatch()
    
    # do activity network
    AN = ActivityNetwork( mixedNetwork, set(PA.collect_hit_Trios() + MA.collect_hit_Trios()) )


    print("\nFinished @ %s\n" %time.asctime())

    #
    #  This is to export data as Python objects
    #
    MCG_JSON = json_export_all(mixedNetwork, PA, MA, AN)

    #print(MCG_JSON)
    print("\n\n~~~~~~~~~~~~~~~~~~~~\n\n")

    # Use Encoder to convert Python to JSON format 
    s = json.JSONEncoder().encode(MCG_JSON )
    with open("mcg_output.json", "w") as O:
        O.write(s)

    #print(json.dumps(s,  indent=4) [:500])
    print("JSON output was written in mcg_output.json.")

  

#
# -----------------------------------------------------------------------------
#

if __name__ == '__main__':

    main()


