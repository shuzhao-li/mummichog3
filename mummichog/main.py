# Copyright (c) 2010-2020 Shuzhao Li.
# All rights reserved.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.


'''
mummichog -
pathway and network analysis for metabolomics

@author: Shuzhao Li
Online documentation: http://mummichog.org



## dev v3, overhaul 


'''

import json

from .functional_analysis import *
from .report.reporting import json_export_all


def main():
    
    print (fishlogo)
    print ( "mummichog version %s \n" %VERSION )
    optdict = dispatcher()

    print_and_loginfo("Started @ %s\n" %time.asctime())
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


    print_and_loginfo("\nFinished @ %s\n" %time.asctime())

    #
    #  This is to export data as Python objects
    #
    MCG_JSON = json_export_all(mixedNetwork, PA, MA, AN)

    #print(MCG_JSON)
    print("\n\n~~~~~~~~~~~~~~~~~~~~\n\n")

    # Use Encoder to convert Python to JSON format 
    s = json.JSONEncoder().encode(MCG_JSON )
    print(json.dumps(s,  indent=4))

  





#
# -----------------------------------------------------------------------------
#

if __name__ == '__main__':

    main()


