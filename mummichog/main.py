# Licensed under the BSD 3-Clause License.
#
# mummichog - pathway and network analysis for metabolomics
# Online documentation: http://mummichog.org
#
#

## dev v3, overhaul 


import json
from mummichog import __version__

from .userData import *
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

import time
import getopt
import sys

#
# Functions to take command line input
#

def cli_options(opts):
    '''
    Ongoing work in version 2, making some options obsolete.
    
    obsolete parameters:
    'analysis': 'total',
    'targeted': False,
    'evidence': 3,
    'visualization': 2,
    
    '''
    time_stamp = str(time.time())
    
    optdict = {
               'cutoff': 0,
               
               'network': 'human_mfn',
               'modeling': None,
               
               'mode': 'pos_default',
               'ppm': 10,
               'instrument': 10,            # same as ppm, phasing out
               'force_primary_ion': True,
               
               'workdir': '',
               'input': '',
               'reference': '',
               'infile': '',
               'output': '',
               'permutation': 100,
               'outdir': 'mcgresult' + time_stamp,
               }
    booleandict = {'T': True, 'F': False, 1: True, 0: False, 
                   'True': True, 'False': False, 'TRUE': True, 'FALSE': False, 'true': True, 'false': False,
                    }
    modedict = {'default': 'pos_default', 'pos': 'pos_default', 'pos_default': 'pos_default',
                'dpj': 'dpj_positive', 'positive': 'generic_positive', 'Positive': 'generic_positive',
                'negative': 'negative', 'Negative': 'negative',
                    }
    # update default from user argument
    for o, a in opts:
        if o in ("-a", "--analysis"): optdict['analysis'] = a
        elif o in ("-c", "--cutoff"): optdict['cutoff'] = float(a)
        elif o in ("-t", "--targeted"): optdict['targeted'] = booleandict.get(a, False)
        elif o in ("-n", "--network"): optdict['network'] = a
        elif o in ("-z", "--force_primary_ion"): optdict['force_primary_ion'] = booleandict.get(a, True)
        elif o in ("-d", "--modeling"): optdict['modeling'] = a
        elif o in ("-e", "--evidence"): optdict['evidence'] = int(a)
        elif o in ("-m", "--mode"): optdict['mode'] = modedict.get(a, a)
        # phasing out `instrument`
        elif o in ("-u", "--instrument"): optdict['ppm'] = a
        elif o in ("-u", "--ppm"): optdict['ppm'] = a
        elif o in ("-v", "--visualization"): optdict['visualization'] = int(a)
        elif o in ("-k", "--workdir"): optdict['workdir'] = a
        elif o in ("-i", "--input"): optdict['input'] = a
        elif o in ("-r", "--reference"): optdict['reference'] = a
        elif o in ("-f", "--infile"): optdict['infile'] = a
        elif o in ("-o", "--output"):
            optdict['output'] = a.replace('.csv', '')
            optdict['outdir'] = '.'.join([time_stamp, a.replace('.csv', '')])
            
        elif o in ("-p", "--permutation"): optdict['permutation'] = int(a)
        else: print ("Unsupported argument ", o)
    
    return optdict



def dispatcher():
    '''
    Dispatch command line arguments to corresponding functions.
    No user supplied id is used in version 1.
    User supplied IDs, str_mz_rtime IDs and targeted metabolites will be supported in version 2.
    

    '''
    helpstr = '''
    Usage example:
    python -m mummichog.main -f mydata.txt -o myoutput
    
        -f, --infile: single file as input, 
              containing all features with tab-delimited columns
              m/z, retention time, p-value, statistic score
        
        -n, --network: network model to use (default human_mfn; models being ported to version 2), 
              [human_mfn, worm]
        
        -o, --output: output file identification string (default 'mcgresult')
        -k, --workdir: directory for all data files.
              Default is current directory.
        
        -m, --mode: analytical mode of mass spec, [positive, negative, pos_defult].
              Default is pos_defult, a short version of positive.
        -u, --ppm: Any integer, treated as ppm of instrument accuracy. Default is 10. 
              
        -p, --permutation: number of permutation to estimate null distributions.
              Default is 100.
        -z,   --force_primary_ion: one of primary ions, 
              ['M+H[1+]', 'M+Na[1+]', 'M-H2O+H[1+]', 'M-H[-]', 'M-2H[2-]', 'M-H2O-H[-]'],  
              must be present for a predicted metabolite, [True, False].
              Default is True.
        
        -c, --cutoff: optional cutoff p-value in user supplied statistics,
              used to select significant list of features. 
        -d, --modeling: modeling permutation data, [no, gamma].
              Default is no.
        '''

    try:
        opts, args = getopt.getopt(sys.argv[1:], "a:c:t:d:e:m:n:u:z:v:k:i:r:f:o:p:", 
                            ["analysis=", "cutoff", "targeted=", "modeling=", "evidence=", "mode=", 
                             "network=", "ppm=", "force_primary_ion",
                             "visualization=", "workdir=", "input=", 
                             "reference=", "infile=", "output=", "permutation="])
        if not opts:
            print (helpstr)
            sys.exit(2)
        
    except getopt.GetoptError as err:
        print (err)
        sys.exit(2)
    
    return cli_options(opts)
    


def main():
    
    print (fishlogo)
    print ( "mummichog version %s \n" %__version__ )
    optdict = dispatcher()

    print("Started @ %s\n" %time.asctime())
    userData = InputUserData(optdict)
    
    theoreticalModel = get_metabolic_model( userData.paradict['network'] )
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


