# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.

'''
v3:

use pandas dataframe for both user input data and metabolic models





'''
import time
import getopt
import base64
import logging
import sys

from io import BytesIO

import pandas as pd

# from .config import *

from .userData import *


def print_and_loginfo(s):
    '''
    Legacy function for logging. This function should retire soon.
    '''
    #print s
    logging.info(s)

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
               'instrument': 'unspecified',
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
        elif o in ("-u", "--instrument"): optdict['instrument'] = a
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
        -u, --instrument: Any integer, treated as ppm of instrument accuracy. Default is 10. 
              
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
                             "network=", "instrument=", "force_primary_ion",
                             "visualization=", "workdir=", "input=", 
                             "reference=", "infile=", "output=", "permutation="])
        if not opts:
            print (helpstr)
            sys.exit(2)
        
    except getopt.GetoptError as err:
        print (err)
        sys.exit(2)
    
    return cli_options(opts)
    




# metabolicNetwork

class DataMeetModel:
    '''
    This returns the tracking map btw massFeatures - EmpiricalCompounds - Compounds.

    Number of EmpiricalCompounds will be used to compute pathway enrichment, and control for module analysis.
    New in v2, but will move to a boutique database in v3.
    
    many to many matches:
    when a Compound matched to multiple MassFeatures, split by retention time to EmpiricalCompounds;
    when a Mass Feature matched to multiple Compounds, no need to do anything.
    
    ??Default primary ion is enforced, so that for an EmpiricalCompound, primary ion needs to exist before other ions.





    # Key change to make
    v3  to separate "annotation"

    Move indexing and query to pandas.dataframe

    Key output:
    empCpd2Features = {empCpd: (), ...,}
    empCpd2Cpds = {empCpd: (), ...,}




    '''
    def __init__(self, theoreticalModel, userData):
        '''
        # from ver 1 to ver 2, major change in .match()
        Trio structure of mapping
        (M.row_number, EmpiricalCompounds, Cpd)
        
        '''
        self.model = theoreticalModel
        self.data = userData
        
        # retention time window for grouping, based on fraction of time or ranks
        self.rtime_tolerance = self.data.max_retention_time * RETENTION_TIME_TOLERANCE_FRAC
        self.rtime_tolerance_rank = len(self.data.ListOfMassFeatures) * RETENTION_TIME_TOLERANCE_FRAC
        
        # major data structures
        # web
        if self.data.web:
            wanted_ions = self.data.paradict['wanted_adduct_list']
        # local
        else:
            wanted_ions = wanted_adduct_list[ self.data.paradict['mode'] ]

        self.IonCpdTree = self.__build_cpdindex__(wanted_ions)

        self.rowDict = self.__build_rowindex__( self.data.ListOfMassFeatures )
        self.ListOfEmpiricalCompounds = self.get_ListOfEmpiricalCompounds()
        
        # this is the reference list
        self.mzrows = [M.row_number for M in self.data.ListOfMassFeatures]
        
        self.rowindex_to_EmpiricalCompounds = self.__make_rowindex_to_EmpiricalCompounds__()
        self.Compounds_to_EmpiricalCompounds = self.__index_Compounds_to_EmpiricalCompounds__()
        
        # this is the sig list
        self.significant_features = self.data.input_featurelist
        self.TrioList = self.batch_rowindex_EmpCpd_Cpd( self.significant_features )

    def __build_cpdindex__(self, wanted_ions):
        '''
        indexed Compound list, to speed up m/z matching.
        Limited to MASS_RANGE (default 50 ~ 2000 dalton).
        
        changing from adduct_function to wanted_adduct_list dictionary
        
        wanted_adduct_list['pos_default'] = ['M[1+]', 'M+H[1+]', 'M+2H[2+]', 'M(C13)+H[1+]', 'M(C13)+2H[2+]', 
                    'M+Na[1+]', 'M+H+Na[2+]', 'M+HCOONa[1+]'
                    ],
        
        # 
        >>> metabolicModels['human_model_mfn']['Compounds'].items()[92]
        ('C00217', {'formula': '', 'mw': 147.0532, 'name': 'D-Glutamate; D-Glutamic acid; D-Glutaminic acid; D-2-Aminoglutaric acid',
         'adducts': {'M+2H[2+]': 74.53387646677, 'M+Br81[-]': 227.9695, 'M-H2O+H[1+]': 130.04987646677, 
         'M-C3H4O2+H[1+]': 76.03937646677, 'M-HCOOH+H[1+]': 102.05507646676999, 'M-HCOONa+H[1+]': 80.07307646677, 
         'M+K[1+]': 186.01597646677, 'M+Cl[-]': 182.0221, 'M+Na-2H[-]': 167.02064706646001, 'M-CO2+H[1+]': 104.07067646677, 
         'M+Na[1+]': 170.04247646677, 'M+Br[-]': 225.9715, 'M(S34)-H[-]': 148.04172353323, 'M+H[1+]': 148.06047646677, 
         'M-H4O2+H[1+]': 112.03927646677, 'M(C13)-H[-]': 147.04932353323, 'M(Cl37)-H[-]': 148.04312353323, 'M+HCOONa[1+]': 216.04787646677, 'M(C13)+2H[2+]': 75.03557646677, 'M+HCOOK[1+]': 232.02177646677, 'M-CO+H[1+]': 120.06547646677, 'M+HCOO[-]': 192.050845, 'M(C13)+3H[3+]': 50.359409800103336, 'M(Cl37)+H[1+]': 150.05767646677, 'M-H[-]': 146.04592353323, 'M+ACN-H[-]': 187.07246853323, 'M+Cl37[-]': 184.0191, 'M-H2O-H[-]': 128.03532353322998, 'M(S34)+H[1+]': 150.05627646677002, 'M-HCOOK+H[1+]': 64.09917646677, 'M+3H[3+]': 50.025009800103334, 'M+CH3COO[-]': 206.066495, 'M(C13)+H[1+]': 149.06387646677, 'M[1+]': 147.0532, 'M-NH3+H[1+]': 131.03397646677, 'M+NaCl[1+]': 206.01907646677, 'M+H+Na[2+]': 85.52487646677, 'M+H2O+H[1+]': 166.07107646677002, 'M-H+O[-]': 162.04083353323, 'M+K-2H[-]': 182.99414706646002, 'M-2H[2-]': 72.51932353323001}})
        >>> len(metabolicModels['human_model_mfn']['Compounds'])
        3560
        '''

        IonCpdTree = []
        
        for ii in range(MASS_RANGE[1]+1): 
            IonCpdTree.append([])       #empty lists for anything below MASS_RANGE
            
        # iteritems vs items is contention of efficiency, but there's change btw Python 2 and Python 3...
        for c,d in self.model.Compounds.items():
            if d['mw']:                 #sanity check; bypass mistake in adducts type
                for ion,mass in d['adducts'].items():
                    if ion in wanted_ions and MASS_RANGE[0] < mass < MASS_RANGE[1]:
                        IonCpdTree[ int(mass) ].append( (c, ion, mass) )
                
        # tree: (compoundID, ion, mass), ion=match form; mass is theoretical
        return IonCpdTree


    def __build_rowindex__(self, ListOfMassFeatures):
        '''
        Index list of MassFeatures by row# in input data
        '''
        rowDict = {}
        for M in ListOfMassFeatures: rowDict[M.row_number] = M
        return rowDict


    def __match_all_to_all__(self):
        '''
        Major change of data structure here in version 2.
        In ver 1, matched m/z is stored in each Compound instance.
        Here, we produce mapping dictionaries for
            * mzFeatures to theoretical ions
            * Compounds to mzFeatures
        Then, 
            * EmpiricalCompounds are determined within Compound matched mzFeatures, considering retention time.
        
        
        '''
        self.__match_to_mzFeatures__()
        self.cpd2mzFeatures = self.index_Compounds_to_mzFeatures()
        return self.compound_to_EmpiricalCompounds()
        

    def __match_to_mzFeatures__(self):
        '''
        Fill mzFeatures with matched ions and compounds
        '''
        for M in self.data.ListOfMassFeatures:
            M.matched_Ions = self.__match_mz_ion__(M.mz, self.IonCpdTree)
        
        
    def index_Compounds_to_mzFeatures(self):
        '''
        compound ID - mzFeatures
        run after self.__match_to_mzFeatures__()
        L: (compoundID, ion, mass)
        cpd2mzFeatures[compoundID] = [(ion, mass, mzFeature), ...]
        '''
        cpd2mzFeatures = {}
        for M in self.data.ListOfMassFeatures:
            for L in M.matched_Ions:
                if L[0] in cpd2mzFeatures:
                    cpd2mzFeatures[L[0]].append( (L[1], L[2], M) )
                else:
                    cpd2mzFeatures[L[0]] = [(L[1], L[2], M)]
        
        print ("Got %d cpd2mzFeatures" %len(cpd2mzFeatures))
        return cpd2mzFeatures
        
        
    def __match_mz_ion__(self, mz, IonCpdTree):
        '''
        L: (compoundID, ion, mass)
        return ions matched to m/z
        '''
        floor = int(mz)
        matched = []
        mztol = mz_tolerance(mz, self.data.paradict['instrument'])
        for ii in [floor-1, floor, floor+1]:
            for L in IonCpdTree[ii]:
                if abs(L[2]-mz) < mztol:
                    matched.append( L )
                    
        return matched

    def compound_to_EmpiricalCompounds(self):
        '''
        EmpiricalCompounds are constructed in this function.
        First splitting features matching to same compound by retention time;
        then merging those matched to same m/z features.
        run after self.index_Compounds_to_mzFeatures()
        '''
        ListOfEmpiricalCompounds = []
        for k,v in self.cpd2mzFeatures.items():
            ListOfEmpiricalCompounds += self.__split_Compound__(k, v)      # getting inital instances of EmpiricalCompound
            
        print ("Got %d ListOfEmpiricalCompounds" %len(ListOfEmpiricalCompounds))
        
        # merge compounds that are not distinguished by analytical platform, e.g. isobaric
        return self.__merge_EmpiricalCompounds__( ListOfEmpiricalCompounds )
        
    def __is_coelution__(self, massFeature1, massFeature2):
        '''
        True if retention times are within a tolerance window in time or ranks.
        Not assuming massFeatures are sorted in this function.
        '''
        if abs(massFeature1.retention_time - massFeature2.retention_time) < self.rtime_tolerance or \
            abs(massFeature1.retention_time_rank - massFeature2.retention_time_rank) < self.rtime_tolerance_rank:
            return True
        else:
            return False

    def __split_Compound__(self, compoundID, list_match_mzFeatures):
        '''
        Determine EmpiricalCompounds among the ions matched to a Compound;
        return list of EmpiricalCompounds (not final, but initiated here).
        
        The retention time is grouped by tolerance value; 
        This method should be updated in the future.
        
        input data format:
        cpd2mzFeatures[compoundID] = list_match_mzFeatures = [(ion, mass, mzFeature), ...]
        
        '''
        # unpacked format: [retention_time, row_number, ion, mass, compoundID]
        all_mzFeatures = [(L[2].retention_time, L[2].row_number, L[0], L[1], compoundID) for L in list_match_mzFeatures]
        all_mzFeatures.sort()
        ECompounds = []
        tmp = [ all_mzFeatures[0] ]
        for ii in range(len(all_mzFeatures)-1):

            if self.__is_coelution__( self.rowDict[all_mzFeatures[ii+1][1]], self.rowDict[all_mzFeatures[ii][1]] ):
                tmp.append(
                            all_mzFeatures[ii+1] )
            else:
                ECompounds.append( EmpiricalCompound( tmp ) )
                tmp = [ all_mzFeatures[ii+1] ]
        
        ECompounds.append( EmpiricalCompound( tmp ) )
        return ECompounds


    def __merge_EmpiricalCompounds__(self, ListOfEmpiricalCompounds):
        '''
        If ion/mzFeatures are the same, merge EmpiricalCompounds
        EmpiricalCompounds.join() adds Compounds
        
        Because EmpiricalCompounds.str_row_ion uses mzFeatures sorted by row_number, this is 
        '''
        mydict = {}
        for L in ListOfEmpiricalCompounds:
            if L.str_row_ion in mydict:
                mydict[ L.str_row_ion ].join(L)
            else:
                mydict[ L.str_row_ion ]= L
        
        print ("Got %d merged ListOfEmpiricalCompounds" %len(mydict))
        return mydict.values()

    def __make_rowindex_to_EmpiricalCompounds__(self):
        mydict = {}
        for E in self.ListOfEmpiricalCompounds:
            for m in E.massfeature_rows:
                if m in mydict:
                    mydict[m].append(E)
                else:
                    mydict[m] = [E]
                    
        return mydict

    def __index_Compounds_to_EmpiricalCompounds__(self):
        '''
        Make dict cpd - EmpiricalCompounds
        '''
        mydict = {}
        for E in self.ListOfEmpiricalCompounds:
            for m in E.compounds:
                if m in mydict:
                    mydict[m].append(E)
                else:
                    mydict[m] = [E]
                    
        return mydict
        

    def batch_rowindex_EmpCpd_Cpd(self, list_features):
        '''
        Batch matching from row feature to Ecpds; Use trio data structure, (M.row_number, EmpiricalCompounds, Cpd).
        Will be used to map for both sig list and permutation lists.
        '''
        new = []
        for f in list_features:
            for E in self.rowindex_to_EmpiricalCompounds.get(f, []):
                for cpd in E.compounds:
                    new.append((f, E, cpd))
            
        return new

            
    def get_ListOfEmpiricalCompounds(self):
        '''
        Collect EmpiricalCompounds.
        Initiate EmpCpd attributes here.
        '''
        ListOfEmpiricalCompounds, ii = [], 1
        for EmpCpd in self.__match_all_to_all__():
            EmpCpd.evaluate()
            EmpCpd.EID = 'E' + str(ii)
            EmpCpd.get_mzFeature_of_highest_statistic( self.rowDict )
            ii += 1
            if self.data.paradict['force_primary_ion']:
                if EmpCpd.primary_ion_present:
                    ListOfEmpiricalCompounds.append(EmpCpd)
            else:
                ListOfEmpiricalCompounds.append(EmpCpd)
        
        print ("Got %d final ListOfEmpiricalCompounds" %len(ListOfEmpiricalCompounds))
        return ListOfEmpiricalCompounds

    
    def to_json(self):
        '''
        JSON export to be consumed by downstream functions

        empCpd2Cpds = {empCpd: (), ...,}

        Will update later in accordance to 
        https://github.com/shuzhao-li/metDataModel

        '''

        empCpd2Features, empCpd2Cpds = {}, {}
        for E in self.ListOfEmpiricalCompounds:
            empCpd2Features[E.EID] = E.massfeature_rows
            empCpd2Cpds[E.EID] = E.compounds

        return {
            'metabolic_model': self.model.version,
            'empCpd2Features': empCpd2Features,
            'empCpd2Cpds': empCpd2Cpds,
        }
