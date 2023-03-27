'''
from jms.modelConvert import convert_json_model
from jms.empiricalCpds import load_epds_from_json


# jms-metabolite-services
from jms.dbStructures import knownCompoundDatabase, ExperimentalEcpdDatabase

from .parameters import adduct_search_patterns, \
                                adduct_search_patterns_neg, \
                                    isotope_search_patterns, \
                                        extended_adducts


        EED = ExperimentalEcpdDatabase(mode=self.mode, mz_tolerance_ppm=self.mz_tolerance_ppm, rt_tolerance=2)
        # passing patterns from .defaul_parameters
        if self.mode == 'pos':
            EED.adduct_patterns = adduct_search_patterns
        else:
            EED.adduct_patterns = adduct_search_patterns_neg
        EED.isotope_search_patterns = isotope_search_patterns
        EED.extended_adducts = extended_adducts

        EED.build_from_list_peaks(self.CMAP.FeatureList)
        EED.extend_empCpd_annotation(self.KCD)
        EED.annotate_singletons(self.KCD)       
        # EED.dict_empCpds misses some features 
        EED.dict_empCpds = self.append_orphans_to_epmCpds(EED.dict_empCpds)

'''

import os
import logging

from metDataModel.mummichog import metabolicNetwork
from mass2chem.adducts import *
from .parameters import *

from .io.JSON_metabolicModels import metabolicModels  


def get_metabolic_model(model):
    '''
    Will change to getting models from JMS
    '''

    return metabolicNetwork(metabolicModels[ 'human_model_mfn' ])


#
# this should be based on Experiment class, and organize by peaks - features - empCpds
#

class InputUserData:
    '''
    
    Changing to JSON list of features and list of epds


    '''
    
    def __init__(self, paradict, web=False):
        '''
        
        '''
        self.web = web
        self.paradict = paradict
        self.header_fields = []
        self.ListOfMassFeatures = []
        self.input_featurelist = []

        # entry point of data input
        self.read()
        self.update()
        
    def update(self):
        '''
        Update retention_time_rank and is_significant to all MassFeatures
        '''
        retention_times = [M['rtime'] for M in self.ListOfMassFeatures]
        self.max_retention_time = max(retention_times)

        self.max_mz = max([M['mz'] for M in self.ListOfMassFeatures])

        self.determine_significant_list(self.ListOfMassFeatures)
        
        


    def text_to_ListOfMassFeatures(self, textValue, delimiter='\t'):
        '''
        Column order is hard coded for now, as mz, retention_time, p_value, statistic, CompoundID_from_user

        use asari style JSON features

        '''
        def _make_id(ii, mz, rt):
            return 'F' + str(ii) + '_' + str(round(mz, 6)) + '@' + str(round(rt, 2))
        #
        lines = self.__check_redundant__( textValue.splitlines() )
        self.header_fields = lines[0].rstrip().split(delimiter)

        excluded_list = []
        for ii in range(len( lines )-1):
            y = lines[ii+1].split('\t')
            
            fid_from_user = ''
            if len(y) > 4: fid_from_user = y[4].strip()

            [mz, rtime, p_value, statistic] = [float(x) for x in y[:4]]
            
            # row_number, mz, retention_time, p_value, statistic, CompoundID_from_user
            if MASS_RANGE[0] < mz < MASS_RANGE[1]:
                # row # human-friendly, numbering from 1
                fid = _make_id(ii+1, mz, rtime)
                peak = {'id_number': fid, 
                        'id': fid,
                        'fid_from_user': fid_from_user,
                        'mz': mz, 
                        'rtime': rtime,
                        'pval': p_value,
                        'statistic': statistic,
                        }
                self.ListOfMassFeatures.append( 
                    peak
                    )
            else:
                excluded_list.append( (ii, mz, rtime) )
        
        if excluded_list:
            print( "Excluding %d features out of m/z range %s." %(len(excluded_list), str(MASS_RANGE)) )

        
    def read_from_file(self, inputFile):
        return open(inputFile).read()
    
    def read_from_webform(self, t):
        return t

    def __check_redundant__(self, L):
        redundant = len(L) - len(set(L))
        if redundant > 0:
            print( "Your input file contains %d redundant features." %(redundant) )
        return L

    def read(self):
        '''
        Read input feature lists to ListOfMassFeatures. 
        Row_numbers (rowii+1) are used as primary ID.
        # not using readlines() to avoid problem in processing some Mac files
        '''
        if self.web:
            self.text_to_ListOfMassFeatures(self.paradict['datatext'])
        else:
            self.text_to_ListOfMassFeatures( 
                open(os.path.join(self.paradict['workdir'], self.paradict['infile'])).read() )

        print("Read %d features as reference list." %len(self.ListOfMassFeatures))
    
    
    # more work?
    def determine_significant_list(self, all_feature_list):
        '''
        For single input file format in ver 2. 
        The significant list, input_mzlist, should be a subset of ref_mzlist,
        determined either by user specificed --cutoff,
        or by automated cutoff close to a p-value hotspot,
        in which case, paradict['cutoff'] is updated accordingly.

        '''
        if not self.paradict['cutoff']:
            # automated cutoff
            new = sorted(all_feature_list, key=lambda x: x['pval'])
            
            p_hotspots = [ 0.2, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0001 ]
            N_hotspots = [ len([x for x in all_feature_list if x['pval'] < pp]) for pp in p_hotspots ]
            
            N_quantile = len(new) / 4
            N_optimum, N_minimum = 300, 30
            chosen = 9999
            for ii in range( len(N_hotspots) ):
                # will get the smallest p as ii increases
                if N_optimum < N_hotspots[ii] < N_quantile:
                    chosen = ii
            
            # if nothing was chosen
            if chosen > 100:
                for ii in range( len(N_hotspots) ):
                    if N_minimum < N_hotspots[ii] < N_quantile:
                        chosen = ii
            
            if chosen > 100:
                N_chosen = int(N_quantile)
                self.paradict['cutoff'] = new[N_chosen+1]['pval']
            else:
                #N_chosen = N_hotspots[chosen]
                
                self.paradict['cutoff'] = p_hotspots[chosen]
        
            print("Automatically choosing (p < %f) as significant cutoff."  %self.paradict['cutoff'])  
        
        # mark MassFeature significant
        for f in self.ListOfMassFeatures:
            if f['pval'] < self.paradict['cutoff']:
                f['is_significant'] = True
        
        self.input_featurelist = [f.row_number for f in self.ListOfMassFeatures if f['is_significant']]
        print("Using %d features (p < %f) as significant list." 
                              %(len(self.input_featurelist), self.paradict['cutoff']))  




# metabolicNetwork

class DataMeetModel:
    '''
    changing in v3

    TrioList is no longer relevant

    '''
    def __init__(self, metabolicModel, userData):
        '''
        # from ver 1 to ver 2, major change in .match()
        Trio structure of mapping
        (M.row_number, EmpiricalCompounds, Cpd)
        
        '''
        self.model = metabolicModel
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



        No more pre-computed adducts




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



        change to JMS/khipu based annotation









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
