import os
import logging

# to import from single place

SIGNIFICANCE_CUTOFF = 0.05
MASS_RANGE = (50, 2000)

# fraction of total retention time, or of ranks of retention time
# used to determine coelution of ions ad hoc
RETENTION_TIME_TOLERANCE_FRAC = 0.02    


# packages maintained separately
# problematic using `pip3 install git+git://github.com/shuzhao-li/mass2chem.git`
# use pypi instead for both
# python3 -m pip install metDataModel
from metDataModel.mummichog import *
from mass2chem.adducts import *

# temporary: getting JSON models
from .models import *


#
# this should be based on Experiment class, and organize by peaks - features - empCpds
#

class InputUserData:
    '''
    
    backward compatibility, 1 or 2-file input formats
    Per Joshua C., there'd be an option to test user designated L_sig, but user specified IDs are required
    
    return ListOfMassFeatures
    self.input_featurelist is "L_sig".


    Changing to pandas.dataframe



    '''
    
    def __init__(self, paradict, web=False):
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
        retention_times = [M.retention_time for M in self.ListOfMassFeatures]
        retention_times.sort()
        self.max_retention_time = retention_times[-1]
        self.max_mz = max([M.mz for M in self.ListOfMassFeatures])

        self.add_retention_time_rank(retention_times)
        self.determine_significant_list(self.ListOfMassFeatures)
        
        
    def add_retention_time_rank(self, retention_times):
        for M in self.ListOfMassFeatures:
            M.retention_time_rank = retention_times.index(M.retention_time)

    def text_to_ListOfMassFeatures(self, textValue, delimiter='\t'):
        '''
        Column order is hard coded for now, as mz, retention_time, p_value, statistic, CompoundID_from_user
        '''
        #
        lines = self.__check_redundant__( textValue.splitlines() )
        self.header_fields = lines[0].rstrip().split(delimiter)
        excluded_list = []
        for ii in range(len( lines )-1):
            y = lines[ii+1].split('\t')
            
            CompoundID_from_user = ''
            if len(y) > 4: CompoundID_from_user = y[4]
            [mz, retention_time, p_value, statistic] = [float(x) for x in y[:4]]
            
            # row_number, mz, retention_time, p_value, statistic, CompoundID_from_user
            if MASS_RANGE[0] < mz < MASS_RANGE[1]:
                # row # human-friendly, numbering from 1
                self.ListOfMassFeatures.append( 
                    MassFeature('row'+str(ii+1), mz, retention_time, p_value, statistic, CompoundID_from_user) 
                    )
            else:
                excluded_list.append( (ii, mz, retention_time) )
        
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
            new = sorted(all_feature_list, key=lambda x: x.p_value)
            
            p_hotspots = [ 0.2, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0001 ]
            N_hotspots = [ len([x for x in all_feature_list if x.p_value < pp]) for pp in p_hotspots ]
            
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
                self.paradict['cutoff'] = new[N_chosen+1].p_value
            else:
                #N_chosen = N_hotspots[chosen]
                
                self.paradict['cutoff'] = p_hotspots[chosen]
        
            print("Automatically choosing (p < %f) as significant cutoff."  %self.paradict['cutoff'])  
        
        # mark MassFeature significant
        for f in self.ListOfMassFeatures:
            if f.p_value < self.paradict['cutoff']:
                f.is_significant = True
        
        self.input_featurelist = [f.row_number for f in self.ListOfMassFeatures if f.is_significant]
        print("Using %d features (p < %f) as significant list." 
                              %(len(self.input_featurelist), self.paradict['cutoff']))  

