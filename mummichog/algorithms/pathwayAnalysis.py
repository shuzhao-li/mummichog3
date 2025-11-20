'''
>>> from _metabolic_models import *
>>> models.keys()
dict_keys(['Staphylococcus_epidermidis_ATCC_12228_massInferred', 'metabolicModel_EBI_20210602_Neisseria_meningitidis_alpha14', 'metabolicModel_AGORA_20210512_Escherichia_coli_O157_H7_str_Sakai', 'metabolicModel_EBI_20210602_Propionibacterium_acnes_SK137', 'metabolicModel_az_HumanGEM_20220302_noCompartmentalization', 'worm_model_icel1273'])

>>> for m in models:
...   print(models[m].keys())
... 
dict_keys(['id', 'version', 'Compounds', 'dict_cpds_def', 'metabolic_rxns', 'cpd_edges', 'edge2rxn', 'edge2enzyme', 'metabolic_pathways', 'cpd2pathways'])
dict_keys(['id', 'version', 'Compounds', 'dict_cpds_def', 'metabolic_rxns', 'cpd_edges', 'edge2rxn', 'edge2enzyme', 'metabolic_pathways', 'cpd2pathways'])
dict_keys(['id', 'version', 'Compounds', 'dict_cpds_def', 'metabolic_rxns', 'cpd_edges', 'edge2rxn', 'edge2enzyme', 'metabolic_pathways', 'cpd2pathways'])
dict_keys(['id', 'version', 'Compounds', 'dict_cpds_def', 'metabolic_rxns', 'cpd_edges', 'edge2rxn', 'edge2enzyme', 'metabolic_pathways', 'cpd2pathways'])
dict_keys(['id', 'version', 'Compounds', 'dict_cpds_def', 'metabolic_rxns', 'cpd_edges', 'edge2rxn', 'edge2enzyme', 'metabolic_pathways', 'cpd2pathways'])
dict_keys(['metabolic_rxns', 'cpd_edges', 'metabolic_pathways', 'Compounds', 'dict_cpds_def', 'cpd2pathways', 'edge2enzyme', 'edge2rxn', 'version', 'dict_cpds_mass'])
>>> 
>>> 
>>> 
>>> list(models['metabolicModel_EBI_20210602_Neisseria_meningitidis_alpha14']['cpd2pathways'].items())[:10]
[]
>>> list(models['metabolicModel_EBI_20210602_Neisseria_meningitidis_alpha14']['edge2enzyme'].items())[:10]
[]
>>> models['worm_model_icel1273']['metabolic_pathways'][:2]
[{'cpds': ['C00267', 'C00221', 'C00031'], 'rxns': ['RCC0138', 'RCC0137'], 'ecs': ['NA'], 'name': '', 'id': 'iCEL1273pathway1'}, {'cpds': ['C00025', 'C00024', 'C00148', 'C00026', 'M00016', 'C00007', 'C00006', 'C00005', 'C00004', 'C00003', 'C00027', 'C00001', 'C03912', 'C02946', 'M00148', 'C00033', 'M01165', 'C00437', 'C05947', 'C05946', 'M01352', 'M03287', 'C00334', 'C00014', 'C00010', 'M00077', 'C02714', 'M00001', 'C00011', 'M00003', 'M00002', 'M00005', 'M00004', 'M00006', 'M00009', 'M00008', 'C00134', 'C00077', 'M00025', 'M00026', 'M03912', 'M05947', 'M05946', 'M05938', 'M04281', 'M01157', 'C04281', 'M00080', 'C00080', 'C00042', 'C01157', 'C05936'], 'rxns': ['RC05050', 'RC01252', 'RC05052', 'RC01251', 'RC03293', 'RC01987', 'RMC0016', 'RC03291', 'RMC0015', 'RM05052', 'RM00239', 'RC00670', 'RM01248', 'RM00667', 'RC01248', 'RC04025', 'RC01154', 'RC00669', 'RM04444', 'RM04445', 'RM03293', 'RM03291', 'RM00245', 'RM05051', 'RM03313', 'RM01251', 'RM01253'], 'ecs': ['1.2.1.88', '1.5.99.8', '2.6.1.13', '1.14.11.2', '3.5.1.14', 'spontaneous', '1.5.1.2', '1.4.3.4', '1.2.1.41', '2.6.1.1', '4.1.1.17', '2.7.2.11', '1.5.1.12', '3.5.1.63', '1.5.-.-', '1.2.1.3', '3.5.1.16', '2.3.1.57'], 'name': 'Arginine and proline metabolism', 'id': 'iCEL1273pathway7'}]
>>> 
>>> models['metabolicModel_EBI_20210602_Neisseria_meningitidis_alpha14']['metabolic_pathways'][:2]
[]
>>> models['metabolicModel_EBI_20210602_Neisseria_meningitidis_alpha14']['metabolic_rxns'][:2]
[{'id': 'MNXR2184_i', 'reactants': ['bigg_fdp_i'], 'products': ['MNXM74_i', 'bigg_dhap_i']}, {'id': 'MNXR648_i', 'reactants': ['bigg_h_i', 'bigg_nadh_i', 'bigg_acald_i'], 'products': ['bigg_etoh_i', 'bigg_nad_i']}]
'''

from scipy import stats


class PathwayAnalysis:
    '''
    Pathway enrichment analysis, considering uncertainty in metabolite annotation.
    
    p-value is from Fisher exact test, 
    adjusted by resampling method in 
    GF Berriz, OD King, B Bryant, C Sander & FP Roth. 
    Characterizing gene sets with FuncAssociate. 
    Bioinformatics 19(18):2502-2504 (2003)
    
    '''
    def __init__(self, pathways, mixedNetwork):
        '''
        mixedNetwork contains both user input data, metabolic model,
        and mapping btw (mzFeature, EmpiricalCompound, cpd)
        
        '''
        self.mixedNetwork = mixedNetwork
        self.network = mixedNetwork.model.network
        self.paradict = mixedNetwork.data.paradict
        
        self.pathways = self.get_pathways(pathways)
        self.resultListOfPathways = []          # will store result of pathway analysis
        
        # to help track wehre sig cpd comes from
        self.TrioList = self.mixedNetwork.TrioList
        self.significant_EmpiricalCompounds = set([x[1] for x in self.TrioList])
        
        self.ListOfEmpiricalCompounds = mixedNetwork.ListOfEmpiricalCompounds
        self.total_number_EmpiricalCompounds = len(self.ListOfEmpiricalCompounds)

        print("\nPathway Analysis...")
        
        
    def get_pathways(self, pathways):
        '''
        convert pathways in JSON formats (import from .py) to list of Pathway class.
        Adding list of EmpiricalCompounds per pathway, which reflects the measured pathway coverage.
        '''
        new = []
        for j in pathways:
            P = metabolicPathway()
            P.json_import(j)
            P.adjusted_p = ''
            P.EmpiricalCompounds = self.__get_empiricalCompounds_by_cpds__(P.cpds)
            new.append(P)
        return new
        

    def __get_empiricalCompounds_by_cpds__(self, cpds):
        '''
        Mapping cpds to empirical_cpds. Also used for counting EmpCpds for each Pathway.
        '''
        cpds_empirical = []
        for c in cpds: cpds_empirical += self.mixedNetwork.Compounds_to_EmpiricalCompounds.get(c, [])
        return set(cpds_empirical)
        
        
    def do_permutations(self, pathways, num_perm):
        '''
        Modified from Berriz et al 2003 method.
        After collecting p-values from resampling, do a Gamma fit.
        
        Permutation is simplified in version 2; no more new TableFeatures instances.
        
        
        May consider fitting Gamma at log scale, to be more accurate --
        
        '''
        self.permutation_record = []
        print("Resampling, %d permutations to estimate background ..." 
                          %num_perm)
        
        # this is feature number not cpd number
        N = len(self.mixedNetwork.significant_features)
        for ii in range(num_perm):
            sys.stdout.write( ' ' + str(ii + 1))
            sys.stdout.flush()
            random_Trios = self.mixedNetwork.batch_rowindex_EmpCpd_Cpd( random.sample(self.mixedNetwork.mzrows, N) )
            query_EmpiricalCompounds = set([x[1] for x in random_Trios])
            self.permutation_record += (self.__calculate_p_ermutation_value__(query_EmpiricalCompounds, pathways))
        
        print("\nPathway background is estimated on %d random pathway values" 
                          %len(self.permutation_record))
        


    def __calculate_p_ermutation_value__(self, query_EmpiricalCompounds, pathways):
        '''
        calculate the FET p-value for all pathways.
        But not save anything to Pathway instances.
        '''
        p_of_pathways = [ ]
        query_set_size = len(query_EmpiricalCompounds)
        total_feature_num = self.total_number_EmpiricalCompounds
        
        for P in pathways:
            overlap_features = query_EmpiricalCompounds.intersection(P.EmpiricalCompounds)
            overlap_size = len(overlap_features)
            ecpd_num = len(P.EmpiricalCompounds)
            if overlap_size > 0:
                negneg = total_feature_num + overlap_size - ecpd_num - query_set_size
                p_val = stats.fisher_exact([[overlap_size, query_set_size - overlap_size],
                                       [ecpd_num - overlap_size, negneg]], 'greater')[1]
                p_of_pathways.append(p_val)
            else: 
                p_of_pathways.append(1)
                
        return p_of_pathways


    def get_adjust_p_by_permutations(self, pathways):
        '''
        EASE score is used as a basis for adjusted p-values,
        as mummichog encourages bias towards more hits/pathway.
        pathways were already updated by first round of Fisher exact test,
        to avoid redundant calculations.
        "Adjusted_p" is not an accurate term. It's rather a permutation based empirical p-value.
        '''
        self.do_permutations(pathways, self.paradict['permutation'])
        
        if self.paradict['modeling'] == 'gamma':
            #vector_to_fit = [-np.log10(x) for x in self.permutation_record if x < 1]
            vector_to_fit = -np.log10(np.array(self.permutation_record))
            self.gamma = stats.gamma.fit(vector_to_fit)
            a, loc, scale = self.gamma
            
            for P in pathways: 
                P.adjusted_p = self.__calculate_gamma_p__(a, loc, scale, P.p_EASE)
        else:
            for P in pathways: P.adjusted_p = self.__calculate_p__(P.p_EASE, self.permutation_record)
        return pathways
        

    def __calculate_p__(self, x, record):
        '''
        calculate p-value based on the rank in record of permutation p-values
        '''
        total_scores = [x] + record
        total_scores.sort()
        D = len(record) + 1.0
        return (total_scores.index(x)+1)/D
    
    def __calculate_gamma_p__(self, a, loc, scale, x):
        '''
        Use -log10 scale for model fitting
        '''
        return 1 - stats.gamma.cdf(-np.log10(x), a, loc, scale)
    
    
    def cpd_enrich_test(self):
        '''
        Fisher Exact Test in cpd space, after correction of detected cpds.
        Fisher exact test is using scipy.stats.fisher_exact
        for right-tail p-value:
        >>> stats.fisher_exact([[12, 5], [29, 2]], 'greater')[1]
        0.99452520602188932
        
        query size is now counted by EmpiricalCompounds.
        adjusted_p should be model p-value, not fdr.
        This returns a list of Pathway instances, with p-values.
        
                        P.p_EASE = stats.fisher_exact([[max(0, overlap_size - 1), query_set_size - overlap_size],
                                   [ecpd_num - overlap_size + 1, negneg]], 'greater')[1]
        '''
        FET_tested_pathways = []
        qset = self.significant_EmpiricalCompounds
        query_set_size = len(qset)
        total_feature_num = self.total_number_EmpiricalCompounds
        
        print("Query number of significant compounds = %d compounds" %query_set_size)
        
        for P in self.pathways:
            # use the measured pathway size
            P.overlap_EmpiricalCompounds = P.overlap_features = qset.intersection(P.EmpiricalCompounds)

            P.overlap_size = overlap_size = len(P.overlap_EmpiricalCompounds)
            P.EmpSize = ecpd_num = len(P.EmpiricalCompounds)
            if overlap_size > 0:
                negneg = total_feature_num + overlap_size - ecpd_num - query_set_size
                # Fisher's exact test
                P.p_FET = stats.fisher_exact([[overlap_size, query_set_size - overlap_size],
                                   [ecpd_num - overlap_size, negneg]], 'greater')[1]
                # EASE score as in Hosack et al 2003
                # taking out EASE, as the new approach of EmpiricalCompound is quite stringent already
                P.p_EASE = P.p_FET
                

            else:
                P.p_FET = P.p_EASE = 1
                
            FET_tested_pathways.append(P)
            #  (enrich_pvalue, overlap_size, overlap_features, P) 
            
        result = self.get_adjust_p_by_permutations(FET_tested_pathways)
        result.sort(key=lambda x: x.adjusted_p, reverse=False)
        self.resultListOfPathways = result

    
    def collect_hit_Trios(self):
        '''
        get [(mzFeature, EmpiricalCompound, cpd),...] for sig pathways.
        Nominate top cpd for EmpCpd here, i.e.
        in an EmpCpd promoted by a significant massFeature, the cpd candidate is chosen from a significant pathway.
        If more than one cpds are chosen, keep multiple.


        ?? where is this used now?

        
        '''
        overlap_EmpiricalCompounds = set([])
        for P in self.resultListOfPathways:
            if P.adjusted_p < SIGNIFICANCE_CUTOFF:
                # print(P.adjusted_p, P.name)
                overlap_EmpiricalCompounds = overlap_EmpiricalCompounds.union(P.overlap_EmpiricalCompounds)
        
        new = []
        for T in self.TrioList:
            # [(mzFeature, EmpiricalCompound, cpd),...]
            if T[1] in overlap_EmpiricalCompounds and T[0] in self.mixedNetwork.significant_features:
                # this does not apply to all sig EmpCpd
                T[1].update_chosen_cpds(T[2])
                T[1].designate_face_cpd()
                new.append(T)
        
        return new
                    
    def to_json(self):
        '''
        Convert result to dataframes, easy JSON export to be consumed by downstream functions

        dict from cpd to empCpd:
        self.Compounds_to_EmpiricalCompounds - needs to clean up for pathway specific
        # 'dicts_cpd2empCpd': self.Compounds_to_EmpiricalCompounds,

        dicts_cpd2empCpd = []
        for P in self.resultListOfPathways:
            for E in P.overlap_EmpiricalCompounds:
                dicts_cpd2empCpd.append ({ E.chosen_compounds: E.EID })

        '''
        self.collect_hit_Trios()    # force update of overlap_EmpiricalCompounds

        L = []
        for P in self.resultListOfPathways:
            L.append({
                'pathway_id': P.id,
                'name': P.name, #.encode("utf-8", "ignore"),   # Pathway name not ascii compliant
                'overlap_size': P.overlap_size,
                'pathway_size': P.EmpSize,
                'p-value': P.adjusted_p ,
                'significant_empCpds': [ E.EID for E in P.overlap_EmpiricalCompounds],
                'significant_compounds': [";".join(E.chosen_compounds) for E in P.overlap_EmpiricalCompounds],
            })
        return L
