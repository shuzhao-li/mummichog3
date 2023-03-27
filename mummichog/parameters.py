'''

not used, but will

Default templates of parameters



'''


SIGNIFICANCE_CUTOFF = 0.05
MASS_RANGE = (50, 2000)

# fraction of total retention time, or of ranks of retention time
# used to determine coelution of ions ad hoc
RETENTION_TIME_TOLERANCE_FRAC = 0.02    




user_parameters = {
               'network': 'human_mfn',
               'mode': 'pos_default',
               'instrument': 'unspecified',
}

algorithm_parameters = {
               'cutoff': 0,
               
               'modeling': None,
               

               'force_primary_ion': True,
               
               'workdir': '',
               'input': '',
               'reference': '',
               'infile': '',
               'output': '',
               'permutation': 100,
               'outdir': 'mcgresult',
}

paradict = {
    **user_parameters, **algorithm_parameters
}



#
# From khipu.utils. Placed here so that people can customize within asari.
#
PROTON = 1.00727646677
electron = 0.000549

# avoid confusing adducts in initial search, e.g. H, H2O
adduct_search_patterns = [  # initial patterns are relative to M+H+
                            (21.9820, 'Na/H'),
                            (41.026549, 'ACN'),     # Acetonitrile
                            (35.9767, 'HCl'),
                            (37.955882, 'K/H'),
                            ]

adduct_search_patterns_neg = [ (35.9767, 'HCl'), 
                            (46.00548, 'HCOOH'),
                            (21.9820, 'Na/H'), 
                            (41.026549, 'ACN'),
                            (37.955882, 'K/H'),
                            ]

isotope_search_patterns = [ (1.003355, '13C/12C', (0, 0.8)),
                            (2.00671, '13C/12C*2', (0, 0.8)),
                            # (3.010065, '13C/12C*3', (0, 0.8)),
                            # (4.01342, '13C/12C*4', (0, 0.8)),
                            # (5.016775, '13C/12C*5', (0, 0.8)),
                            # (6.02013, '13C/12C*6', (0, 0.8)),
                            # (7.023485, '13C/12C*7', (0, 0.8)),
                            # (8.02684, '13C/12C*8', (0, 0.8)),
                            # (9.030195, '13C/12C*9', (0, 0.8)),
                            # (10.03355, '13C/12C*10', (0, 0.8)),
                            # (11.036905, '13C/12C*11', (0, 0.8)),
                            # (12.04026, '13C/12C*12', (0, 0.8)),
                            ]

extended_adducts = [(1.0078, 'H'),
                            (-1.0078, '-H'),
                            (10.991, 'Na/H, double charged'),
                            (0.5017, '13C/12C, double charged'),
                            (117.02655, '-NH3'),
                            (17.02655, 'NH3'),
                            (-18.0106, '-H2O'),
                            (18.0106, 'H2O'),      # easy to confuse with bio reactions
                            (18.033823, 'NH4'),
                            (27.01089904, 'HCN'),
                            (37.94694, 'Ca/H2'),
                            (32.026215, 'MeOH'),
                            (43.96389, 'Na2/H2'),
                            (67.987424, 'NaCOOH'),
                            (83.961361, 'KCOOH'),
                            (97.96737927, 'H2SO4'),
                            (97.97689507, 'H3PO4'),
]


