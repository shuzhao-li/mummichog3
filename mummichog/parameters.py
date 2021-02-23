'''

not used, but will

Default templates of parameters



'''

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