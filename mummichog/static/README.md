# current data exchange formats in mummichog


## user input data, part 1 - parameters, defined in 

`io.userData.InputUserData.paradict`, 
which is updated from either cli_options or web form:

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

Not all parameters are valid. Working on cleanup. 
But a loosely defined dictionary is good for almost everything.

## user input data, part 2 - data matrix

People often use tabulated data, but more standardized formats should be used programmatically.

This can be either 
- list of Features
or
- list of EmpericalCompounds

if the former,
an annotation function is applied to convert to EmpericalCompounds.

The definition of 
    peak
    feature
    empirical compound
    experiment
    compound
    reaction
    pathway
    network
is given in the metDataModel repo (https://github.com/shuzhao-li/metDataModel),
and should be imported from there whenever feasible.

## Metabolic model

Defined by a list of reactions, and a mid-layer of organization (pathways, or other kind of groups).

Reactions are consisted of compounds (both reactants and products).

Compounds are looked up in a table for attributes (name, neutral_base_mass, formula, etc).

Pathways are defined by a set of reactions.

There's only one global network, connected by all reactions. 
Modules are sub-networks defined by a property, e.g. enrichment in an experiment.


## Mummichog output

including significant/enriched pathways, network modules and one activity network.

Will box in a new class mummichogResult.


- pathway result

refer to current pathway result table.

- network result, including module result and activity network

refer to current module result table.
refer to current activity network result table.

- mapping of empCpds to Cpds, with probablistic update

To be updated.
A list from cpd->empCpd, 
and updated empCpd annotation by mummichog.











