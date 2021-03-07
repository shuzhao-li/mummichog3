'''
# metabolic models
# to do: models to be pulled from the Azimuth DB

'''

# to use 
# from metDataModel.derived import metabolicModel

from metDataModel.mummichog import metabolicNetwork

# from metDataModel.__metabolicModels__ import metabolicModels
from .JSON_metabolicModels import metabolicModels


def get_metabolic_model(model):
    '''
    #specify which metabolic model 
    if model in ['human', 'hsa', 'Human', 'human_mfn', 'hsa_mfn', '']:
        theoreticalModel = metabolicNetwork(metabolicModels[ 'human_model_mfn' ])
    elif model in ['worm', 'C. elegans', 'icel1273', 'Caenorhabditis elegans']:
        theoreticalModel = metabolicNetwork(metabolicModels[ 'worm_model_icel1273' ])
        
    else:
        raise KeyError( "Unsupported species/model. Pls contact author." )
    '''

    return metabolicNetwork(metabolicModels[ 'human_model_mfn' ])


def get_metabolic_model_azimuth(model):
    '''
    pull model from the Azimuth DB
    Doc on Google Firestore client: https://cloud.google.com/firestore/docs/quickstart-servers

    Used for Azimuth: set credentials $ export GOOGLE_APPLICATION_CREDENTIALS="/***/***.json"
    Need install google.cloud and likely more clients
    >>> from google.cloud import firestore
    >>> db = firestore.Client()

    Use Firestore API
    '''

    pass
