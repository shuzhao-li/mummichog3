'''
metabolic models from
JMS
or mummichogDB


from jms.modelConvert import convert_json_model
from jms.empiricalCpds import load_epds_from_json

'''




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
