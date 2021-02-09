import logging
import sys
import random

from scipy import stats

from metDataModel.mummichog import *

USE_DEBUG = False

SEARCH_STEPS = 4
MODULE_SIZE_LIMIT = 100
SIGNIFICANCE_CUTOFF = 0.05
MASS_RANGE = (50, 2000)