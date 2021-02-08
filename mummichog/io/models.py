# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.

'''
Data models in mummichog

Import from serialized Python objects.
(A server version will use or a central database)

metabolicModels = {
    human_model_default: [
        ListOfCompounds: [Compound: {ID, name, mw, formula, AdductsAndDerivatives,
                                    }, ...],
        Pathway
        MetabolicNetwork
    
        ],
    
    ...,
}

@author: Shuzhao Li
'''

# the azimuth module will be replaced by its own separate package
from .JSON_metabolicModels import metabolicModels

# from .config import primary_ions, dict_weight_adduct



primary_ions = ['M+H[1+]', 'M+Na[1+]', 'M-H2O+H[1+]', 'M-H[-]', 'M-2H[2-]', 'M-H2O-H[-]']




# weighting function of isotopic derivatives and adducts
dict_weight_adduct = {
            'M[1+]': 5, 
            'M+H[1+]': 5,
            'M+2H[2+]': 3,
            'M+3H[3+]': 1,
            'M(C13)+H[1+]': 2,
            'M(C13)+2H[2+]': 1,
            'M(C13)+3H[3+]': 1,
            'M(S34)+H[1+]': 1,
            'M(Cl37)+H[1+]': 1,
            'M+Na[1+]': 3, 
            'M+H+Na[2+]': 2,
            'M+K[1+]': 2, 
            'M+H2O+H[1+]': 1, 
            'M-H2O+H[1+]': 1, 
            'M-H4O2+H[1+]': 1,
            'M-NH3+H[1+]': 1,
            'M-CO+H[1+]': 1,
            'M-CO2+H[1+]': 1,
            'M-HCOOH+H[1+]': 1,
            'M+HCOONa[1+]': 1,
            'M-HCOONa+H[1+]': 1,
            'M+NaCl[1+]': 1, 
            'M-C3H4O2+H[1+]': 1,
            'M+HCOOK[1+]': 1,
            'M-HCOOK+H[1+]': 1,
            # negative
            'M-H[-]': 5,
            'M-2H[2-]': 3,
            'M(C13)-H[-]': 2,
            'M(S34)-H[-]': 1,
            'M(Cl37)-H[-]': 1,
            'M+Na-2H[-]': 2,
            'M+K-2H[-]': 1,
            'M-H2O-H[-]': 1,
            'M+Cl[-]': 1,
            'M+Cl37[-]': 1,
            'M+Br[-]': 1,
            'M+Br81[-]': 1,
            'M+ACN-H[-]': 1,
            'M+HCOO[-]': 1,
            'M+CH3COO[-]': 1,
            'M-H+O[-]': 1,
                    }




