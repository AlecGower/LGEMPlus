import sys
sys.path.append("/usr/src/logmet")

from lgemcore.io import logical_model_from_sbml

logical_model_from_sbml('model-files/yeastGEM.xml')
logical_model_from_sbml('model-files/iMM904.xml')

print("Both models loaded successfully.")