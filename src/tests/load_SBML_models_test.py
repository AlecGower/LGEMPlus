import sys
sys.path.append("/usr/src/logmet")

from genesis.io import logical_model_from_sbml

logical_model_from_sbml('ModelFiles/yeastGEM.xml')
logical_model_from_sbml('ModelFiles/iMM904.xml')

print("Both models loaded successfully.")