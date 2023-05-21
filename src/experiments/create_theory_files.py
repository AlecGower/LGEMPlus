# Import packages
from genesis.io import logical_model_from_sbml

# Set parameters
model_xml = "ModelFiles/yeastGEM.xml"

GENE_ACTIVATIONS = True
SPLIT_ESSENTIAL = False
COMPARTMENTLESS = False

# Load model
lm = logical_model_from_sbml(model_file=model_xml)

# Load lists of essential and ubiquitous compounds
compounds = {"essential": [], "ubiquitous": [], "media": []}
media_name = ""

