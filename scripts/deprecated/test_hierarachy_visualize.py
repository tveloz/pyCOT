import os  # Imports os to interact with the operating system
import sys # Imports sys to manipulate interpreter parameters and functions


import warnings  # Imports warnings to manage code warnings
warnings.filterwarnings("ignore", category=UserWarning, message="Pick support for FancyArrowPatch is missing.")

# Add the root directory to PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pyCOT.Hierarchy_visualize import * 
from pyCOT.RN_visualize import *
from pyCOT.reaction_network import * 
from pyCOT.closure_structure import * 
from pyCOT.file_manipulation import * 

#####################################################################################
# Call the function with the input data
# file_path = 'Txt/autopoietic.txt'
# file_path = 'Txt/2019fig1.txt'
# file_path = 'Txt/2019fig2.txt'
# file_path = 'Txt/non_connected_example.txt'
file_path = 'Txt/Farm.txt'  
file_path = r'C:\Users\tvelo\Dropbox\Public\AcademicWork\Europe\CLEA\Postdocs\TempletonPostdoc\sftw\pyCOT\networks\testing\Stentor_1.txt'

testRN = load_pyCOT_from_file(file_path)

reactive_sorgs = reactive_semi_orgs(testRN)
print("Semi-organisations:", reactive_sorgs)

input_data = convert_to_input_data(reactive_sorgs)
print("Label, Semi-Organisation, level:", input_data)

# Execute Hierarchy by Levels Code
Hierarchy_visualize(input_data)
