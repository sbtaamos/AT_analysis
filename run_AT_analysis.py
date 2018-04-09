# import generic python modules

import argparse
import operator
from operator import itemgetter
import sys, os, shutil, itertools
import os.path
import datetime

version_nb = "0.1.0"
parser = argparse.ArgumentParser(prog='run_AT_analysis', usage='', add_help=False, formatter_class=argparse.RawDescriptionHelpFormatter, description=\
'''
**********************************************
v''' + version_nb + '''
author: Sarah-Beth Amos (sarah-beth.amos@bioch.ox.ac.uk)
git: 
**********************************************

	
[ DESCRIPTION ]

[ REQUIREMENTS ]
The following MD-specific python modules are needed :
 - MDAnalysis
 - MDTraj

 Other generic python modules are required. If dependencies are not met, the script will exit with an appropriate error message.

[ USAGE ]

Option	      Default  	Description                    
-----------------------------------------------------

-f			: structure file [.gro] (required)
-x			: trajectory file [.xtc] (required)
-o			: name of output folder

Other options
-----------------------------------------------------
--version		: show version number and exit
-h, --help		: show this menu and exit
  
''')


# options

#data options
parser.add_argument('-f', nargs=1, dest='grofilename', default=['no'], help=argparse.SUPPRESS, required=True)
parser.add_argument('-x', nargs=1, dest='xtcfilename', default=['no'], help=argparse.SUPPRESS, required=True)
#parser.add_argument('-s', nargs=1, dest='tprfilename', default=['no'], help=argparse.SUPPRESS, required=True)
parser.add_argument('-o', nargs=1, dest='output_folder', default=['no'], help=argparse.SUPPRESS, required=True)
parser.add_argument('--version', action='version', version='%(prog)s v' + version_nb, help=argparse.SUPPRESS)
parser.add_argument('-h','--help', action='help', help=argparse.SUPPRESS)


print "========================================"
print "Atomistic Simulation Analysis - S-B. T. A. Amos 2018"
print "========================================"
print "Running preparation:"


#=========================================================================================
# Store inputs
#=========================================================================================

print "Parsing..."
# parse user inputs
#------------------------
args = parser.parse_args()
args.grofilename = args.grofilename[0]
args.xtcfilename = args.xtcfilename[0]
args.output_folder = args.output_folder[0]
print "Parsed OK."
#=========================================================================================
# import modules 
#=========================================================================================
#generic science modules

print "Checking modules..."

try:
	import math
except:
	print "Error: you need to install the maths module."
	sys.exit(1)
try:
	import numpy as np
except:
	print "Error: you need to install the np module."
	sys.exit(1)
try:
	import scipy
except:
	print "Error: you need to install the scipy module."
	sys.exit(1)
try:
	import matplotlib as mpl
	mpl.use('Agg')
	import matplotlib.colors as mcolors
	mcolorconv = mcolors.ColorConverter()
	import matplotlib.cm as cm				#colours library
	import matplotlib.ticker
	from matplotlib.ticker import MaxNLocator, NullFormatter, ScalarFormatter
	from matplotlib.font_manager import FontProperties
	fontP=FontProperties()
except:
	print "Error: you need to install the matplotlib module."
	sys.exit(1)
try:
	import pylab as plt
except:
	print "Error: you need to install the pylab module."
	sys.exit(1)
#MDTraj module
try:
	import mdtraj as mdt
except:
	print "Error: you need to install the MDTraj molecule first. See http://mdtraj.org/1.9.0/"
	sys.exit(1)
#MDAnalysis module
try:
	import MDAnalysis
	from MDAnalysis import *
	import MDAnalysis.analysis
	import MDAnalysis.analysis.leaflet
	import MDAnalysis.analysis.distances
	#set MDAnalysis to use periodic boundary conditions
	MDAnalysis.core.flags['use_periodic_selections'] = True
	MDAnalysis.core.flags['use_KDTree_routines'] = False
except:
	print "Error: you need to install the MDAnalysis module first. See http://mdanalysis.googlecode.com"
	sys.exit(1)

print "Modules OK."
#=========================================================================================
# Checks
#=========================================================================================

print "Checking files..."
if not os.path.isfile(args.grofilename):
	print "Error: file " + str(args.grofilename) + "not found."
	sys.exit(1)

if not os.path.isfile(args.xtcfilename):
	print "Error: file " + str(args.xtcfilename) + "not found."
	sys.exit(1)
print "Files OK."
#=========================================================================================
# create folders and log file
#=========================================================================================

# create folders
print "Creating Folders..."
if os.path.isdir(args.output_folder):
	print "Error: folder " + str(args.output_folder) + "already exists. Choose a different output name using -o."

else:

	os.mkdir(args.output_folder)

	#Bilayer interactions - incl. penetration
	os.mkdir(args.output_folder + "/bilayer_interactions")
	os.mkdir(args.output_folder + "/bilayer_interactions/contacts")
	os.mkdir(args.output_folder + "/bilayer_interactions/penetration")
	# Protein characteristics - RMSD, RMSF, DSSP
	os.mkdir(args.output_folder + "/protein_properties")
	os.mkdir(args.output_folder + "/protein_properties/RMSD")
	os.mkdir(args.output_folder + "/protein_properties/RMSF")
	os.mkdir(args.output_folder + "/protein_properties/DSSP")
	os.mkdir(args.output_folder + "/protein_properties/ramachandran")
	# PCs
	os.mkdir(args.output_folder + "/PCA")
	# Lipid properties
	os.mkdir(args.output_folder + "/lipid_properties")
	os.mkdir(args.output_folder + "/lipid_properties/RDF")
	os.mkdir(args.output_folder + "/lipid_properties/density_maps")

	#create log
	#------------------------
	filename_log=os.getcwd() + '/' + str(args.output_folder) + '/AT_analysis.log'
	log_start_time = datetime.datetime.now()
	output_log = open(filename_log, 'w')
	output_log.write("[AT analysis v" + str(version_nb) + "]\n")
	output_log.write("Date/Time: " + str(log_start_time) + "\n")
	output_log.write("=========================================\n")
	output_log.write("\nThis folder and its content were created using the following command:\n\n")
	tmp_log = "python run_AT_analysis.py"
	for c in sys.argv[1:]:
		tmp_log+=" " + c
	output_log.write(tmp_log + "\n")
	output_log.close()

print "Folders OK."
###############################################################################
# FUNCTIONS
###############################################################################

def load_MDA_universe():
	global U
	print "\nLoading trajectory..."
	U = Universe(args.grofilename, args.xtcfilename)
	### might need some more globals here
	test_prot = U.select_atoms("protein")
	if test_prot.numberOfAtoms() == 0:
		print "Error: no protein found in the system."
		sys.exit(1)
	else:
		print "MDA Trajectory loaded successfully."


def calculate_BilayerPenetration():
	OUTPUT_BP = open("(args.output_folder + "/bilayer_interactions/penetration/bilayer_penetration.dat", 'w')")
	CA_res = U.select_atoms("protein and name CA")
	bilayer = U.select_atoms("resname POPC")    
	for residue in CA_res:
    	penetration_list = []
    	for frame in u.trajectory:
        	bilayer_z = bilayer.center_of_mass()[2]
        	CA_res_z = residue.position[2]
        	penetration = (CA_res_z - bilayer_z)
        	penetration_list.append(penetration)
    	np_penetration_list = np.asarray(penetration_list)
    	resid = np.asarray(residue.resid)
    	avg_penetration = np.mean(np_penetration_list)
    	std_penetration = np.std(np_penetration_list)
    	DAT = np.column_stack((resid, avg_penetration, std_penetration))
    	np.savetxt(OUTPUT_BP, (DAT), delimiter="   ", fmt="%s")
    OUTPUT_BP.close()

def calculate_LipidContacts():
	OUTPUT_contacts = open("(args.output_folder + "/bilayer_interactions/contacts/contacts_per_frame.dat", 'w')")
	for residue in CA_res:
    	lipid_list = []
		resi = residue.resid
		for frame in u.trajectory:
			lipids = u.select_atoms("resname PI4P and around 4 resid %i"%(resi))
			lipidNumber = lipids.n_residues
			lipid_list.append(lipidNumber)
		np_lipid_list = np.asarray(lipid_list)
		DAT = np.column_stack(np_lipid_list)
		np.savetxt(OUTPUT_contacts, (DAT), delimiter=" ", fmt="%s")
    OUTPUT_contacts.close()

# load files in MDTraj

print "Loading files (this may take a while)..."

#traj_mdt = mdt.load(args.xtcfilename, top=args.grofilename)
#print "Simulation length (ps):"
#print traj_mdt.time[-1]
#print "Done."


# Select proteins and lipids
#mdt_topology = traj_mdt.topology

#protein_mdt = topology.select("protein")
#CA_protein_mdt = topology.select("protein and name CA")
#residue_list_mdt = list(enumerate(CA_protein_mdt,1))
#lipids_set = ["POPC", "POPS", "PI4P", "PIP2"]

#bilayer = topology.select

# load files in MDAnalysis

load_MDA_universe()

# Bilayer interactions

print "Analysing bilayer interactions..."

calculate_BilayerPenetration()


print "Done."

# Protein properties
#print "Analysing protein properties..."
#print "Done."


# PCA
#print "Analysing principal components..."
#print "Done."
print "========================================"
print "Analysis complete! Check output folder for results."