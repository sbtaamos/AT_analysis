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
-s			: topology file [.tpr] (required)
-o			: name of output folder (required)

Other options
-----------------------------------------------------
--version		: show version number and exit
-h, --help		: show this menu and exit
  
''')


# options

#data options
parser.add_argument('-f', nargs=1, dest='grofilename', default=['no'], help=argparse.SUPPRESS, required=True)
parser.add_argument('-x', nargs=1, dest='xtcfilename', default=['no'], help=argparse.SUPPRESS, required=True)
parser.add_argument('-s', nargs=1, dest='tprfilename', default=['no'], help=argparse.SUPPRESS, required=True)
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
try:
	import pandas as pd
except:
	print "Error: you need to install the pandas module."

try:
	import seaborn as sns
except:
	print "Error: you need to install the seaborn module."	
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
	import MDAnalysis.analysis.rms
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
	print "Error: file " + str(args.grofilename) + " not found."
	sys.exit(1)

if not os.path.isfile(args.xtcfilename):
	print "Error: file " + str(args.xtcfilename) + " not found."
	sys.exit(1)
print "Files OK."
#=========================================================================================
# create folders and log file
#=========================================================================================

# create folders
print "Creating Folders..."
if os.path.isdir(args.output_folder):
	print "Error: folder " + str(args.output_folder) + " already exists. Choose a different output name using -o."
	sys.exit(1)
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
	print "\nLoading MDA trajectory..."
	U = Universe(args.grofilename, args.xtcfilename)
	### might need some more globals here
	test_prot = U.select_atoms("protein")	
	if test_prot.atoms.n_atoms == 0:
		print "Error: no protein found in the system."
		sys.exit(1)
	else:
		print "MDA Trajectory loaded successfully."

def load_MDT_trajectory():
	print "\n Loading MDT trajectory..."
	traj_mdt = mdt.load(args.xtcfilename, top=args.grofilename)
	traj_mdt_protein = traj_mdt.topology.select("protein")
	print "Simulation length (ps):"
	print traj_mdt.time[-1]
	print "MDT Trajectory loaded successfully."

###############################################################################
#
# Analysis functions
# 
###############################################################################

def calculate_BilayerPenetration():
	# file for averaged data
	filename_details_1 = os.getcwd() + '/' + str(args.output_folder) + '/bilayer_interactions/penetration/bilayer_penetration.dat'
	OUTPUT_BP = open(filename_details_1,'w')
	# file for time-resolved data
	filename_details_2 = os.getcwd() + '/' + str(args.output_folder) + '/bilayer_interactions/penetration/bilayer_penetration_over_time.dat'
	OUTPUT_BP_TIME = open(filename_details_2,'w')
	CA_res = U.select_atoms("protein and name CA")
	bilayer = U.select_atoms("resname POPC")    
	for residue in CA_res:
    		penetration_list = []
    		for frame in U.trajectory:
        		bilayer_z = bilayer.center_of_mass()[2]
        		CA_res_z = residue.position[2]
        		penetration = (CA_res_z - bilayer_z)
        		penetration_list.append(penetration)
    		np_penetration_list = np.asarray(penetration_list)
    		resid = np.asarray(residue.resid)
    		avg_penetration = np.mean(np_penetration_list)
    		std_penetration = np.std(np_penetration_list)
    		DAT1 = np.column_stack((resid, avg_penetration, std_penetration))
    		np.savetxt(OUTPUT_BP, (DAT1), delimiter="   ", fmt="%s")
    		DAT2 = np_penetration_list
    		np.savetxt(OUTPUT_BP_TIME, (DAT2), delimiter=" ", fmt="%s")
    	OUTPUT_BP.close()
    	OUTPUT_BP_TIME.close()
    # plot the data
    # full (DAT1)


    # create filename
    #----------------
    filename_svg = os.getcwd() + '/' + str(args.output_folder) + '/bilayer_interactions/penetration/bilayer_penetration.svg'

    # create figure
    #--------------

    fig, ax = plt.subplots()
    fig.suptitle("Bilayer Penetration")


    # plot the data
    #--------------

    data = np.genfromtxt(filename_details_1)
    df = pd.DataFrame(data, columns=('residue','avg','std'))
	sns.barplot(x="residue", y="avg", error="std", data=df)
    fig.save(filename_svg)


def calculate_LipidContacts():
	filename_details_1 = os.getcwd() + '/' + str(args.output_folder) + '/bilayer_interactions/contacts/contacts_full.dat'
	OUTPUT_contacts = open(filename_details_1, 'w')
	filename_details_2 = os.getcwd() + '/' + str(args.output_folder) + '/bilayer_interactions/contacts/contacts_per_frame.dat'
	OUTPUT_contacts_TIME = open(filename_details_2, 'w')
	CA_res = U.select_atoms("protein and name CA")
	for residue in CA_res:
			lipid_list = []
			resi = residue.resid
			for frame in U.trajectory:
				lipids = U.select_atoms("resname PI4P and around 4 resid %i"%(resi))
				lipidNumber = lipids.n_residues
				lipid_list.append(lipidNumber)
			np_lipid_list = np.asarray(lipid_list)
			resid = np.asarray(residue.resid)
			SUM = np.sum(lipid_list)
			DAT1 = np.column_stack((resid, SUM))
			np.savetxt(OUTPUT_contacts, (DAT1), delimiter=" ", fmt="%s")
			DAT2 = np.column_stack(np_lipid_list)
			np.savetxt(OUTPUT_contacts_TIME, (DAT2), delimiter=" ", fmt="%s")
	OUTPUT_contacts.close()
	OUTPUT_contacts_TIME.close()


def calculate_RMSD():
	filename_details = os.getcwd() + '/' + str(args.output_folder) + '/protein_properties/RMSD/RMSD.dat'
	OUTPUT_RMSD = open(filename_details, 'w')
	ref = U.trajectory.rewind()
	R = MDAnalysis.analysis.rms.RMSD(U, ref, select="backbone", filename="RMSD_test.dat")
	R.run()
	R.save(OUTPUT_RMSD)
	print "RMSD calculation complete."
	OUTPUT_RMSD.close()

def calculate_RMSF():
	filename_details = os.getcwd() + '/' + str(args.output_folder) + '/protein_properties/RMSF/RMSF.dat'
	OUTPUT_RMSF = open(filename_details, 'w')
	u_rmsf = U.select_atoms("protein")
	protein = U.select_atoms("protein")
	R = MDAnalysis.analysis.rms.RMSF(u_rmsf, verbose=True)
	R.run()
	DAT = R.rmsf
	np.savetxt(OUTPUT_RMSF, (DAT), fmt="%s")
	OUTPUT_RMSF.close()
#	protein.bfactors = R.rmsf
#	U.trajectory[-1]
#	protein.write(os.getcwd() + '/' + str(args.output_folder) + '/protein_properties/RMSF/' + "protein_with_bfactors.pdb")
	print "RMSF calculation complete."

def calculate_DSSP():
	filename_details = os.getcwd() + '/' + str(args.output_folder) + '/protein_properties/DSSP/DSSP.dat'
	OUTPUT_DSSP = open(filename_details, 'w')
	dssp = mdt.compute_dssp(traj_mdt.atom_slice(atom_indices=traj.topology.select('protein')), simplified=True)
	dssp_over_time = dssp.T
	np.savetxt(OUTPUT_DSSP, (dssp_over_time), delimiter=" ")
	OUTPUT_DSSP.close()
	print "DSSP calculation complete"

def calculate_RDF():
	filename_details = os.getcwd() + '/' + str(args.output_folder) + '/lipid_properties/RDF/RDF.dat'
	OUTPUT_RDF = open(filename_details, 'w')
	u_rdf = Universe(args.tprfilename, args.xtcfilename, in_memory=True)
	rdf_prot = u_rdf.select_atoms("protein")
	rdf_lipids = u.rdf.select_atoms("resname PI4P")
	rdf = InderRDF(rdf_prot,rdf_lipids)
	rdf.run()
	rdf.save(OUTPUT_RDF)
	OUTPUT_RDF.close()
	print "RDF calculation complete."


###############################################################################

### END OF FUNCTIONS

###############################################################################


# load files in MDTraj

load_MDT_trajectory()

print "Loading files (this may take a while)..."

# load files in MDAnalysis

load_MDA_universe()

# Bilayer interactions

print "Analysing bilayer interactions..."

calculate_BilayerPenetration()

calculate_LipidContacts()

print "Done."

# Protein properties

print "Analysing protein properties..."

calculate_RMSD()

calculate_RMSF()

calculate_DSSP()

print "Done."

print "Analysing lipid properties..."

calculate_RDF()

print "========================================"
print "Analysis complete! Check output folder for results."
