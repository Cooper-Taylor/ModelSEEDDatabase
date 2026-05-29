#!/usr/bin/env python
import os,sys
sys.path.append('../../Libs/Python/')
from BiochemPy import Compounds

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()
structures_dict = compounds_helper.loadStructures(["SMILE","InChIKey"],["ModelSEED"])

############################################################################
##
## First we apply Group Contribution Energies
##
############################################################################
#In the case where there was originally conflicting structures, we only want
#The energy for the structure that was curated, and we're using the provenance here
#To delete the aliases from structures_dict that are from conflicting structures
structures_root=os.path.dirname(__file__)+"/../../Biochemistry/Structures/"
file_name=structures_root+'All_ModelSEED_Structures.txt'
all_structures_dict=dict()
with open(file_name) as file_handle:
    for line in file_handle.readlines():
        line=line.strip()
        array=line.split('\t')
        if(array[0] not in all_structures_dict):
            all_structures_dict[array[0]]=dict()
        if(array[7] not in all_structures_dict[array[0]]):
            all_structures_dict[array[0]][array[7]]=dict()
        all_structures_dict[array[0]][array[7]][array[3]]=1

for cpd in structures_dict:
    structure_type = 'InChIKey' if 'InChIKey' in structures_dict[cpd] else 'SMILE'
    struct_entry = structures_dict[cpd][structure_type]
    structure = next(iter(struct_entry))
    curated = all_structures_dict[cpd][structure]
    struct_entry[structure]['alias'] = [a for a in struct_entry[structure]['alias'] if a in curated]
############################################################################

thermodynamics_root=os.path.dirname(__file__)+"/../../Biochemistry/Thermodynamics/"
thermodynamics_dict=dict()
for source in ["KEGG","MetaCyc"]:
    for process in ["Charged","Original"]:
        file_name=thermodynamics_root+'ModelSEED/'+source+'_'+process+'_MolAnalysis.tbl'
        with open(file_name) as file_handle:
            for line in file_handle.readlines():
                line=line.strip()
                array=line.split('\t')
                if(array[0] not in thermodynamics_dict):
                    thermodynamics_dict[array[0]]={'dg':"{0:.2f}".format(float(array[7])),'dge':"{0:.2f}".format(float(array[8]))}
                else:
                    #There's a few (~20) cases where the protonated mol file had a 'NoGroup' cue added by MFAToolkit
                    #So using Original energy
                    if(thermodynamics_dict[array[0]]['dg'] == "10000000.00" and array[7] != "1e+07"):
                        thermodynamics_dict[array[0]]={'dg':"{0:.2f}".format(float(array[7])),'dge':"{0:.2f}".format(float(array[8]))}

for cpd in sorted (compounds_dict.keys()):

    #Default energy and error
    lowest_dg=10000000.0
    lowest_dge=10000000.0

    # Condition 1, no structure, use default
    # Condition 2, structure is InChIKey or SMILE

    structure = None
    if cpd in structures_dict:
        for structure_type in ('InChIKey', 'SMILE'):
            if structure_type in structures_dict[cpd]:
                structure = next(iter(structures_dict[cpd][structure_type]))
                break

    if(structure is not None):
        energies_dict=dict()
        for alias in structures_dict[cpd][structure_type][structure]['alias']:
            if(alias not in thermodynamics_dict):
                continue
            energies_dict[float(thermodynamics_dict[alias]['dg'])]=float(thermodynamics_dict[alias]['dge'])

        #In case where multiple energies because of distribution of bonds
        #Take lowest energy as most likely result of equilibrium
        #If the lowest energy is the default energy (i.e. 10000000)
        #We will still save it
        for energy in energies_dict:
            if(energy < lowest_dg):
                lowest_dg=energy
                lowest_dge=energies_dict[energy]

    # values always saved as list of energy and error
    thermo = compounds_dict[cpd].get('thermodynamics')
    if not isinstance(thermo, dict):
        thermo = dict()
        compounds_dict[cpd]['thermodynamics'] = thermo
    thermo['Group contribution'] = [lowest_dg, lowest_dge]

print("Saving compounds")
compounds_helper.saveCompounds(compounds_dict)
