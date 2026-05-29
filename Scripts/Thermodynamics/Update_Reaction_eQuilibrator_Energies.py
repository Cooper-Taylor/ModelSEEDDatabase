#!/usr/bin/env python
import os,sys
sys.path.append('../../Libs/Python/')
from BiochemPy import Reactions

label="eQuilibrator"
reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()

############################################################################
##
## We apply/overwrite with eQuilibrator Energies
##
############################################################################

thermodynamics_root=os.path.dirname(__file__)+"/../../Biochemistry/Thermodynamics/"
file_name=thermodynamics_root+'eQuilibrator/MetaNetX_Reaction_Energies.tbl'
eq_reactions=dict()
with open(file_name) as file_handle:
    for line in file_handle.readlines():
        line = line.strip()
        array= line.split('\t')

        if('energy' in array[1] or array[1] == 'nan'):
            continue

        eq_reactions[array[0]]=[float("{0:.2f}".format(float(array[1]))),float("{0:.2f}".format(float(array[2])))]

# print(len(eq_reactions))
# 13, 874 ModelSEED Reactions

for rxn in sorted (reactions_dict.keys()):

    if(rxn not in eq_reactions):
        #NB There are a number of reactions for which there are structures available
        #for every eQuilibrator record (and labeled EQC in reaction notes)
        #but, for whatever reason, we couldn't retrieve an estimated energy for
        #the reaction. This is likely because we couldn't retrieve an estimated
        #energy for every compound with a structure in eQuilibrator
        continue

    # values always saved as list of energy and error
    rxn_thermo = reactions_dict[rxn].get('thermodynamics')
    if(not isinstance(rxn_thermo, dict)):
        rxn_thermo = dict()
        reactions_dict[rxn]['thermodynamics'] = rxn_thermo
    rxn_thermo[label] = eq_reactions[rxn]

print("Saving reactions")
reactions_helper.saveReactions(reactions_dict)