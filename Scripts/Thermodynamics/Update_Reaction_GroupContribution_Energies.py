#!/usr/bin/env python
import sys
sys.path.append('../../Libs/Python/')
from BiochemPy import Compounds, Reactions

label = 'Group contribution'

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()

gc_cpds_dict=dict()
for cpd in compounds_dict:
    cpd_obj = compounds_dict[cpd]
    cpd_thermo = cpd_obj.get('thermodynamics')
    if(isinstance(cpd_thermo, dict) and label in cpd_thermo and cpd_thermo[label][0] != 10000000):
        gc_cpds_dict[cpd]=1

        if(cpd_obj['is_obsolete']):
            for link in cpd_obj['linked_compound'].split(';'):
                if(link in gc_cpds_dict):
                    continue

                link_obj = compounds_dict[link]
                link_thermo = link_obj.get('thermodynamics')
                if(isinstance(link_thermo, dict) and label in link_thermo and link_thermo[label][0] != 10000000):
                    gc_cpds_dict[link]=1

reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()

for rxn in reactions_dict:
    if(reactions_dict[rxn]['status']=='EMPTY'):
        continue

    rxn_cpds_array = reactions_dict[rxn]["stoichiometry"]
    complete = all(rgt['compound'] in gc_cpds_dict for rgt in rxn_cpds_array)

    dg_dge_list = [10000000.0, 10000000.0]

    if(complete):

        # build deltaG of reaction
        dg_sum = 0.0
        dge_sum = 0.0
        for rgt in rxn_cpds_array:
            if(rgt['compound'] not in gc_cpds_dict):
                print("Warning: wrong reaction: "+rxn)

            (dg, dge) = compounds_dict[rgt['compound']]['thermodynamics'][label]

            dg_sum += dg * rgt['coefficient']
            dge_sum += (dge * rgt['coefficient']) ** 2

        dg_dge_list = [float("{0:.2f}".format(dg_sum)),
                       float("{0:.2f}".format(dge_sum ** 0.5))]

    # values always saved as list of energy and error
    rxn_thermo = reactions_dict[rxn].get('thermodynamics')
    if(not isinstance(rxn_thermo, dict)):
        rxn_thermo = dict()
        reactions_dict[rxn]['thermodynamics'] = rxn_thermo
    rxn_thermo[label] = dg_dge_list

print("Saving reactions")
reactions_helper.saveReactions(reactions_dict)
