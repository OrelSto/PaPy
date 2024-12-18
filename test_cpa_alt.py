# run in your local directory where the github repository is copied
#$ python3 -m test_cpa_alt

import time
import math
import json
import numpy as np
import matplotlib.pyplot as plt

from matplotlib import gridspec
from matplotlib import rcParams
rcParams['text.usetex'] = True

from src.chem_pathway_analysis import cpa

altitude = np.linspace(45,80,36)

BP_species = 'SO2'

total_start = time.time()
for z in altitude:
    print("Starting CPA for altitude:",z)
    # record start time
    start = time.time()
    cpa.run_cpa(timestep=120.0,rate_threshold=1.0e-12,t_min=1,
                BP_species=BP_species,
                filename_model='reactions_VenusPCM_vmr_'+str(int(z))+'.txt',
                filename_concentration='concentrations_VenusPCM_vmr_'+str(int(z))+'.txt',
                final_AP_file='ActP_VenusPCM_'+BP_species+'_'+str(int(z))+'.json',
                final_DP_file='DelP_VenusPCM_'+BP_species+'_'+str(int(z))+'.json',
                final_CS_file='ChemS_VenusPCM_'+BP_species+'_'+str(int(z))+'.json',
                final_SL_file='SpecL_VenusPCM_'+BP_species+'_'+str(int(z))+'.json',
                chronicle_writing=False,
                steps_save=False,
                rate_threshold_BP_auto=0.1)
    # record end time
    end = time.time()
    # print the difference between start 
    # and end time in milli. secs
    print("The time of execution of above program is :",(end-start) * 10**3, "ms")

total_end = time.time()
print("The total time of execution of above program is :",(total_end-total_start), "s")


# for z in altitude:
#     print(z)
#     cpa.out.pie_output(target_species=['SO2'],act_P_json='ActP_VenusPCM_'+BP_species+'_'+str(int(z))+'.json',del_P_json='DelP_VenusPCM_'+BP_species+'_'+str(int(z))+'.json',chem_R_json='ChemS_VenusPCM_'+BP_species+'_'+str(int(z))+'.json',spec_L_json='SpecL_VenusPCM_'+BP_species+'_'+str(int(z))+'.json',slow_percent=1.0)


# A new for loop over altitude where we exploit the saved pathways state for each altitude
pathways_contrib = []
pathways_used = []
meaningful_pathways_used = []
mask_pathways_used = []
p_slow = []
deltas_C = []
slow_percent = 5.0

print('ADDING THE DELTA C FOR',BP_species)
for z in altitude:
    t_species,_ = cpa.d_tools.get_compound_dict_from_results(compound=BP_species,SpecL='SpecL_VenusPCM_'+BP_species+'_'+str(int(z))+'.json')
    # Adding the delta C
    # print('At Altitude '+str(z)+' km the '+BP_species+' have a delta C of '+str(t_species["delta"]))
    deltas_C.append(t_species["delta"])

print('ADDING THE PATHWAYS USED FOR',BP_species)
for z in altitude:
    with open('ActP_VenusPCM_'+BP_species+'_'+str(int(z))+'.json', 'r') as active_pathways_file:
        # Parse the JSON data and store it in a variable
        active_pathways_data = json.load(active_pathways_file)
    
    for p in active_pathways_data:
        if not cpa.d_tools.is_pathway_in_list(pathway_to_be_checked=p,list_of_pathways=pathways_used):
            # print('Pathway at altitude '+str(z)+' km not in pathways_used')
            if cpa.d_tools.find_compound_in_merged_list(p["branching points"],BP_species):
                # print('Pathway at altitude '+str(z)+' km has '+BP_species)
                # print('adding pathway in pathway used')
                # Adding a new pathway!
                pathways_used.append(p)

print('SO we have a total of ',len(pathways_used),' in pathways_used')

print('ADDING THE CONTRIBUTIONS FOR',BP_species)
for p in pathways_used:
    p_percent = []
    p_slow_contrib = 0.0
    for z in altitude:
        t_species,_ = cpa.d_tools.get_compound_dict_from_results(compound=BP_species,SpecL='SpecL_VenusPCM_'+BP_species+'_'+str(int(z))+'.json')

        with open('ActP_VenusPCM_'+BP_species+'_'+str(int(z))+'.json', 'r') as active_pathways_file:
            # Parse the JSON data and store it in a variable
            active_pathways_data = json.load(active_pathways_file)
        
        # The pathways is present at this altitude
        # So we add it contribution
        cond_is_it_in = cpa.d_tools.is_pathway_in_list(pathway_to_be_checked=p,list_of_pathways=active_pathways_data)
        # print('cond_is_it_in',cond_is_it_in)
        if cond_is_it_in:
            # print('We have pathway in active_pathways_data at ',z,' km')
            ind_p = cpa.d_tools.find_pathway_in_list(pathway_to_be_found=p,list_of_pathways=active_pathways_data)
            pathway = active_pathways_data[ind_p]
            # print('BPs: ',pathway["branching points"])
            ind = cpa.d_tools.find_compound_in_merged_list(pathway["branching points"],BP_species)[0]
            # We must have the stoich coeff
            k = abs(pathway["branching points"][ind]["stoichiometry"])
            # And we add the percent contribution
            p_contrib = 100.0 * k * pathway["rate"] / abs(t_species["delta"])
            # print('We have the contribution: ',p_contrib)
            if p_contrib >= slow_percent :
                p_percent.append(p_contrib)
            else:
                p_percent.append(0.0)
                p_slow_contrib += p_contrib
        # if not, adding 0.0
        else:
            # print('We DONT have pathway in active_pathways_data at ',z,' km')
            p_percent.append(0.0)
    # Adding this pathway to the array of pathways
    pathways_contrib.append(p_percent)
    # Adding this pathway to the array of pathways
    p_slow.append(p_slow_contrib)

# So, now we have all the pathways but not sorted
# We can sort them by their maximum contribution
max_p_contrib = []
for p in pathways_contrib:
    # we add the maximum contrib from each pathway
    max_p_contrib.append(max(p))
    if max(p) < slow_percent:
        mask_pathways_used.append(False)
    else:
        mask_pathways_used.append(True)


# Now we sort max_p_contrib and pathways_contrib from the same order
pathways_contrib = [x for _, x in sorted(zip(max_p_contrib, pathways_contrib), key=lambda pair: pair[0], reverse=True)]
pathways_used = [x for _, x in sorted(zip(max_p_contrib, pathways_used), key=lambda pair: pair[0], reverse=True)]
mask_pathways_used = [x for _, x in sorted(zip(max_p_contrib, mask_pathways_used), key=lambda pair: pair[0], reverse=True)]

# We save the meaningful pathways
for cond,p in zip(mask_pathways_used,pathways_used):
    if cond:
        meaningful_pathways_used.append(p)
# NOW we have a sorted list of meaningful pathways. It can be saved as a JSON file
cpa.d_tools.save_pathways_to_JSON(pathways=meaningful_pathways_used,filename='MeaningfulPathways_'+BP_species+'.json')

print('We have ',len(meaningful_pathways_used),' meaningful pathways with more than ',slow_percent,' % contribution over ',len(pathways_used),'pathways used in the range of altitudes')

# Now the plot
fig = plt.figure()
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])
ax0 = plt.subplot(gs[0])
ax1 = plt.subplot(gs[1])
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage[version=4]{mhchem}')

ind = 0
# Using Numpy arrays to use Masks
pathways_contrib = np.array(pathways_contrib)
mask_pathways_used = np.array(mask_pathways_used)
print(len(pathways_contrib[mask_pathways_used]))
for p in pathways_contrib[mask_pathways_used]:
    ax1.plot(p,altitude,label=r'P'+str(ind))
    ind += 1


ax1.plot(p,altitude,'k',label=r'P_slow')

ax1.set_xlabel(r'Contribution \%')
ax1.set_ylabel(r'altitude (km)')

ax1.grid()
ax1.legend(loc='upper right', bbox_to_anchor=(1.2, 1.0))

ax0.plot(deltas_C,altitude)
ax0.grid()
# minor grid on too
ax0.xaxis.grid(which='minor')
lim_up = 10.0 ** math.ceil(math.log10(max(deltas_C)))
lim_down = - 10.0 ** math.ceil(math.log10(abs(min(deltas_C))))
ax0.set_xlim(lim_down,lim_up)
# ax0.set_xscale('symlog',linthresh=max(lim_up,abs(lim_down))/100.0)
lin_threshold = min(np.abs(np.array(deltas_C)))
lthresh_cond = max(lim_up,abs(lim_down))/1000.0
if lin_threshold > lthresh_cond:
    ax0.set_xscale('symlog',linthresh=lin_threshold)
else:
    # We force the symlog to have 3 orders of mag around 0 max
    ax0.set_xscale('symlog',linthresh=lthresh_cond)

ax0.set_xlabel(r'VMR.hr$^{-1}$')
ax0.set_ylabel(r'altitude (km)')

fig.suptitle(r'\ce{'+BP_species+r'} results')

plt.show()

# Saving the meaningful pathways into a tex table
cpa.out.list_pathways_Tex(list_P_json='MeaningfulPathways_'+BP_species+'.json',chem_R_json='ChemS_VenusPCM_'+BP_species+'_50.json',filename_sav='SO2_meaningful_pathways.tex')

