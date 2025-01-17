# run in your local directory where the github repository is copied
#$ python3 -m test_tools

import tools.convert_VenusPCMoutputs_to_inputsCPA as VPCM
import numpy as np

altitude = np.linspace(45.0e3,80.0e3,36)
for z in altitude:
    VPCM.convert(lat=0.0,lon=180.0,alt=z,filename='/home/asto/Science/VENUS/VenusPCM_dump/Data/Xave-CPAP-48perVd.021_A.nc',unit='vmr',
             file_reactions='reactions_VenusPCM_vmr_'+str(int(z/1e3))+'.txt',file_concentrations='concentrations_VenusPCM_vmr_'+str(int(z/1e3))+'.txt')
# VPCM.convert(lat=0.0,lon=0.0,alt=62.0e3,filename='/home/asto/Science/VENUS/VenusPCM_dump/Data/Xave-CPAP-48perVd.021_A.nc',unit='molec',
#              file_reactions='reactions_VenusPCM_molec.txt',file_concentrations='concentrations_VenusPCM_molec_76.txt')
