import tools.convert_VenusPCMoutputs_to_inputsCPA as VPCM

VPCM.convert(lat=0.0,lon=0.0,alt=60.0e3,filename='/home/asto/IAA_nodes_sshfs/VenusPCM_Rundir/test_para_056/Xave-CPAP-48perVd.021_A.nc',unit='vmr',
             file_reactions='reactions_VenusPCM_vmr.txt',file_concentrations='concentrations_VenusPCM_vmr.txt')
VPCM.convert(lat=0.0,lon=0.0,alt=60.0e3,filename='/home/asto/IAA_nodes_sshfs/VenusPCM_Rundir/test_para_056/Xave-CPAP-48perVd.021_A.nc',unit='molec',
             file_reactions='reactions_VenusPCM_molec.txt',file_concentrations='concentrations_VenusPCM_molec.txt')