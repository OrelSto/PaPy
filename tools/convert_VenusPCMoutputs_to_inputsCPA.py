import numpy as np
import scipy.constants as scc
import matplotlib.pyplot as plt
from copy import deepcopy
from scipy.io import netcdf
from scipy.interpolate import RegularGridInterpolator
from typing import TypeAlias

def convert(lat:float,lon:float,alt:float,filename:str,unit:str,file_reactions:str,file_concentrations:str) -> None:
    """convert: converting Venus PCM output as input for CPA

    convert() uses the NetCDF output file from the Venus PCM
    to edit the concentrations.txt and reactions.txt files 
    for CPA to run. convert() will interpolate the data of 
    the Venus PCM outputs ta the point (lat,lon,alt) given 
    by the user.
    The file eactions_VenusPCM.txt is the reaction list in 
    the CPA format for the Venus PCM
 
    Parameters
    ----------
    lat : float
        latitude of the desired point
    lon : float
        longitude of the desired point
    alt : float
        altitude of the desired point
    filename : string
        name of the output file of the Venus PCM
    """
    # Opening the Venus PCM output file
    # f=netcdf.netcdf_file("Xave-CPAP-48perVd.021.nc", 'r')
    f=netcdf.netcdf_file(filename, 'r')

    temp=deepcopy(f.variables["temp"].data)
    pres=deepcopy(f.variables["pres"].data)
    vitu=deepcopy(f.variables["vitu"].data)
    vitv=deepcopy(f.variables["vitv"].data)
    vitw=deepcopy(f.variables["vitw"].data)

    # VMR as inputs
    if unit == 'vmr':
        # convert to ppmv
        # convert_to_unit = 1.0e6
        # convert to vmr
        convert_to_unit = 1.0
    elif unit == 'molec':
        # convert to #.cm-3
        convert_to_unit = pres*1e-6*scc.N_A/(scc.R*temp)
        print('convert_to_unit:',convert_to_unit[0,30,48,0])
    else:
        convert_to_unit = 1.0
    
    tracers = []

    co2=deepcopy(f.variables["co2"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(co2,lon,lat,alt,f))
    h2=deepcopy(f.variables["h2"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(h2,lon,lat,alt,f))
    co=deepcopy(f.variables["co"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(co,lon,lat,alt,f))
    h2o=deepcopy(f.variables["h2o"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(h2o,lon,lat,alt,f))
    o1d=deepcopy(f.variables["o1d"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(o1d,lon,lat,alt,f))
    o=deepcopy(f.variables["o"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(o,lon,lat,alt,f))
    o2=deepcopy(f.variables["o2"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(o2,lon,lat,alt,f))
    o2dg=deepcopy(f.variables["o2dg"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(o2dg,lon,lat,alt,f))
    o3=deepcopy(f.variables["o3"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(o3,lon,lat,alt,f))
    h=deepcopy(f.variables["h"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(h,lon,lat,alt,f))
    oh=deepcopy(f.variables["oh"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(oh,lon,lat,alt,f))
    ho2=deepcopy(f.variables["ho2"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(ho2,lon,lat,alt,f))
    h2o2=deepcopy(f.variables["h2o2"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(h2o2,lon,lat,alt,f))
    cl=deepcopy(f.variables["cl"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(cl,lon,lat,alt,f))
    clo=deepcopy(f.variables["clo"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(clo,lon,lat,alt,f))
    cl2=deepcopy(f.variables["cl2"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(cl2,lon,lat,alt,f))
    hcl=deepcopy(f.variables["hcl"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(hcl,lon,lat,alt,f))
    hocl=deepcopy(f.variables["hocl"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(hocl,lon,lat,alt,f))
    clco=deepcopy(f.variables["clco"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(clco,lon,lat,alt,f))
    clco3=deepcopy(f.variables["clco3"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(clco3,lon,lat,alt,f))
    cocl2=deepcopy(f.variables["cocl2"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(cocl2,lon,lat,alt,f))
    s=deepcopy(f.variables["s"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(s,lon,lat,alt,f))
    so=deepcopy(f.variables["so"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(so,lon,lat,alt,f))
    so2=deepcopy(f.variables["so2"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(so2,lon,lat,alt,f))
    so3=deepcopy(f.variables["so3"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(so3,lon,lat,alt,f))
    s2o2=deepcopy(f.variables["s2o2"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(s2o2,lon,lat,alt,f))
    ocs=deepcopy(f.variables["ocs"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(ocs,lon,lat,alt,f))
    hso3=deepcopy(f.variables["hso3"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(hso3,lon,lat,alt,f))
    h2so4=deepcopy(f.variables["h2so4"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(h2so4,lon,lat,alt,f))
    s2=deepcopy(f.variables["s2"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(s2,lon,lat,alt,f))
    clso2=deepcopy(f.variables["clso2"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(clso2,lon,lat,alt,f))
    n2=deepcopy(f.variables["n2"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(n2,lon,lat,alt,f))
    n=deepcopy(f.variables["n"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(n,lon,lat,alt,f))
    n2d=deepcopy(f.variables["n2d"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(n2d,lon,lat,alt,f))
    no=deepcopy(f.variables["no"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(no,lon,lat,alt,f))
    no2=deepcopy(f.variables["no2"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(no2,lon,lat,alt,f))
    no3=deepcopy(f.variables["no3"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(no3,lon,lat,alt,f))
    clno=deepcopy(f.variables["clno"].data)*convert_to_unit
    tracers.append(interp_3D_Venus(clno,lon,lat,alt,f))
    # h2oliq=deepcopy(f.variables["h2oliq"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    # h2so4liq=deepcopy(f.variables["h2so4liq"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    # d_tr_micro_H2SO4liq=deepcopy(f.variables["d_tr_micro H2SO4 liq"].data)
    # d_tr_micro_H2O_liq=deepcopy(f.variables["d_tr_micro H2O liq"].data)
    # oscl=deepcopy(f.variables["oscl"].data)*pres*1e-6*scc.N_A/(scc.R*temp)

    # Inputs rates are in molec.cm-3.s-1
    if unit == 'vmr':
        # ppmv.hr-1
        convert_time = 3600.0
        # ppmv.s-1
        # convert_time = 1.0
        convert_to_rates_unit = convert_time * (pres*1e-6*scc.N_A) / (scc.R*temp)
        print('convert_to_rates_unit:',convert_to_rates_unit[0,30,48,0])
    elif unit == 'molec':
        convert_time = 3600.0
        convert_to_rates_unit = convert_time * 1.0
    else:
        convert_time = 1.0
        convert_to_rates_unit = 1.0
    
    # rates in #.cm-3.s-1 or vmr.hr-1
    # photodissociations
    # No need to convert to units of concentrations since it's in s-1
    k_phot_o2_o=deepcopy(f.variables["k_phot_o2_o"].data)*o2*convert_time
    k_phot_o2_o1d=deepcopy(f.variables["k_phot_o2_o1d"].data)*o2*convert_time
    k_phot_co2_o=deepcopy(f.variables["k_phot_co2_o"].data)*co2*convert_time
    k_phot_co2_o1d=deepcopy(f.variables["k_phot_co2_o1d"].data)*co2*convert_time
    k_phot_o3_o1d=deepcopy(f.variables["k_phot_o3_o1d"].data)*o3*convert_time
    k_phot_o3_o=deepcopy(f.variables["k_phot_o3_o"].data)*o3*convert_time
    k_phot_h2=deepcopy(f.variables["k_phot_h2"].data)*h2*convert_time
    k_phot_h2o=deepcopy(f.variables["k_phot_h2o"].data)*h2o*convert_time
    k_phot_ho2=deepcopy(f.variables["k_phot_ho2"].data)*ho2*convert_time
    k_phot_h2o2=deepcopy(f.variables["k_phot_h2o2"].data)*h2o2*convert_time
    k_phot_hcl=deepcopy(f.variables["k_phot_hcl"].data)*hcl*convert_time
    k_phot_cl2=deepcopy(f.variables["k_phot_cl2"].data)*cl2*convert_time
    k_phot_hocl=deepcopy(f.variables["k_phot_hocl"].data)*hocl*convert_time
    k_phot_so2=deepcopy(f.variables["k_phot_so2"].data)*so2*convert_time
    k_phot_so=deepcopy(f.variables["k_phot_so"].data)*so*convert_time
    k_phot_so3=deepcopy(f.variables["k_phot_so3"].data)*so3*convert_time
    k_phot_s2=deepcopy(f.variables["k_phot_s2"].data)*s2*convert_time
    k_phot_clo=deepcopy(f.variables["k_phot_clo"].data)*clo*convert_time
    k_phot_ocs=deepcopy(f.variables["k_phot_ocs"].data)*ocs*convert_time
    k_phot_cocl2=deepcopy(f.variables["k_phot_cocl2"].data)*cocl2*convert_time
    k_phot_h2so4=deepcopy(f.variables["k_phot_h2so4"].data)*h2so4*convert_time
    k_phot_no2=deepcopy(f.variables["k_phot_no2"].data)*no2*convert_time
    k_phot_no=deepcopy(f.variables["k_phot_no"].data)*no*convert_time
    k_phot_n2=deepcopy(f.variables["k_phot_n2"].data)*n2*convert_time
    # Neutral chemistry
    k_4_a001=deepcopy(f.variables["k_4_a001"].data)*o*o2*convert_to_rates_unit
    k_3_a002=deepcopy(f.variables["k_3_a002"].data)*o**2*convert_to_rates_unit
    k_4_a003=deepcopy(f.variables["k_4_a003"].data)*o*o3*convert_to_rates_unit
    k_phot_b001=deepcopy(f.variables["k_phot_b001"].data)*o1d*convert_time
    k_4_b002=deepcopy(f.variables["k_4_b002"].data)*o1d*h2o*convert_to_rates_unit
    k_4_b003=deepcopy(f.variables["k_4_b003"].data)*o1d*h2*convert_to_rates_unit
    k_phot_b004=deepcopy(f.variables["k_phot_b004"].data)*o1d*convert_time
    k_4_b005=deepcopy(f.variables["k_4_b005"].data)*o1d*o3*convert_to_rates_unit
    k_4_b006=deepcopy(f.variables["k_4_b006"].data)*o1d*o3*convert_to_rates_unit
    k_4_c001=deepcopy(f.variables["k_4_c001"].data)*o*ho2*convert_to_rates_unit
    k_3_c002=deepcopy(f.variables["k_3_c002"].data)*o*oh*convert_to_rates_unit
    k_4_c003=deepcopy(f.variables["k_4_c003"].data)*h*o3*convert_to_rates_unit
    k_4_c004=deepcopy(f.variables["k_4_c004"].data)*h*ho2*convert_to_rates_unit
    k_4_c005=deepcopy(f.variables["k_4_c005"].data)*h*ho2*convert_to_rates_unit
    k_4_c006=deepcopy(f.variables["k_4_c006"].data)*h*ho2*convert_to_rates_unit
    k_4_c007=deepcopy(f.variables["k_4_c007"].data)*oh*ho2*convert_to_rates_unit
    k_3_c008=deepcopy(f.variables["k_3_c008"].data)*ho2**2*convert_to_rates_unit
    k_4_c009=deepcopy(f.variables["k_4_c009"].data)*oh*h2o2*convert_to_rates_unit
    k_4_c010=deepcopy(f.variables["k_4_c010"].data)*oh*h2*convert_to_rates_unit
    k_4_c011=deepcopy(f.variables["k_4_c011"].data)*h*o2*convert_to_rates_unit
    k_4_c012=deepcopy(f.variables["k_4_c012"].data)*o*h2o2*convert_to_rates_unit
    k_3_c013=deepcopy(f.variables["k_3_c013"].data)*oh**2*convert_to_rates_unit
    k_4_c014=deepcopy(f.variables["k_4_c014"].data)*oh*o3*convert_to_rates_unit
    k_4_c015=deepcopy(f.variables["k_4_c015"].data)*ho2*o3*convert_to_rates_unit
    k_3_c016=deepcopy(f.variables["k_3_c016"].data)*ho2**2*convert_to_rates_unit
    k_3_c017=deepcopy(f.variables["k_3_c017"].data)*oh**2*convert_to_rates_unit
    k_3_c018=deepcopy(f.variables["k_3_c018"].data)*h**2*convert_to_rates_unit
    k_4_d001=deepcopy(f.variables["k_4_d001"].data)*no2*o*convert_to_rates_unit
    k_4_d002=deepcopy(f.variables["k_4_d002"].data)*no*o3*convert_to_rates_unit
    k_4_d003=deepcopy(f.variables["k_4_d003"].data)*no*ho2*convert_to_rates_unit
    k_4_d004=deepcopy(f.variables["k_4_d004"].data)*n*no*convert_to_rates_unit
    k_4_d005=deepcopy(f.variables["k_4_d005"].data)*n*o2*convert_to_rates_unit
    k_4_d006=deepcopy(f.variables["k_4_d006"].data)*no2*h*convert_to_rates_unit
    k_4_d007=deepcopy(f.variables["k_4_d007"].data)*n*o*convert_to_rates_unit
    k_4_d008=deepcopy(f.variables["k_4_d008"].data)*n*ho2*convert_to_rates_unit
    k_4_d009=deepcopy(f.variables["k_4_d009"].data)*n*oh*convert_to_rates_unit
    k_phot_d010=deepcopy(f.variables["k_phot_d010"].data)*n2d*convert_time
    k_phot_d011=deepcopy(f.variables["k_phot_d011"].data)*n2d*convert_time
    k_4_d012=deepcopy(f.variables["k_4_d012"].data)*n2d*co2*convert_to_rates_unit
    k_4_d013=deepcopy(f.variables["k_4_d013"].data)*n*o*convert_to_rates_unit
    k_phot_d014=deepcopy(f.variables["k_phot_d014"].data)*n2d*convert_time
    k_4_d015=deepcopy(f.variables["k_4_d015"].data)*no*o*convert_to_rates_unit
    k_4_d016=deepcopy(f.variables["k_4_d016"].data)*no2*o*convert_to_rates_unit
    k_3_d017=deepcopy(f.variables["k_3_d017"].data)*no3*no*convert_to_rates_unit
    k_4_d018=deepcopy(f.variables["k_4_d018"].data)*no3*o*convert_to_rates_unit
    k_4_d019=deepcopy(f.variables["k_4_d019"].data)*no*cl*convert_to_rates_unit
    k_4_d020=deepcopy(f.variables["k_4_d020"].data)*clno*cl*convert_to_rates_unit
    k_4_d021=deepcopy(f.variables["k_4_d021"].data)*clno*o*convert_to_rates_unit
    k_4_d022=deepcopy(f.variables["k_4_d022"].data)*no*clo*convert_to_rates_unit
    k_phot_d023=deepcopy(f.variables["k_phot_d023"].data)*no3*convert_time
    k_phot_d024=deepcopy(f.variables["k_phot_d024"].data)*no3*convert_time
    k_phot_d025=deepcopy(f.variables["k_phot_d025"].data)*clno*convert_time
    k_4_d026=deepcopy(f.variables["k_4_d026"].data)*so*no2*convert_to_rates_unit
    k_4_e001=deepcopy(f.variables["k_4_e001"].data)*co*oh*convert_to_rates_unit
    k_4_e002=deepcopy(f.variables["k_4_e002"].data)*co*o*convert_to_rates_unit
    k_4_f001=deepcopy(f.variables["k_4_f001"].data)*hcl*o1d*convert_to_rates_unit
    k_4_f002=deepcopy(f.variables["k_4_f002"].data)*hcl*o1d*convert_to_rates_unit
    k_4_f003=deepcopy(f.variables["k_4_f003"].data)*hcl*o*convert_to_rates_unit
    k_4_f004=deepcopy(f.variables["k_4_f004"].data)*hcl*oh*convert_to_rates_unit
    k_4_f005=deepcopy(f.variables["k_4_f005"].data)*clo*o*convert_to_rates_unit
    k_4_f006=deepcopy(f.variables["k_4_f006"].data)*clo*oh*convert_to_rates_unit
    k_4_f007=deepcopy(f.variables["k_4_f007"].data)*clo*oh*convert_to_rates_unit
    k_4_f008=deepcopy(f.variables["k_4_f008"].data)*clo*h2*convert_to_rates_unit
    k_4_f009=deepcopy(f.variables["k_4_f009"].data)*clo*o3*convert_to_rates_unit
    k_4_f010=deepcopy(f.variables["k_4_f010"].data)*clo*ho2*convert_to_rates_unit
    k_4_f011=deepcopy(f.variables["k_4_f011"].data)*clo*ho2*convert_to_rates_unit
    k_4_f012=deepcopy(f.variables["k_4_f012"].data)*clo*h2o2*convert_to_rates_unit
    k_4_f013=deepcopy(f.variables["k_4_f013"].data)*clo*co*convert_to_rates_unit
    k_phot_f014=deepcopy(f.variables["k_phot_f014"].data)*clco*convert_time
    k_4_f015=deepcopy(f.variables["k_4_f015"].data)*clco*o2*convert_to_rates_unit
    k_4_f016a=deepcopy(f.variables["k_4_f016a"].data)*clco3*0.5*cl*0.5*convert_to_rates_unit
    k_4_f016b=deepcopy(f.variables["k_4_f016b"].data)*clco3*0.5*cl*0.5*convert_to_rates_unit
    k_4_f017a=deepcopy(f.variables["k_4_f017a"].data)*clco3*0.5*o*0.5*convert_to_rates_unit
    k_4_f017b=deepcopy(f.variables["k_4_f017b"].data)*clco3*0.5*o*0.5*convert_to_rates_unit
    k_4_f018=deepcopy(f.variables["k_4_f018"].data)*clo*ho2*convert_to_rates_unit
    k_4_f019=deepcopy(f.variables["k_4_f019"].data)*oh*hocl*convert_to_rates_unit
    k_4_f020=deepcopy(f.variables["k_4_f020"].data)*o*hocl*convert_to_rates_unit
    k_3_f021=deepcopy(f.variables["k_3_f021"].data)*cl*cl*convert_to_rates_unit
    k_4_f022=deepcopy(f.variables["k_4_f022"].data)*clco*o*convert_to_rates_unit
    k_4_f023=deepcopy(f.variables["k_4_f023"].data)*cl2*o1d*convert_to_rates_unit
    k_4_f024=deepcopy(f.variables["k_4_f024"].data)*cl2*h*convert_to_rates_unit
    k_4_f025=deepcopy(f.variables["k_4_f025"].data)*cl*clco*convert_to_rates_unit
    k_3_f026=deepcopy(f.variables["k_3_f026"].data)*clco*clco*convert_to_rates_unit
    k_4_f027=deepcopy(f.variables["k_4_f027"].data)*cl*so2*convert_to_rates_unit
    k_4_f028=deepcopy(f.variables["k_4_f028"].data)*clso2*o*convert_to_rates_unit
    k_4_f029=deepcopy(f.variables["k_4_f029"].data)*clso2*h*convert_to_rates_unit
    k_3_f030=deepcopy(f.variables["k_3_f030"].data)*clso2*clso2*convert_to_rates_unit
    k_4_f031=deepcopy(f.variables["k_4_f031"].data)*cl*o*convert_to_rates_unit
    k_4_f032=deepcopy(f.variables["k_4_f032"].data)*cl2*o*convert_to_rates_unit
    k_4_f033=deepcopy(f.variables["k_4_f033"].data)*clco*oh*convert_to_rates_unit
    k_4_f034=deepcopy(f.variables["k_4_f034"].data)*cl2*oh*convert_to_rates_unit
    k_4_f035=deepcopy(f.variables["k_4_f035"].data)*clco*o*convert_to_rates_unit
    k_4_f036=deepcopy(f.variables["k_4_f036"].data)*clco*cl2*convert_to_rates_unit
    k_4_f037=deepcopy(f.variables["k_4_f037"].data)*hcl*h*convert_to_rates_unit
    k_4_f038=deepcopy(f.variables["k_4_f038"].data)*clco*h*convert_to_rates_unit
    k_4_f039=deepcopy(f.variables["k_4_f039"].data)*cl*h*convert_to_rates_unit
    k_4_g001=deepcopy(f.variables["k_4_g001"].data)*s*o2*convert_to_rates_unit
    k_4_g002=deepcopy(f.variables["k_4_g002"].data)*s*o3*convert_to_rates_unit
    k_4_g003=deepcopy(f.variables["k_4_g003"].data)*so*o2*convert_to_rates_unit
    k_4_g004=deepcopy(f.variables["k_4_g004"].data)*so*o3*convert_to_rates_unit
    k_4_g005=deepcopy(f.variables["k_4_g005"].data)*so*oh*convert_to_rates_unit
    k_4_g006=deepcopy(f.variables["k_4_g006"].data)*s*oh*convert_to_rates_unit
    k_4_g007=deepcopy(f.variables["k_4_g007"].data)*so*o*convert_to_rates_unit
    k_4_g008=deepcopy(f.variables["k_4_g008"].data)*so*ho2*convert_to_rates_unit
    k_4_g009=deepcopy(f.variables["k_4_g009"].data)*so2*o*convert_to_rates_unit
    k_4_g010=deepcopy(f.variables["k_4_g010"].data)*s*o*convert_to_rates_unit
    k_4_g011=deepcopy(f.variables["k_4_g011"].data)*so3*h2o*convert_to_rates_unit
    k_4_g012=deepcopy(f.variables["k_4_g012"].data)*so*clo*convert_to_rates_unit
    k_4_g013=deepcopy(f.variables["k_4_g013"].data)*so*so3*convert_to_rates_unit
    k_4_g014=deepcopy(f.variables["k_4_g014"].data)*so3*o*convert_to_rates_unit
    k_3_g015=deepcopy(f.variables["k_3_g015"].data)*so*so*convert_to_rates_unit
    k_phot_g016=deepcopy(f.variables["k_phot_g016"].data)*s2o2*convert_time
    k_4_g017a=deepcopy(f.variables["k_4_g017a"].data)*clco3*0.5*so*0.5*convert_to_rates_unit
    k_4_g017b=deepcopy(f.variables["k_4_g017b"].data)*clco3*0.5*so*0.5*convert_to_rates_unit
    k_4_g018=deepcopy(f.variables["k_4_g018"].data)*s*co*convert_to_rates_unit
    k_4_g019=deepcopy(f.variables["k_4_g019"].data)*clco*so*convert_to_rates_unit
    k_4_g020=deepcopy(f.variables["k_4_g020"].data)*so2*oh*convert_to_rates_unit
    k_4_g021=deepcopy(f.variables["k_4_g021"].data)*hso3*o2*convert_to_rates_unit
    k_3_g022=deepcopy(f.variables["k_3_g022"].data)*s*s*convert_to_rates_unit
    k_4_g023=deepcopy(f.variables["k_4_g023"].data)*s2*o*convert_to_rates_unit
    k_4_g024=deepcopy(f.variables["k_4_g024"].data)*s*ocs*convert_to_rates_unit
    k_4_g025=deepcopy(f.variables["k_4_g025"].data)*ocs*o*convert_to_rates_unit
    k_4_g026=deepcopy(f.variables["k_4_g026"].data)*s*so3*convert_to_rates_unit
    k_4_g027=deepcopy(f.variables["k_4_g027"].data)*s*ho2*convert_to_rates_unit
    k_4_g028=deepcopy(f.variables["k_4_g028"].data)*s*clo*convert_to_rates_unit
    k_phot_g029=deepcopy(f.variables["k_phot_g029"].data)*h2so4*convert_time
    k_4_g030=deepcopy(f.variables["k_4_g030"].data)*so3*ocs*convert_to_rates_unit
    k_4_g031a=deepcopy(f.variables["k_4_g031a"].data)*s2o2*0.5*ocs*0.5*convert_to_rates_unit
    k_4_g031b=deepcopy(f.variables["k_4_g031b"].data)*s2o2*0.5*ocs*0.5*convert_to_rates_unit
    k_3_g032=deepcopy(f.variables["k_3_g032"].data)*so*so*convert_to_rates_unit
    k_phot_j001=deepcopy(f.variables["k_phot_j001"].data)*o2dg*convert_time
    k_phot_j002=deepcopy(f.variables["k_phot_j002"].data)*o2dg*convert_time

    rates=[]
    k_phot_o2_o = interp_3D_Venus(k_phot_o2_o,lon,lat,alt,f)
    rates.append(k_phot_o2_o)
    k_phot_o2_o1d = interp_3D_Venus(k_phot_o2_o1d,lon,lat,alt,f)
    rates.append(k_phot_o2_o1d)
    k_phot_co2_o = interp_3D_Venus(k_phot_co2_o,lon,lat,alt,f)
    rates.append(k_phot_co2_o)
    k_phot_co2_o1d = interp_3D_Venus(k_phot_co2_o1d,lon,lat,alt,f)
    rates.append(k_phot_co2_o1d)
    k_phot_o3_o1d = interp_3D_Venus(k_phot_o3_o1d,lon,lat,alt,f)
    rates.append(k_phot_o3_o1d)
    k_phot_o3_o = interp_3D_Venus(k_phot_o3_o,lon,lat,alt,f)
    rates.append(k_phot_o3_o)
    k_phot_h2 = interp_3D_Venus(k_phot_h2,lon,lat,alt,f)
    rates.append(k_phot_h2)
    k_phot_h2o = interp_3D_Venus(k_phot_h2o,lon,lat,alt,f)
    rates.append(k_phot_h2o)
    k_phot_ho2 = interp_3D_Venus(k_phot_ho2,lon,lat,alt,f)
    rates.append(k_phot_ho2)
    k_phot_h2o2 = interp_3D_Venus(k_phot_h2o2,lon,lat,alt,f)
    rates.append(k_phot_h2o2)
    k_phot_hcl = interp_3D_Venus(k_phot_hcl,lon,lat,alt,f)
    rates.append(k_phot_hcl)
    k_phot_cl2 = interp_3D_Venus(k_phot_cl2,lon,lat,alt,f)
    rates.append(k_phot_cl2)
    k_phot_hocl = interp_3D_Venus(k_phot_hocl,lon,lat,alt,f)
    rates.append(k_phot_hocl)
    k_phot_so2 = interp_3D_Venus(k_phot_so2,lon,lat,alt,f)
    rates.append(k_phot_so2)
    k_phot_so = interp_3D_Venus(k_phot_so,lon,lat,alt,f)
    rates.append(k_phot_so)
    k_phot_so3 = interp_3D_Venus(k_phot_so3,lon,lat,alt,f)
    rates.append(k_phot_so3)
    k_phot_s2 = interp_3D_Venus(k_phot_s2,lon,lat,alt,f)
    rates.append(k_phot_s2)
    k_phot_clo = interp_3D_Venus(k_phot_clo,lon,lat,alt,f)
    rates.append(k_phot_clo)
    k_phot_ocs = interp_3D_Venus(k_phot_ocs,lon,lat,alt,f)
    rates.append(k_phot_ocs)
    k_phot_cocl2 = interp_3D_Venus(k_phot_cocl2,lon,lat,alt,f)
    rates.append(k_phot_cocl2)
    k_phot_h2so4 = interp_3D_Venus(k_phot_h2so4,lon,lat,alt,f)
    rates.append(k_phot_h2so4)
    k_phot_no2 = interp_3D_Venus(k_phot_no2,lon,lat,alt,f)
    rates.append(k_phot_no2)
    k_phot_no = interp_3D_Venus(k_phot_no,lon,lat,alt,f)
    rates.append(k_phot_no)
    k_phot_n2 = interp_3D_Venus(k_phot_n2,lon,lat,alt,f)
    rates.append(k_phot_n2)
    k_4_a001 = interp_3D_Venus(k_4_a001,lon,lat,alt,f)
    rates.append(k_4_a001)
    k_3_a002 = interp_3D_Venus(k_3_a002,lon,lat,alt,f)
    rates.append(k_3_a002)
    k_4_a003 = interp_3D_Venus(k_4_a003,lon,lat,alt,f)
    rates.append(k_4_a003)
    k_phot_b001 = interp_3D_Venus(k_phot_b001,lon,lat,alt,f)
    rates.append(k_phot_b001)
    k_4_b002 = interp_3D_Venus(k_4_b002,lon,lat,alt,f)
    rates.append(k_4_b002)
    k_4_b003 = interp_3D_Venus(k_4_b003,lon,lat,alt,f)
    rates.append(k_4_b003)
    k_phot_b004 = interp_3D_Venus(k_phot_b004,lon,lat,alt,f)
    rates.append(k_phot_b004)
    k_4_b005 = interp_3D_Venus(k_4_b005,lon,lat,alt,f)
    rates.append(k_4_b005)
    k_4_b006 = interp_3D_Venus(k_4_b006,lon,lat,alt,f)
    rates.append(k_4_b006)
    k_4_c001 = interp_3D_Venus(k_4_c001,lon,lat,alt,f)
    rates.append(k_4_c001)
    k_3_c002 = interp_3D_Venus(k_3_c002,lon,lat,alt,f)
    rates.append(k_3_c002)
    k_4_c003 = interp_3D_Venus(k_4_c003,lon,lat,alt,f)
    rates.append(k_4_c003)
    k_4_c004 = interp_3D_Venus(k_4_c004,lon,lat,alt,f)
    rates.append(k_4_c004)
    k_4_c005 = interp_3D_Venus(k_4_c005,lon,lat,alt,f)
    rates.append(k_4_c005)
    k_4_c006 = interp_3D_Venus(k_4_c006,lon,lat,alt,f)
    rates.append(k_4_c006)
    k_4_c007 = interp_3D_Venus(k_4_c007,lon,lat,alt,f)
    rates.append(k_4_c007)
    k_3_c008 = interp_3D_Venus(k_3_c008,lon,lat,alt,f)
    rates.append(k_3_c008)
    k_4_c009 = interp_3D_Venus(k_4_c009,lon,lat,alt,f)
    rates.append(k_4_c009)
    k_4_c010 = interp_3D_Venus(k_4_c010,lon,lat,alt,f)
    rates.append(k_4_c010)
    k_4_c011 = interp_3D_Venus(k_4_c011,lon,lat,alt,f)
    rates.append(k_4_c011)
    k_4_c012 = interp_3D_Venus(k_4_c012,lon,lat,alt,f)
    rates.append(k_4_c012)
    k_3_c013 = interp_3D_Venus(k_3_c013,lon,lat,alt,f)
    rates.append(k_3_c013)
    k_4_c014 = interp_3D_Venus(k_4_c014,lon,lat,alt,f)
    rates.append(k_4_c014)
    k_4_c015 = interp_3D_Venus(k_4_c015,lon,lat,alt,f)
    rates.append(k_4_c015)
    k_3_c016 = interp_3D_Venus(k_3_c016,lon,lat,alt,f)
    rates.append(k_3_c016)
    k_3_c017 = interp_3D_Venus(k_3_c017,lon,lat,alt,f)
    rates.append(k_3_c017)
    k_3_c018 = interp_3D_Venus(k_3_c018,lon,lat,alt,f)
    rates.append(k_3_c018)
    k_4_d001 = interp_3D_Venus(k_4_d001,lon,lat,alt,f)
    rates.append(k_4_d001)
    k_4_d002 = interp_3D_Venus(k_4_d002,lon,lat,alt,f)
    rates.append(k_4_d002)
    k_4_d003 = interp_3D_Venus(k_4_d003,lon,lat,alt,f)
    rates.append(k_4_d003)
    k_4_d004 = interp_3D_Venus(k_4_d004,lon,lat,alt,f)
    rates.append(k_4_d004)
    k_4_d005 = interp_3D_Venus(k_4_d005,lon,lat,alt,f)
    rates.append(k_4_d005)
    k_4_d006 = interp_3D_Venus(k_4_d006,lon,lat,alt,f)
    rates.append(k_4_d006)
    k_4_d007 = interp_3D_Venus(k_4_d007,lon,lat,alt,f)
    rates.append(k_4_d007)
    k_4_d008 = interp_3D_Venus(k_4_d008,lon,lat,alt,f)
    rates.append(k_4_d008)
    k_4_d009 = interp_3D_Venus(k_4_d009,lon,lat,alt,f)
    rates.append(k_4_d009)
    k_phot_d010 = interp_3D_Venus(k_phot_d010,lon,lat,alt,f)
    rates.append(k_phot_d010)
    k_phot_d011 = interp_3D_Venus(k_phot_d011,lon,lat,alt,f)
    rates.append(k_phot_d011)
    k_4_d012 = interp_3D_Venus(k_4_d012,lon,lat,alt,f)
    rates.append(k_4_d012)
    k_4_d013 = interp_3D_Venus(k_4_d013,lon,lat,alt,f)
    rates.append(k_4_d013)
    k_phot_d014 = interp_3D_Venus(k_phot_d014,lon,lat,alt,f)
    rates.append(k_phot_d014)
    k_4_d015 = interp_3D_Venus(k_4_d015,lon,lat,alt,f)
    rates.append(k_4_d015)
    k_4_d016 = interp_3D_Venus(k_4_d016,lon,lat,alt,f)
    rates.append(k_4_d016)
    k_3_d017 = interp_3D_Venus(k_3_d017,lon,lat,alt,f)
    rates.append(k_3_d017)
    k_4_d018 = interp_3D_Venus(k_4_d018,lon,lat,alt,f)
    rates.append(k_4_d018)
    k_4_d019 = interp_3D_Venus(k_4_d019,lon,lat,alt,f)
    rates.append(k_4_d019)
    k_4_d020 = interp_3D_Venus(k_4_d020,lon,lat,alt,f)
    rates.append(k_4_d020)
    k_4_d021 = interp_3D_Venus(k_4_d021,lon,lat,alt,f)
    rates.append(k_4_d021)
    k_4_d022 = interp_3D_Venus(k_4_d022,lon,lat,alt,f)
    rates.append(k_4_d022)
    k_phot_d023 = interp_3D_Venus(k_phot_d023,lon,lat,alt,f)
    rates.append(k_phot_d023)
    k_phot_d024 = interp_3D_Venus(k_phot_d024,lon,lat,alt,f)
    rates.append(k_phot_d024)
    k_phot_d025 = interp_3D_Venus(k_phot_d025,lon,lat,alt,f)
    rates.append(k_phot_d025)
    k_4_d026 = interp_3D_Venus(k_4_d026,lon,lat,alt,f)
    rates.append(k_4_d026)
    k_4_e001 = interp_3D_Venus(k_4_e001,lon,lat,alt,f)
    rates.append(k_4_e001)
    k_4_e002 = interp_3D_Venus(k_4_e002,lon,lat,alt,f)
    rates.append(k_4_e002)
    k_4_f001 = interp_3D_Venus(k_4_f001,lon,lat,alt,f)
    rates.append(k_4_f001)
    k_4_f002 = interp_3D_Venus(k_4_f002,lon,lat,alt,f)
    rates.append(k_4_f002)
    k_4_f003 = interp_3D_Venus(k_4_f003,lon,lat,alt,f)
    rates.append(k_4_f003)
    k_4_f004 = interp_3D_Venus(k_4_f004,lon,lat,alt,f)
    rates.append(k_4_f004)
    k_4_f005 = interp_3D_Venus(k_4_f005,lon,lat,alt,f)
    rates.append(k_4_f005)
    k_4_f006 = interp_3D_Venus(k_4_f006,lon,lat,alt,f)
    rates.append(k_4_f006)
    k_4_f007 = interp_3D_Venus(k_4_f007,lon,lat,alt,f)
    rates.append(k_4_f007)
    k_4_f008 = interp_3D_Venus(k_4_f008,lon,lat,alt,f)
    rates.append(k_4_f008)
    k_4_f009 = interp_3D_Venus(k_4_f009,lon,lat,alt,f)
    rates.append(k_4_f009)
    k_4_f010 = interp_3D_Venus(k_4_f010,lon,lat,alt,f)
    rates.append(k_4_f010)
    k_4_f011 = interp_3D_Venus(k_4_f011,lon,lat,alt,f)
    rates.append(k_4_f011)
    k_4_f012 = interp_3D_Venus(k_4_f012,lon,lat,alt,f)
    rates.append(k_4_f012)
    k_4_f013 = interp_3D_Venus(k_4_f013,lon,lat,alt,f)
    rates.append(k_4_f013)
    k_phot_f014 = interp_3D_Venus(k_phot_f014,lon,lat,alt,f)
    rates.append(k_phot_f014)
    k_4_f015 = interp_3D_Venus(k_4_f015,lon,lat,alt,f)
    rates.append(k_4_f015)
    k_4_f016a = interp_3D_Venus(k_4_f016a,lon,lat,alt,f)
    k_4_f016b = interp_3D_Venus(k_4_f016b,lon,lat,alt,f)
    k_4_f016 = k_4_f016a + k_4_f016b
    rates.append(k_4_f016)
    k_4_f017a = interp_3D_Venus(k_4_f017a,lon,lat,alt,f)
    k_4_f017b = interp_3D_Venus(k_4_f017b,lon,lat,alt,f)
    k_4_f017 = k_4_f017a + k_4_f017a
    rates.append(k_4_f017)
    k_4_f018 = interp_3D_Venus(k_4_f018,lon,lat,alt,f)
    rates.append(k_4_f018)
    k_4_f019 = interp_3D_Venus(k_4_f019,lon,lat,alt,f)
    rates.append(k_4_f019)
    k_4_f020 = interp_3D_Venus(k_4_f020,lon,lat,alt,f)
    rates.append(k_4_f020)
    k_3_f021 = interp_3D_Venus(k_3_f021,lon,lat,alt,f)
    rates.append(k_3_f021)
    k_4_f022 = interp_3D_Venus(k_4_f022,lon,lat,alt,f)
    rates.append(k_4_f022)
    k_4_f023 = interp_3D_Venus(k_4_f023,lon,lat,alt,f)
    rates.append(k_4_f023)
    k_4_f024 = interp_3D_Venus(k_4_f024,lon,lat,alt,f)
    rates.append(k_4_f024)
    k_4_f025 = interp_3D_Venus(k_4_f025,lon,lat,alt,f)
    rates.append(k_4_f025)
    k_3_f026 = interp_3D_Venus(k_3_f026,lon,lat,alt,f)
    rates.append(k_3_f026)
    k_4_f027 = interp_3D_Venus(k_4_f027,lon,lat,alt,f)
    rates.append(k_4_f027)
    k_4_f028 = interp_3D_Venus(k_4_f028,lon,lat,alt,f)
    rates.append(k_4_f028)
    k_4_f029 = interp_3D_Venus(k_4_f029,lon,lat,alt,f)
    rates.append(k_4_f029)
    k_3_f030 = interp_3D_Venus(k_3_f030,lon,lat,alt,f)
    rates.append(k_3_f030)
    k_4_f031 = interp_3D_Venus(k_4_f031,lon,lat,alt,f)
    rates.append(k_4_f031)
    k_4_f032 = interp_3D_Venus(k_4_f032,lon,lat,alt,f)
    rates.append(k_4_f032)
    k_4_f033 = interp_3D_Venus(k_4_f033,lon,lat,alt,f)
    rates.append(k_4_f033)
    k_4_f034 = interp_3D_Venus(k_4_f034,lon,lat,alt,f)
    rates.append(k_4_f034)
    k_4_f035 = interp_3D_Venus(k_4_f035,lon,lat,alt,f)
    rates.append(k_4_f035)
    k_4_f036 = interp_3D_Venus(k_4_f036,lon,lat,alt,f)
    rates.append(k_4_f036)
    k_4_f037 = interp_3D_Venus(k_4_f037,lon,lat,alt,f)
    rates.append(k_4_f037)
    k_4_f038 = interp_3D_Venus(k_4_f038,lon,lat,alt,f)
    rates.append(k_4_f038)
    k_4_f039 = interp_3D_Venus(k_4_f039,lon,lat,alt,f)
    rates.append(k_4_f039)
    k_4_g001 = interp_3D_Venus(k_4_g001,lon,lat,alt,f)
    rates.append(k_4_g001)
    k_4_g002 = interp_3D_Venus(k_4_g002,lon,lat,alt,f)
    rates.append(k_4_g002)
    k_4_g003 = interp_3D_Venus(k_4_g003,lon,lat,alt,f)
    rates.append(k_4_g003)
    k_4_g004 = interp_3D_Venus(k_4_g004,lon,lat,alt,f)
    rates.append(k_4_g004)
    k_4_g005 = interp_3D_Venus(k_4_g005,lon,lat,alt,f)
    rates.append(k_4_g005)
    k_4_g006 = interp_3D_Venus(k_4_g006,lon,lat,alt,f)
    rates.append(k_4_g006)
    k_4_g007 = interp_3D_Venus(k_4_g007,lon,lat,alt,f)
    rates.append(k_4_g007)
    k_4_g008 = interp_3D_Venus(k_4_g008,lon,lat,alt,f)
    rates.append(k_4_g008)
    k_4_g009 = interp_3D_Venus(k_4_g009,lon,lat,alt,f)
    rates.append(k_4_g009)
    k_4_g010 = interp_3D_Venus(k_4_g010,lon,lat,alt,f)
    rates.append(k_4_g010)
    k_4_g011 = interp_3D_Venus(k_4_g011,lon,lat,alt,f)
    rates.append(k_4_g011)
    k_4_g012 = interp_3D_Venus(k_4_g012,lon,lat,alt,f)
    rates.append(k_4_g012)
    k_4_g013 = interp_3D_Venus(k_4_g013,lon,lat,alt,f)
    rates.append(k_4_g013)
    k_4_g014 = interp_3D_Venus(k_4_g014,lon,lat,alt,f)
    rates.append(k_4_g014)
    k_3_g015 = interp_3D_Venus(k_3_g015,lon,lat,alt,f)
    rates.append(k_3_g015)
    k_phot_g016 = interp_3D_Venus(k_phot_g016,lon,lat,alt,f)
    rates.append(k_phot_g016)
    k_4_g017a = interp_3D_Venus(k_4_g017a,lon,lat,alt,f)
    k_4_g017b = interp_3D_Venus(k_4_g017b,lon,lat,alt,f)
    k_4_g017 = k_4_g017a + k_4_g017b
    rates.append(k_4_g017)
    k_4_g018 = interp_3D_Venus(k_4_g018,lon,lat,alt,f)
    rates.append(k_4_g018)
    k_4_g019 = interp_3D_Venus(k_4_g019,lon,lat,alt,f)
    rates.append(k_4_g019)
    k_4_g020 = interp_3D_Venus(k_4_g020,lon,lat,alt,f)
    rates.append(k_4_g020)
    k_4_g021 = interp_3D_Venus(k_4_g021,lon,lat,alt,f)
    rates.append(k_4_g021)
    k_3_g022 = interp_3D_Venus(k_3_g022,lon,lat,alt,f)
    rates.append(k_3_g022)
    k_4_g023 = interp_3D_Venus(k_4_g023,lon,lat,alt,f)
    rates.append(k_4_g023)
    k_4_g024 = interp_3D_Venus(k_4_g024,lon,lat,alt,f)
    rates.append(k_4_g024)
    k_4_g025 = interp_3D_Venus(k_4_g025,lon,lat,alt,f)
    rates.append(k_4_g025)
    k_4_g026 = interp_3D_Venus(k_4_g026,lon,lat,alt,f)
    rates.append(k_4_g026)
    k_4_g027 = interp_3D_Venus(k_4_g027,lon,lat,alt,f)
    rates.append(k_4_g027)
    k_4_g028 = interp_3D_Venus(k_4_g028,lon,lat,alt,f)
    rates.append(k_4_g028)
    k_phot_g029 = interp_3D_Venus(k_phot_g029,lon,lat,alt,f)
    rates.append(k_phot_g029)
    k_4_g030 = interp_3D_Venus(k_4_g030,lon,lat,alt,f)
    rates.append(k_4_g030)
    k_4_g031a = interp_3D_Venus(k_4_g031a,lon,lat,alt,f)
    k_4_g031b = interp_3D_Venus(k_4_g031b,lon,lat,alt,f)
    k_4_g031 = k_4_g031a + k_4_g031b
    rates.append(k_4_g031)
    k_3_g032 = interp_3D_Venus(k_3_g032,lon,lat,alt,f)
    rates.append(k_3_g032)
    # k_phot_h001 = interp_3D_Venus(k_phot_h001,lon,lat,alt,f)
    # rates.append(k_phot_h001)
    # k_phot_h002 = interp_3D_Venus(k_phot_h002,lon,lat,alt,f)
    # rates.append(k_phot_h002)
    # k_phot_h003 = interp_3D_Venus(k_phot_h003,lon,lat,alt,f)
    # rates.append(k_phot_h003)
    k_phot_j001 = interp_3D_Venus(k_phot_j001,lon,lat,alt,f)
    rates.append(k_phot_j001)
    k_phot_j002 = interp_3D_Venus(k_phot_j002,lon,lat,alt,f)
    rates.append(k_phot_j002)

    # closing the NetCDF file
    f.close()

    # Open the file 'reactions_VenusPCM.txt' in read mode to read lines
    with open('./tools/reactions_VenusPCM.txt', 'r') as file:
        # Read all lines from the file
        lines = file.readlines()

    print('lines',len(lines))
    print('rates',len(rates))

    with open(file_reactions, 'w') as file:
        # For each line, add the rate of the reaction and write back to the file
        for line,rate in zip(lines,rates):
            # print('line: ',line)
            # print('rate: ',rate)
            # Strip any trailing newline characters from the original line, then append "END OF LINE"
            file.write(line.strip() + " " + str(rate) + "\n")
    

    # Open the file 'reactions_VenusPCM.txt' in read mode to read lines
    with open('./tools/concentrations_VenusPCM.txt', 'r') as file:
        # Read all lines from the file
        lines = file.readlines()

    print('lines',len(lines))
    print('tracers',len(tracers))

    with open(file_concentrations, 'w') as file:
        # For each line, add the rate of the reaction and write back to the file
        for line,tracer in zip(lines,tracers):
            # print('line: ',line)
            # print('rate: ',tracer)
            # Strip any trailing newline characters from the original line, then append "END OF LINE"
            file.write(line.strip() + " " + str(tracer) + "\n")


def interp_3D_Venus(variable:list,lon:float,lat:float,alt:float,file:object) -> float:
    """interp_3D_Venus: 3D interpolation of data from the Venus PCM to a user
    specific given location

    _extended_summary_

    Parameters
    ----------
    variable : list
        _description_
    lon : float
        _description_
    lat : float
        _description_
    alt : float
        _description_
    filename : str
        _description_

    Returns
    -------
    list
        _description_
    """

    # Opening the Venus PCM output file
    tme_dim = deepcopy(file.variables["Time"].data)
    # print(tme_dim)
    alt_dim = deepcopy(file.variables["altitude"].data)
    # print(alt_dim)
    # print(len(alt_dim))
    lat_dim = deepcopy(file.variables["latitude"].data)
    # print(lat_dim)
    # print(len(lat_dim))
    lon_dim = deepcopy(file.variables["longitude"].data)
    # lon_dim = np.append(lon_dim,[-180.0])
    # print(lon_dim)
    # print(len(lon_dim))

    fit_points = [np.array(alt_dim), np.array(lat_dim), np.array(lon_dim)]
    interp = RegularGridInterpolator(fit_points, variable[0])
    
    # we return the interpolated value of variable at lon,lat,alt and time=first timestep
    # return interp(np.array([alt,lat,lon]),method='cubic')

    print(str(variable[0,34,48,0]))
    return variable[0,34,48,0]
# def F(u, v):
#     return u * np.cos(u * v) + v * np.sin(u * v)


# fit_points = [np.linspace(0, 3, 8), np.linspace(0, 3, 11)]
# values = F(*np.meshgrid(*fit_points, indexing='ij'))
# ut, vt = np.meshgrid(np.linspace(0, 3, 80), np.linspace(0, 3, 80), indexing='ij')
# true_values = F(ut, vt)
# test_points = np.array([ut.ravel(), vt.ravel()]).T
# print(test_points)
# interp = RegularGridInterpolator(fit_points, values)
# print(interp(np.array([0.,0.]),method='cubic'))
# print(values.shape)
# fig, axes = plt.subplots(2, 3, figsize=(10, 6))
# axes = axes.ravel()
# fig_index = 0
# for method in ['linear', 'nearest', 'slinear', 'cubic', 'quintic']:
#     im = interp(test_points, method=method).reshape(80, 80)
#     axes[fig_index].imshow(im)
#     axes[fig_index].set_title(method)
#     axes[fig_index].axis("off")
#     fig_index += 1
# axes[fig_index].imshow(true_values)
# axes[fig_index].set_title("True values")
# fig.tight_layout()

# # fig.show() does not show anything and you egt stuck at error
# # qt.qpa.plugin: Could not find the Qt platform plugin "wayland" in ""
# # plt.show() does show the plots but you still have
# # qt.qpa.plugin: Could not find the Qt platform plugin "wayland" in ""
# # in the terminal
# plt.show()