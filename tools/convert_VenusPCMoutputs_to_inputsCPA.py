import numpy as np
import scipy.constants as scc
import matplotlib.pyplot as plt
from copy import copy
from scipy.io import netcdf
from scipy.interpolate import RegularGridInterpolator
from typing import TypeAlias

def convert(lat:float,lon:float,alt:float,filename:str) -> None:
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

    temp=copy(f.variables["temp"].data)
    pres=copy(f.variables["pres"].data)
    vitu=copy(f.variables["vitu"].data)
    vitv=copy(f.variables["vitv"].data)
    vitw=copy(f.variables["vitw"].data)

    # VMR to number densities (#.cm-3) for chemical species
    co2=copy(f.variables["co2"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    h2=copy(f.variables["h2"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    co=copy(f.variables["co"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    h2o=copy(f.variables["h2o"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    o1d=copy(f.variables["o1d"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    o=copy(f.variables["o"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    o2=copy(f.variables["o2"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    o2dg=copy(f.variables["o2dg"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    o3=copy(f.variables["o3"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    h=copy(f.variables["h"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    oh=copy(f.variables["oh"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    ho2=copy(f.variables["ho2"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    h2o2=copy(f.variables["h2o2"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    cl=copy(f.variables["cl"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    clo=copy(f.variables["clo"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    cl2=copy(f.variables["cl2"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    hcl=copy(f.variables["hcl"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    hocl=copy(f.variables["hocl"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    clco=copy(f.variables["clco"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    clco3=copy(f.variables["clco3"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    cocl2=copy(f.variables["cocl2"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    s=copy(f.variables["s"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    so=copy(f.variables["so"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    so2=copy(f.variables["so2"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    so3=copy(f.variables["so3"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    s2o2=copy(f.variables["s2o2"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    ocs=copy(f.variables["ocs"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    hso3=copy(f.variables["hso3"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    h2so4=copy(f.variables["h2so4"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    s2=copy(f.variables["s2"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    clso2=copy(f.variables["clso2"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    oscl=copy(f.variables["oscl"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    n2=copy(f.variables["n2"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    n=copy(f.variables["n"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    n2d=copy(f.variables["n2d"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    no=copy(f.variables["no"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    no2=copy(f.variables["no2"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    no3=copy(f.variables["no3"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    clno=copy(f.variables["clno"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    h2oliq=copy(f.variables["h2oliq"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    h2so4liq=copy(f.variables["h2so4liq"].data)*pres*1e-6*scc.N_A/(scc.R*temp)
    d_tr_micro_H2SO4liq=copy(f.variables["d_tr_micro H2SO4 liq"].data)
    d_tr_micro_H2O_liq=copy(f.variables["d_tr_micro H2O liq"].data)

    # rates in #.cm-3.s-1
    # photodissociations
    k_phot_o2_o=copy(f.variables["k_phot_o2_o"].data)*o2
    k_phot_o2_o1d=copy(f.variables["k_phot_o2_o1d"].data)*o2
    k_phot_co2_o=copy(f.variables["k_phot_co2_o"].data)*co2
    k_phot_co2_o1=copy(f.variables["k_phot_co2_o1"].data)*co2
    k_phot_o3_o1d=copy(f.variables["k_phot_o3_o1d"].data)*o3
    k_phot_o3_o=copy(f.variables["k_phot_o3_o"].data)*o3
    k_phot_h2=copy(f.variables["k_phot_h2"].data)*h2
    k_phot_h2o=copy(f.variables["k_phot_h2o"].data)*h2o
    k_phot_ho2=copy(f.variables["k_phot_ho2"].data)*ho2
    k_phot_h2o2=copy(f.variables["k_phot_h2o2"].data)*h2o2
    k_phot_hcl=copy(f.variables["k_phot_hcl"].data)*hcl
    k_phot_cl2=copy(f.variables["k_phot_cl2"].data)*cl2
    k_phot_hocl=copy(f.variables["k_phot_hocl"].data)*hocl
    k_phot_so2=copy(f.variables["k_phot_so2"].data)*so2
    k_phot_so=copy(f.variables["k_phot_so"].data)*so
    k_phot_so3=copy(f.variables["k_phot_so3"].data)*so3
    k_phot_s2=copy(f.variables["k_phot_s2"].data)*s2
    k_phot_clo=copy(f.variables["k_phot_clo"].data)*clo
    k_phot_ocs=copy(f.variables["k_phot_ocs"].data)*ocs
    k_phot_cocl2=copy(f.variables["k_phot_cocl2"].data)*cocl2
    k_phot_h2so4=copy(f.variables["k_phot_h2so4"].data)*h2so4
    k_phot_no2=copy(f.variables["k_phot_no2"].data)*no2
    k_phot_no=copy(f.variables["k_phot_no"].data)*no
    k_phot_n2=copy(f.variables["k_phot_n2"].data)*n2
    # Neutral chemistry
    k_4_a001=copy(f.variables["k_4_a001"].data)*o*o2*co2
    k_3_a002=copy(f.variables["k_3_a002"].data)*o**2*co2
    k_4_a003=copy(f.variables["k_4_a003"].data)*o*o3
    k_phot_b001=copy(f.variables["k_phot_b001"].data)*o1d*co2
    k_4_b002=copy(f.variables["k_4_b002"].data)*o1d*h2o
    k_4_b003=copy(f.variables["k_4_b003"].data)*o1d*h2
    k_phot_b004=copy(f.variables["k_phot_b004"].data)*o1d*o2
    k_4_b005=copy(f.variables["k_4_b005"].data)*o1d*o3
    k_4_b006=copy(f.variables["k_4_b006"].data)*o1d*o3
    k_4_c001=copy(f.variables["k_4_c001"].data)*o*ho2
    k_3_c002=copy(f.variables["k_3_c002"].data)*o*oh
    k_4_c003=copy(f.variables["k_4_c003"].data)*h*o3
    k_4_c004=copy(f.variables["k_4_c004"].data)*h*ho2
    k_4_c005=copy(f.variables["k_4_c005"].data)*h*ho2
    k_4_c006=copy(f.variables["k_4_c006"].data)*h*ho2
    k_4_c007=copy(f.variables["k_4_c007"].data)*oh*ho2
    k_3_c008=copy(f.variables["k_3_c008"].data)*ho2**2
    k_4_c009=copy(f.variables["k_4_c009"].data)*oh*h2o2
    k_4_c010=copy(f.variables["k_4_c010"].data)*oh*h2
    k_4_c011=copy(f.variables["k_4_c011"].data)*h*o2*co2
    k_4_c012=copy(f.variables["k_4_c012"].data)*o*h2o2
    k_3_c013=copy(f.variables["k_3_c013"].data)*oh**2
    k_4_c014=copy(f.variables["k_4_c014"].data)*oh*o3
    k_4_c015=copy(f.variables["k_4_c015"].data)*ho2*o3
    k_3_c016=copy(f.variables["k_3_c016"].data)*ho2**2*co2
    k_3_c017=copy(f.variables["k_3_c017"].data)*oh**2*co2
    k_3_c018=copy(f.variables["k_3_c018"].data)*h**2*co2
    k_4_d001=copy(f.variables["k_4_d001"].data)*no2*o
    k_4_d002=copy(f.variables["k_4_d002"].data)*no*o3
    k_4_d003=copy(f.variables["k_4_d003"].data)*no*ho2
    k_4_d004=copy(f.variables["k_4_d004"].data)*n*no
    k_4_d005=copy(f.variables["k_4_d005"].data)*n*o2
    k_4_d006=copy(f.variables["k_4_d006"].data)*no2*h
    k_4_d007=copy(f.variables["k_4_d007"].data)*n*o
    k_4_d008=copy(f.variables["k_4_d008"].data)*n*ho2
    k_4_d009=copy(f.variables["k_4_d009"].data)*n*oh
    k_phot_d010=copy(f.variables["k_phot_d010"].data)*n2d*o
    k_phot_d011=copy(f.variables["k_phot_d011"].data)*n2d*n2
    k_4_d012=copy(f.variables["k_4_d012"].data)*n2d*co2
    k_4_d013=copy(f.variables["k_4_d013"].data)*n*o*co2
    k_phot_d014=copy(f.variables["k_phot_d014"].data)*n2d*co
    k_4_d015=copy(f.variables["k_4_d015"].data)*no*o*co2
    k_4_d016=copy(f.variables["k_4_d016"].data)*no2*o*co2
    k_3_d017=copy(f.variables["k_3_d017"].data)*no3*no
    k_4_d018=copy(f.variables["k_4_d018"].data)*no3*o
    k_4_d019=copy(f.variables["k_4_d019"].data)*no*cl*co2
    k_4_d020=copy(f.variables["k_4_d020"].data)*clno*cl
    k_4_d021=copy(f.variables["k_4_d021"].data)*clno*o
    k_4_d022=copy(f.variables["k_4_d022"].data)*no*clo
    k_phot_d023=copy(f.variables["k_phot_d023"].data)*no3
    k_phot_d024=copy(f.variables["k_phot_d024"].data)*no3
    k_phot_d025=copy(f.variables["k_phot_d025"].data)*clno
    k_4_d026=copy(f.variables["k_4_d026"].data)*so*no2
    k_4_e001=copy(f.variables["k_4_e001"].data)*co*oh
    k_4_e002=copy(f.variables["k_4_e002"].data)*co*o
    k_4_f001=copy(f.variables["k_4_f001"].data)*hcl*o1d
    k_4_f002=copy(f.variables["k_4_f002"].data)*hcl*o1d
    k_4_f003=copy(f.variables["k_4_f003"].data)*hcl*o
    k_4_f004=copy(f.variables["k_4_f004"].data)*hcl*oh
    k_4_f005=copy(f.variables["k_4_f005"].data)*clo*o
    k_4_f006=copy(f.variables["k_4_f006"].data)*clo*oh
    k_4_f007=copy(f.variables["k_4_f007"].data)*clo*oh
    k_4_f008=copy(f.variables["k_4_f008"].data)*clo*h2
    k_4_f009=copy(f.variables["k_4_f009"].data)*clo*o3
    k_4_f010=copy(f.variables["k_4_f010"].data)*clo*ho2
    k_4_f011=copy(f.variables["k_4_f011"].data)*clo*ho2
    k_4_f012=copy(f.variables["k_4_f012"].data)*clo*h2o2
    k_4_f013=copy(f.variables["k_4_f013"].data)*clo*co
    k_phot_f014=copy(f.variables["k_phot_f014"].data)*clco
    k_4_f015=copy(f.variables["k_4_f015"].data)*clco*o2
    k_4_f016a=copy(f.variables["k_4_f016a"].data)*0.5*clco3**0.5*cl
    k_4_f016b=copy(f.variables["k_4_f016b"].data)*0.5*clco3**0.5*cl
    k_4_f017a=copy(f.variables["k_4_f017a"].data)*0.5*clco3**0.5*o
    k_4_f017b=copy(f.variables["k_4_f017b"].data)*0.5*clco3**0.5*o
    k_4_f018=copy(f.variables["k_4_f018"].data)*clo*ho2
    k_4_f019=copy(f.variables["k_4_f019"].data)*oh*hocl
    k_4_f020=copy(f.variables["k_4_f020"].data)*o*hocl
    k_3_f021=copy(f.variables["k_3_f021"].data)*cl*cl
    k_4_f022=copy(f.variables["k_4_f022"].data)*clco*o
    k_4_f023=copy(f.variables["k_4_f023"].data)*cl2*o1d
    k_4_f024=copy(f.variables["k_4_f024"].data)*cl2*h
    k_4_f025=copy(f.variables["k_4_f025"].data)*cl*clco
    k_3_f026=copy(f.variables["k_3_f026"].data)*clco*clco
    k_4_f027=copy(f.variables["k_4_f027"].data)*clso2
    k_4_f028=copy(f.variables["k_4_f028"].data)*clso2*o
    k_4_f029=copy(f.variables["k_4_f029"].data)*clso2*h
    k_3_f030=copy(f.variables["k_3_f030"].data)*clso2*clso2
    k_4_f031=copy(f.variables["k_4_f031"].data)*cl*o
    k_4_f032=copy(f.variables["k_4_f032"].data)*cl2*o
    k_4_f033=copy(f.variables["k_4_f033"].data)*clco*oh
    k_4_f034=copy(f.variables["k_4_f034"].data)*cl2*oh
    k_4_f035=copy(f.variables["k_4_f035"].data)*clco*o
    k_4_f036=copy(f.variables["k_4_f036"].data)*clco*cl2
    k_4_f037=copy(f.variables["k_4_f037"].data)*hcl*h
    k_4_f038=copy(f.variables["k_4_f038"].data)*clco*h
    k_4_f039=copy(f.variables["k_4_f039"].data)*cl*h
    k_4_g001=copy(f.variables["k_4_g001"].data)*s*o2
    k_4_g002=copy(f.variables["k_4_g002"].data)*s*o3
    k_4_g003=copy(f.variables["k_4_g003"].data)*so*o2
    k_4_g004=copy(f.variables["k_4_g004"].data)*so*o3
    k_4_g005=copy(f.variables["k_4_g005"].data)*so*oh
    k_4_g006=copy(f.variables["k_4_g006"].data)*s*oh
    k_4_g007=copy(f.variables["k_4_g007"].data)*so*o
    k_4_g008=copy(f.variables["k_4_g008"].data)*so*ho2
    k_4_g009=copy(f.variables["k_4_g009"].data)*so2*o
    k_4_g010=copy(f.variables["k_4_g010"].data)*s*o
    k_4_g011=copy(f.variables["k_4_g011"].data)*so3*h2o
    k_4_g012=copy(f.variables["k_4_g012"].data)*so*clo
    k_4_g013=copy(f.variables["k_4_g013"].data)*so*so3
    k_4_g014=copy(f.variables["k_4_g014"].data)*so3*o
    k_3_g015=copy(f.variables["k_3_g015"].data)*so*so
    k_phot_g016=copy(f.variables["k_phot_g016"].data)*s2o2
    k_4_g017a=copy(f.variables["k_4_g017a"].data)*clco3*0.5*so*0.5
    k_4_g017b=copy(f.variables["k_4_g017b"].data)*clco3*0.5*so*0.5
    k_4_g018=copy(f.variables["k_4_g018"].data)*s*co
    k_4_g019=copy(f.variables["k_4_g019"].data)*clco*so
    k_4_g020=copy(f.variables["k_4_g020"].data)*so2*oh
    k_4_g021=copy(f.variables["k_4_g021"].data)*hso3*o2
    k_3_g022=copy(f.variables["k_3_g022"].data)*s*s
    k_4_g023=copy(f.variables["k_4_g023"].data)*s2*o
    k_4_g024=copy(f.variables["k_4_g024"].data)*s*ocs
    k_4_g025=copy(f.variables["k_4_g025"].data)*ocs*o
    k_4_g026=copy(f.variables["k_4_g026"].data)*s*so3
    k_4_g027=copy(f.variables["k_4_g027"].data)*s*ho2
    k_4_g028=copy(f.variables["k_4_g028"].data)*s*clo
    k_phot_g029=copy(f.variables["k_phot_g029"].data)*h2so4
    k_4_g030=copy(f.variables["k_4_g030"].data)*so3*ocs
    k_4_g031a=copy(f.variables["k_4_g031a"].data)*s2o2*0.5*ocs*0.5
    k_4_g031b=copy(f.variables["k_4_g031b"].data)*s2o2*0.5*ocs*0.5
    k_3_g032=copy(f.variables["k_3_g032"].data)*so*so
    k_phot_h001=copy(f.variables["k_phot_h001"].data)*0.0
    k_phot_h002=copy(f.variables["k_phot_h002"].data)*0.0
    k_phot_h003=copy(f.variables["k_phot_h003"].data)*0.0
    k_phot_j001=copy(f.variables["k_phot_j001"].data)*o2dg
    k_phot_j002=copy(f.variables["k_phot_j002"].data)*o2dg

    rates=[]
    k_phot_o2_o = interp_3D_Venus(k_phot_o2_o,lon,lat,alt)
    rates.append(k_phot_o2_o)
    k_phot_o2_o1d = interp_3D_Venus(k_phot_o2_o1d,lon,lat,alt)
    rates.append(k_phot_o2_o1d)
    k_phot_co2_o = interp_3D_Venus(k_phot_co2_o,lon,lat,alt)
    rates.append(k_phot_co2_o)
    k_phot_co2_o1 = interp_3D_Venus(k_phot_co2_o1,lon,lat,alt)
    rates.append(k_phot_co2_o1)
    k_phot_o3_o1d = interp_3D_Venus(k_phot_o3_o1d,lon,lat,alt)
    rates.append(k_phot_o3_o1d)
    k_phot_o3_o = interp_3D_Venus(k_phot_o3_o,lon,lat,alt)
    rates.append(k_phot_o3_o)
    k_phot_h2 = interp_3D_Venus(k_phot_h2,lon,lat,alt)
    rates.append(k_phot_h2)
    k_phot_h2o = interp_3D_Venus(k_phot_h2o,lon,lat,alt)
    rates.append(k_phot_h2o)
    k_phot_ho2 = interp_3D_Venus(k_phot_ho2,lon,lat,alt)
    rates.append(k_phot_ho2)
    k_phot_h2o2 = interp_3D_Venus(k_phot_h2o2,lon,lat,alt)
    rates.append(k_phot_h2o2)
    k_phot_hcl = interp_3D_Venus(k_phot_hcl,lon,lat,alt)
    rates.append(k_phot_hcl)
    k_phot_cl2 = interp_3D_Venus(k_phot_cl2,lon,lat,alt)
    rates.append(k_phot_cl2)
    k_phot_hocl = interp_3D_Venus(k_phot_hocl,lon,lat,alt)
    rates.append(k_phot_hocl)
    k_phot_so2 = interp_3D_Venus(k_phot_so2,lon,lat,alt)
    rates.append(k_phot_so2)
    k_phot_so = interp_3D_Venus(k_phot_so,lon,lat,alt)
    rates.append(k_phot_so)
    k_phot_so3 = interp_3D_Venus(k_phot_so3,lon,lat,alt)
    rates.append(k_phot_so3)
    k_phot_s2 = interp_3D_Venus(k_phot_s2,lon,lat,alt)
    rates.append(k_phot_s2)
    k_phot_clo = interp_3D_Venus(k_phot_clo,lon,lat,alt)
    rates.append(k_phot_clo)
    k_phot_ocs = interp_3D_Venus(k_phot_ocs,lon,lat,alt)
    rates.append(k_phot_ocs)
    k_phot_cocl2 = interp_3D_Venus(k_phot_cocl2,lon,lat,alt)
    rates.append(k_phot_cocl2)
    k_phot_h2so4 = interp_3D_Venus(k_phot_h2so4,lon,lat,alt)
    rates.append(k_phot_h2so4)
    k_phot_no2 = interp_3D_Venus(k_phot_no2,lon,lat,alt)
    rates.append(k_phot_no2)
    k_phot_no = interp_3D_Venus(k_phot_no,lon,lat,alt)
    rates.append(k_phot_no)
    k_phot_n2 = interp_3D_Venus(k_phot_n2,lon,lat,alt)
    rates.append(k_phot_n2)
    k_4_a001 = interp_3D_Venus(k_4_a001,lon,lat,alt)
    rates.append(k_4_a001)
    k_3_a002 = interp_3D_Venus(k_3_a002,lon,lat,alt)
    rates.append(k_3_a002)
    k_4_a003 = interp_3D_Venus(k_4_a003,lon,lat,alt)
    rates.append(k_4_a003)
    k_phot_b001 = interp_3D_Venus(k_phot_b001,lon,lat,alt)
    rates.append(k_phot_b001)
    k_4_b002 = interp_3D_Venus(k_4_b002,lon,lat,alt)
    rates.append(k_4_b002)
    k_4_b003 = interp_3D_Venus(k_4_b003,lon,lat,alt)
    rates.append(k_4_b003)
    k_phot_b004 = interp_3D_Venus(k_phot_b004,lon,lat,alt)
    rates.append(k_phot_b004)
    k_4_b005 = interp_3D_Venus(k_4_b005,lon,lat,alt)
    rates.append(k_4_b005)
    k_4_b006 = interp_3D_Venus(k_4_b006,lon,lat,alt)
    rates.append(k_4_b006)
    k_4_c001 = interp_3D_Venus(k_4_c001,lon,lat,alt)
    rates.append(k_4_c001)
    k_3_c002 = interp_3D_Venus(k_3_c002,lon,lat,alt)
    rates.append(k_3_c002)
    k_4_c003 = interp_3D_Venus(k_4_c003,lon,lat,alt)
    rates.append(k_4_c003)
    k_4_c004 = interp_3D_Venus(k_4_c004,lon,lat,alt)
    rates.append(k_4_c004)
    k_4_c005 = interp_3D_Venus(k_4_c005,lon,lat,alt)
    rates.append(k_4_c005)
    k_4_c006 = interp_3D_Venus(k_4_c006,lon,lat,alt)
    rates.append(k_4_c006)
    k_4_c007 = interp_3D_Venus(k_4_c007,lon,lat,alt)
    rates.append(k_4_c007)
    k_3_c008 = interp_3D_Venus(k_3_c008,lon,lat,alt)
    rates.append(k_3_c008)
    k_4_c009 = interp_3D_Venus(k_4_c009,lon,lat,alt)
    rates.append(k_4_c009)
    k_4_c010 = interp_3D_Venus(k_4_c010,lon,lat,alt)
    rates.append(k_4_c010)
    k_4_c011 = interp_3D_Venus(k_4_c011,lon,lat,alt)
    rates.append(k_4_c011)
    k_4_c012 = interp_3D_Venus(k_4_c012,lon,lat,alt)
    rates.append(k_4_c012)
    k_3_c013 = interp_3D_Venus(k_3_c013,lon,lat,alt)
    rates.append(k_3_c013)
    k_4_c014 = interp_3D_Venus(k_4_c014,lon,lat,alt)
    rates.append(k_4_c014)
    k_4_c015 = interp_3D_Venus(k_4_c015,lon,lat,alt)
    rates.append(k_4_c015)
    k_3_c016 = interp_3D_Venus(k_3_c016,lon,lat,alt)
    rates.append(k_3_c016)
    k_3_c017 = interp_3D_Venus(k_3_c017,lon,lat,alt)
    rates.append(k_3_c017)
    k_3_c018 = interp_3D_Venus(k_3_c018,lon,lat,alt)
    rates.append(k_3_c018)
    k_4_d001 = interp_3D_Venus(k_4_d001,lon,lat,alt)
    rates.append(k_4_d001)
    k_4_d002 = interp_3D_Venus(k_4_d002,lon,lat,alt)
    rates.append(k_4_d002)
    k_4_d003 = interp_3D_Venus(k_4_d003,lon,lat,alt)
    rates.append(k_4_d003)
    k_4_d004 = interp_3D_Venus(k_4_d004,lon,lat,alt)
    rates.append(k_4_d004)
    k_4_d005 = interp_3D_Venus(k_4_d005,lon,lat,alt)
    rates.append(k_4_d005)
    k_4_d006 = interp_3D_Venus(k_4_d006,lon,lat,alt)
    rates.append(k_4_d006)
    k_4_d007 = interp_3D_Venus(k_4_d007,lon,lat,alt)
    rates.append(k_4_d007)
    k_4_d008 = interp_3D_Venus(k_4_d008,lon,lat,alt)
    rates.append(k_4_d008)
    k_4_d009 = interp_3D_Venus(k_4_d009,lon,lat,alt)
    rates.append(k_4_d009)
    k_phot_d010 = interp_3D_Venus(k_phot_d010,lon,lat,alt)
    rates.append(k_phot_d010)
    k_phot_d011 = interp_3D_Venus(k_phot_d011,lon,lat,alt)
    rates.append(k_phot_d011)
    k_4_d012 = interp_3D_Venus(k_4_d012,lon,lat,alt)
    rates.append(k_4_d012)
    k_4_d013 = interp_3D_Venus(k_4_d013,lon,lat,alt)
    rates.append(k_4_d013)
    k_phot_d014 = interp_3D_Venus(k_phot_d014,lon,lat,alt)
    rates.append(k_phot_d014)
    k_4_d015 = interp_3D_Venus(k_4_d015,lon,lat,alt)
    rates.append(k_4_d015)
    k_4_d016 = interp_3D_Venus(k_4_d016,lon,lat,alt)
    rates.append(k_4_d016)
    k_3_d017 = interp_3D_Venus(k_3_d017,lon,lat,alt)
    rates.append(k_3_d017)
    k_4_d018 = interp_3D_Venus(k_4_d018,lon,lat,alt)
    rates.append(k_4_d018)
    k_4_d019 = interp_3D_Venus(k_4_d019,lon,lat,alt)
    rates.append(k_4_d019)
    k_4_d020 = interp_3D_Venus(k_4_d020,lon,lat,alt)
    rates.append(k_4_d020)
    k_4_d021 = interp_3D_Venus(k_4_d021,lon,lat,alt)
    rates.append(k_4_d021)
    k_4_d022 = interp_3D_Venus(k_4_d022,lon,lat,alt)
    rates.append(k_4_d022)
    k_phot_d023 = interp_3D_Venus(k_phot_d023,lon,lat,alt)
    rates.append(k_phot_d023)
    k_phot_d024 = interp_3D_Venus(k_phot_d024,lon,lat,alt)
    rates.append(k_phot_d024)
    k_phot_d025 = interp_3D_Venus(k_phot_d025,lon,lat,alt)
    rates.append(k_phot_d025)
    k_4_d026 = interp_3D_Venus(k_4_d026,lon,lat,alt)
    rates.append(k_4_d026)
    k_4_e001 = interp_3D_Venus(k_4_e001,lon,lat,alt)
    rates.append(k_4_e001)
    k_4_e002 = interp_3D_Venus(k_4_e002,lon,lat,alt)
    rates.append(k_4_e002)
    k_4_f001 = interp_3D_Venus(k_4_f001,lon,lat,alt)
    rates.append(k_4_f001)
    k_4_f002 = interp_3D_Venus(k_4_f002,lon,lat,alt)
    rates.append(k_4_f002)
    k_4_f003 = interp_3D_Venus(k_4_f003,lon,lat,alt)
    rates.append(k_4_f003)
    k_4_f004 = interp_3D_Venus(k_4_f004,lon,lat,alt)
    rates.append(k_4_f004)
    k_4_f005 = interp_3D_Venus(k_4_f005,lon,lat,alt)
    rates.append(k_4_f005)
    k_4_f006 = interp_3D_Venus(k_4_f006,lon,lat,alt)
    rates.append(k_4_f006)
    k_4_f007 = interp_3D_Venus(k_4_f007,lon,lat,alt)
    rates.append(k_4_f007)
    k_4_f008 = interp_3D_Venus(k_4_f008,lon,lat,alt)
    rates.append(k_4_f008)
    k_4_f009 = interp_3D_Venus(k_4_f009,lon,lat,alt)
    rates.append(k_4_f009)
    k_4_f010 = interp_3D_Venus(k_4_f010,lon,lat,alt)
    rates.append(k_4_f010)
    k_4_f011 = interp_3D_Venus(k_4_f011,lon,lat,alt)
    rates.append(k_4_f011)
    k_4_f012 = interp_3D_Venus(k_4_f012,lon,lat,alt)
    rates.append(k_4_f012)
    k_4_f013 = interp_3D_Venus(k_4_f013,lon,lat,alt)
    rates.append(k_4_f013)
    k_phot_f014 = interp_3D_Venus(k_phot_f014,lon,lat,alt)
    rates.append(k_phot_f014)
    k_4_f015 = interp_3D_Venus(k_4_f015,lon,lat,alt)
    rates.append(k_4_f015)
    k_4_f016a = interp_3D_Venus(k_4_f016a,lon,lat,alt)
    rates.append(k_4_f016a)
    k_4_f016b = interp_3D_Venus(k_4_f016b,lon,lat,alt)
    rates.append(k_4_f016b)
    k_4_f017a = interp_3D_Venus(k_4_f017a,lon,lat,alt)
    rates.append(k_4_f017a)
    k_4_f017b = interp_3D_Venus(k_4_f017b,lon,lat,alt)
    rates.append(k_4_f017b)
    k_4_f018 = interp_3D_Venus(k_4_f018,lon,lat,alt)
    rates.append(k_4_f018)
    k_4_f019 = interp_3D_Venus(k_4_f019,lon,lat,alt)
    rates.append(k_4_f019)
    k_4_f020 = interp_3D_Venus(k_4_f020,lon,lat,alt)
    rates.append(k_4_f020)
    k_3_f021 = interp_3D_Venus(k_3_f021,lon,lat,alt)
    rates.append(k_3_f021)
    k_4_f022 = interp_3D_Venus(k_4_f022,lon,lat,alt)
    rates.append(k_4_f022)
    k_4_f023 = interp_3D_Venus(k_4_f023,lon,lat,alt)
    rates.append(k_4_f023)
    k_4_f024 = interp_3D_Venus(k_4_f024,lon,lat,alt)
    rates.append(k_4_f024)
    k_4_f025 = interp_3D_Venus(k_4_f025,lon,lat,alt)
    rates.append(k_4_f025)
    k_3_f026 = interp_3D_Venus(k_3_f026,lon,lat,alt)
    rates.append(k_3_f026)
    k_4_f027 = interp_3D_Venus(k_4_f027,lon,lat,alt)
    rates.append(k_4_f027)
    k_4_f028 = interp_3D_Venus(k_4_f028,lon,lat,alt)
    rates.append(k_4_f028)
    k_4_f029 = interp_3D_Venus(k_4_f029,lon,lat,alt)
    rates.append(k_4_f029)
    k_3_f030 = interp_3D_Venus(k_3_f030,lon,lat,alt)
    rates.append(k_3_f030)
    k_4_f031 = interp_3D_Venus(k_4_f031,lon,lat,alt)
    rates.append(k_4_f031)
    k_4_f032 = interp_3D_Venus(k_4_f032,lon,lat,alt)
    rates.append(k_4_f032)
    k_4_f033 = interp_3D_Venus(k_4_f033,lon,lat,alt)
    rates.append(k_4_f033)
    k_4_f034 = interp_3D_Venus(k_4_f034,lon,lat,alt)
    rates.append(k_4_f034)
    k_4_f035 = interp_3D_Venus(k_4_f035,lon,lat,alt)
    rates.append(k_4_f035)
    k_4_f036 = interp_3D_Venus(k_4_f036,lon,lat,alt)
    rates.append(k_4_f036)
    k_4_f037 = interp_3D_Venus(k_4_f037,lon,lat,alt)
    rates.append(k_4_f037)
    k_4_f038 = interp_3D_Venus(k_4_f038,lon,lat,alt)
    rates.append(k_4_f038)
    k_4_f039 = interp_3D_Venus(k_4_f039,lon,lat,alt)
    rates.append(k_4_f039)
    k_4_g001 = interp_3D_Venus(k_4_g001,lon,lat,alt)
    rates.append(k_4_g001)
    k_4_g002 = interp_3D_Venus(k_4_g002,lon,lat,alt)
    rates.append(k_4_g002)
    k_4_g003 = interp_3D_Venus(k_4_g003,lon,lat,alt)
    rates.append(k_4_g003)
    k_4_g004 = interp_3D_Venus(k_4_g004,lon,lat,alt)
    rates.append(k_4_g004)
    k_4_g005 = interp_3D_Venus(k_4_g005,lon,lat,alt)
    rates.append(k_4_g005)
    k_4_g006 = interp_3D_Venus(k_4_g006,lon,lat,alt)
    rates.append(k_4_g006)
    k_4_g007 = interp_3D_Venus(k_4_g007,lon,lat,alt)
    rates.append(k_4_g007)
    k_4_g008 = interp_3D_Venus(k_4_g008,lon,lat,alt)
    rates.append(k_4_g008)
    k_4_g009 = interp_3D_Venus(k_4_g009,lon,lat,alt)
    rates.append(k_4_g009)
    k_4_g010 = interp_3D_Venus(k_4_g010,lon,lat,alt)
    rates.append(k_4_g010)
    k_4_g011 = interp_3D_Venus(k_4_g011,lon,lat,alt)
    rates.append(k_4_g011)
    k_4_g012 = interp_3D_Venus(k_4_g012,lon,lat,alt)
    rates.append(k_4_g012)
    k_4_g013 = interp_3D_Venus(k_4_g013,lon,lat,alt)
    rates.append(k_4_g013)
    k_4_g014 = interp_3D_Venus(k_4_g014,lon,lat,alt)
    rates.append(k_4_g014)
    k_3_g015 = interp_3D_Venus(k_3_g015,lon,lat,alt)
    rates.append(k_3_g015)
    k_phot_g016 = interp_3D_Venus(k_phot_g016,lon,lat,alt)
    rates.append(k_phot_g016)
    k_4_g017a = interp_3D_Venus(k_4_g017a,lon,lat,alt)
    rates.append(k_4_g017a)
    k_4_g017b = interp_3D_Venus(k_4_g017b,lon,lat,alt)
    rates.append(k_4_g017b)
    k_4_g018 = interp_3D_Venus(k_4_g018,lon,lat,alt)
    rates.append(k_4_g018)
    k_4_g019 = interp_3D_Venus(k_4_g019,lon,lat,alt)
    rates.append(k_4_g019)
    k_4_g020 = interp_3D_Venus(k_4_g020,lon,lat,alt)
    rates.append(k_4_g020)
    k_4_g021 = interp_3D_Venus(k_4_g021,lon,lat,alt)
    rates.append(k_4_g021)
    k_3_g022 = interp_3D_Venus(k_3_g022,lon,lat,alt)
    rates.append(k_3_g022)
    k_4_g023 = interp_3D_Venus(k_4_g023,lon,lat,alt)
    rates.append(k_4_g023)
    k_4_g024 = interp_3D_Venus(k_4_g024,lon,lat,alt)
    rates.append(k_4_g024)
    k_4_g025 = interp_3D_Venus(k_4_g025,lon,lat,alt)
    rates.append(k_4_g025)
    k_4_g026 = interp_3D_Venus(k_4_g026,lon,lat,alt)
    rates.append(k_4_g026)
    k_4_g027 = interp_3D_Venus(k_4_g027,lon,lat,alt)
    rates.append(k_4_g027)
    k_4_g028 = interp_3D_Venus(k_4_g028,lon,lat,alt)
    rates.append(k_4_g028)
    k_phot_g029 = interp_3D_Venus(k_phot_g029,lon,lat,alt)
    rates.append(k_phot_g029)
    k_4_g030 = interp_3D_Venus(k_4_g030,lon,lat,alt)
    rates.append(k_4_g030)
    k_4_g031a = interp_3D_Venus(k_4_g031a,lon,lat,alt)
    rates.append(k_4_g031a)
    k_4_g031b = interp_3D_Venus(k_4_g031b,lon,lat,alt)
    rates.append(k_4_g031b)
    k_3_g032 = interp_3D_Venus(k_3_g032,lon,lat,alt)
    rates.append(k_3_g032)
    k_phot_h001 = interp_3D_Venus(k_phot_h001,lon,lat,alt)
    rates.append(k_phot_h001)
    k_phot_h002 = interp_3D_Venus(k_phot_h002,lon,lat,alt)
    rates.append(k_phot_h002)
    k_phot_h003 = interp_3D_Venus(k_phot_h003,lon,lat,alt)
    rates.append(k_phot_h003)
    k_phot_j001 = interp_3D_Venus(k_phot_j001,lon,lat,alt)
    rates.append(k_phot_j001)
    k_phot_j002 = interp_3D_Venus(k_phot_j002,lon,lat,alt)
    rates.append(k_phot_j002)

    # closing the NetCDF file
    f.close()

    # Open the file 'reactions_VenusPCM.txt' in read mode to read lines
    with open('reactions_VenusPCM.txt', 'r') as file:
        # Read all lines from the file
        lines = file.readlines()

        # For each line, add the rate of the reaction and write back to the file
        for line,rate in zip(lines,rates):
            # Strip any trailing newline characters from the original line, then append "END OF LINE"
            file.write(line.strip() + " " + rate + "\n")


def interp_3D_Venus(rates:list,lon:float,lat:float,alt:float,file) -> float:
    """interp_3D_Venus: 3D interpolation of data from the Venus PCM to a user
    specific given location

    _extended_summary_

    Parameters
    ----------
    rates : list
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
    lon_dim = copy(file.variables["longitude"].data)
    lat_dim = copy(file.variables["latitude"].data)
    alt_dim = copy(file.variables["altitude"].data)
    tme_dim = copy(file.variables["Time"].data)

    fit_points = [np.array(tme_dim), np.array(alt_dim), np.array(lon_dim), np.array(lat_dim)]
    interp = RegularGridInterpolator(fit_points, rates)
    
    # we return the interpolated value of rates at lon,lat,alt and time=first timestep
    return interp(np.array([tme_dim[0],alt,lon,lat]),method='cubic')

def F(u, v):
    return u * np.cos(u * v) + v * np.sin(u * v)


fit_points = [np.linspace(0, 3, 8), np.linspace(0, 3, 11)]
values = F(*np.meshgrid(*fit_points, indexing='ij'))
ut, vt = np.meshgrid(np.linspace(0, 3, 80), np.linspace(0, 3, 80), indexing='ij')
true_values = F(ut, vt)
test_points = np.array([ut.ravel(), vt.ravel()]).T
print(test_points)
interp = RegularGridInterpolator(fit_points, values)
print(interp(np.array([0.,0.]),method='cubic'))
print(values.shape)
fig, axes = plt.subplots(2, 3, figsize=(10, 6))
axes = axes.ravel()
fig_index = 0
for method in ['linear', 'nearest', 'slinear', 'cubic', 'quintic']:
    im = interp(test_points, method=method).reshape(80, 80)
    axes[fig_index].imshow(im)
    axes[fig_index].set_title(method)
    axes[fig_index].axis("off")
    fig_index += 1
axes[fig_index].imshow(true_values)
axes[fig_index].set_title("True values")
fig.tight_layout()

# fig.show() does not show anything and you egt stuck at error
# qt.qpa.plugin: Could not find the Qt platform plugin "wayland" in ""
# plt.show() does show the plots but you still have
# qt.qpa.plugin: Could not find the Qt platform plugin "wayland" in ""
# in the terminal
plt.show()