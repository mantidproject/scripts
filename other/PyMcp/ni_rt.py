from pymcp import *
aa=mcp()
#aa.prefix='/Users/jon/Desktop/HE3579/id15/ni_250K_1T_normalstick_a'
#prefix2='/Users/jon/Desktop/HE3579/id15/Gd_235K_c_axis_1T_a_'
#prefix='/Users/jon/Desktop/HE3579/id15/EuFeCoAs_2K_2T_b_'
aa.prefix='/Users/jon/Desktop/HE3579/id15/ni_320k_2_5T_a_'
aa.loadspecfile()
#dat.calibrate_energy(25.004,72.804,420,1210,10)
aa.replaceEnergyScale(energyScale)
aa.setXUnits('energy')
aa.ConvertToPz(118,219,5,fixei=True)
mcpout=aa.rebinMCP([-20,.09,20])