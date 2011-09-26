aa=mcp()
aa.prefix='/Users/jon/Desktop/HE3579/id15/ni_250K_1T_normalstick_a'
aa.loadspecfile()
aa.replaceEnergyScale(energyScale)
aa.setXUnits('energy')
aa.ConvertToPz(118,219,5,fixei=True)
aa.rebinMCP([-20,.09,20])