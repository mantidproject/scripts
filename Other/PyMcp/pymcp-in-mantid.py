#from pymcp import *
execfile('pymcp-mantid.py')
test=mcp()
test.prefix='/Users/jon/Desktop/HE3579/id15/Gd_235K_c_axis_1T_a_'
test.loadspecfile()
test.calibrate_energy(25.004,72.804,420,1210,10)
test.ConvertToPz(118,219,5)
test.RawDataToMantid('data')
out=test.rebinMCP([-20,.1,20])
out.McpToMantid('Gd')
out.RawMcpToMantid('Gd')