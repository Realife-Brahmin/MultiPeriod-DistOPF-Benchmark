import matplotlib.pyplot as pyplt
import py_dss_interface
import pandas as pd
import numpy as np
import os
import pathlib
import math
import openpyxl
data=pd.read_csv(r"C:\Users\anup.parajuli\Desktop\pythonProject\1440_similar.csv")
print(data)
# dss=py_dss_interface.DSSDLL(r"C:\Program Files\OpenDSS")
# dss_file=r"C:\Users\anup.parajuli\Desktop\pythonProject\prev_123JPT.dss"
script_path=os.path.dirname(os.path.abspath(__file__))
dss_file=pathlib.Path(script_path).joinpath("prev_123JPT.dss")
dss=py_dss_interface.DSSDLL()
interval=data.iloc[:,0]
print(interval)
loadshape=data.iloc[:,1]
dss.text(f"compile [{dss_file}]")
dss.text("New LoadShape.similar npts=1440 minterval=0 csvfile=1440_similar.csv")
# dss.text("New energymeter.meter element=Line.1_2 terminal=1")
dss.text("Batchedit Load..* daily=similar")
print(dss.circuit_all_bus_names())
# dss.text("set mode=daily")
# dss.text("set number=1440")
# dss.text("set stepsize=1m")
# dss.solution_solve()
# dss.text("show voltages LN")
voltage_data = []
for interval in range(1440):
    dss.text(f"set mode=daily")
    dss.text(f"set number={interval+1}")
    dss.text(f"set stepsize=1m")
    dss.solution_solve()
    bus_voltages = dss.circuit_all_bus_vmag_pu()
    voltage_data.append(bus_voltages)
voltage_df = pd.DataFrame(voltage_data)
voltage_df.to_excel("all_bus_voltages_similar.xlsx", index=True)
dss.text("show voltages LN")










