## This file is meant for experimentation with PDSim features and may not
## work for you.  Caveat emptor!!
#
from __future__ import division, print_function

# If being run from the folder that contains the PDSim source tree, 
# remove the current location from the python path and use the 
# site-packages version of PDSim
import sys, os
from math import pi

from PDSim.flow.flow_models import IsentropicNozzleWrapper
from PDSim.flow.flow import FlowPath
from PDSim.scroll import scroll_geo
from PDSim.core.core import struct
from PDSim.scroll.core import Scroll
from PDSim.core.containers import ControlVolume, Tube
from PDSim.core.motor import Motor
from PDSim.scroll.plots import plotScrollSet
from PDSim.plot.plots import debug_plots # (Uncomment if you want to do the debug_plots)
from CoolProp import State
from CoolProp import CoolProp as CP

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import time



def Compressor(Te = 273, Tc = 300,Tamb = 25, f = None, OneCycle = False, Ref = 'R410A'):

    ScrollComp=Scroll()
    #This runs if the module code is run directly 
    # geo=ScrollComp.geoVals(rb=0.003522,phi_i0=0.19829,phi_is=4.7,phi_ie=15.5,phi_o0=-1.1248,phi_os=1.8,phi_oe=15.5,h=0.03289)
    # setDiscGeo(geo)
    #ScrollComp.set_scroll_geo(2.0556e-4, 2.63, 0.00459582, 0.00565,0.0,1.074,3.861477) #Set the scroll wrap geometry
    ScrollComp.set_scroll_geo(2.5323e-4, 2.6, 0.00519, 0.0058,0.0,1.57,3.639)
    ScrollComp.set_disc_geo('2Arc',r2 = 0)
    ScrollComp.geo.delta_flank = 1.9e-5
    ScrollComp.geo.delta_radial = 1.9e-5
    X_d=0.7
    ScrollComp.geo.delta_suction_offset = 0.0e-3
    #ScrollComp.geo.phi_ie_offset = 0.0
    ScrollComp.Tamb = Tamb + 273.15 #[K]
    ScrollComp.omega = 3000/60*2*pi
    #ScrollComp.Tamb = 298.0
    
    #Temporarily set the bearing dimensions
    ScrollComp.mech = struct()
    ScrollComp.mech.D_upper_bearing = 0.04
    ScrollComp.mech.L_upper_bearing = 0.04
    ScrollComp.mech.c_upper_bearing = 20e-6
    ScrollComp.mech.D_crank_bearing = 0.04
    ScrollComp.mech.L_crank_bearing = 0.04
    ScrollComp.mech.c_crank_bearing = 20e-6
    ScrollComp.mech.D_lower_bearing = 0.025
    ScrollComp.mech.L_lower_bearing = 0.025
    ScrollComp.mech.c_lower_bearing = 20e-6
    ScrollComp.mech.thrust_ID = 0.05
    ScrollComp.mech.thrust_friction_coefficient = 0.028 #From Chen thesis
    ScrollComp.mech.orbiting_scroll_mass = 2.5
    ScrollComp.mech.L_ratio_bearings = 5
    ScrollComp.mech.mu_oil = 0.02
    
    ScrollComp.mech.oldham_ring_radius = 0.06 #m
    ScrollComp.mech.oldham_mass = 0.1 #kg
    ScrollComp.mech.oldham_thickness = 0.008 #m
    ScrollComp.mech.oldham_key_height = 0.006 #m
    ScrollComp.mech.oldham_key_width = 0.006 #m
    ScrollComp.mech.oldham_key_friction_coefficient = 0.01 #-
    ScrollComp.mech.oldham_rotation_beta = 0 #rad
    
    # ScrollComp.h_shell = 0.02
    # ScrollComp.A_shell = 0.05
    # ScrollComp.HTC = 3
    
    
    ScrollComp.h_shell = 0.02
    ScrollComp.A_shell = 450.75e-3*((246.126e-3)**2*pi/4)
    #ScrollComp.HTC = 0.0
    ScrollComp.HT_corr = 'Pereira-Deschamps' #'Jang-Jeong'
    
    # Temperature Lumps
    ScrollComp.OEB_type = 'single-lump' #'single-lump'
    ScrollComp.OEB_solver = 'MDNR'
    ScrollComp.Rshell_oil = 190 #K/kW  from Chen (2000) - PhD thesis
    
    
    ScrollComp.motor = Motor()
    ScrollComp.motor.set_eta(0.9)
    ScrollComp.motor.suction_fraction = 1.0
    
    Te = Te+ 273.15
    Tc = Tc+ 273.15
    Tin = Te + 11.1
    DT_sc = 7
    temp = State.State(Ref,{'T':Te,'Q':1})
    pe = temp.p
    #pe=2950
    temp.update(dict(T=Tc, Q=1))
    pc = temp.p
    #pc=4500
    inletState = State.State(Ref,{'T':Tin,'P':pe})

    T2s = ScrollComp.guess_outlet_temp(inletState,pc)
    outletState = State.State(Ref,{'T':T2s,'P':pc})
    
    mdot_guess = inletState.rho*ScrollComp.Vdisp*ScrollComp.omega/(2*pi)
    
    ScrollComp.add_tube(Tube(key1='inlet.1',
                             key2='inlet.2',
                             L=0.3,
                             ID=0.16,
                             mdot=mdot_guess, 
                             State1=inletState.copy(),
                             fixed=1,
                             TubeFcn=ScrollComp.TubeCode))
    ScrollComp.add_tube(Tube(key1='outlet.1',
                             key2='outlet.2',
                             L=0.3,
                             ID=0.0206,
                             mdot=mdot_guess, 
                             State2=outletState.copy(),
                             fixed=2,
                             TubeFcn=ScrollComp.TubeCode))
    
    ScrollComp.auto_add_CVs(inletState, outletState)
    
    ScrollComp.auto_add_leakage(flankFunc = ScrollComp.FlankLeakage, 
                                radialFunc = ScrollComp.RadialLeakage)
    
    FP = FlowPath(key1='inlet.2', 
                  key2='sa', 
                  MdotFcn=IsentropicNozzleWrapper(),
                  )
    FP.A = pi*0.16**2/4
    ScrollComp.add_flow(FP)
    
    ScrollComp.add_flow(FlowPath(key1='sa', 
                                 key2='s1',
                                 MdotFcn=ScrollComp.SA_S1,
                                 MdotFcn_kwargs = dict(X_d = X_d)
                                 )
                        )
    ScrollComp.add_flow(FlowPath(key1 = 'sa',
                                 key2 = 's2',
                                 MdotFcn = ScrollComp.SA_S2,
                                 MdotFcn_kwargs = dict(X_d = X_d)
                                 )
                        )
    
    ScrollComp.add_flow(FlowPath(key1 = 'outlet.1',
                                 key2 = 'dd',
                                 MdotFcn = ScrollComp.DISC_DD,
                                 MdotFcn_kwargs = dict(X_d = X_d)
                                 )
                        )
       
    ScrollComp.add_flow(FlowPath(key1 = 'outlet.1',
                                 key2 = 'ddd',
                                 MdotFcn = ScrollComp.DISC_DD,
                                 MdotFcn_kwargs = dict(X_d = X_d)
                                 )
                        )
#     ScrollComp.add_flow(FlowPath(key1 = 'outlet.1',
#                                  key2 = 'd1',
#                                  MdotFcn = ScrollComp.DISC_D1,
#                                  MdotFcn_kwargs = dict(X_d = 0.7)
#                                  )
#                         )
#     
#     FP = FlowPath(key1='outlet.1', 
#                   key2='dd', 
#                   MdotFcn=IsentropicNozzleWrapper(),
#                   )
#     FP.A = pi*0.006**2/4
#     ScrollComp.add_flow(FP)
#       
#     FP = FlowPath(key1='outlet.1', 
#                   key2='ddd', 
#                   MdotFcn=IsentropicNozzleWrapper(),
#                   )
#     FP.A = pi*0.006**2/4
#     ScrollComp.add_flow(FP)
    
    ScrollComp.add_flow(FlowPath(key1='d1',
                                 key2='dd',
                                 MdotFcn=ScrollComp.D_to_DD, MdotFcn_kwargs = dict(X_d = 1)))
    ScrollComp.add_flow(FlowPath(key1='d2',
                                 key2='dd',
                                 MdotFcn=ScrollComp.D_to_DD, MdotFcn_kwargs = dict(X_d = 1)))
    
    #Connect the callbacks for the step, endcycle, heat transfer and lump energy balance
    ScrollComp.connect_callbacks(step_callback = ScrollComp.step_callback,
                                 endcycle_callback = ScrollComp.endcycle_callback,
                                 heat_transfer_callback = ScrollComp.heat_transfer_callback,
                                 lumps_energy_balance_callback = ScrollComp.lump_energy_balance_callback
                                 )
    
    from time import clock
    t1=clock()
    ScrollComp.RK45_eps = 1e-6
    ScrollComp.eps_cycle = 3e-4
    try:
        ScrollComp.precond_solve(key_inlet='inlet.1',
                                 key_outlet='outlet.2',
                                 solver_method='RK45',
                                 OneCycle = OneCycle,
                                 plot_every_cycle= False,
                                 x0 = [320,320], #Guesses [Td,Tlump[0],Tlump[1]]
                                 #hmin = 1e-3
                                 eps_cycle = 1e-3,
                                 eps_energy_balance = 0.1 #relaxed multi-lump convergence
                                 )
    except BaseException as E:
        print(E)
        raise

    print('time taken', clock()-t1)
    
    #debug_plots(ScrollComp)
    
    del ScrollComp.FlowStorage
    from PDSim.misc.hdf5 import HDF5Writer
    h5 = HDF5Writer()
    import CoolProp
    Tes=str(Te)
    Tcs=str(Tc)


    Tes=str(Te)
    Tcs=str(Tc)


    Filepath="/"
    FileName= Filepath+'Results'+'_'+Tes+'_'+Tcs+'.h5'
    h5.write_to_file(ScrollComp, FileName)
    
    return ScrollComp
    
## Running the model     
if __name__=='__main__':
    #ScrollComp=Compressor(Te = 7, Tc = 68,Tamb = 35, f = None, OneCycle = False, Ref = 'R410A')
    
    Te=[7,20]
    Tc=[54,40]
    mdot=np.zeros((len(Te),len(Tc)))
    for i in range(len(Te)):
        for j in range(len(Tc)):
            ScrollComp=Compressor(Te = Te[i], Tc = Tc[j], f = None, OneCycle = False, Ref = 'R410A')
            mdot[i,j]=ScrollComp.mdot
    
    mdot=mdot*1000           # unit conversion from kg/s to g/s
    Te=[str(i) for i in Te]
    Tc=[str(i) for i in Tc]
    df1=pd.DataFrame(mdot,columns=Tc,index=Te)  
    df1.to_excel(excel_writer = "", sheet_name='Sheet1',startcol=3)  


# if __name__=='__main__':
#     Te=[-23,-10,0,7,20]
#     Tc=[10,25,49,68]
#     Wdot=np.zeros((len(Te),len(Tc)))
#     mdot=np.zeros((len(Te),len(Tc)))
#     for i in range(len(Te)):
#         for j in range(len(Tc)):
#             if Te[i]==-23 and (Tc[j]==49 or Tc[j]==68):
#                 Wdot[i,j]=np.NaN
#                 mdot[i,j]=np.NaN
#         
#             elif Te[i]==-10 and (Tc[j]==68):
#                 Wdot[i,j]=np.NaN
#                 mdot[i,j]=np.NaN
#             elif Te[i]==0 and (Tc[j]==10):
#                 Wdot[i,j]=np.NaN
#                 mdot[i,j]=np.NaN
#             elif Te[i]==7 and (Tc[j]==10):
#                 Wdot[i,j]=np.NaN
#                 mdot[i,j]=np.NaN
#             elif Te[i]==20 and (Tc[j]==10 or Tc[j]==25 or Tc[j]==68):
#                 Wdot[i,j]=np.NaN
#                 mdot[i,j]=np.NaN
#             else:
#                 ScrollComp=Compressor(Te = Te[i], Tc = Tc[j], Tamb = 35, f = None, OneCycle = False, Ref = 'R410A')
#                 Wdot[i,j]=ScrollComp.Wdot
#                 mdot[i,j]=ScrollComp.mdot
#     
# mdot=mdot*1000
    
Te=[str(i) for i in Te]
Tc=[str(i) for i in Tc]
df=pd.DataFrame(mdot,columns=Tc,index=Te)
df1=pd.DataFrame(Wdot,columns=Tc,index=Te)  
df.to_excel(excel_writer = "", sheet_name='Sheet1',startcol=3)  
df1.to_excel(excel_writer = "", sheet_name='Sheet1',startcol=3)  
    
    
