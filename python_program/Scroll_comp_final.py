"""
This file is a modified version fo the scroll compressor example available in the 
PDSim. It is the part of the electronic ennex for the work presenterd in Tanveer and Bradshaw (2021) 
"Performance evaluation of low-GWP refrigerants in 1-100 ton scroll compressors" 
presented in Int. J. of Ref.


It requires an excel file as an input that contains the scroll geometry. The required file
can be generated using the file named  "Scroll_comp_scalling".

The output is .h5 filefor each cooling capacity that can be read by the file provided
named "Post processing".

"""
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



def Compressor(hs,Vdisp,Vratio,Thickness,OrbitingRadius,phi_i0=0.0,phi_os=0.3, phi_is = pi,phi_ie=21.76458,phi_o0=-1.06276,rb=0.004319,f = None, OneCycle = False, Ref = 'R410A',Qct=1,do=0.123,r_upper_bearing=0.02,r_lower_bearing=0.02,c_ub=1e-5):

#def Compressor(Te = 273, Tc = 300,Tamb = 25, f = None, OneCycle = False, Ref = 'R410A'):
    
    Te = 280.37
    Tc = 327.594
    Tamb = 35

    ScrollComp=Scroll()
    #This runs if the module code is run directly 
    #ScrollComp.set_scroll_geo(1.53204e-5, 2.63, 0.00459, 0.004099,0.0,1.442,3.6)
    print(Vdisp,Vratio,Thickness,OrbitingRadius,phi_os,phi_is)
    ScrollComp.set_scroll_geo(Vdisp,Vratio,Thickness,OrbitingRadius,phi_i0=0.0,phi_os=phi_os, phi_is = phi_is,phi_ie=phi_ie,phi_o0=phi_o0,hs=hs,rb=rb)
    
    #ScrollComp.set_scroll_geo(2.1723e-4, 2.63, 0.00459582, 0.00565,0.0,1.442,3.639) #Set the scroll wrap geometry
    ScrollComp.set_disc_geo('2Arc',r2 = 0)
    ScrollComp.geo.delta_flank = 1.9e-5
    ScrollComp.geo.delta_radial = 1.9e-5
    #X_d=1
    ScrollComp.geo.delta_suction_offset = 0.0e-3
    #ScrollComp.geo.phi_ie_offset = 0.0
    ScrollComp.Tamb = Tamb + 273.15 #[K]
    ScrollComp.omega = 3600/60*2*pi
    f=3600/60             # frequency for port size calculaton
    #ScrollComp.Tamb = 298.0
    #ScrollComp.journal_tune_factor=1.3
    
    
    
    
    # if Qct<2.6:
    #     d_inlet=0.019
    # elif Qct> 2.6 and Qct<7.6:
    #     d_inlet=0.022
    # elif Qct> 7.6 and Qct<15:
    #     d_inlet=0.02857
    # elif Qct> 15 and Qct<20:
    #     d_inlet=0.0412
    # elif Qct> 20 and Qct<29:
    #     d_inlet=0.0539
    # elif Qct> 29 and Qct<41:
    #     d_inlet=0.0457
    # elif Qct> 41 and Qct<88:
    #     d_inlet=0.0666
    # elif Qct> 88 and Qct<105:
    #     d_inlet=0.092
    # else:
    #     d_inlet=0.092
    # 
    # 
    if Qct<4.5:
        d_outlet=0.0127
    elif Qct> 4.5 and Qct<13.5:
        d_outlet=0.020
    elif Qct> 13.5 and Qct<16.5:
        d_outlet=0.02857
    elif Qct> 16.5 and Qct<29:
        d_outlet=0.0349
    elif Qct> 29 and Qct<87:
        d_outlet=0.041
    elif Qct> 87 and Qct<105:
        d_outlet=0.0539
    else:
        d_outlet=0.0539
    
    
    # d_inlet=0.01
    # d_outlet=0.003
    
    #Temporarily set the bearing dimensions
    ScrollComp.mech = struct()
    
    
    ScrollComp.mech.D_upper_bearing = 2*r_upper_bearing
    ScrollComp.mech.L_upper_bearing = 2*r_upper_bearing
    ScrollComp.mech.c_upper_bearing = c_ub
    ScrollComp.mech.D_crank_bearing = 2*r_upper_bearing
    ScrollComp.mech.L_crank_bearing = 2*r_upper_bearing
    ScrollComp.mech.c_crank_bearing = c_ub
    ScrollComp.mech.D_lower_bearing = 2*r_lower_bearing
    ScrollComp.mech.L_lower_bearing = 2*r_lower_bearing
    ScrollComp.mech.c_lower_bearing = 20e-6
    ScrollComp.mech.thrust_ID = do
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
    
    
    ScrollComp.h_shell = 0.01
    ScrollComp.A_shell = (2*do+hs)*((do)**2*pi/4)
    #ScrollComp.HTC = 0.0
    ScrollComp.HT_corr = 'Dittus-Boelter' #'Jang-Jeong'
    
    # Temperature Lumps
    ScrollComp.OEB_type = 'single-lump' #'single-lump'
    ScrollComp.OEB_solver = 'MDNR'
    ScrollComp.Rshell_oil = 190 #K/kW  from Chen (2000) - PhD thesis
    
    
    ScrollComp.motor = Motor()
    ScrollComp.motor.set_eta(0.9)
    ScrollComp.motor.suction_fraction = 1.0
    
    # Te = Te+ 273.15
    # Tc = Tc+ 273.15
    Tin = Te + 11.1
    DT_sc = 8.3
    temp = State.State(Ref,{'T':Te,'Q':1})
    pe = temp.p
    #pe=2950
    temp.update(dict(T=Tc, Q=1))
    pc = temp.p
    
    PR=pc/pe
    
    X_d=0.7
    #pc=4500
    inletState = State.State(Ref,{'T':Tin,'P':pe})
    
    # defining the port area w.r.t displacement volume
    V=10
    A=Vdisp*f/V
    d_inlet=(4*A/3.14)**0.5

    T2s = ScrollComp.guess_outlet_temp(inletState,pc)
    outletState = State.State(Ref,{'T':T2s,'P':pc})
    
    mdot_guess = inletState.rho*ScrollComp.Vdisp*ScrollComp.omega/(2*pi)
    
    ScrollComp.add_tube(Tube(key1='inlet.1',
                             key2='inlet.2',
                             L=0.3,
                             ID=d_inlet,
                             mdot=mdot_guess, 
                             State1=inletState.copy(),
                             fixed=1,
                             TubeFcn=ScrollComp.TubeCode))
    ScrollComp.add_tube(Tube(key1='outlet.1',
                             key2='outlet.2',
                             L=0.3,
                             ID=d_outlet,
                             mdot=mdot_guess, 
                             State2=outletState.copy(),
                             fixed=2,
                             TubeFcn=ScrollComp.TubeCode))
    
    ScrollComp.auto_add_CVs(inletState, outletState)
    # 
    ScrollComp.auto_add_leakage(flankFunc = ScrollComp.FlankLeakage, 
                                radialFunc = ScrollComp.RadialLeakage)
    
    FP = FlowPath(key1='inlet.2', 
                  key2='sa', 
                  MdotFcn=IsentropicNozzleWrapper(),
                  )
    FP.A = pi*d_inlet**2/4
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
                                 MdotFcn=ScrollComp.D_to_DD, MdotFcn_kwargs = dict(X_d =1)))
    
    #Connect the callbacks for the step, endcycle, heat transfer and lump energy balance
    ScrollComp.connect_callbacks(step_callback = ScrollComp.step_callback,
                                 endcycle_callback = ScrollComp.endcycle_callback,
                                 heat_transfer_callback = ScrollComp.heat_transfer_callback,
                                 lumps_energy_balance_callback = ScrollComp.lump_energy_balance_callback
                                 )
    
    from time import clock
    t1=clock()
    ScrollComp.RK45_eps = 1e-4
    ScrollComp.eps_cycle = 1e-3
    try:
        ScrollComp.precond_solve(key_inlet='inlet.1',
                                 key_outlet='outlet.2',
                                 solver_method='RK45',
                                 OneCycle = OneCycle,
                                 plot_every_cycle= False,
                                 x0 = [320,320], #Guesses [Td,Tlump[0],Tlump[1]]
                                 hmin = 1e-4,
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


    Filepath1='C:\Users\Mohsin\OneDrive - Oklahoma A and M System\Documents\Phd\Scroll comp scalling\cases_final\case 8/'

    Filepath2='E:\OneDrive - Oklahoma A and M System\Documents\Phd\Scroll comp scalling\cases_final\case 8/'
        
    if os.path.isfile(Filepath1)==True:
        Filepath=Filepath1
    else:
        Filepath=Filepath2
    FileName= Filepath+'Qc_t'+'_'+str(Qct)+'.h5'
    h5.write_to_file(ScrollComp, FileName)
    
    return ScrollComp

def ReadingVar(Filepath='',Filename='Wdot',sheet_name='Sheet1'):
    """ A function to read .xlsx files and removing unnamed columns. The output form of Var is a list while column key is a numpy array  """
    Filename=Filepath+Filename +'.xlsx'
    Var = pd.read_excel(Filename, sheet_name=sheet_name)
    Var = Var.loc[:, ~Var.columns.str.contains('^Unnamed')] # removing unnamed columns from the list
    column_key=np.array(Var.columns)
    row_key=np.array(Var.index)

    return Var,column_key


    
## Running the model     
 
# if __name__=='__main__':
#      ScrollComp=Compressor(Te = 7.37, Tc =54.59 ,Tamb = 35, f = None, OneCycle = False, Ref = 'R410A')
# 
if __name__=='__main__':
    
    

    Filepath1='C:\Users\Mohsin\OneDrive - Oklahoma A and M System\Documents\Phd\Scroll comp scalling\cases_final\case 8'

    Filepath2='E:\OneDrive - Oklahoma A and M System\Documents\Phd\Scroll comp scalling\cases_final\case 8'

        
    if os.path.isfile(Filepath1)==True:
        Filepath=Filepath1
    else:
        Filepath=Filepath2
    Filename='/geom1'
    
    data,index_list=ReadingVar(Filepath=Filepath,Filename='/bearing_size',sheet_name='Sheet1')
    for i in range(len(index_list)):
        vars()[index_list[i]]=data[index_list[i]].to_numpy().tolist()
    
    data,index_list=ReadingVar(Filepath=Filepath,Filename=Filename,sheet_name='Sheet1')
    for i in range(len(index_list)):
        vars()[index_list[i]]=data[index_list[i]].to_numpy().tolist()
    
    
    
    # ScrollComp=Compressor(hs=0.0052,Vdisp=2.17e-4,
    #                       Vratio=2.63,
    #                       Thickness=0.00459,
    #                       OrbitingRadius=0.0058,
    #                       phi_i0=0.0,
    #                       phi_os=1.47,
    #                       phi_is =3.8,
    #                       f = None,
    #                       OneCycle = False,
    #                       Ref = 'R410A')
    
    
    bearingLosses_total=np.zeros((len(ro)))
    bearing_lower=np.zeros((len(ro)))
    bearing_crank=np.zeros((len(ro)))
    bearing_thrust=np.zeros((len(ro)))
    
    for i in range(len(ro)):
        print('Vdisp',Vdisp[i],'ro',ro[i],'Qc_t',Qc_t[i])
        ScrollComp=Compressor(hs=hs[i],Vdisp=round(Vdisp[i],8),
                            Vratio=Vr[i],
                            Thickness=Thickness[i],
                            OrbitingRadius=round(ro[i],5),
                            phi_i0=0.0,
                            phi_os=phi_os[i],
                            phi_is =phi_is[i],
                            phi_ie=phi_ie[i],
                            phi_o0=phi_o0[i],
                            rb=rb[i],
                            f = None,
                            OneCycle = False,
                            Ref = 'R410A',Qct=Qc_t[i],
                            do=do_optimal[i],
                            r_upper_bearing=r_ub[i],
                            r_lower_bearing=r_lb[i],
                            c_ub=c_ub[i])
    
        bearingLosses_total[i]=ScrollComp.losses.bearings
        bearing_lower[i]=ScrollComp.losses.lower_bearing[-1]
        bearing_crank[i]=ScrollComp.losses.crank_bearing[-1]
        bearing_thrust[i]=ScrollComp.losses.thrust_bearing[-1]

    data={'Qc_t':Qc_t,'bearingLosses_total':bearingLosses_total,'bearing_lower':bearing_lower,'bearing_crank':bearing_crank,'bearing_thrust':bearing_thrust}
# 
    df=pd.DataFrame(data)
    df.to_excel(excel_writer = Filepath+"/inletState.xlsx", sheet_name='Sheet1',startcol=0)
