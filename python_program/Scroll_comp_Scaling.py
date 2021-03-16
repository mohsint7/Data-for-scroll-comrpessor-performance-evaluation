"""
This file contains method to calculate the scroll geometry based on various manufacturing
and perfroamcne constraints. It is the part of the electronic ennex
for the work presenterd in Tanveer and Bradshaw (2021) "Performance evaluation of low-GWP
refrigerants in 1-100 ton scroll compressors" presented in Int. J. of Ref.


It requires the path to the file named "Scroll_scaling_functions".
The evaporator and condenser temperatures need to be specified.
The list of the capacities should be provide for which the geometry needs to be calculated.
Volume ratio and some fixed parameters of the scroll geometry as described in the article can be
modified if required.

The file produce an excel file containing the geometry for all specified cooling capacities.

"""










import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
from math import pi
import h5py
import itertools
from scipy.optimize import fsolve
import textwrap
import warnings
import os


import wx
from PDSim.scroll.plots import ScrollAnimForm

from CoolProp import State
from CoolProp import CoolProp as CP

from scroll_scalling_functions import *

## Execution

if __name__=='__main__':
    import json

    # add following command to provide path to refprop if mixtures are required
    #CP.set_config_string(CP.ALTERNATIVE_REFPROP_PATH, 'C:\\Program Files (x86)\\REFPROP\\')


    # evaporator and condensert temperatures
    Te=7
    Tc=54
    Qc_t=[1,3,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100] # list of cooling capacitites  [ton]
    Qc=[x*3.51685 for x in Qc_t] # ton to kW conversion

    # C to Kelvin conversion
    Te=C_to_k(Te)
    Tc=C_to_k(Tc)
    Ref='R1234ZEE' # refrigerant

    Vr=2.63    # volume ratio

    # scroll profile paramters that are kept fixed
    phi_i0=0.0
    phi_os=np.array([1.442])
    phi_is=np.array([3.6])
    Thickness=np.array([0.00459])
    hs1=np.linspace(0.002,0.09,40)
    key='Vdisp'

    # Converting capacity to displacement volume
    Vdis=[0]*len(Qc)     # list initialization
    for i in range(len(Qc)):
        Vdis[i]=Vdisp(Q_c=Qc[i],Te=Te,Tc=Tc,Ref=Ref,f=60)

    # calculation of remaning scroll profile parameters
    phi_ie1,ro1,rb1,phi_o01,Aleakage_total1,Vdisp_Calc=parametric(Vdisp=Vdis,
                                                Vr=Vr,phi_i0=phi_i0,
                                                phi_os=phi_os,phi_is =phi_is,
                                                Thickness=Thickness,hs=hs1,
                                                F=3,delta=19e-6,key=key)

    # calculation of dimensionless parameters for scaling
    gw,gc,do=param_dimen(hs1,rb1,Thickness,phi_ie1,phi_is,Vdis,key=key)


##########################
## stacked plots       ###
##########################
    """ This portion of the code is not compulsory but can be used for visualization of the various
    dimesnionless parameters w.r.t to scroll height"""

    hs_plot=hs1*1000  # convert to mm
    do_plot=do*1000  # convert to mm

    D_limit=1
    Gw_limit=8.5
    Gc_limit=2.5

    limits=[Gw_limit,Gc_limit,D_limit,min([min(element) for element in Aleakage_total1])]

    ylabel=['$G_{w}$','$G_{c}$',r'$D_r$','$A_{leak}$ [$mm^2$]']
    xlabel='$h_{s} [mm]$'

    n_plots=4

    d_dim=D_r(do_plot,Qc,Vdis,Te,Tc,Ref,f,hs1,eta_a=1)
    n_curves=len(Vdis)

    params = {'legend.fontsize': 'large',
            'axes.labelsize': 'large',
            'axes.titlesize':'large',
            'xtick.labelsize':'large',
            'ytick.labelsize':'large'}
    plt.rcParams.update(params)
    stacked_plot(n_plots,n_curves,[gw,gc,d_dim,Aleakage_total1],hs_plot,ylabel,xlabel,limits,legend_param=Qc_t,key=key,frameon=True,legend=True)

########################################
## Selection of optimal geometry     ###
########################################

    # defining maximum limit for the scroll regidity and cuting tool constraint
    Gw_limit=8.5
    Gc_limit=2.5

    # initialization of variables
    hs=np.zeros(len(Qc))
    phi_ie=np.zeros(len(Qc))
    ro=np.zeros(len(Qc))
    rb=np.zeros(len(Qc))
    phi_o0=np.zeros(len(Qc))
    Aleakage_total=np.zeros(len(Qc))
    do_optimal=np.zeros(len(Qc))
    Vdisp_verify=np.zeros(len(Qc))

    # Applying constraints

    P=[0]*len(Qc)
    P1=[0]*len(Qc)
    cons_Do=[0]*len(Qc)
    cons_Gw=[0]*len(Qc)
    cons_Gc=[0]*len(Qc)
    index_n=[0]*len(Qc)
    index_optimal=[0]*len(Qc)
    D_limit=[0]*len(Qc)
    for i in range(len(Qc)):
        P[i]=guess_comp_P(Vdis[i],Te=Te,Tc=Tc,Ref=Ref,f=60, eta_a=0.7)
        D_limit[i]=motor_size(P[i])
        print('d', D_limit[i])

        # return the indices for which the condition satisfies
        cons_Do[i]=find_indices(do[:,i],lambda e: e<D_limit[i])
        cons_Gw[i]=find_indices(gw[:,i],lambda e: e<Gw_limit)
        cons_Gc[i]=find_indices(gc[:,i],lambda e: e<Gc_limit)

        # intersection of all three limits. returns list of indices within limit
        index_n[i]=intersection1(cons_Do[i],cons_Gw[i],cons_Gc[i])


        if len(index_n[i])>1:
            index=index_n[i]
            Aleakage_total2=Aleakage_total1[:,i]
            Aleakage_total3=Aleakage_total2[index[0]:index[-1]] # Slicing the list to within constraint limits

            # finding the minimum leakage area within constraint limit
            index_optimal[i]=find_indices(Aleakage_total2,lambda e: e==min(Aleakage_total3))
            index_optimal[i]=index_optimal[i][0]

        else:
            index=cons_Gw[i]
            index_optimal[i]=index[-2]
            print('gc',gc[index_optimal[i],i])

            if gc[index_optimal[i],i]>Gc_limit:
                index=cons_Gc[i]
                index_optimal[i]=index[-2]
                print(index_optimal)
            else:
                pass



        # Creating the vectors of the optimal geometry
        hs[i]=hs1[index_optimal[i]]
        phi_ie[i]=phi_ie1[index_optimal[i],i]
        phi_o0[i]=phi_o01[index_optimal[i],i]
        ro[i]=ro1[index_optimal[i],i]
        rb[i]=rb1[index_optimal[i],i]
        Aleakage_total[i]=Aleakage_total1[index_optimal[i],i]
        do_optimal[i]=do[index_optimal[i],i]


        Vdisp_verify[i]=-2*pi*hs[i]*rb[i]*ro[i]*(3*pi-2*phi_ie[i]+phi_i0+phi_o0[i])
        print(Vdisp_verify[i])
        # leakage ratio

    plt.figure()
    Vdis1=[0]*len(Vdis)
    for i in range(len(Vdis)):
        Vdis1[i]=Vdis[i]*1e6

    leakage_ratio=Aleakage_total/(Vdis1)

    plt.plot(Qc_t,leakage_ratio)
    plt.xlabel('Capacity [ton]')
    plt.ylabel('leakage area/displacement volume')
    plt.show()
##########################################
## Writing data to excell file        ###
#########################################
    path=''       # provide the path to save the file

    Vr=Vr*np.ones(len(hs))
    Thickness=Thickness*np.ones(len(hs))
    phi_os= phi_os*np.ones(len(hs))
    phi_is= phi_is*np.ones(len(hs))

    data={'Qc':Qc,'Qc_t':Qc_t,'Vdisp':Vdis,'Vr':Vr,'Thickness':Thickness,'ro':ro,'rb':rb,'phi_os':phi_os,'phi_is':phi_is,'phi_o0':phi_o0,'phi_ie':phi_ie,'hs':hs,'A_leakage':Aleakage_total,'do_optimal':do_optimal,'d_limit':D_limit}

    hs2=[round(i,2) for i in hs]
    hs2=[str(i) for i in hs2]


    df=pd.DataFrame(data,index=hs2)
    df.to_excel(excel_writer = path+"/geometry.xlsx", sheet_name='Sheet1',startcol=0)




###############################################
###  Scroll geometry plot                  ####
###############################################

""" This portion of the code can be utilized to visulaize the scroll geometry"""

    # Vdis=Vdis.tolist()
    # Thickness=Thickness.tolist()
    # ro=ro.tolist()
    # phi_os=phi_os.tolist()
    # phi_is=phi_is.tolist()
    #
    #
    #
    # from PDSim.scroll.core import Scroll
    # ScrollComp=Scroll()
    # #This runs if the module code is run directly
    # i=4
    # phi_oe=phi_ie[i]
    # #ScrollComp.set_scroll_geo(2.1723e-4, 2.63, 0.00459582, 0.00565,0.0,1.442,3.639)
    # ScrollComp.set_scroll_geo(Vdisp=Vdis[i],Vratio=Vr,Thickness=Thickness[0],OrbitingRadius=ro[i],phi_i0=phi_i0,phi_os=phi_os[0], phi_is = phi_is[0])
    # #ScrollComp.set_scroll_geo(Vdis[i], Vr, Thickness[i], ro[i],phi_i0,phi_os[i], phi_is[i]) #Set the scroll wrap geometry
    # ScrollComp.set_disc_geo('2Arc',r2 = 0)
    # ScrollComp.geo.delta_flank = 10e-6
    # ScrollComp.geo.delta_radial = 10e-6
    #
    # ScrollComp.geo.delta_suction_offset = 0.0e-3
    #
    # app = wx.App()
    # frame = ScrollAnimForm(ScrollComp.geo)
    # frame.Show()
    # app.MainLoop()
    #
    # pass

    #