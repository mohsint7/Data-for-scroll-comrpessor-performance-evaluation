"""
This file contains functions that will be used for scaling the scroll compressor
for various cooling capacities and refrigerants. It is the part of the electronic ennex
for the work presenterd in Tanveer and Bradshaw (2021) "Performance evaluation of low-GWP
refrigerants in 1-100 ton scroll compressors" presented in Int. J. of Ref.


This file only contains the required function for the scroll scaling and do not run independently

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


def C_to_k(T):

    """ Centigrade to kelvin conversion """
    T=T+273
    return T

def in_mm(inch):
   mm=inch*25.4
   return mm

def kW_hp(kW):
    hp=kW*1.34102
    return hp

def P_sat(T,Q,Ref):

    """ saturation pressure at specific condition """
    temp = State.State(Ref,{'T':T,'Q':1})
    p = temp.p
    return p

def guess_comp_P(V_dis,Te=280,Tc=300,dt_sh=11.1,dt_sc=8.3,Ref='R410A',f=60, eta_a=0.7):

    # Calculation of the compressor inlet and condenser outlet temperature
    T_in=Te+dt_sh
    T_out_cond=Tc-dt_sc

    # Evaporator and condenser pressure
    pe=P_sat(T=Te,Q=1,Ref=Ref)
    if Ref=='R744' or 'carbondioxide':
        # Evaporator and condenser pressure
        pe=P_sat(T=Te,Q=1,Ref=Ref)
        # calculation of optimum high side pressure for carbon dioxide

        p_hp=(2.778-0.0157*(Te-273))*(T_out_cond-273)+0.381*(Te-273)-9.34     #calculated pressure is in bar and temperature used is in C
        p_hp=p_hp*100                #bar to kPa
        pc=p_hp
    else:
        pc=P_sat(T=Tc,Q=1,Ref=Ref)

    # Compressor inlet and Condenser outlet state
    inlet_state = State.State(Ref,{'T':T_in,'P':pe})

    h1 = inlet_state.h
    out_state = inlet_state.copy()
    out_state.update(dict(S = inlet_state.s, P = pc))
    h2s = out_state.h
    if pc > inlet_state.p:
        # Compressor Mode
        h2 = h1 + (h2s-h1)/eta_a
    else:
        # Expander Mode
        h2 = h1 + (h2s-h1)*eta_a
    out_state.update(dict(H = h2, P = pc))

    rho_s=inlet_state.rho
    mdot=rho_s*V_dis*f
    P=mdot*(h2-h1)
    P=kW_hp(P)

    return P

def Vdisp(Q_c=49,Te=280,Tc=300,dt_sh=11.1,dt_sc=8.3,Ref='R410A',eta_v=0.92,f=50):

    """ This function calculates the displacement volume of the compressor from the required cooling capacity of the evaporator.
    Q_c=cooling capacity
    Te=evaporator temperature
    Tc=condenser temperature
    dt_sh=super heat
    dt_Sc=sub cooling
    Ref=refrigerant
    eta_V=volumatric efficiency
    f= frequency"""



    # Calculation of the compressor inlet and condenser outlet temperature
    T_in=Te+dt_sh
    T_out_cond=Tc-dt_sc
    # Evaporator and condenser pressure
    pe=P_sat(T=Te,Q=1,Ref=Ref)

    if Ref=='R744' or 'carbondioxide':
        # Evaporator and condenser pressure
        pe=P_sat(T=Te,Q=1,Ref=Ref)
        # calculation of optimum high side pressure for carbon dioxide

        p_hp=(2.778-0.0157*(Te-273))*(T_out_cond-273)+0.381*(Te-273)-9.34     #calculated pressure is in bar and temperature used is in C
        p_hp=p_hp*100                #bar to kPa
        pc=p_hp
    else:
        pc=P_sat(T=Tc,Q=1,Ref=Ref)

    # Compressor inlet and Condenser outlet state
    inletState = State.State(Ref,{'T':T_in,'P':pe})
    outletState_cond=State.State(Ref,{'P':pc,'T':T_out_cond})

    # Enthalpies at both states mentioned above
    h4=outletState_cond.h  # Evaporator inlet
    h1=inletState.h

    rho_s=inletState.rho

    # required mass flow rate to achieve required cooling capacity
    mdot_evap=Q_c/(h1-h4)
    Vdot=mdot_evap/rho_s

    # Compressor displced volume
    Vdisp=Vdot/(eta_v*f)

    return Vdisp

def guess_Vratio(Te=280,Tc=300,n=1.26):

    """ This function guess the the volume ratio, assuming polytropic compression process
    Te= Evaporator temperature
    Tc=condenser temp.
    n=polytropic constant """

    pe=P_sat(T=Te,Q=1,Ref='R410A')
    pc=P_sat(T=Tc,Q=1,Ref='R410A')

    Vr=(pc/pe)**(1/n)

    return Vr

def f(x,hs,phi_os,Vdisp_goal,Vratio_goal,t_goal,phi_i0_goal):

    """ Function containing scroll compressor dimention.
    This function will be used in solve to calculate the
    remaining parameters of the scroll profile """
    pi=math.pi
    phi_ie=x[0]
    ro=x[1]
    rb=x[2]
    phi_o0=x[3]

    t=(rb*pi)-ro
    phi_i0=(t/rb)+phi_o0
    #t=rb*(phi_i0-phi_o0)
    #ro=rb*pi-t
    Vdisp=-2*pi*hs*rb*ro*(3*pi-2*phi_ie+phi_i0+phi_o0)
    Vratio=(3*pi-2*phi_ie+phi_i0+phi_o0)/(-2*phi_os-3*pi+phi_i0+phi_o0)

    r1=Vdisp-Vdisp_goal
    r2=Vratio-Vratio_goal
    r3=t-t_goal
    r4=phi_i0-phi_i0_goal
    #bnprint([r1,r2,r3,r4])
    return [r1,r2,r3,r4]

def pitch(rb):
    p=2*pi*rb
    return p

def G_w(h,t):

    """ Dimensionless parameter representing the rigidness of the scroll """
    G_w=h/t
    return G_w

def G_c(h,rb,t):

    """ Dimensionless paramter describing the  cutting toll constraints """
    p=pitch(rb)
    G_c=h/(p-t)
    return G_c

def D_o(rb,phi_oe):

    """ Minimum outer diamater for the scroll with set compressor profile """
    d=2*rb*math.sqrt(1+phi_oe**2)
    return d

def parametric(Vdisp=2.1723e-4,Vr=2.63,phi_i0=0.0,phi_os=0.3,phi_is =math.pi,Thickness=0.0045,hs=[0.01,0.02],F=3,delta=16e-6,key='phi_is'):

    """ This function takes the input scroll dimensional paramters to calculate the scroll profiles for different scroll height. Also, it calculates the total leakage area (radial+ flank) fromt the calculated profile """



    if key=='phi_os':
        n_column=len(phi_os)

        # initializing arrays
        phi_ie=np.zeros((len(hs),n_column))
        ro=np.zeros((len(hs),n_column))
        rb=np.zeros((len(hs),n_column))
        phi_o0=np.zeros((len(hs),n_column))
        Aleakage_total=np.zeros((len(hs),n_column))
        A_radial=np.zeros((len(hs),n_column))
        A_flank=np.zeros((len(hs),n_column))

        for j in range(n_column):

            for i in range(len(hs)):
                phi_ie[i,j],ro[i,j],rb[i,j],phi_o0[i,j] = fsolve(f,[20,1.3,0.03,0.003],args=(hs[i],phi_os[j],Vdisp,Vr,Thickness[0],phi_i0))
                A_radial[i,j],A_flank[i,j]=leakage_area(F,hs[i],delta,rb[i,j],phi_ie[i,j],phi_is,phi_i0)
                Aleakage_total[i,j]=A_radial[i,j]+A_flank[i,j]


    elif key=='Thickness':
        n_column=len(Thickness)

        # initializing arrays
        phi_ie=np.zeros((len(hs),n_column))
        ro=np.zeros((len(hs),n_column))
        rb=np.zeros((len(hs),n_column))
        phi_o0=np.zeros((len(hs),n_column))
        Aleakage_total=np.zeros((len(hs),n_column))
        A_radial=np.zeros((len(hs),n_column))
        A_flank=np.zeros((len(hs),n_column))

        for j in range(n_column):

            for i in range(len(hs)):
                phi_ie[i,j],ro[i,j],rb[i,j],phi_o0[i,j] = fsolve(f,[20,0.03,0.03,0.003],args=(hs[i],phi_os[0],Vdisp,Vr,Thickness[j],phi_i0))
                A_radial[i,j],A_flank[i,j]=leakage_area(F,hs[i],delta,rb[i,j],phi_ie[i,j],phi_is,phi_i0)
                Aleakage_total[i,j]=A_radial[i,j]+A_flank[i,j]

    elif key=='Vdisp':
        n_column=len(Vdisp)

        # initializing arrays
        phi_ie=np.zeros((len(hs),n_column))
        ro=np.zeros((len(hs),n_column))
        rb=np.zeros((len(hs),n_column))
        phi_o0=np.zeros((len(hs),n_column))
        Aleakage_total=np.zeros((len(hs),n_column))
        A_radial=np.zeros((len(hs),n_column))
        A_flank=np.zeros((len(hs),n_column))
        Vdisp_Calc=np.zeros((len(hs),n_column))

        for j in range(n_column):

            for i in range(len(hs)):
                phi_ie[i,j],ro[i,j],rb[i,j],phi_o0[i,j] = fsolve(f,[20,1.3,0.03,0.003],args=(hs[i],phi_os[0],Vdisp[j],Vr,Thickness[0],phi_i0),xtol=1e-12,maxfev=10000)

                A_radial[i,j],A_flank[i,j]=leakage_area(F,hs[i],delta,rb[i,j],phi_ie[i,j],phi_is,phi_i0)
                Aleakage_total[i,j]=A_radial[i,j]+A_flank[i,j]

                Vdisp_Calc[i,j]=-2*pi*hs[i]*rb[i,j]*ro[i,j]*(3*pi-2*phi_ie[i,j]+phi_i0+phi_o0[i,j])
    else:
        raise ValueError('wrong key input. acceptable keys are phi_is, Vdisp and Thickness')

    return phi_ie,ro,rb,phi_o0,Aleakage_total,Vdisp_Calc



    # loop to calculate the compressor prifile and total leakage area for various scroll heights



def param_dimen(hs,rb,Thickness,phi_ie,phi_is,Vdisp,key='phi_is'):

    if key=='phi_os':
        n_column=len(phi_os)
    elif key=='Thickness':
        n_column=len(Thickness)
    elif key=='Vdisp':
        n_column=len(Vdisp)
    else:
        raise ValueError('wrong key input. acceptable keys are phi_os and Thickness')

    """ This function calculates the dimensionless paramters for each scroll profile """
    # initiazlizing the arays
    p=np.zeros((len(hs),n_column))
    gw=np.zeros((len(hs),n_column))
    gc=np.zeros((len(hs),n_column))
    do=np.zeros((len(hs),n_column))

    if key=='phi_os':
        n_column=len(phi_os)

        for j in range(n_column):
            for i in range(len(hs)):
                p[i,j]=pitch(rb[i,j])
                gw[i,j]=G_w(hs[i],Thickness)
                gc[i,j]=G_c(hs[i],rb[i,j],Thickness)
                do[i,j]=D_o(rb[i,j],phi_ie[i,j])
    elif key=='Thickness':

        for j in range(n_column):
            for i in range(len(hs)):
                p[i,j]=pitch(rb[i,j])
                gw[i,j]=G_w(hs[i],Thickness[j])
                gc[i,j]=G_c(hs[i],rb[i,j],Thickness[j])
                do[i,j]=D_o(rb[i,j],phi_ie[i,j])

    elif key=='Vdisp':

        for j in range(n_column):
            for i in range(len(hs)):
                p[i,j]=pitch(rb[i,j])
                gw[i,j]=G_w(hs[i],Thickness)
                gc[i,j]=G_c(hs[i],rb[i,j],Thickness)
                do[i,j]=D_o(rb[i,j],phi_ie[i,j])
    else:
        raise ValueError('wrong key input. acceptable keys are phi_os and Thickness')


    return gw,gc,do


def stacked_plot(n_plots,n_curves,ydata,xdata,ylabel,xlabel,limits,legend_param,key,frameon=True,legend=True):

    """ This function can be used to plot stacked pltos in same figure with same x axis and multiple yaxis
    n_plots= number of plots to be stacked over each other
    n_curve= number of curves to be added in each plot
    ydata=a list array to be plotted on the y axis
    x_data= x azis data. common for all plots
    y_label= list of y-labels
    xlabel=x axis label
    limits=points to highlight the min/max point of a paramter
    legend_param= the ariable changing within a single plot for multiple lines"""



    fig, (ax) = plt.subplots(nrows=n_plots, sharex=True, subplot_kw=dict(frameon=True)) # frameon=False removes frames

    plt.subplots_adjust(hspace=.0) # gap between subplots
    # plot specifications
    #marker = itertools.cycle((',', '+', '.', 'o', '*'))
    #color = itertools.cycle(('b', 'r', 'm', 'k', 'c'))
    #linestyle = itertools.cycle(('-', '--', '-.',':'))

    marker = itertools.cycle((',', '+', '.'))
    color = itertools.cycle(('b', 'r', 'm'))
    linestyle = itertools.cycle(('-', '--', '-.'))

    for i in range(n_plots):
        ydata1=ydata[i]

        for j in range(n_curves):
            ydata2=[row[j] for row in ydata1]

            if key=='phi_os':
                label=key+'='+str(round(legend_param[j],2))
            elif key=='Thickness':
                label=key+'='+str(round(legend_param[j]*1000,2))+'mm'
            elif key=='Vdisp':
                label=str(legend_param[j])+'ton'
            else:
                raise ValueError('wrong key input. acceptable keys are phi_os and Thickness')


            #ax[i].grid()
            ax[i].plot(xdata, ydata2, color=next(color), linestyle=next(linestyle),label=label)


        wrap_label_string = '\n'.join(textwrap.wrap(ylabel[i],15))
        ax[i].set_ylabel(wrap_label_string)
        if i<3:
            ax[i].axhline(y=limits[i],color='c',linestyle='-.')
        else:
            pass

        if legend==True:

            ax[i].legend()
        else:
            pass

    plt.xlabel(xlabel)
    plt.show()

def leakage_area(F,hs,delta,rb,phi_ie,phi_is,phi_i0):

    """ This function calculates the radial and flank leakage areas """
    A_radial=2*delta*rb*((phi_ie-pi)**2/2-(phi_is+pi)**2/2-phi_i0*(phi_ie-phi_is-2*pi))*1e6
    Nfl=(phi_ie-phi_is)/(2*pi)
    A_flank=2*Nfl*delta*hs*1e6*F
    return A_radial,A_flank

def D_r(do,Qc,Vdis,Te,Tc,Ref,f,hs,eta_a):
    P=[0]*len(Qc)
    P1=[0]*len(Qc)
    D_limit=[0]*len(Qc)
    for i in range(len(Qc)):
        P[i]=guess_comp_P(Vdis[i],Te=Te,Tc=Tc,Ref=Ref,f=60, eta_a=1)
        #P1[i]=P_motor(Qc[i],COP)
        D_limit[i]=motor_size(P[i])*1000
    do1=np.zeros((len(hs),len(Qc)))
    for j in range(len(Qc)):
        for i in range(len(hs)):
            do1[i,j]=do[i,j]/D_limit[j]

    return do1

def find_indices(lst, condition):
    return [i for i, elem in enumerate(lst) if condition(elem)]


def intersection(*arg):
    lst=[0]*len(arg)
    for i in range(len(arg)):
        lst[i]=set(arg[i])
        if i>0:
            inters=inters & lst[i]
        else:
            inters=lst[i]
    return list(inters)

def intersection1(*arg):
    lst=[0]*len(arg)
    for i in range(len(arg)):
        lst[i]=arg[i]
        if i>0:
            inters=[j for j in inters if j in lst[i]]
        else:
            inters=lst[i]
    return list(inters)

def P_ideal(Qc,COP):
    P=Qc/COP
    return P

def P_motor(Qc,COP):
    """ This function returns an estimate of the motar power based on cycle COP.
    This can be used for selcting the motor size"""
    P=P_ideal(Qc,COP)
    P=kW_hp(P)
    return P

def Analyze_Vdisp(hs,rb,ro,phi_ie,phi_i0,phi_o0):
    Vdisp=-2*pi*hs*rb*ro*(3*pi-2*phi_ie+phi_i0+phi_o0)
    return Vdisp

def motor_size(P):
    if P<6:
        d=0.149
    elif P<10 and P>6:
        d=0.17018
    elif P>10 and P<23:
        d=0.179
    elif P>13.33 and P<26.66:
        d=0.19
    elif P>10 and P<40:
        d=0.199
    elif P>16.668 and P<93.34:
        d=0.223
    elif P>60 and P<120:
        d=0.25717
    else:
        d=0.25717
        warnings.warn('Motor size not found for the required motor power')
    return d



