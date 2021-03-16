"""
This file contains method to design the journal bearing based on the maximum oil
film pressure and minimum oil film thickness. It is the part of the electronic ennex
for the work presenterd in Tanveer and Bradshaw (2021) "Performance evaluation of low-GWP
refrigerants in 1-100 ton scroll compressors" presented in Int. J. of Ref.


It requires an excel file as an input that contains the scroll geometry. The required file
can be generated using the file named  "Scroll_comp_scalling".

It requires the .h5 files named in specific formate for each cooling capacity listed the geometric
file. The required files can be generated using the file named "scroll_comp_final"

It will produce the two excel files that contains the bearing forces and the bearing sizes.

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


from PDSim.core.bearings import journal_bearing_design
def h5toxlsx(path="/"):
    """ This function reads h5 files, removes data points that are not available in the data sheet and writes the results to excell shett. It need path of the files where .h5 files are located for a particular case"""

    data,data_key=ReadingVar(Filepath=path,Filename='geom',sheet_name='Sheet1')
    Qc_t=data['Qc_t']


    #Vdis=[1.743,2.178,2.614,3.049,3.485,3.921,4.356,4.792,5.228,5.663,6.97,6.099,6.534]

    # Initialzing the arrays
    Fr=np.zeros((len(Qc_t)))  # compressor power
    inertial=np.zeros((len(Qc_t))) # mass flow rate
    Ft=np.zeros((len(Qc_t)))  # overall isentropic efficiency
    Fz=np.zeros((len(Qc_t)))   # isentropic work
    Qc_t1=np.zeros((len(Qc_t)))
    F_bearing=np.zeros((len(Qc_t)))

    for i in range(len(Qc_t)):

            Qc_t1[i]=str(Qc_t[i])
            FileName=path+'Qc_t'+'_'+str(int(Qc_t1[i]))+'.h5'
            print(FileName)
            hdf=h5py.File(FileName, 'r')   # reading the .h5 file
            G_forces = hdf.get('forces')
            # ls=list(hdf.keys())
            # print('list of keys in this file: \n', ls)

            # Fr[i]=max(abs(np.array(G_forces.get('summed_Fr'))))
            # inertial[i]=np.array(G_forces.get('inertial'))
            # Ft[i]=max(abs(np.array(G_forces.get('summed_Ft'))))
            # Fz[i]=max(abs(np.array(G_forces.get('summed_Fz'))))


            Fr[i]=np.array(G_forces.get('mean_Fr'))
            inertial[i]=np.array(G_forces.get('inertial'))
            Ft[i]=np.array(G_forces.get('mean_Ft'))
            Fz[i]=np.array(G_forces.get('mean_Fz'))
            F_bearing[i]=np.sqrt((Fr[i]+inertial[i])**2+Ft[i]**2)
            F_bearing[i]=F_bearing[i]*1000                    #kN to N


    # Te and Tc for column and row indexing in excell sheet

    data={'Qc_t':Qc_t,'Fr_max':Fr,'inertial':inertial,'Ft_max':Ft,'Fz_max':Fz,'F_bearing':F_bearing}

    # writing data to excll sheets
    df=pd.DataFrame(data)
    df.to_excel(excel_writer = path+"/forces.xlsx", sheet_name='Sheet1',startcol=0)


def ReadingVar(Filepath='',Filename='Wdot',sheet_name='Sheet1'):
    """ A function to read .xlsx files and removing unnamed columns. The output form of Var is a list while column key is a numpy array  """
    Filename=Filepath+Filename +'.xlsx'
    Var = pd.read_excel(Filename, sheet_name=sheet_name)
    Var = Var.loc[:, ~Var.columns.str.contains('^Unnamed')] # removing unnamed columns from the list
    column_key=np.array(Var.columns)
    row_key=np.array(Var.index)

    return Var,column_key


print(journal_bearing_design(r_b = 0.02,
                    L = 0.04,
                    design = 'friction',
                    W = 2200,
                    eta_0 = 0.17,
                    omega = 3600/60.0*2*pi
                    ))


def f(x,design,W,eta_0,omega,goal,key):
    r_b=x[0]
    L=2*r_b
    results=journal_bearing_design(r_b = r_b,
                        L = L,
                        design = 'friction',
                        W = W,
                        eta_0 = eta_0,
                        omega = omega
                        )
    h_min=results['h_min']
    Pmax=results['Pmax']
    c=results['c']

    if key=='hmin':
        r1=h_min-goal
    elif key=='Pmax':
        r1=Pmax-goal
    else:
        raise ValueError('only acceptable keys are hmin and Pmax')
    print(r1)
    return [r1]

######################################
###          Execution            ###
#####################################
if __name__=='__main__':

    Filepath=''

    Filename='/forces'

    # reading variables from the excel files and assigning them variables names based on column header
    data,index_list=ReadingVar(Filepath=Filepath,Filename=Filename,sheet_name='Sheet1')
    for i in range(len(index_list)):
        vars()[index_list[i]]=data[index_list[i]].to_numpy().tolist()

    # Initializing variables
    r_ub=np.zeros(len(F_bearing))        # upper bearing radius
    r_lb=np.zeros(len(F_bearing))        # upper bearing radius
    h_min_ub=np.zeros(len(F_bearing))
    Pmax_ub=np.zeros(len(F_bearing))
    c_ub=np.zeros(len(F_bearing))
    h_min_lb=np.zeros(len(F_bearing))
    Pmax_lb=np.zeros(len(F_bearing))

    FOS=1.2       # factor of safety

    ## Designing upper bearing based on maximum oil film pressure
    for i in range(len(F_bearing)):

        design=0.5
        W=F_bearing[i]*FOS
        eta_0=0.02
        omega=3000/60.0*2*pi
        key='Pmax'
        if key=='hmin':
            goal=0.000006
        elif key=='Pmax':
            goal=9000000
        else:
            raise ValueError('only acceptable keys are hmin and Pmax')



        r_ub[i] = fsolve(f,[0.02],args=(design,W,eta_0,omega,goal,key))
        results_ub=journal_bearing_design(r_b = r_ub[i],
                            L = 2*r_ub[i],
                            design = 'friction',
                            W = W,
                            eta_0 = eta_0,
                            omega = omega
                            )
        h_min_ub[i]=results_ub['h_min']
        Pmax_ub[i]=results_ub['Pmax']
        c_ub=results_ub['c']

    ## Designing lower bearing based on minimum oil film thickness
    for i in range(len(F_bearing)):

        design=0.5
        W=F_bearing[i]*FOS
        eta_0=0.02
        omega=3000/60.0*2*pi
        key='hmin'
        if key=='hmin':
            goal=0.000006
        elif key=='Pmax':
            goal=9000000
        else:
            raise ValueError('only acceptable keys are hmin and Pmax')


        r_lb[i] = fsolve(f,[0.02],args=(design,W,eta_0,omega,goal,key))
        results_lb=journal_bearing_design(r_b = r_lb[i],
                            L = 2*r_lb[i],
                            design = 'friction',
                            W = W,
                            eta_0 = eta_0,
                            omega = omega
                            )
        h_min_lb[i]=results_lb['h_min']
        Pmax_lb[i]=results_lb['Pmax']

    data={'Qc_t':Qc_t,'r_ub':r_ub,'r_lb':r_lb,'h_min_ub':h_min_ub,'Pmax_ub':Pmax_ub,'h_min_lb':h_min_lb,'Pmax_lb':Pmax_lb,'c_ub':c_ub}

    df=pd.DataFrame(data)
    df.to_excel(excel_writer = Filepath+"/bearing_size.xlsx", sheet_name='Sheet1',startcol=0)