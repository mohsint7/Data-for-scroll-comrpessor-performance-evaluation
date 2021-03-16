"""
This file is reads the .h5 files named in a specific formate and plot some basic
compressor perfroamcne graphs. It is the part of the electronic ennex for the work
presenterd in Tanveer and Bradshaw (2021) "Performance evaluation of low-GWP refrigerants
in 1-100 ton scroll compressors"
presented in Int. J. of Ref.


It requires an excel file as an input that contains the scroll geometry. The required file
can be generated using the file named  "Scroll_comp_scalling".

"""






import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import h5py
import itertools
import os

def Te_Tc():
    """ Define the compressor operating conditions here"""
    Te=[-20,-10,0,7,20]   # evaporator temperature
    Tc=[25,40,54,68]   # Condenser temperature
    return Te,Tc




def h5toxlsx(path="/"):
    """ This function reads h5 files, removes data points that are not available in the data sheet and writes the results to excell shett. It need path of the files where .h5 files are located for a particular case"""

    data,data_key=ReadingVar(Filepath=path,Filename='geom1',sheet_name='Sheet1')
    Qc_t=data['Qc_t']


    #Vdis=[1.743,2.178,2.614,3.049,3.485,3.921,4.356,4.792,5.228,5.663,6.97,6.099,6.534]

    # Initialzing the arrays
    Wdot=np.zeros((len(Qc_t)))  # compressor power
    mdot=np.zeros((len(Qc_t))) # mass flow rate
    eta_oi=np.zeros((len(Qc_t)))  # overall isentropic efficiency
    Wdoti=np.zeros((len(Qc_t)))   # isentropic work
    dh=np.zeros((len(Qc_t)))     # isentropic enthalpy chage
    eta_v=np.zeros((len(Qc_t)))    # volumetric efficiency
    Vdis1=np.zeros((len(Qc_t)))
    Vdis2=np.zeros((len(Qc_t)))
    Wdot_mech=np.zeros((len(Qc_t)))
    Wdot_pv=np.zeros((len(Qc_t)))
    Wdoti=np.zeros((len(Qc_t)))
    dh=np.zeros((len(Qc_t)))
    T_outlet=np.zeros((len(Qc_t)))
    P_outlet=np.zeros((len(Qc_t)))
    losses_t=np.zeros((len(Qc_t)))
    eta_m=np.zeros((len(Qc_t)))

    for i in range(len(Qc_t)):

            Vdis1[i]=str(Qc_t[i])
            FileName=path+'Qc_t'+'_'+str(int(Vdis1[i]))+'.h5'
            print(FileName)
            hdf=h5py.File(FileName, 'r')   # reading the .h5 file

            #ls=list(hdf.keys())
            #print('list of keys in this file: \n', ls)

            Wdot[i]=np.array(hdf.get('Wdot_electrical'))
            mdot[i]=np.array(hdf.get('mdot'))
            eta_oi[i]=np.array(hdf.get('eta_oi'))
            Wdoti[i]=np.array(hdf.get('Wdot_i'))
            dh[i]=Wdoti[i]/(mdot[i])
            eta_v[i]=np.array(hdf.get('eta_v'))
            Vdis2[i]=np.array(hdf.get('Vdisp'))
            Wdot_mech[i]=np.array(hdf.get('Wdot_mechanical'))
            Wdot_pv[i]=np.array(hdf.get('Wdot_pv'))
            Wdoti[i]=np.array(hdf.get('Wdot_i'))
            dh[i]=Wdoti[i]/(mdot[i])

            outlet_state= hdf.get('outlet_state')
            #ls=list(outlet_state.keys())
            #print('list of keys in this file: \n', ls)
            T_outlet[i]=np.array(outlet_state.get('T'))
            P_outlet[i]=np.array(outlet_state.get('p'))

            losses= hdf.get('losses')
            #ls=list(losses.keys())
            #print('list of keys in this file: \n', ls)
            losses_t[i]=np.array(losses.get('bearings'))
            eta_m[i]=(Wdot_mech[i]-losses_t[i])/Wdot_mech[i]

    mdot=mdot*1000      #kg/s to g/s
    eta_oi=eta_oi*100   #%
    eta_v=eta_v*100     #%

    # Te and Tc for column and row indexing in excell sheet

    data={'Qc_t':Qc_t,'Vdisp':Vdis2,'mdot':mdot,'Wdot':Wdot,'eta_oi':eta_oi,'eta_v':eta_v,'Wdot_mech':Wdot_mech,'Wdot_pv':Wdot_pv,'T_outlet':T_outlet,'P_outlet':P_outlet,'dh_is':dh,'eta_m':eta_m}

    # writing data to excll sheets
    df=pd.DataFrame(data)
    df.to_excel(excel_writer = path+"/results.xlsx", sheet_name='Sheet1',startcol=0)




def ReadingVar(Filepath='',Filename='Wdot',sheet_name='Sheet1'):
    """ A function to read .xlsx files and removing unnamed columns. The output form of Var is a list while column key is a numpy array  """
    Filename=Filepath+Filename +'.xlsx'
    Var = pd.read_excel(Filename, sheet_name=sheet_name)
    Var = Var.loc[:, ~Var.columns.str.contains('^Unnamed')] # removing unnamed columns from the list
    column_key=np.array(Var.columns)
    row_key=np.array(Var.index)

    return Var,column_key





def perity_plot(PDSim,PDSim_list,Data,Data_list,key):

    """ This function can be called to generate perity plot """

    # data for error lines
    no_error=np.linspace(0,np.max(np.max(Data))+10,5)
    error_p=no_error+0.1*no_error
    error_n=no_error-0.1*no_error


    label=['']*len(PDSim_list) # initializing array for labels
    marker = itertools.cycle((',', '+', '.', 'o', '*'))

    for i in range(len(PDSim_list)):
        label[i]="$T_{evap}$="+PDSim_list[i]+'C'
        plt.scatter(Data[Data_list[i]].dropna(),PDSim[PDSim_list[i]].dropna(),marker = marker.next(),label=label[i])
    plt.hold

    #axis labels
    if key=='Wdot':
        plt.xlabel('Power(Data Sheet) [kW]')
        plt.ylabel('Power(PDSim) [kW]')
    elif key=='mdot':
        plt.xlabel('Mass flow rate(Data Sheet) [g/s]')
        plt.ylabel('Mass flow rate(PDSim) [g/s]')
    elif key=='eta_is':
        plt.xlabel('Isentropic Enfficiency(Data Sheet) [%]')
        plt.ylabel('Isentropic Enfficiency(PDSim) [%]')
    else:
        raise ValueError('wrong input for key. valid inputs are Wdot, mdot,eta_is')

    # error lines
    plt.plot(no_error,error_n,'--',label='-10%')
    plt.plot(no_error,error_p,'--',label='+10%')
    plt.plot(no_error,no_error)
    plt.legend()
    plt.show()

def Analysis(Filepath='',Filename_PDSim='Wdot',Filename_Data='/Comp2_data',sheet_Data='Wdot',key='Wdot'):

    """ This function calculates MAE and generate perity plot """
    # generating variables and lists
    PDSim,PDSim_list=ReadingVar(Filepath=Filepath,Filename=Filename_PDSim,sheet_name='Sheet1')
    Data,Data_list=ReadingVar(Filepath=Filepath,Filename=Filename_Data,sheet_name=sheet_Data)
    # Calculation of MAE
    error=abs(PDSim-Data)
    error=error/Data
    N=sum(error.count())
    MAE=error.sum(numeric_only=True).sum()*100/N
    print('MAE for'+ '  '+ key+'='+ str(MAE))
    # Perity plot
    perity_plot(PDSim,PDSim_list,Data,Data_list,key)
    plt.savefig(key)



def eta_isen_Data(Filepath='',Filename_Data='Comp1_data',Filename_PDSim='/dh'):

    """ This function generate an excell sheet containing data for actual isentropic efficiency i.e. calcuated from data sheet """
    dh,dh_list=ReadingVar(Filepath=Filepath,Filename=Filename_PDSim,sheet_name='Sheet1')

    Wdot,Wdot_list=ReadingVar(Filepath=Filepath,Filename=Filename_Data,sheet_name='Wdot')
    mdot,mdot_list=ReadingVar(Filepath=Filepath,Filename=Filename_Data,sheet_name='mdot')
    mdot=mdot/1000      #g/s to kg/s
    eta_is=mdot*dh*100/Wdot

    df=pd.DataFrame(eta_is)
    df.to_excel(excel_writer = Filepath+"/eta_is_data.xlsx", sheet_name='Sheet1',startcol=0)

def eta_vol_plot(Filepath='',Filename='eta_v'):
    PDSim,PDSim_list=ReadingVar(Filepath=Filepath,Filename=Filename,sheet_name='Sheet1')
    Tc=[[25,40],[25,40,54],[25,40,54,68],[25,40,54,68],[40,54]]

    label=['']*len(PDSim_list) # initializing array for labels
    marker = itertools.cycle((',', '+', '.', 'o', '*'))

    for i in range(len(PDSim_list)):
        label[i]="$T_{evap}$="+PDSim_list[i]+'C'
        plt.scatter(Tc[i],PDSim[PDSim_list[i]].dropna(),marker = marker.next(),label=label[i])

    plt.legend()
    plt.xlabel('$T_{cond}$ [C]')
    plt.ylabel('Volumetric efficiency [%]')
    plt.show()
    plt.savefig('eta_v')






## Execution
# generating excell files
Filepath1='C:\Users\Mohsin\OneDrive - Oklahoma A and M System\Documents\Phd\Scroll comp scalling\cases_final\case 8'

Filepath2='E:\OneDrive - Oklahoma A and M System\Documents\Phd\Scroll comp scalling\cases_final\case 8'

if os.path.isfile(Filepath1)==True:
    Filepath=Filepath1
else:
    Filepath=Filepath2

# Filepath=Filepath1
h5toxlsx(path=Filepath+'/')



Filename='/results'


data,index_list=ReadingVar(Filepath=Filepath,Filename=Filename,sheet_name='Sheet1')
for i in range(len(index_list)):
    vars()[index_list[i]]=data[index_list[i]].to_numpy().tolist()

plt.figure()
plt.plot(Qc_t,mdot)
plt.xlabel('Capacity [ton]')
plt.ylabel('Mass flow rate [g/s]')
plt.show()
plt.savefig('Mdot')

plt.figure()
plt.plot(Qc_t,Wdot)
plt.xlabel('Capacity [ton]')
plt.ylabel('Power [kW]')
plt.show()
plt.savefig('Wdot')

plt.figure()
plt.plot(Qc_t,eta_v,label='Volumetric')
plt.plot(Qc_t,eta_oi,label='Isentropic')
plt.xlabel('Capacity [ton]')
plt.ylabel('Efficeincy [%]')
plt.legend()
plt.show()
plt.savefig('eta')

plt.figure()
plt.plot(Qc_t,T_outlet)
plt.xlabel('Capacity [ton]')
plt.ylabel('Outlet temperature [K]')
plt.legend()
plt.show()
plt.savefig('T_out')

plt.figure()
plt.plot(Qc_t,P_outlet)
plt.xlabel('Capacity [ton]')
plt.ylabel('Outlet pressure [kPa]')
plt.legend()
plt.show()
plt.savefig('P_out')






