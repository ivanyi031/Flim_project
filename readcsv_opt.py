import pandas as pd
import numpy as np
def Readcsv_opt(filenm):
    #print(filenm)
    df=pd.read_csv(filenm,nrows=3,header=None)
    #print(df)
    data=df.values
    
    data=list(map(list,zip(*data)))
    
   
    nphase=float(data[1][0])
    fmod=float(data[1][1])
    tauRef=float(data[1][2])

    #print(type(nphase),type(fmod),type(tauRef))
    Data=pd.read_csv(filenm,skiprows=4,low_memory=False)
    coef=np.zeros((1008,1008,4))
    #print(Data)
    #print(type(Data.iloc[0,6]),float(Data.iloc[0,6]))


    #coef[Data.iloc(0,0)-1,Data.iloc(0,1)-1,0]=Data.iloc(0,2)
    #print(coef[Data.iloc(0,0)-1,Data.iloc(0,1)-1,0])

    # for idata in range(1008*1008):
    #     coef[Data.iloc[idata,0]-1,Data.iloc[idata,1]-1,0]=float(Data.iloc[idata,2])
    #     coef[Data.iloc[idata,0]-1,Data.iloc[idata,1]-1,1]=float(Data.iloc[idata,3])
    #     coef[Data.iloc[idata,0]-1,Data.iloc[idata,1]-1,2]=float(Data.iloc[idata,4])
    #     if (Data.iloc[idata,6])=='-1.#INF00':
    #         coef[Data.iloc[idata,0]-1,Data.iloc[idata,1]-1,3]=0
    #     else:
    #         coef[Data.iloc[idata,0]-1,Data.iloc[idata,1]-1,3]=float(Data.iloc[idata,6])
    Data[Data['R-square'].str.contains("INF")]=0
    coef[:,:,0]=np.transpose(np.reshape(np.array(Data['c1']),(1008,1008)))
    coef[:,:,1]=np.transpose(np.reshape(np.array(Data['c2']),(1008,1008)))
    coef[:,:,2]=np.transpose(np.reshape(np.array(Data['phi']),(1008,1008)))
    coef[:,:,3]=np.transpose(np.reshape(np.array(Data['R-square']),(1008,1008)))   
    datas=[nphase,fmod,tauRef,coef]
    return datas

