import pandas as pd
import numpy as np
def Readcsv_opt(filenm):
    print(filenm)
    df=pd.read_csv(filenm,nrows=3,header=None)
    #print(df)
    data=df.values
    
    data=list(map(list,zip(*data)))
    
   
    nphase=data[1][0]
    fmod=data[1][1]
    tauRef=data[1][2]

    #print(nphase,fmod,tauRef)
    Data=pd.read_csv(filenm,skiprows=4,low_memory=False)
    coef=np.zeros((1008,1008,4))
    #print(Data)
    #print(type(Data.iloc[0,6]),float(Data.iloc[0,6]))


    #coef[Data.iloc(0,0)-1,Data.iloc(0,1)-1,0]=Data.iloc(0,2)
    #print(coef[Data.iloc(0,0)-1,Data.iloc(0,1)-1,0])

    for idata in range(1008*1008):
        coef[Data.iloc[idata,0]-1,Data.iloc[idata,1]-1,0]=float(Data.iloc[idata,2])
        coef[Data.iloc[idata,0]-1,Data.iloc[idata,1]-1,1]=float(Data.iloc[idata,3])
        coef[Data.iloc[idata,0]-1,Data.iloc[idata,1]-1,2]=float(Data.iloc[idata,4])
        if (Data.iloc[idata,6])=='-1.#INF00':
            coef[Data.iloc[idata,0]-1,Data.iloc[idata,1]-1,3]=0
        else:
            coef[Data.iloc[idata,0]-1,Data.iloc[idata,1]-1,3]=float(Data.iloc[idata,6])
        
    datas=[nphase,fmod,tauRef,coef]
    return datas

