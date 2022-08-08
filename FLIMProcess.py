
import tkinter as tk
from tkinter  import filedialog
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import imagesc as imagesc 
import math
sys.path.append('E:\V03-04_EXE\V03-04_EXE\readcsv_opt.py')
import readcsv_opt
sys.path.append('E:\V03-04_EXE\V03-04_EXE\Li_2dSmooth.py')
import Li_2dSmooth
window=tk.Tk
file_path=filedialog.askopenfilename()
print(file_path)
f=open(file_path,'r')

Para_data=[]
count=1
for line in f:
    if count%2==0:
        Para_data.append(line.split())
    count+=1
f.close
#print(Para_data)
iRef=int(Para_data[0][0])
iExp_start=int(Para_data[1][0])
iExp_end=int(Para_data[2][0])
crop_x=[int(Para_data[3][0]),int(Para_data[3][1])]
crop_y=[int(Para_data[4][0]),int(Para_data[4][1])]
Rsq_criterion =float(Para_data[5][0])
flag_Rsq=int(Para_data[6][0])
flag_tauPhi=int(Para_data[7][0])
flag_Phasor=int(Para_data[8][0])
flag_O2=int(Para_data[9][0])
flag_O2grad=int(Para_data[10][0])
#print(iRef,iExp_start,iExp_end,crop_x,crop_y,Rsq_criterion)
typeProbe=int(Para_data[11][0])
if typeProbe==1:
    tauN2 = 532.84e-9;     # lifetime of N2 of RTDP (unit: sec)
    kquench = 2.0658;      # reaction constant of RTDP
elif typeProbe==2:    
    tauN2 = 2387e-9;       # lifetime of N2 of TRC (unit: sec)
    kquench = 5.2835;      # reaction constant of TRC
elif typeProbe== 3:#other type
        tauN2 = Para_data[12][0]     # lifetime of N2 of other dye (unit: sec)
        kquench = Para_data[13][0]   # reaction constant of other dye
Lscalebar = 20;            # unit: um
scalebar_edge = 3;         # unit: um
spansmooth = 25;    #the span for the moving average (smooth the O2 gradient data)
#Output sata
OptFile='E:\V03-04_EXE\V03-04_EXE\Output'
optfilelist=os.walk(OptFile)
OptFile_Name=[]
for dirname,subdir,files in optfilelist:
    OptFile_Name.append(subdir)
print('OptFile_Name:',OptFile_Name)   
#reference
OptRef_file=OptFile+"\\"+OptFile_Name[0][0]# E:\V03-04_EXE\V03-04_EXE\Output\00_REF
print('OptRef_file:',OptRef_file)
OptRef_filelist=os.walk(OptRef_file)
OptRef_filename=[]
for dirname,subdir,files in OptRef_filelist:
    OptRef_filename.append(files)

#experience data
OptExp_file=OptFile+"\\"+ OptFile_Name[0][1]# E:\V03-04_EXE\V03-04_EXE\Output\01_EXP
print('OptExp_file:',OptExp_file)
OptExp_filelist=os.walk(OptExp_file)
OptExp_filename=[]
for dirname,subdir,files in OptExp_filelist:
    OptExp_filename.append(files)
print('OptExp_filename:',OptExp_filename)

OptPost_file=OptFile_Name[3]

OptPlot_file=OptFile_Name[5]
#print("OptRef_file:",OptFile_Name)
ref_post_filenm=OptRef_file+'\\'+OptRef_filename[0][iRef-1]
print(ref_post_filenm)

ref_data=readcsv_opt.Readcsv_opt(ref_post_filenm)
print (ref_data)
nphase=ref_data[0]
fmod=ref_data[1]
tauRef=ref_data[2]
coef_ref=ref_data[3][crop_y[0]-1:crop_y[1],crop_x[0]-1:crop_x[1],:]
Mref=coef_ref[:,:,0]/coef_ref[:,:,1]
print(Mref)

exp_post_filenm= OptExp_file +'\\'+ OptExp_filename[0][iExp_start-1]

flag_file=exp_post_filenm[-13:-11]
print(flag_file,type(flag_file))

if flag_file=='05':
    pixel_size=1.5759
elif flag_file=='25':
    pixel_size=3.1834
elif flag_file=='10':
    pixel_size=0.7857
elif flag_file =='20':
    pixel_size=0.3978
elif flag_file=='40':
    pixel_size=0.1971
    print('a')
elif flag_file=='63':
    pixel_size=0.1256
print(pixel_size)
nrow=len(coef_ref)
ncol=len(coef_ref[0])  
print (nrow,ncol)
# question
xaxis=np.arange(1,nrow+1,1)*pixel_size
yaxis=np.arange(1,ncol+1,1)*pixel_size
max_xaxis=math.floor(max(xaxis))
for iExp in range(iExp_start,iExp_end+1):
    print("The program is going to run the experiment No.%d"%(iExp))
    flag_temp=0;

    exp_post_filenm= OptExp_file +'\\'+ OptExp_filename[0][iExp-1]

    exp_data=readcsv_opt.Readcsv_opt(exp_post_filenm)
    coef_exp=exp_data[3][crop_y[0]-1:crop_y[1],crop_x[0]-1:crop_x[1],:]
    Mexp=coef_exp[:,:,0]/coef_exp[:,:,1]
    print(Mexp)
    #check R square fitting results
    IM=np.ones((nrow,ncol))
    Rsq_ref=coef_ref[:,:,3]
    nongoodfit_ind=np.nonzero(Rsq_ref<Rsq_criterion)
    IM[nongoodfit_ind]=0
    Rsq_exp = coef_exp[:,:,3]
    nongoodfit_ind=np.nonzero(Rsq_exp<Rsq_criterion)
    IM[nongoodfit_ind]=0
    
    if flag_Rsq==1:
        flag_temp=flag_temp+flag_Rsq;
        #if flag_temp==1
            #mkdir(OptPlot_file,OptExp_filelist(iExp+2).name(1:end-4))
        #reference
        Rsq_ref= np.rot90(Rsq_ref,1)
        plt.imshow(Rsq_ref, extent=[crop_x[0], crop_x[1], crop_y[0],crop_y[1]],cmap=plt.get_cmap('jet'),vmax=1,vmin=0.95)
        plt.xlabel('pixel') 
        plt.ylabel('pixel')
        plt.title('R Square Value(Ref)')
        plt.colorbar()
        plt.show()
        
        plt.hist(Rsq_ref.compress((Rsq_ref!=0).flat),range=(0.9,1),bins=20)
        plt.xlabel('R Square Value'); 
        plt.ylabel('Counts'); 
        plt.title('Histogram of R Square Value (Ref)');
        plt.show()
        #experiment
        Rsq_exp= np.rot90(Rsq_ref,1)
        plt.imshow(Rsq_exp, extent=[crop_x[0], crop_x[1], crop_y[0],crop_y[1]],cmap=plt.get_cmap('jet'),vmax=1,vmin=0.95)
        plt.xlabel('pixel') 
        plt.ylabel('pixel')
        plt.title('R Square Value(Exp)')
        plt.colorbar()
        plt.show()
        
        plt.hist(Rsq_exp.compress((Rsq_ref!=0).flat),range=(0.9,1),bins=20)
        plt.xlabel('R Square Value'); 
        plt.ylabel('Counts'); 
        plt.title('Histogram of R Square Value (Ref)');
        plt.show()

        #Life time calculated from phase angle (tauPhi)
        tauPhi=np.tan((-coef_exp[:,:,2]+coef_ref[:,:,2])+math.atan(2*math.pi*fmod*tauRef))/(2*math.pi*fmod)
        tauPhi=tauPhi*IM
        tauPhi=np.rot90(tauPhi,-1)
        tauPhi_plot=Li_2dSmooth.Li2dSmooth(Li_2dSmooth.Li2dSmooth(tauPhi))

        if flag_tauPhi == 1:
            flag_temp = flag_temp + flag_tauPhi;
            #if flag_temp == 1:
                #mkdir(OptPlot_file,OptExp_filelist(iExp+2).name(1:end-4))
            if typeProbe==1:
                plt.imshow(tauPhi_plot*1e9,extent=[crop_x[0], crop_x[1], crop_y[0],crop_y[1]],cmap=plt.get_cmap('jet'),vmax=500,vmin=100)
            elif typeProbe==2:
                plt.imshow(tauPhi_plot*1e6,extent=[crop_x[0], crop_x[1], crop_y[0],crop_y[1]],cmap=plt.get_cmap('jet'),vmax=2.5,vmin=1)
            elif typeProbe==3:
                if math.floor(math.log10(tauN2*1e9))+1 >3:
                    plt.imshow(tauPhi_plot*1e6,extent=[crop_x[0], crop_x[1], crop_y[0],crop_y[1]],cmap=plt.get_cmap('jet'))
                else:
                    plt.imshow(tauPhi_plot*1e9,extent=[crop_x[0], crop_x[1], crop_y[0],crop_y[1]],cmap=plt.get_cmap('jet'))
        plt.xlabel('pixel')
        plt.ylabel('pixel')
        plt.title('Phase Lifetime(ns)')
        plt.colorbar()
        plt.show()






        
         

        
