# %%

import tkinter as tk
from tkinter  import filedialog
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import math
sys.path.append('E:\V03-04_EXE\V03-04_EXE\python\readcsv_opt.py')
import readcsv_opt
sys.path.append('E:\V03-04_EXE\V03-04_EXE\python\Li_2dSmooth.py')
import Li_2dSmooth
sys.path.append('E:\V03-04_EXE\V03-04_EXE\python\DrawCircle.py')
import DrawCircle
import scipy.io
PATH = 'E:\V03-04_EXE\V03-04_EXE\FLIM_colormap.mat'
FLIM_colormap = scipy.io.loadmat(PATH)#load matlab's colormap
def smooth(a,cnt):
    length=len(a)
    yy=np.zeros(length)
    N=int((cnt-1)/2)
    for i in range(0,length):
        if i in range(0,N):
            yy[i]=np.mean(np.array(a[0:2*i+1]))
            yy[-i-1]=np.mean(np.array(a[length-1-2*i:length]))
        else:
            yy[i]=np.mean(np.array(a[i-N:i+N+1]))
    return yy
def accumarray_mean(mat):
    # calculate the mean of nonzero idex in each row
     a=np.nonzero(mat)
     rows=np.shape(mat)[0]
     c=np.bincount(a[0])
     #print(c)
     mean=np.zeros(rows)
     pnt_start=0
     pnt_end=-1
     for i in range(len(c)):
         pnt_end+=c[i]
         mean[i]=np.mean(np.array (mat[a][pnt_start:pnt_end+1]))
         pnt_start=pnt_end+1
     return mean

window=tk.Tk()
file_path=filedialog.askopenfilename()
#print(file_path)
window.destroy()
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
#%%
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
#%%
#filename
#Output data
OptFile='E:\V03-04_EXE\V03-04_EXE\Output'
optfilelist=os.walk(OptFile)
OptFile_Name=[]
for dirname,subdir,files in optfilelist:
    OptFile_Name.append(subdir)
#print('OptFile_Name:',OptFile_Name)   
#reference
OptRef_file=OptFile+"\\"+OptFile_Name[0][0]
# E:\V03-04_EXE\V03-04_EXE\Output\00_REF
print('OptRef_file:',OptRef_file)
OptRef_filelist=os.walk(OptRef_file)
OptRef_filename=[]
for dirname,subdir,files in OptRef_filelist:
    OptRef_filename.append(files)

#experience data
OptExp_file=OptFile+"\\"+ OptFile_Name[0][1]
# E:\V03-04_EXE\V03-04_EXE\Output\01_EXP
#print('OptExp_file:',OptExp_file)
OptExp_filelist=os.walk(OptExp_file)
OptExp_filename=[]
for dirname,subdir,files in OptExp_filelist:
    OptExp_filename.append(files)
#print('OptExp_filename:',OptExp_filename)

#output data
OptPost_file=OptFile+"\\"+ OptFile_Name[0][2]
#OPtPost_file: E:\V03-04_EXE\V03-04_EXE\Output\02_POST
#print('OPtPost_file:',OptPost_file)
OptPlot_file=OptFile+"\\"+ OptFile_Name[0][3]
#OptPlot_file: E:\V03-04_EXE\V03-04_EXE\Output\03_PLOT
path_python='python'
OptPlot_file=OptPlot_file+'\\'+path_python
OptPost_file=OptPost_file+'\\'+path_python


#%%
ref_post_filenm=OptRef_file+'\\'+OptRef_filename[0][iRef-1]

ref_data=readcsv_opt.Readcsv_opt(ref_post_filenm)

nphase=ref_data[0]
fmod=ref_data[1]
tauRef=ref_data[2]
coef_ref=ref_data[3][crop_y[0]-1:crop_y[1],crop_x[0]-1:crop_x[1],:]
Mref=coef_ref[:,:,0]/coef_ref[:,:,1]


exp_post_filenm= OptExp_file +'\\'+ OptExp_filename[0][iExp_start-1]
#%%
#倍率
flag_file=exp_post_filenm[-13:-11]


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
elif flag_file=='63':
    pixel_size=0.1256
#print(pixel_size)
nrow=len(coef_ref)
ncol=len(coef_ref[0])  
#print (nrow,ncol)

xaxis=np.arange(1,nrow+1,1)*pixel_size
yaxis=np.arange(1,ncol+1,1)*pixel_size
max_xaxis=math.floor(max(xaxis))
#%%
for iExp in range(iExp_start,iExp_end+1):
    print("The program is going to run the experiment No.%d"%(iExp))
    flag_temp=0;

    exp_post_filenm= OptExp_file +'\\'+ OptExp_filename[0][iExp-1]

    exp_data=readcsv_opt.Readcsv_opt(exp_post_filenm)
    coef_exp=exp_data[3][crop_y[0]-1:crop_y[1],crop_x[0]-1:crop_x[1],:]
    Mexp=coef_exp[:,:,0]/coef_exp[:,:,1]
    
    #%%
    #0.check R square fitting results
    IM=np.ones((nrow,ncol))
    Rsq_ref=coef_ref[:,:,3]
    nongoodfit_ind=np.nonzero(Rsq_ref<Rsq_criterion)
    IM[nongoodfit_ind]=0
    Rsq_exp = coef_exp[:,:,3]
    nongoodfit_ind=np.nonzero(Rsq_exp<Rsq_criterion)
    IM[nongoodfit_ind]=0
    #%%
    if flag_Rsq==1:
        flag_temp=flag_temp+flag_Rsq;
        if flag_temp==1:
            path=OptPlot_file+'\\'+OptExp_filename[0][iExp-1][0:-4]
            if not os.path.exists(path):
                os.mkdir(path)
        #reference
        Rsq_ref= np.flipud(np.rot90(Rsq_ref,1))
        
        Rsq_ref_10=Rsq_ref[0:10,0:10]
        Rsq_ref_100=Rsq_ref[0:100,0:100]
        Rsq_ref_500=Rsq_ref[0:500,0:500]
        #Rsq_ref_1000=Rsq_ref[0:1000,0:1000]
        
        plt.imshow(Rsq_ref, origin='upper',extent=[xaxis[0], xaxis[-1], yaxis[-1],yaxis[0]],
                   cmap=ListedColormap(FLIM_colormap['flimcolormap_black']),vmax=1,vmin=0.95)
        #plt.xlim([0,max(xaxis)])
        plt.xlabel('pixel') 
        plt.ylabel('pixel')
        plt.title('R Square Value(Ref)')
        plt.colorbar()
        plotfilenm=path+'\\'+OptExp_filename[0][iExp-1][0:-4]+'_Rsq_ref'
        plt.savefig(plotfilenm+'.tif',format='tif')
        plt.savefig(plotfilenm+'.pdf',format='pdf')
        plt.show()
        
        plt.imshow(Rsq_ref_10, origin='upper',extent=[xaxis[0], xaxis[-1], yaxis[-1],yaxis[0]],cmap=ListedColormap(FLIM_colormap['flimcolormap_black']),vmax=1,vmin=0.95)
        #plt.xlim([0,max(xaxis)])
        plt.xlabel('pixel') 
        plt.ylabel('pixel')
        plt.title('R Square Value(Ref_10)')
        plt.colorbar()
        plotfilenm=path+'\\'+OptExp_filename[0][iExp-1][0:-4]+'_Rsq_ref_10'
        plt.savefig(plotfilenm+'.tif',format='tif')
        plt.savefig(plotfilenm+'.pdf',format='pdf')
        plt.show()

        plt.imshow(Rsq_ref_100, origin='upper',extent=[xaxis[0], xaxis[-1], yaxis[-1],yaxis[0]],cmap=ListedColormap(FLIM_colormap['flimcolormap_black']),vmax=1,vmin=0.95)
        #plt.xlim([0,max(xaxis)])
        plt.xlabel('pixel') 
        plt.ylabel('pixel')
        plt.title('R Square Value(Ref_100)')
        plt.colorbar()
        plotfilenm=path+'\\'+OptExp_filename[0][iExp-1][0:-4]+'_Rsq_ref_100'
        plt.savefig(plotfilenm+'.tif',format='tif')
        plt.savefig(plotfilenm+'.pdf',format='pdf')
        plt.show()

        plt.imshow(Rsq_ref_500, origin='upper',extent=[xaxis[0], xaxis[-1], yaxis[-1],yaxis[0]],cmap=ListedColormap(FLIM_colormap['flimcolormap_black']),vmax=1,vmin=0.95)
        #plt.xlim([0,max(xaxis)])
        plt.xlabel('pixel') 
        plt.ylabel('pixel')
        plt.title('R Square Value(Ref_500)')
        plt.colorbar()
        plotfilenm=path+'\\'+OptExp_filename[0][iExp-1][0:-4]+'_Rsq_ref_500'
        plt.savefig(plotfilenm+'.tif',format='tif')
        plt.savefig(plotfilenm+'.pdf',format='pdf')
        plt.show()
        
        #%%
        plt.hist(Rsq_ref[np.nonzero(Rsq_ref>=0.8)],bins=np.arange(0.8,1,0.005))
        plt.xlabel('R Square Value'); 
        plt.ylabel('Counts'); 
        plt.title('Histogram of R Square Value (Ref)');
        plotfilenm=path+'\\'+OptExp_filename[0][iExp-1][0:-4]+'_Rsqhis_ref'
        plt.savefig(plotfilenm+'.tif',format='tif')
        plt.savefig(plotfilenm+'.pdf',format='pdf')
        plt.show()
        #experiment
        #%%
        Rsq_exp= np.flipud(np.rot90(Rsq_exp,1))
        #Rsq_exp_small=Rsq_exp[0:99,0:99]
        plt.imshow(Rsq_exp,origin='upper' ,extent=[xaxis[0], xaxis[-1], yaxis[-1],yaxis[0]],cmap=ListedColormap(FLIM_colormap['flimcolormap_black']),vmax=1,vmin=0.95)
        plt.xlim([0,max(xaxis)])
        plt.xlabel('pixel') 
        plt.ylabel('pixel')
        plt.title('R Square Value(Exp)')
        plt.colorbar()
        plotfilenm=path+'\\'+OptExp_filename[0][iExp-1][0:-4]+'_Rsq_exp'
        plt.savefig(plotfilenm+'.tif',format='tif')
        plt.savefig(plotfilenm+'.pdf',format='pdf')
        plt.show()
        
        plt.hist(Rsq_exp[np.nonzero(Rsq_exp>=0.8)],bins=np.arange(0.8,1,0.005))
        plt.xlabel('R Square Value'); 
        plt.ylabel('Counts'); 
        plt.title('Histogram of R Square Value (EXP)');
        plotfilenm=path+'\\'+OptExp_filename[0][iExp-1][0:-4]+'_Rsqhis_exp'
        plt.savefig(plotfilenm+'.tif',format='tif')
        plt.savefig(plotfilenm+'.pdf',format='pdf')
        plt.show()
        #%%
        #1.Life time calculated from phase angle (tauPhi)
        tauPhi=np.tan((-coef_exp[:,:,2]+coef_ref[:,:,2])+math.atan(2*math.pi*fmod*tauRef))/(2*math.pi*fmod)
        tauPhi=tauPhi*IM
        tauPhi=np.flipud(np.rot90(tauPhi,-1))
        tauPhi_plot=Li_2dSmooth.Li2dSmooth(Li_2dSmooth.Li2dSmooth(tauPhi))

        if flag_tauPhi == 1:
            flag_temp = flag_temp + flag_tauPhi;
            if flag_temp==1:
                path=OptPlot_file+'\\'+OptExp_filename[0][iExp-1][0:-4]
                if not os.path.exists(path):
                    os.mkdir(path)
            if typeProbe==1:
                plt.imshow(tauPhi_plot*1e9,origin='upper',extent=[xaxis[0], xaxis[-1], yaxis[0],yaxis[-1]],
                           cmap=plt.get_cmap('jet'),vmax=500,vmin=100)
            elif typeProbe==2:
                plt.imshow(tauPhi_plot*1e6,origin='upper',extent=[xaxis[0], xaxis[-1], yaxis[0],yaxis[-1]],
                           cmap=plt.get_cmap('jet'),vmax=2.5,vmin=1)
            elif typeProbe==3:
                if math.floor(math.log10(tauN2*1e9))+1 >3:
                    plt.imshow(tauPhi_plot*1e6,origin='upper',extent=[xaxis[0], xaxis[-1], yaxis[0],yaxis[-1]],
                               cmap=plt.get_cmap('jet'))
                else:
                    plt.imshow(tauPhi_plot*1e9,origin='upper',extent=[xaxis[0], xaxis[-1], yaxis[0],yaxis[-1]],
                               cmap=plt.get_cmap('jet'))
        plt.xlabel('pixel')
        plt.ylabel('pixel')
        plt.title('Phase Lifetime(ns)')
        plt.colorbar()
        plt.plot([max_xaxis-scalebar_edge-Lscalebar+1, max_xaxis-scalebar_edge],[Lscalebar*0.9, Lscalebar*0.9],'-k',lw=5)
        textscalebar=str(Lscalebar)+r'$\mu m$'
        plt.text(max_xaxis-scalebar_edge-0.5*Lscalebar, Lscalebar*0.1,textscalebar, horizontalalignment='center',
                 color='k',fontsize=12)
        plotfilenm=path+'\\'+OptExp_filename[0][iExp-1][0:-4]+'-tauPhi_Scale'
        plt.savefig(plotfilenm+'.tif',format='tif')
        plt.savefig(plotfilenm+'.pdf',format='pdf')
        plt.show()
        #%%
        #2.Phasor plot
        M=(Mexp/Mref)*((1+(2*math.pi*fmod*tauRef)**2)**(-1/2))
        M=M*IM
        
        Phi=(-coef_exp[:,:,2]+coef_ref[:,:,2])+math.atan(2*math.pi*fmod*tauRef)
        Phi=Phi*IM#PHI correct
        phasor_x=np.abs(M*np.cos(Phi))
        phasor_x2=np.reshape(phasor_x,(1,1008*1008))
        phasor_y=np.abs(M*np.sin(Phi))
        phasor_y2=np.reshape(phasor_y,(1,1008*1008))
    
        #plot figure
        if flag_Phasor==1:
            flag_temp=flag_temp+flag_Phasor
            if flag_temp==1:
                path=OptPlot_file+'\\'+OptExp_filename[0][iExp-1][0:-4]
                if not os.path.exists(path):
                    os.mkdir(path)
            if typeProbe==1:
                tauO2 = 173.8*1e-9;
                tauAir = tauRef;
                tau = np.array([tauN2, tauO2, tauAir]);
                tmp = (1+(2*math.pi*fmod*tau)**2);
                phasor_x_ref = 1/tmp;
                phasor_y_ref = (2*math.pi*fmod*tau)/tmp;
            else:
                tauAir = tauRef;
                tau = np.array([tauN2, tauAir]);
                tmp = (1+(2*math.pi*fmod*tau)**2);
                phasor_x_ref = 1/tmp;
                phasor_y_ref = (2*math.pi*fmod*tau)/tmp;
            cnt=[np.arange(0,1,7e-3),np.arange(0,0.6,7e-3)]
        
            plt.hist2d(phasor_x2[0],phasor_y2[0],bins=[cnt[0],cnt[1]],cmap =ListedColormap(FLIM_colormap['flimcolormap']))
            DrawCircle.DrawCircle(0.5,0,0.5,128,'k-')
            plt.plot(phasor_x_ref[0],phasor_y_ref[0],'k-x',ms=10,lw=2)
            plt.plot(phasor_x_ref[1],phasor_y_ref[1],'k-x',ms=10,lw=2)
            if typeProbe==1:
                plt.plot(phasor_x_ref[2],phasor_y_ref[2],'k-x',ms=10,lw=2)
            plt.xlim([0,1.0])
            plt.ylim([0,0.6])
            plt.title('Phasor Plot')
            plotfilenm=path+'\\'+OptExp_filename[0][iExp-1][0:-4]+'-Phasor'
            plt.savefig(plotfilenm+'.tif',format='tif')
            plt.savefig(plotfilenm+'.pdf',format='pdf')
            plt.show()
        #%%
        # 3. Oxygen
        sm_tauPhi= Li_2dSmooth.Li2dSmooth(Li_2dSmooth.Li2dSmooth(tauPhi))
        O2exp=(tauN2/sm_tauPhi-1)/kquench
        O2exp_plot=O2exp*100
        O2_zero=np.nonzero(np.rot90(Rsq_ref,-2)<Rsq_criterion)
        O2exp_plot[O2_zero]=0
        O2_zero=np.nonzero(np.rot90(Rsq_exp,-2)<Rsq_criterion)
        O2exp_plot[O2_zero]=0
        if flag_O2==1:
            flag_temp=flag_temp+flag_O2
            if flag_temp==1:
                path=OptPlot_file+'\\'+OptExp_filename[0][iExp-1][0:-4]
                if not os.path.exists(path):
                    os.mkdir(path)
            plt.imshow(O2exp_plot,origin='upper',extent=[xaxis[0], xaxis[-1], yaxis[-1],yaxis[0]],cmap=plt.get_cmap('jet'),vmax=25,vmin=0)
            plt.xlabel('pixel')
            plt.ylabel('pixel')
            plt.title('Oxygen Map(%)')
            plt.colorbar()
            plotfilenm=path+'\\'+OptExp_filename[0][iExp-1][0:-4]+'-O2Map'
            plt.savefig(plotfilenm+'.tif',format='tif')
            plt.savefig(plotfilenm+'.pdf',format='pdf')
            plt.plot([max_xaxis-scalebar_edge-Lscalebar+1, max_xaxis-scalebar_edge],[Lscalebar*0.9, Lscalebar*0.9],'-k',lw=5)
            textscalebar=str(Lscalebar)+r'$\mu m$ '
            plt.text(max_xaxis-scalebar_edge-0.5*Lscalebar, Lscalebar*0.1,textscalebar, horizontalalignment='center',color='k',fontsize=12)
            plotfilenm=path+'\\'+OptExp_filename[0][iExp-1][0:-4]+'-O2Map_Scale'
            plt.savefig(plotfilenm+'.tif',format='tif')
            plt.savefig(plotfilenm+'.pdf',format='pdf')
            plt.show()
        #%%
        #4. Oxygen gradient plot
        # ''
        if flag_O2grad == 1:         
            flag_temp = flag_temp + flag_O2grad;
            if flag_temp==1:
                path=OptPlot_file+'\\'+OptExp_filename[0][iExp-1][0:-4]
                if not os.path.exists(path):
                    os.mkdir(path)
            #O2exp=np.flipud(O2exp)
            O2exp_up=smooth(O2exp[7,:],spansmooth)*100
            O2exp_low=smooth(O2exp[1008-8,:],spansmooth)*100
            O2exp_mean=accumarray_mean(np.transpose(O2exp))
            O2exp_mean=smooth(O2exp_mean,spansmooth)*100
            print(O2exp_mean)
            plt.plot(xaxis,O2exp_up,'b.',ms=3,label='upstream')
            plt.plot(xaxis,O2exp_low,'r.',ms=3,label='downstream')
            plt.plot(xaxis,O2exp_mean,'k.',ms=3,label='average')
            plt.ylim(0,25)
            plt.yticks([0, 5, 10, 15, 20, 25])
            plt.xlabel(r'$\mu m$ ') 
            plt.ylabel('oxygen level(%)')
            plt.title('Oxygen Gradient')
            plt.legend()
            plotfilenm=path+'\\'+OptExp_filename[0][iExp-1][0:-4]+'-O2Gradient'
            plt.savefig(plotfilenm+'.tif',format='tif')
            plt.savefig(plotfilenm+'.pdf',format='pdf')
            plt.show()

            #Save Data
            postfilenm=OptPost_file+'\\'+OptExp_filename[0][iExp-1][0:-4]+'.mat'
            scipy.io.savemat(postfilenm,{'tauPhi':tauPhi,'M':M,'Phi':Phi,'phasor_x':phasor_x,
                                         'phasor_y':phasor_y,'O2exp':O2exp,'Rsq_ref':Rsq_ref,
                                         'Rsq_exp':Rsq_exp,'xaxis':xaxis,'yaxis':yaxis,
                                         'tauPhi_plot':tauPhi_plot,'O2exp_plot':O2exp_plot})

            



        
        


            









        
         

        
