import numpy as np

def Li2dSmooth(img):
    nrow=len(img[0])
    ncol=len(img)
    smimg=np.zeros((ncol,nrow))
    #middle section
    for irow in range(1,nrow-1):
        for icol in range(1,ncol-1):
            smimg[icol,irow]=\
            (img[icol-1,irow-1]+img[icol-1,irow+1]+img[icol+1,irow-1]+img[icol+1,irow+1]\
            +2*(img[icol-1,irow]+img[icol+1,irow]+img[icol,irow-1]+img[icol,irow+1])\
            +4*img[icol,irow])/16
    # edge section
    for icol in range(1,ncol-1):
        smimg[icol,0]=\
            (img[icol-1,1]+img[icol+1,1]+2*img[icol,1]+\
                3*(img[icol-1,0]+img[icol+1,0]+6*img[icol,0]))/16
        smimg[icol,nrow-1]=\
            (img[icol-1,nrow-2]+img[icol+1,nrow-2]+2*img[icol,nrow-2]+\
                3*(img[icol-1,nrow-1]+img[icol+1,nrow-1])+6*img[icol,nrow-1])/16
        
    for irow in range(1,nrow-1):
        smimg[0,irow]=\
            (img[1,irow-1]+img[1,irow+1]+2*img[1,irow]+\
                3*(img[0,irow-1]+img[0,irow+1])+6*img[0,irow])/16
        smimg[ncol-1,irow]=\
            (img[ncol-2,irow-1]+img[ncol-2,irow+1]+2*img[ncol-2,irow]+\
                3*(img[ncol-1,irow-1]+img[ncol-1,irow+1])+6*img[ncol-1,irow])/16
    #corner section
    smimg[0,0]=(9*img[0,0]+3*(img[1,0]+img[0,1])+img[1,1])/16
    smimg[ncol-1,0]=(9*img[ncol-1,0]+3*(img[ncol-2,0]+img[ncol-1,1]+img[ncol-2,1]))/16
    smimg[0,nrow-1]=(9*img[0,nrow-1]+3*(img[0,nrow-2]+img[1,nrow-1]+img[1,nrow-2]))/16
    smimg[ncol-1,nrow-1]=(9*img[ncol-1,nrow-1]+3*(img[ncol-2,nrow-1]+img[ncol-1,nrow-2])+img[ncol-2,nrow-2])/16

    return smimg
