# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 18:37:01 2022

@author: vulpine
"""

import os

#from tensorflow import keras
import numpy as np
#from tensorflow.keras.preprocessing.image import load_img
#from IPython.display import Image, display
import PIL
from skimage.morphology import thin
from matplotlib import pyplot as plt
import statistics
 
# pos=range(6)
it = 6
# for it in pos:  
# the code is for looping through "Pos" folders, 
# I removed this look to make easier to deploy, but 
# this has to be reinstated in the final version
    
"""
 1. vectorize images' paths'
"""
    
input_images = "C:/Users/vulpine/Desktop/OAM_180515_FULLLIFECYCLE_A/" + "Pos" + str(it) + "/"
mask_paths = sorted(
    [
        os.path.join(input_images, fname)
        for fname in os.listdir(input_images)
        if fname.endswith("_cp_masks.tif")
    ]
 )
# it = 6

"""
 2. load first mask and allocate 
"""   

m0 = PIL.Image.open(mask_paths[0]).convert("L") # zero start #
IS1 = np.array(m0) # image to array # plt.imshow(IS1)
IblankG=np.zeros(IS1.shape, dtype="uint16") # allocate
# masks[:,:,0]=IS1
masks = np.zeros((IS1.shape[0], IS1.shape[1], len(mask_paths))) # allocate

masks[:,:,0]= IS1

"""
 3. loop through time points so that each cell receives the same label 
"""   
    
for it0 in range(1,len(mask_paths)):
######## it0 = 1
        m0b = PIL.Image.open(mask_paths[it0]).convert("L") # plt.imshow(m0b)
        IS2 = np.array(m0b) # plt.imshow(IS2) #
        IS2C = IS2 # plt.imshow(IS2C)
        IS1B = IS1 # plt.imshow(IS1B)
        IS1B= (IS1 > 0.5).astype(np.uint16) # binarize # plt.imshow(IS1B)
        
        IS3 = np.multiply(IS1B,IS2) # plt.imshow(IS3)
        
        tr_cells=(np.unique(IS1.astype(np.int0)))
        tr_cells=tr_cells[1:len(tr_cells)]
        gap_cells=(np.unique(IblankG.astype(np.int0)))
        gap_cells=gap_cells[1:len(gap_cells)]
        
        # gap_cells=np.array((10,122)).astype(np.uint64) 
        # cells_tr =  cells_tr + gap_cells
        if gap_cells.size==0 :
            cells_tr = tr_cells;
        else:
           # cells_tr = np.array([tr_cells, gap_cells]) # !
           # cells_tr =  (tr_cells) + (gap_cells)
            cells_tr = np.concatenate((tr_cells,gap_cells),axis=0) # right concatenation
            cells_tr = np.sort(cells_tr)
            
        Iblank0=np.zeros(IS1.shape, dtype="uint16") # plt.imshow(Iblank0)
     
        for it1 in cells_tr: # ascending, default
            # it1=1
            IS5=(IS1==it1) # plt.imshow(IS5) # sum(sum(IS5))
            IS5=IS5.astype(np.uint16) # plt.imshow(IS5)
            IS56 = thin(IS5,1) # plt.imshow(IS56)  #  sum(sum(IS56))
            IS6A = np.multiply(IS56,IS3) # plt.imshow(IS6A)  #   sum(sum(IS6A))
           
            if sum(sum(IS5))==0 :
               IS5=(IblankG==it1)
               IS6A = np.multiply(IS56,IS2C)
               IblankG[IblankG==it1]=0
              # IblankG=np.where(IblankG==it1,IblankG,0)
               
            if sum(sum(IS6A))>0 :
                IS2ind=(statistics.mode(IS6A[np.nonzero(IS6A)])) # nonzero gives the indexes
                Iblank0[IS2==IS2ind]=it1
                IS3[IS3==IS2ind]=0;
                IS2C[IS2==IS2ind]=0;
                


        bl_cells = np.unique(Iblank0) # plt.imshow(Iblank0)
        bl_cells = bl_cells[1:len(bl_cells)]
        seg_gap = np.setdiff1d(tr_cells,bl_cells) 
        
                   
        if seg_gap.size>0 :
             for itG in seg_gap :
              IblankG[IS1==itG]=itG;
           
         
      
        Iblank0B = np.copy(Iblank0)      # plt.imshow(Iblank0)
        Iblank0B[np.nonzero(Iblank0B)] = 1;  
        Iblank0B = (Iblank0B < 0.5).astype(np.uint16) # plt.imshow(Iblank0B) # invert Iblank0B
        
        ISB = np.multiply(IS2,Iblank0B) # plt.imshow(ISB)
        
        newcells=np.unique(ISB)
        newcells = newcells[1:len(bl_cells)]
        Iblank=Iblank0  # plt.imshow(Iblank)
        A=1;
        
        if newcells.size>0 :
            for it2 in newcells :
             Iblank[IS2==it2]=max(cells_tr)+A; 
             A=A+1;
          # plt.imshow(Iblank)       
                
        masks[:,:,it0]=Iblank # plt.imshow(masks[:,:,it0-1])
        IS1=masks[:,:,it0] # plt.imshow(IS1)
         
       
        
       
        
"""
 4.  Re-assign labels in linear sequence 
"""   
       

ccell2=np.unique(masks[:,:,len(masks[0][0])-1]);
ccell2 = ccell2[1:len(ccell2)] # remove zero since its background


Mask2 = np.zeros((masks.shape[0], masks.shape[1], len(masks[0][0])))

            
for itt4 in range(0,len(masks[0][0])) : # masks
    mask3=np.zeros((masks.shape[0], masks.shape[1]))
            
            
    for itt3 in range(0,len(ccell2)) : # % itt3=23
        pix2 = np.nonzero(masks[:,:,itt4]==ccell2[itt3])
        mask3[pix2]=itt3 + 1 # zero correction    
            
    Mask2[:,:,itt4] = mask3   

 

"""
 5.  Calculate size
"""   

        
no_obj=len(ccell2);
all_ob=np.zeros((no_obj,len(Mask2[0][0])))

for ccell in range(0,no_obj): # ccell=0
           
        for itt in range(0,len(masks[0][0])): # itt = 0
            Mas=Mask2[:,:,itt]==ccell + 1 # % imagesc(Mask2)
            pix=sum(sum(Mas))
            all_ob[ccell,itt]=pix #  figure;imagesc(all_ob)

 
"""
 6.  Find birth times
"""  

cell_exists=np.zeros((1,len(all_ob[:,0]))).astype("uint16")

for itt2 in range(0,len(all_ob[:,0])) : # itt2 = 20

    etime = sum(all_ob[itt2,:]==0) #
    
    if etime == 0 :

       cell_exists[0,itt2]=1

    else:
        
      cell_exists[0,itt2]=np.max(np.nonzero(all_ob[itt2,:]==0)) + 1  


"""
 7.  Find birth times
"""  

     #    Save Mask2,cell_exists,all_obj  
     
     
     
            
            
            
            
        
