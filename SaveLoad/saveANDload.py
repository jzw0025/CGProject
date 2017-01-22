####### save #########

import scipy.io as sio

directoryFileName = '/Users/junchaowei/desktop/realdata/iter2'

# save the image system
sio.savemat(directoryFileName+"/image1", {'par1':image1}) # image1 
sio.savemat(directoryFileName+"/image2", {'par1':image2}) # image2 
sio.savemat(directoryFileName+"/image2_3", {'par1':image2_3}) # image2_3

# save the particle system
sio.savemat(directoryFileName+"/point_arr1", {'par1':point_arr1}) # point_arr1 
sio.savemat(directoryFileName+"/point_arr2", {'par1':point_arr2}) # point_arr2 
sio.savemat(directoryFileName+"/referennce_point_arr1", {'par1':referennce_point_arr1}) # referennce_point_arr1 
sio.savemat(directoryFileName+"/referennce_point_arr2", {'par1':referennce_point_arr2}) # referennce_point_arr2 

# save the equal-histogram plot (optional)
sio.savemat(directoryFileName+"/im2", {'par1':im2}) # point_arr1 
sio.savemat(directoryFileName+"/im3", {'par1':im3}) # point_arr1 

# save the mesh system (optional)

# save the computation system
sio.savemat(directoryFileName+"/dispx", {'par1':dispx}) # point_arr1 
sio.savemat(directoryFileName+"/dispxx", {'par1':dispxx}) # point_arr1 
sio.savemat(directoryFileName+"/dispy", {'par1':dispy}) # point_arr1 
sio.savemat(directoryFileName+"/dispyy", {'par1':dispyy}) # point_arr1 
sio.savemat(directoryFileName+"/dispz", {'par1':dispz}) # point_arr1 
sio.savemat(directoryFileName+"/dispzz", {'par1':dispzz}) # point_arr1 

# save the strain computation
sio.savemat(directoryFileName+"/s11", {'par1':s11}) # point_arr1 
sio.savemat(directoryFileName+"/s12", {'par1':s12}) # point_arr1 
sio.savemat(directoryFileName+"/s13", {'par1':s13}) # point_arr1 
sio.savemat(directoryFileName+"/s22", {'par1':s22}) # point_arr1 
sio.savemat(directoryFileName+"/s23", {'par1':s23}) # point_arr1 
sio.savemat(directoryFileName+"/s33", {'par1':s33}) # point_arr1 
sio.savemat(directoryFileName+"/s_effective", {'par1':s_effective}) # point_arr1 

#################################################################################
################################    Load   ######################################
#################################################################################

directoryFileName = '/Users/junchaowei/desktop/realdata/iter2'

mat_contents1 = sio.loadmat(directoryFileName+"/image1")
image1 = mat_contents1['par1'] # default saved data name.
mat_contents1 = sio.loadmat(directoryFileName+"/image2")
image2 = mat_contents1['par1'] # default saved data name.
mat_contents1 = sio.loadmat(directoryFileName+"/image2_3")
image2_3 = mat_contents1['par1'] # default saved data name.

mat_contents1 = sio.loadmat(directoryFileName+"/im2")
im2 = mat_contents1['par1'] # default saved data name.
mat_contents1 = sio.loadmat(directoryFileName+"/im3")
im3 = mat_contents1['par1'] # default saved data name.

mat_contents1 = sio.loadmat(directoryFileName+"/point_arr1")
point_arr1 = mat_contents1['par1'] # default saved data name.
mat_contents1 = sio.loadmat(directoryFileName+"/point_arr2")
point_arr2 = mat_contents1['par1'] # default saved data name.
mat_contents1 = sio.loadmat(directoryFileName+"/referennce_point_arr1")
referennce_point_arr1 = mat_contents1['par1'] # default saved data name.
mat_contents1 = sio.loadmat(directoryFileName+"/referennce_point_arr2")
referennce_point_arr2 = mat_contents1['par1'] # default saved data name.

mat_contents1 = sio.loadmat(directoryFileName+"/dispx")
dispx = mat_contents1['par1'] # default saved data name.
mat_contents1 = sio.loadmat(directoryFileName+"/dispxx")
dispxx = mat_contents1['par1'] # default saved data name.
mat_contents1 = sio.loadmat(directoryFileName+"/dispy")
dispy = mat_contents1['par1'] # default saved data name.
mat_contents1 = sio.loadmat(directoryFileName+"/dispyy")
dispyy = mat_contents1['par1'] # default saved data name.
mat_contents1 = sio.loadmat(directoryFileName+"/dispz")
dispz = mat_contents1['par1'] # default saved data name.
mat_contents1 = sio.loadmat(directoryFileName+"/dispzz")
dispzz = mat_contents1['par1'] # default saved data name.

from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as ndimage


################### save the displacement ##############
inter_x = interpolate.LinearNDInterpolator(point_arr1, dispx, fill_value=0.0, rescale=True) 
inter_y = interpolate.LinearNDInterpolator(point_arr1, dispy, fill_value=0.0, rescale=True) 
inter_z = interpolate.LinearNDInterpolator(point_arr1, dispz, fill_value=0.0, rescale=True) 

#startz1 = 67
#level_z = startz1
sx, sy = image1[:,:,0].shape
mx, my = np.meshgrid(np.linspace(1,sx,sx),np.linspace(1,sy,sy))
# meshgrid change the orientation of x,y: x becomes the horizontal axis, y becomes the vertical
# this change would affect the return statement of reshape.
Cord_xy = np.vstack([mx.flatten(), my.flatten()]) # column major vector
directoryFileName21 = '/Users/junchaowei/desktop/realdata/iter2/xdisp'
directoryFileName22 = '/Users/junchaowei/desktop/realdata/iter2/ydisp'
directoryFileName23 = '/Users/junchaowei/desktop/realdata/iter2/zdisp'

for i in range(im2.shape[2]):
    selected = i
    xy = np.concatenate((Cord_xy, selected*np.ones([1,Cord_xy.shape[1]])), axis=0)
    extrapo_x_layer = inter_x.__call__(xy.T)
    extrapo_y_layer = inter_y.__call__(xy.T)
    extrapo_z_layer = inter_z.__call__(xy.T)
    reshape_extrapo_x_layer = ndimage.gaussian_filter(extrapo_x_layer.reshape(sy,sx),3) # gaussian smoother
    reshape_extrapo_y_layer = ndimage.gaussian_filter(extrapo_y_layer.reshape(sy,sx),3) # gaussian smoother
    reshape_extrapo_z_layer = ndimage.gaussian_filter(extrapo_z_layer.reshape(sy,sx),3) # gaussian smoother  
    ######################  x  #####################
    plt.figure(i)
    levels =[dispx.min(),dispx.min()/2, dispx.min()/4, 0, dispx.max()/4,dispx.max()/2, dispx.max()]
    plt.imshow(image1[:,:,selected].T,interpolation='nearest', cmap='Greys_r')
    plt.contourf(mx, my, reshape_extrapo_x_layer, levels,
                    alpha=0.2)
    plt.colorbar()
    plt.title("displacement x")
    plt.xlabel("X Direction")
    plt.ylabel("Y Direction")
    plt.savefig(directoryFileName21+'/M6_OD_125_C'+str(i)+'.png', bbox_inches='tight')
    ######################  y  #####################
    plt.figure(i+im2.shape[2])
    levels =[dispy.min(),dispy.min()/2, dispy.min()/4, 0, dispy.max()/4,dispy.max()/2, dispy.max()]
    plt.imshow(image1[:,:,selected].T,interpolation='nearest', cmap='Greys_r')
    plt.contourf(mx, my, reshape_extrapo_y_layer, levels,
                    alpha=0.2)
    plt.colorbar()
    plt.title("displacement y")
    plt.xlabel("X Direction")
    plt.ylabel("Y Direction")
    plt.savefig(directoryFileName22+'/M6_OD_125_C'+str(i)+'.png', bbox_inches='tight')
    ######################  z  #####################
    plt.figure(i+2*im2.shape[2])
    levels =[dispz.min(),dispz.min()/2, dispz.min()/4, 0, dispz.max()/4,dispz.max()/2, dispz.max()]
    plt.imshow(image1[:,:,selected].T,interpolation='nearest', cmap='Greys_r')
    plt.contourf(mx, my, reshape_extrapo_z_layer, levels,
                    alpha=0.2)
    plt.colorbar()
    plt.title("displacement z")
    plt.xlabel("X Direction")
    plt.ylabel("Y Direction")
    plt.savefig(directoryFileName23+'/M6_OD_125_C'+str(i)+'.png', bbox_inches='tight')
    
plt.close("all")
 
###################  save the image ###############       
directoryFileName3 = '/Users/junchaowei/desktop/realdata/iter2/img'
for i in range(im2.shape[2]):
    selected = i
    plt.imshow(image1[:,:,selected].T,interpolation='nearest', cmap='Greys_r')
    plt.xlabel("X Direction")
    plt.ylabel("Y Direction")
    plt.savefig(directoryFileName3+'/M6_OD_125_C'+str(i)+'.png', bbox_inches='tight')
    
plt.close("all")

###################  save the strain ################
inter_s11 = interpolate.LinearNDInterpolator(point_arr1, s11, fill_value=0.0, rescale=True) 
inter_s12 = interpolate.LinearNDInterpolator(point_arr1, s12, fill_value=0.0, rescale=True) 
inter_s13 = interpolate.LinearNDInterpolator(point_arr1, s13, fill_value=0.0, rescale=True) 
inter_s22 = interpolate.LinearNDInterpolator(point_arr1, s22, fill_value=0.0, rescale=True) 
inter_s23 = interpolate.LinearNDInterpolator(point_arr1, s23, fill_value=0.0, rescale=True) 
inter_s33 = interpolate.LinearNDInterpolator(point_arr1, s33, fill_value=0.0, rescale=True)
inter_effective = interpolate.LinearNDInterpolator(point_arr1, s_effective, fill_value=0.0, rescale=True) 

sx, sy = image1[:,:,0].shape
mx, my = np.meshgrid(np.linspace(1,sx,sx),np.linspace(1,sy,sy))
# meshgrid change the orientation of x,y: x becomes the horizontal axis, y becomes the vertical
# this change would affect the return statement of reshape.
Cord_xy = np.vstack([mx.flatten(), my.flatten()]) # column major vector
directoryFileName31 = '/Users/junchaowei/desktop/realdata/iter2/s11'
directoryFileName32 = '/Users/junchaowei/desktop/realdata/iter2/s12'
directoryFileName33 = '/Users/junchaowei/desktop/realdata/iter2/s13'
directoryFileName34 = '/Users/junchaowei/desktop/realdata/iter2/s22'
directoryFileName35 = '/Users/junchaowei/desktop/realdata/iter2/s23'
directoryFileName36 = '/Users/junchaowei/desktop/realdata/iter2/s33'
directoryFileName37 = '/Users/junchaowei/desktop/realdata/iter2/effectivestrain'

for i in range(im2.shape[2]):
    selected = i
    xy = np.concatenate((Cord_xy, selected*np.ones([1,Cord_xy.shape[1]])), axis=0)
    extrapo_s11_layer = inter_s11.__call__(xy.T)
    extrapo_s12_layer = inter_s12.__call__(xy.T)
    extrapo_s13_layer = inter_s13.__call__(xy.T)
    extrapo_s22_layer = inter_s22.__call__(xy.T)
    extrapo_s23_layer = inter_s23.__call__(xy.T)
    extrapo_s33_layer = inter_s33.__call__(xy.T)
    extrapo_effective_layer = inter_s33.__call__(xy.T)
    reshape_extrapo_s11_layer = ndimage.gaussian_filter(extrapo_s11_layer.reshape(sy,sx),3) # gaussian smoother
    reshape_extrapo_s12_layer = ndimage.gaussian_filter(extrapo_s12_layer.reshape(sy,sx),3) # gaussian smoother
    reshape_extrapo_s13_layer = ndimage.gaussian_filter(extrapo_s13_layer.reshape(sy,sx),3) # gaussian smoother
    reshape_extrapo_s22_layer = ndimage.gaussian_filter(extrapo_s22_layer.reshape(sy,sx),3) # gaussian smoother
    reshape_extrapo_s23_layer = ndimage.gaussian_filter(extrapo_s23_layer.reshape(sy,sx),3) # gaussian smoother
    reshape_extrapo_s33_layer = ndimage.gaussian_filter(extrapo_s33_layer.reshape(sy,sx),3) # gaussian smoother
    reshape_extrapo_effective_layer = ndimage.gaussian_filter(extrapo_effective_layer.reshape(sy,sx),3) # gaussian smoother
    
    ######################  s11  #####################
    plt.figure(i)
    levels =[s11.min(),s11.min()/2, s11.min()/4, 0, s11.max()/4,s11.max()/2, s11.max()]
    plt.imshow(image1[:,:,selected].T,interpolation='nearest', cmap='Greys_r')
    plt.contourf(mx, my, reshape_extrapo_s11_layer, levels,
                    alpha=0.2)
    plt.colorbar()
    plt.title("Strain 11")
    plt.xlabel("X Direction")
    plt.ylabel("Y Direction")
    plt.savefig(directoryFileName31+'/M6_OD_125_C'+str(i)+'.png', bbox_inches='tight')
    plt.close("all")
    ######################  s12  #####################
    plt.figure(i+im2.shape[2])
    levels =[s12.min(),s12.min()/2, s12.min()/4, 0, s12.max()/4,s12.max()/2, s12.max()]
    plt.imshow(image1[:,:,selected].T,interpolation='nearest', cmap='Greys_r')
    plt.contourf(mx, my, reshape_extrapo_s12_layer, levels,
                    alpha=0.2)
    plt.colorbar()
    plt.title("Strain 12")
    plt.xlabel("X Direction")
    plt.ylabel("Y Direction")
    plt.savefig(directoryFileName32+'/M6_OD_125_C'+str(i)+'.png', bbox_inches='tight')
    plt.close("all")
    ######################  s13  #####################
    plt.figure(i+2*im2.shape[2])
    levels =[s13.min(),s13.min()/2, s13.min()/4, 0, s13.max()/4,s13.max()/2, s13.max()]
    plt.imshow(image1[:,:,selected].T,interpolation='nearest', cmap='Greys_r')
    plt.contourf(mx, my, reshape_extrapo_s13_layer, levels,
                    alpha=0.2)
    plt.colorbar()
    plt.title("Strain 13")
    plt.xlabel("X Direction")
    plt.ylabel("Y Direction")
    plt.savefig(directoryFileName33+'/M6_OD_125_C'+str(i)+'.png', bbox_inches='tight')
    plt.close("all")
    ######################  s22  #####################
    plt.figure(i+3*im2.shape[2])
    levels =[s22.min(),s22.min()/2, s22.min()/4, 0, s22.max()/4,s22.max()/2, s22.max()]
    plt.imshow(image1[:,:,selected].T,interpolation='nearest', cmap='Greys_r')
    plt.contourf(mx, my, reshape_extrapo_s22_layer, levels,
                    alpha=0.2)
    plt.colorbar()
    plt.title("Strain 22")
    plt.xlabel("X Direction")
    plt.ylabel("Y Direction")
    plt.savefig(directoryFileName34+'/M6_OD_125_C'+str(i)+'.png', bbox_inches='tight')
    plt.close("all")
    ######################  s23  #####################
    plt.figure(i+4*im2.shape[2])
    levels =[s23.min(),s23.min()/2, s23.min()/4, 0, s23.max()/4,s23.max()/2, s23.max()]
    plt.imshow(image1[:,:,selected].T,interpolation='nearest', cmap='Greys_r')
    plt.contourf(mx, my, reshape_extrapo_s23_layer, levels,
                    alpha=0.2)
    plt.colorbar()
    plt.title("Strain 23")
    plt.xlabel("X Direction")
    plt.ylabel("Y Direction")
    plt.savefig(directoryFileName35+'/M6_OD_125_C'+str(i)+'.png', bbox_inches='tight')
    plt.close("all")
    ######################  s33  #####################
    plt.figure(i+5*im2.shape[2])
    levels =[s33.min(),s33.min()/2, s33.min()/4, 0, s33.max()/4,s33.max()/2, s33.max()]
    plt.imshow(image1[:,:,selected].T,interpolation='nearest', cmap='Greys_r')
    plt.contourf(mx, my, reshape_extrapo_s33_layer, levels,
                    alpha=0.2)
    plt.colorbar()
    plt.title("Strain 33")
    plt.xlabel("X Direction")
    plt.ylabel("Y Direction")
    plt.savefig(directoryFileName36+'/M6_OD_125_C'+str(i)+'.png', bbox_inches='tight')
    plt.close("all")
    ######################  effective strain  #####################
    plt.figure(i+6*im2.shape[2])
    levels = [1e6*s_effective.min(),  1e6*s_effective.mean()/4, 1e6*s_effective.mean()/2, 1e6*s_effective.mean(), 1e6*s_effective.mean()*1.5,
                1e6*s_effective.mean()*2,1e6*s_effective.mean()*3,1e6*s_effective.mean()*4]
    plt.imshow(image1[:,:,selected].T,interpolation='nearest', cmap='Greys_r')
    plt.contourf(mx, my, 1e6*reshape_extrapo_effective_layer, levels,
                    alpha=0.2)
    plt.colorbar()
    plt.title("Effective Micro Strain ")
    plt.xlabel("X Direction")
    plt.ylabel("Y Direction")
    plt.savefig(directoryFileName37+'/M6_OD_125_C'+str(i)+'.png', bbox_inches='tight')
    plt.close("all")
    
