#J.-H. Sim

#How to use: qsband.plot.py [bandPlot: lda, qs, both]  [band_num::Int]  [xtic_label1]  [2] ...

import matplotlib.pyplot as plt
import numpy as np
import sys
import argparse

parser = argparse.ArgumentParser()


parser.add_argument("-p", "--plot_type", type=str,  help ="\"lda\", \"qsband\", or \"both\"", default="qsband")
parser.add_argument("-n", "--nkpath", type=int, help="number of kpath")
parser.add_argument("-g", "--kgrid", type=int, help="k-points grid for each k-lines")
parser.add_argument("-kp", "--kpoint_label", type=str, default='faluse')





fig=plt.figure(num=None, figsize=(3,3), dpi=100, facecolor='w', edgecolor='k') #figsize in inches



args = parser.parse_args()


data_raw=0
numOf_K_point=0
numOf_w_point=0
k_grid=[]


grid_per_seg = 0
plot_type = args.plot_type
xtic_label = (args.kpoint_label).split()
num_of_seg = len(xtic_label)-1
print("num_of_seg: ", num_of_seg)



if(plot_type=="qsband" or plot_type=="both"):
  data_raw = np.loadtxt("qsband.dat")
  print("Read qsband.dat")
  numOf_K_point = int( int(data_raw[-1,0])+1)
  numOf_w_point = int(data_raw.shape[0]/ numOf_K_point)


  if(num_of_seg>0):
    grid_per_seg = int(numOf_K_point/num_of_seg);
  
  print("k-grid:", numOf_K_point)
  print("w-grid:", numOf_w_point)
  
  
  
  k_grid = np.zeros(numOf_K_point)
  w_grid = np.zeros(numOf_w_point)
  
  Akw = np.zeros((numOf_w_point,numOf_K_point))
  
  
  
  ######DMFT spectral ftn
  for k in range(numOf_K_point):
      k_grid[k] = data_raw[k*numOf_w_point][1]
      
  for w in range(numOf_w_point):
      w_grid[w] = data_raw[w][2]
      
  for k in range(numOf_K_point):
      for w in range(numOf_w_point):
          Akw[w][k] = data_raw[k*numOf_w_point+w][3]
  
  
  
  max_intensity = np.amax(Akw)
#  k_origin = k_grid[0]
#  k_last = k_grid[numOf_K_point-1]
#  omega_origin = w_grid[0]
#  omega_last = w_grid[numOf_w_point-1]
  
  X,Y =  np.meshgrid(k_grid,w_grid)
  
  
  #plt.pcolormesh(X,Y, Akw, cmap='PuBu_r', vmin=0, vmax=max_intensity, shading='gouraud')
  plt.pcolormesh(X,Y, Akw, cmap='YlOrRd', vmin=0, vmax=max_intensity, shading='gouraud')
  plt.colorbar()



##################################################################
####LDA BAND

elif(plot_type=="lda" or plot_type=="both"):
  lda_band_raw = np.loadtxt("band.dat")
  numOf_K_point    = int( np.amax(lda_band_raw[:,0]+1))
  numOf_bands = int(lda_band_raw.shape[0] / numOf_K_point)

  if(num_of_seg>0):
    grid_per_seg = int(numOf_K_point/num_of_seg);
  
  print("k-grid_lda  :", numOf_K_point)
  print("num_of_Bands:", numOf_bands)
  
  
  ax=fig.add_subplot(111)
  lda_band_disper = np.zeros((numOf_K_point, numOf_bands))
  k_grid = np.zeros(numOf_K_point)
  for k in range(numOf_K_point):
      k_grid[k] = lda_band_raw[k][1]
      for w in range(numOf_bands):
          lda_band_disper[k][w] = lda_band_raw[w*numOf_K_point+k][2]
  
  ax.plot(k_grid,lda_band_disper,c='k',ls='-',linewidth="0.5")
##################################################################
#



##PLOT option


plt.ylabel("Energy (eV)")

if(num_of_seg>0):
  xtic_points =[k_grid[k] for k in np.arange(0,numOf_K_point,grid_per_seg)]
  plt.xticks(xtic_points, xtic_label)
  print(xtic_points, xtic_label)






plt.ylim(-20,20)
#plt.ylim(-10,5)
plt.axhline(0, c='w',ls='--',lw='0.2') #fermi level
plt.grid(b='TRUE',which='major',axis='x',c='black',ls='-',lw='0.1')


plt.tight_layout()
fig.set_size_inches(3.0,2.1)
plt.show()
fig.savefig('qsband.png', dpi=600, bbox_inches='tight')
