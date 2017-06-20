#simulating visibilities
from astropy.io import fits
from astropy import units as u
import numpy as np
import math
import anta_pos as ant
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import healpy as hp
from astropy.io import fits
#import generate_gaussian_data as gdata
# read map
#map_I = hp.read_map('GLEAM_EGC.fits')
#hp.mollview(map_I, coord=['G','E'], title='Histogram equalized Ecliptic', unit='mK', norm='hist', min=-1,max=1, xsize=2000)

RA_src, Dec_src = np.load('src_pos.npy')[0], np.load('src_pos.npy')[1]
src_flux = np.load('src_flux.npy')[1]

xdim = int(sys.argv[1])
ydim = int(sys.argv[2])
long_0= 21.4278 # deg
lalt_0=  -30.7224 # deg
ant_positions= ant.get_antenna_pos(xdim,ydim,long_0,lalt_0)
#x,y = ant.get_antenna_pos(xdim,ydim,long_0,lalt_0)
zen_vec = [0.0,1.0,0.0] #np.round(ant_positions[0]/np.linalg.norm(ant_positions[0]))





x, y, z = np.cos(np.radians(Dec_src))*np.cos(np.radians(RA_src)), np.cos(np.radians(Dec_src))*np.sin(np.radians(RA_src)), np.sin(np.radians(Dec_src))
xyz = np.concatenate((x,y,z))
mydot = x*zen_vec[0]+ y*zen_vec[1] + z*zen_vec[2]

theta = np.arccos(mydot)
sigma = 0.5*(1.22*(0.5/6.0))
beam = np.exp(-0.5*(theta/sigma)**2)
ii =  beam > 0.1

RA_src= RA_src[ii]
Dec_src= Dec_src[ii]
src_flux = src_flux[ii]

# ants epsilon

Ant_eps = np.load('Ant_epsolon0.05.npy')
blVectors, blPairs = ant.get_blvectors(ant_positions)
	#np.save('redund_blVectors.npy',[blVectors,blPairs])
	

ublDict = {}
for b in range(len(blVectors)):
     	# grouping togather all blPairs with same blVector
    	if ublDict.has_key(blVectors[b]): ublDict[blVectors[b]].append(blPairs[b])
    	else: ublDict[blVectors[b]] = [blPairs[b]]

ublIndexDict = {antPair: i for i,antPairs in enumerate(ublDict.values()) for antPair in antPairs }
ublVectors = np.array([ant_positions[antList[0][0]]-ant_positions[antList[0][1]] for antList in ublDict.values()])



print "With", len(ant_positions), "antennas there are", len(ublDict.items()), "unique baselines."
#return [ublVectors,ublIndexDict]

	
zen = zen_vec 
ant1, ant2 =[],[]
for k in range(len(blPairs)):
		ant1.append(np.array(blPairs)[k][0])
		ant2.append(np.array(blPairs)[k][1])

vis_map =[]
for i in range(len(ant1)):
		vis_map.append(ublIndexDict[(ant1[i],ant2[i])])

ant_eps_i = Ant_eps[0][0]
corr_matrix= np.zeros((len(blVectors),len(blVectors)),dtype='complex')
l_wave = 0.5
for bl_i in range(len(blVectors)):
		for bl_j in range(len(blVectors)):

	
		       	 	x, y, z = np.cos(np.radians(Dec_src))*np.cos(np.radians(RA_src)), np.cos(np.radians(Dec_src))*np.sin(np.radians(RA_src)), np.sin(np.radians(Dec_src))
				r=np.sqrt(x**2 + y**2 + z**2)
				mydot = (x/r)*zen[0] + (y/r)*zen[1] + (z/r)*zen[2]
  	       			theta = np.arccos(mydot)
				#print np.linalg.norm(blVectors[bl_i]) -theta ,ant_eps_i[ant1[bl_i]]
				beam1 = np.exp(-(0.5*(theta)**2)*(1.0/(sigma*(1.0 + ant_eps_i[ant1[bl_i]]))**2))
				beam2 = np.exp(-0.5*((theta)**2)*(1.0/(sigma*(1.0 +ant_eps_i[ant2[bl_j]]))**2))
				#plt.plot(beam1,'*')
				#plt.plot(beam2,'.')
				#print blVectors[bl_i][0]
				mydotxyz =  (x/r)*blVectors[bl_i][0] + (y/r)*blVectors[bl_i][1] + (z/r)*blVectors[bl_i][2]
				exp_phase_shift = np.exp(-2j*np.pi*(mydotxyz/l_wave))
				
				I_u_i = np.sum(src_flux*exp_phase_shift)
				mydotxyz_j =  (x/r)*blVectors[bl_j][0] + (y/r)*blVectors[bl_j][1] + (z/r)*blVectors[bl_j][2]
				exp_phase_shift_j = np.exp(-2j*np.pi*(mydotxyz_j/l_wave))
				I_u_j = np.sum(src_flux*exp_phase_shift_j)
				c_u_ij = np.conj(I_u_i)*I_u_j
				corr_matrix[bl_i][bl_j]= np.sum(beam1*beam2*c_u_ij)
				#plt.plot(blVectors[bl_i][0],corr,'*')



bl_len = np.sqrt(np.array(blVectors)[:,0]**2 + np.array(blVectors)[:,1]**2 + np.array(blVectors)[:,2]**2)
for k in range(len(corr_matrix[0])):
	plt.plot(vis_map,corr_matrix[k],'*')
				
			

			 
			
			 




