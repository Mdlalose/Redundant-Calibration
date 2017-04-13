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
zen_vec = ant_positions[0]/np.linalg.norm(ant_positions[0])
"""
plt.figure()
plt.plot(x,y,'o')
plt.xlabel('E-W [m]')
plt.ylabel('N-S [m]')
plt.grid()

Nside = 256




t,p =hp.pix2ang(Nside,np.arange(hp.nside2npix(Nside)))

Dec=np.pi/2-t
H=p
Dec[0:len(Dec_src)] = Dec_src
H[0:len(RA_src)] = RA_src
x, y, z = np.cos(np.radians(Dec))*np.cos(np.radians(H)), np.cos(np.radians(Dec))*np.sin(np.radians(H)), np.sin(np.radians(Dec))

mydot = x*zen_vec[0] + y*zen_vec[1] + z*zen_vec[2]
theta = np.arccos(mydot)
def sigma(l_wave,dia_m):
    return (0.5*(1.22*l_wave))/dia_m
sigma_ = sigma(0.5,6.0)
beam = np.exp(-0.5*(theta/sigma_)**2)
Beam = np.asarray(beam) 
hp.mollview(Beam)
plt.show()

"""
"""
nu = np.arange(170.0,231,2.0)
nu_0 = (231.0-170.0)/2.0
flux_nu =[]
for nu_k in range(len(nu)):
                x, y, z = np.cos(np.radians(Dec_src))*np.cos(np.radians(RA_src)), np.cos(np.radians(Dec_src))*np.sin(np.radians(RA_src)), np.sin(np.radians(Dec_src))	
		mydot = x*zen_vec[0] + y*zen_vec[1] + z*zen_vec[2]
		theta = np.arccos(mydot)
		theta= theta[theta <= 2.0*np.pi]
		theta2 = theta*theta
		sigma_ = (0.5*(1.22*0.5))/6.0
		beam = np.exp(-theta2/(0.5*0.05))		
		tmp =beam*src_flux*(nu[nu_k]/nu_0)**-2.5
		flux_nu.append(tmp)
		

x, y, z = np.cos(np.radians(Dec_src))*np.cos(np.radians(RA_src)), np.cos(np.radians(Dec_src))*np.sin(np.radians(RA_src)), np.sin(np.radians(Dec_src))	
mydot = x*zen_vec[0] + y*zen_vec[1] + z*zen_vec[2]
theta = np.arccos(mydot)
theta= theta[theta <= 2.0*np.pi]
sigma_ = 0.1
beam = np.exp(-0.5*(theta/sigma_)**2)		
tmp =beam #*src_flux*(nu[nu_k]/nu_0)**-2.5
flux = beam*src_flux
flux = flux.sort()
flux_10b = flux[len(flux)-10:-1]
"""

def get_V_uvw(blVectors,RA_src,Dec_src,src_flux,l_wave,dia_m,zen_vec,nu_min,nu_max):
	nu = np.arange(nu_min,nu_max,2.0)
	nu_0 = (nu_max-nu_min)/2.0
	V_nu = []
	for nu_i in range(len(nu)):
		#print (nu[nu_i]/nu_0)**-2.5
		V_uvw = np.zeros(len(blVectors), dtype = "complex")
       		x, y, z = np.cos(np.radians(Dec_src))*np.cos(np.radians(RA_src)), np.cos(np.radians(Dec_src))*np.sin(np.radians(RA_src)), np.sin(np.radians(Dec_src))	
		for bl in range(len(blVectors)):
			mydotxyz = x*blVectors[bl][0] + y*blVectors[bl][1] + z*blVectors[bl][2]
                        phase_shift = 2.0*np.pi*(mydotxyz/l_wave)
			exp_phase_shift = np.cos(phase_shift) -1j*np.sin(phase_shift)
               		#print exp_phase_shift
			sigma_ = (0.5*(1.22*l_wave))/dia_m

			mydot = x*zen_vec[0] + y*zen_vec[1] + z*zen_vec[2]
			theta = np.arccos(mydot)
			theta = np.array(theta)
		        theta= theta[theta <= 2.0*np.pi]
			theta2 = theta*theta
			beam = np.exp(-theta2/(2.0*sigma_**2))
			#dmydot = x*dblVectors[bl][0] + y*dblVectors[bl][1] + z*dblVectors[bl][2]
			#tmp0 = -1j*2.0*np.sum(beam*src_flux*(dmydot/l_wave))
			tmp= np.sum(beam*src_flux*np.power(nu[nu_i]/nu_0,-2.5)*exp_phase_shift) 
			V_uvw[bl]=tmp
			
		V_nu.append(V_uvw)
		
	return np.array(V_nu)

ant_positions_err=[]
dx = 0.05
dy = dx
dz = dx
for k in range(len(ant_positions)):
	ant_positions_err.append(ant_positions[k]+ [np.random.normal(dx),np.random.normal(dy),np.random.normal(dz)])

blVectors, blPairs = ant.get_blvectors(ant_positions_err)
#ublVectors,ublIndexDick = ant.get_ublvectors(ant_positions,blVectors,blPairs)
#def get_ublvectors(Antenna_positions,blVectors,blPairs):
ublDict = {}
for b in range(len(blVectors)):
     		# grouping togather all blPairs with same blVector
    		if ublDict.has_key(blVectors[b]): ublDict[blVectors[b]].append(blPairs[b])
    		else: ublDict[blVectors[b]] = [blPairs[b]]

ublIndexDict = {antPair: i for i,antPairs in enumerate(ublDict.values()) for antPair in antPairs }
ublVectors = np.array([ant_positions_err[antList[0][0]]-ant_positions_err[antList[0][1]] for antList in ublDict.values()])

"""
bl_i_j =[]
db_i_j = []
b_k = []
dp = 1.0
for antList in ublDict.values():
	for i in range(len(antList)):
		dtmp0 =[np.random.normal(0.0,dp),np.random.normal(0.0,dp),np.random.normal(0.0,dp)]
		bl_i_j.append(ant_positions[antList[0][0]]-ant_positions[antList[0][1]] + dtmp0)
		db_i_j.append(dtmp0)
		b_k.append(ant_positions[antList[0][0]]-ant_positions[antList[0][1]])
"""

print "With", len(ant_positions_err), "antennas there are", len(ublDict.items()), "unique baselines."
#return [ublVectors,ublIndexDict]


 
ant1, ant2 =[],[]
for k in range(len(blPairs)):
	ant2.append(np.array(blPairs)[k][0])
	ant1.append(np.array(blPairs)[k][1])

vis_map =[]
for (key,value) in ublIndexDict.iteritems():
	vis_map.append(value)   
	 

	
visTrue = get_V_uvw(ublVectors,RA_src,Dec_src,src_flux,0.5,6.0,zen_vec,170.0,231.0)
np.save('obstrue_vis_8x8_0.05.npy',[visTrue,ant1,ant2,vis_map])
print " successfully calculated visibilities"

"""
u= []
v= []
uk= []
vk= []

for i in range(len(bl_i_j)):
	u.append(bl_i_j[i][0])
	v.append(bl_i_j[i][1])

for i in range(len(ublVectors)):
	uk.append(ublVectors[i][0])
	vk.append(ublVectors[i][1])	
plt.title('Near-Redundant N(0.0,0.05m)')
plt.plot(np.array(u)*(1.0/0.5),np.array(v)*(1.0/0.5),'.')
plt.plot(np.array(uk)*(1.0/0.5),np.array(vk)*(1.0/0.5),'*r')
plt.xlabel('u [m/$\lambda$] ')
plt.ylabel('v [m/$\lambda$]')
#plt.legend(loc='best')
plt.grid()
"""

