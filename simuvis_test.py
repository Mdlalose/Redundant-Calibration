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
"""
plt.figure()
plt.plot(x,y,'o')
plt.xlabel('E-W [m]')
plt.ylabel('N-S [m]')
plt.grid()
"""





x, y, z = np.cos(np.radians(Dec_src))*np.cos(np.radians(RA_src)), np.cos(np.radians(Dec_src))*np.sin(np.radians(RA_src)), np.sin(np.radians(Dec_src))
xyz = np.concatenate((x,y,z))
mydot = x*zen_vec[0]+ y*zen_vec[1] + z*zen_vec[2]

theta = np.arccos(mydot)
sigma = 0.5*(1.22*(0.5/6.0))
beam = np.exp(-0.5*(theta/sigma)**2)
ii =  beam > 0.1
#beam = beam[beam>0.1]

RA_src= RA_src[ii]
Dec_src= Dec_src[ii]
src_flux = src_flux[ii]


x, y, z = np.cos(np.radians(Dec_src))*np.cos(np.radians(RA_src)), np.cos(np.radians(Dec_src))*np.sin(np.radians(RA_src)), np.sin(np.radians(Dec_src))
mydot = x*zen_vec[0] + y*zen_vec[1] + z*zen_vec[2]
theta = np.arccos(mydot)
theta= theta[theta <= 2.0*np.pi]
sigma_ = 0.05
beam = np.exp(-0.5*(theta/0.01*sigma_)**2)
tmp =beam #*src_flux*(nu[nu_k]/nu_0)**-2.5
flux = beam*src_flux
flux_ = beam*src_flux
flux.sort()
flux_10b = flux[len(x)-10:]

RA_b, Dec_b, src_b =[], [], []
for src in range(len(flux_10b)):
	for i in range(len(src_flux)):
		if np.round(flux_[i],5) == np.round(flux_10b[src],5):
			# np.round(flux_10b[src],5), np.round(flux_[i],5)
			RA_b.append(RA_src[i])
			Dec_b.append(Dec_src[i])
			src_b.append(src_flux[i])

		else:
			pass


np.save('10_brightest_src.npy',[RA_b,Dec_b, src_b])



def get_V_uvw(blVectors,RA_src,Dec_src,src_flux,l_wave,dia_m,zen,nu_min,nu_max,ant1,ant2,ant_eps):
	nu = np.arange(nu_min,nu_max,1.0)
	nu_0 = (nu_max-nu_min)/2.0
	V_nu = []
	myvis=[]

	expec_err =[]

	sigma = 0.5*(1.22*(l_wave/dia_m))
	for nu_i in range(len(nu)):
		#print (nu[nu_i]/nu_0)**-2.5
		V_uvw =   np.zeros(len(blVectors), dtype = "complex")
		#V_uvw_error =  np.zeros(len(blVectors), dtype = "complex")
		
		nsrc = len(RA_src)
		#x,y, z = np.random.randn(nsrc), np.random.randn(nsrc), np.random.randn(nsrc)
		x, y, z = np.cos(np.radians(Dec_src))*np.cos(np.radians(RA_src)), np.cos(np.radians(Dec_src))*np.sin(np.radians(RA_src)), np.sin(np.radians(Dec_src))
		r=np.sqrt(x**2 + y**2 + z**2)
		mydot = (x/r)*zen[0] + (y/r)*zen[1] + (z/r)*zen[2]
  	        theta = np.arccos(mydot)
  	        ant_eps_i = ant_eps[nu_i]
		#print ant_eps_i.size, nu[nu_i]

	
		for k in range(len(blVectors)):
			 
			
		         
			 beam12 = np.exp(-(0.5*theta**2)*(1.0/(sigma*(1.0 + ant_eps_i[ant1[k]]))**2 + 1.0/(sigma*(1.0 +ant_eps_i[ant2[k]]))**2))
			 beam = np.exp(-(theta/sigma)**2)
			 #plt.title('0.001$\epsilon$')
			 #plt.plot(theta,beam,'.',label='identical Beam')
			 #plt.plot(theta,beam12,'.',label='Imperfection Beams')
			 #plt.legend(loc='best')
	
			 
                         
                         
  			 r=np.sqrt(x**2 + y**2 + z**2)
  			 
			 mydotxyz =  (x/r)*blVectors[k][0] + (y/r)*blVectors[k][1] + (z/r)*blVectors[k][2]
			 exp_phase_shift = np.exp(-2j*np.pi*(mydotxyz/l_wave))
			 tmp= np.sum(beam*src_flux*np.power(nu[nu_i]/nu_0,-2.5)*exp_phase_shift)
			 V_uvw[k]=tmp
			 # first order error in visibility due to imperfections in antenna locations
			 #delta_blVectors= np.array(blVectors_err[bl]) - np.array(blVectors[bl])
			 #mydotxyz_err =  (x/r)*delta_blVectors[0] + (y/r)*delta_blVectors[1] + (z/r)*delta_blVectors[2]
			 #tmp_err = -1j*np.sum(beam*src_flux*((2.0*np.pi*mydotxyz_err)/l_wave))

			 #V_uvw_error[bl] = tmp_err
		
			 
                         
		V_nu.append(V_uvw)
	

	return np.array(V_nu)
		        	


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

	
       
ant1, ant2 =[],[]
for k in range(len(blPairs)):
		ant1.append(np.array(blPairs)[k][0])
		ant2.append(np.array(blPairs)[k][1])

vis_map =[]
for i in range(len(ant1)):
		vis_map.append(ublIndexDict[(ant1[i],ant2[i])])


#RA_src, Dec_src, src_flux = np.load('10_brightest_src.npy')[0], np.load('10_brightest_src.npy')[1], np.load('10_brightest_src.npy')[2]
sim_vis_n_sample =[]
for n_i in range(1):

	visTrue = get_V_uvw(blVectors,RA_src,Dec_src,src_flux,0.5,6.0,zen_vec,170,230,ant1,ant2,Ant_eps[n_i])
	sim_vis_n_sample.append(visTrue)
	




#np.save('obstrue_vis_' + repr(xdim) + 'x' + repr(ydim)  + '20_data_redundant'+ '.npy',[sim_vis_n_sample,ant1,ant2,vis_map])

"""
error = float(sys.argv[3])
n_sample = int(sys.argv[4])
ant_1=[]
sg=[]
bl_=[]
v1 = np.load('obstrue_vis_8x80.0.npy')[0]
#bl_red = np.load('redund_blVectors.npy')[0]
for n in range(n_sample):

	ant_positions_err= np.random.randn(ant_positions.shape[0],ant_positions.shape[1])
	#print ant_positions_err[0][2]

	ant_positions_err[:,1] = 0.0
	ant_positions_err = ant_positions  + error*ant_positions_err
        ant_p = np.array(ant_positions_err)
	nAntennas = len(ant_p)

	
       

	blVectors, blPairs = ant.get_blvectors(ant_positions_err)
	#np.save('redund_blVectors.npy',[blVectors,blPairs])
	

	ublDict = {}
	for b in range(len(blVectors)):
     		# grouping togather all blPairs with same blVector
    		if ublDict.has_key(blVectors[b]): ublDict[blVectors[b]].append(blPairs[b])
    		else: ublDict[blVectors[b]] = [blPairs[b]]

	ublIndexDict = {antPair: i for i,antPairs in enumerate(ublDict.values()) for antPair in antPairs }
	ublVectors = np.array([ant_positions_err[antList[0][0]]-ant_positions_err[antList[0][1]] for antList in ublDict.values()])



	print "With", len(ant_positions_err), "antennas there are", len(ublDict.items()), "unique baselines."
	#return [ublVectors,ublIndexDict]

	
        
	ant1, ant2 =[],[]
	for k in range(len(blPairs)):
		ant1.append(np.array(blPairs)[k][0])
		ant2.append(np.array(blPairs)[k][1])

	vis_map =[]
	for i in range(len(ant1)):
		vis_map.append(ublIndexDict[(ant1[i],ant2[i])])
	
	
	visTrue = get_V_uvw(blVectors,RA_src,Dec_src,src_flux,0.5,6.0,zen_vec,170,230,ant1,ant2,ant_eps)
	
	print np.std(visTrue-v1, axis=1)/np.std(v1, axis=1), np.exp(2j*np.pi*error*(0.05/0.5)).imag
	np.save('obstrue_vis_' + repr(xdim) + 'x' + repr(ydim) + repr(error) + 'beam' +'.npy',[visTrue,ant1,ant2,vis_map])
	
	#print np.linalg.norm(v1-visTrue[0]),  np.exp(2j*np.pi*error*(0.05/0.5)).imag
	if  np.mean(np.std(visTrue-v1, axis=1)/np.std(v1, axis=1)) <= np.exp(2j*np.pi*error*(0.05/0.5)).imag*1.3:
		break

	#print blVectors[2]

        #sg.append(np.std(visTrue-v1)/np.std(v1))
	#bl_.append(blVectors[6][0])
        
	print " successfully calculated visibilities"
"""
	
"""
for pos in range(len(ant_positions)):
	plt.title('$\mathbf{r} + \delta \mathbf{r}(\sigma=0.05)$')
	plt.plot(ant_positions[pos][0],ant_positions[pos][1], '*', label = 'Perfect Redundant')
	#plt.plot(ant_positions_err[pos][0],ant_positions_err[pos][1], '.', label = 'Near-Redundant')
	plt.xlabel('E-W [m]')
	plt.ylabel('N-S [m]')
	plt.grid()

"""
