import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy.polynomial.polynomial as poly

### BEAM IMPERFECTIONS
gain_true= np.load('GAIN_SOL_beamvar_0.05_100.npy')[0]
#sky, vis_map, best_fit_sky_beam = np.load('obstrue_vis_8x8.npy')[0],np.load('obstrue_vis_8x8beamvardata_100.npy')[3], 

Gain_Sol = np.load('lincal_results_100_sim88.npy')


# PERFECT CASE
gains = np.load('GAIN_SOL_100_redu.npy')[0]
#best_fit_sky_red =  np.load('SKY_SOL_perf_red.npy')
Gain_Sol_red = np.load('lincal_results_100_sim_redu88.npy')
xdim =8
ydim =8
nu = np.arange(170.0,230,1.0)
Gains_Sol =[]
for q in range(len(Gain_Sol)):
		
		#print len(Gain_Sol[q][n][0]),q,n
		Gain={}
		ant = np.arange(xdim*ydim)
		for ii in range(len(ant)):
			tmp=[]
			for nu_i in range(len(nu)):
				for i in range(len(Gain_Sol[q][nu_i][0])):
					if ii == i:
						#print q,ii,i,nu_i
						tmp.append(Gain_Sol[q][nu_i][0][i])
			Gain[ii]= tmp



		Gains_Sol.append(Gain)



GAIN_SOL ={}

ant = np.arange(xdim*ydim)

for ii in range(len(ant)):
	tmp = []
	for n in range(len(Gains_Sol)):
		for key in Gains_Sol[n].keys():
			if ii == key:
				tmp.append(Gains_Sol[n][key])
	ant_ii = np.column_stack(tmp)
	#print ant_ii.shape, ii

	GAIN_SOL[ii]=ant_ii



GAIN_SOL_AVG={}
for ant_i in range(len(ant)):
	tmp =[]
	for q in range(len(nu)):
		#print np.mean(GAIN_SOL[ant_i][q]), ant_i,q
		tmp.append(np.mean(GAIN_SOL[ant_i][q]))
	GAIN_SOL_AVG[ant_i]=tmp/np.mean(tmp)
	"""	
	plt.plot(nu,np.arctan2(np.real(tmp/np.mean(tmp)), np.imag(tmp/np.mean(tmp)))- np.arctan2(np.real(gain_true[ant_i]/np.mean(gain_true[ant_i])), np.imag(gain_true[ant_i]/np.mean(gain_tr[ant_i]))),'o')
	plt.xlabel('Frequency (MHz)')
	plt.ylabel('Gain Phase Residual (radians)')
	"""
	
	


"""
plt.title('230MHz ')
plt.plot(np.real(sky[-1]), np.imag(sky[-1]),'ro',label='True Visibilities')
plt.plot(np.real(best_fit_sky_red[0][-1][1]), np.imag(best_fit_sky_red[0][-1][1]),'*',label='Visibility Solution' )
plt.legend(loc='best')
plt.xlabel('Real Part')
plt.ylabel('Imaginary Part')
plt.grid(True)
"""
"""
sky_sample =[]

for s_i in range(len(best_fit_sky_beam)):
	tmp =[]
	for nu_i in range(len(best_fit_sky_beam[s_i])):
	        #plt.title('50% offset at 170MHz ')
		plt.plot(np.real(sky[nu_i]), np.imag(sky[nu_i]),'r.',label='True Visibilities')
		#plt.plot(np.real(best_fit_sky), np.imag(best_fit_sky),'*',label='Visibility Solution' )
		plt.legend(loc='best')
		plt.xlabel('Real Part')
		plt.ylabel('Imaginary Part')
		plt.grid(True)
	        #tmp.append(np.std(sky[nu_i]-best_fit_sky_beam[s_i][nu_i][1][vis_map])/np.std(sky[nu_i]))
		#sky_sample.append(tmp)
"""	

"""
sky_sample = np.column_stack(sky_sample)
		
for i in range(len(nu)):
	#print nu[i], np.std(sky[i]),len(best_fit_sky_red[i][1]) #, np.std(sky[i]-best_fit_sky_red[i])/np.std(sky[i])
	
	
	plt.plot(nu[i], np.mean(sky_sample[i]),'o')
	plt.xlabel('Frequency (MHz)')
	plt.ylabel('$\sigma$(True Visibility - Visibility Solution)/$\sigma$(True Visibility)') 
	

"""






	

n_data =  np.array(GAIN_SOL_AVG[0]).size
phase_autcorr = np.zeros(n_data-1)
amp_autcorr = np.zeros(n_data-1)
Phase_autocorr = []
Amp_autocorr =[]

for ant in range(len(GAIN_SOL_AVG)):

	phase_autocorr ={}
	amp_autocorr ={}
	for k in range(n_data-1):
		#phase = np.arctan2(np.real(GAIN_SOL_AVG[ant]),np.imag(GAIN_SOL_AVG[ant]))
		phase = np.angle(GAIN_SOL_AVG[ant]) 
		phase = phase - np.mean(phase)
		#amp = np.sqrt(np.real(GAIN_SOL_AVG[ant])**2 + np.imag(GAIN_SOL_AVG[ant])**2)
		amp = np.absolute(GAIN_SOL_AVG[ant])
		amp = amp - np.mean(amp)
		tmp0=[]
		tmp1 =[]
		for n in range(n_data-k):
			tmp0.append(phase[n]*phase[n+k])
			tmp1.append(amp[n]*amp[n+k])
		
				
		phase_autocorr[k]=np.sum(tmp0)/(n_data-k)
		amp_autocorr[k]=np.sum(tmp1)/(n_data-k)
	Phase_autocorr.append(phase_autocorr.values())
	Amp_autocorr.append(amp_autocorr.values())


		


Phase_autocorr = np.column_stack(Phase_autocorr)
Amp_autocorr =  np.column_stack(Amp_autocorr)

Phase_autocorr_avg = []
Amp_autocorr_avg = []

for row in range(n_data-1):
	Phase_autocorr_avg.append(np.mean(Phase_autocorr[row]))
	Amp_autocorr_avg.append(np.mean(Amp_autocorr[row]))
	

def sincfunc(x,A,sigma):
	return A*np.sin(np.pi*x*sigma)
def gauss_function(x, a, sigma):
    return a*np.exp(-(x)**2/(2*sigma**2))
xx = np.arange(len(Phase_autocorr_avg))
xx[0] = 0.01
k_phase = np.arange(len(Phase_autocorr_avg))
k_amp = np.arange(len(Phase_autocorr_avg))
pp, cov = curve_fit(gauss_function,k_phase[0:22],Phase_autocorr_avg[0:22])
pp3, cov3 = curve_fit(sincfunc,k_phase[22:],Phase_autocorr_avg[22:])
pp2, cov2 = curve_fit(gauss_function,k_amp[0:22],Amp_autocorr_avg[0:22])
k_phase = np.arange(len(Phase_autocorr_avg))
k_amp = np.arange(len(Phase_autocorr_avg))
k_amp_poly =k_amp[22:]
amp_coeff= poly.polyfit(k_amp_poly,Amp_autocorr_avg[22:],12)
new_k_amp = np.linspace(k_amp_poly[0],k_amp_poly[-1], num = len(k_amp_poly)*5)
amp_ffit = poly.polyval(new_k_amp,amp_coeff)

k_phase_poly =k_phase[22:]
phase_coeff= poly.polyfit(k_phase_poly,Phase_autocorr_avg[22:],14)
new_k_phase = np.linspace(k_phase_poly[0],k_phase_poly[-1], num = len(k_phase_poly)*5)
phase_ffit = poly.polyval(new_k_phase,phase_coeff)


#poly_phase= np.polyfit(k_amp[7:],Amp_autocorr_avg[7:],10)

#plt.title('Average over all antennas')	
plt.plot(k_phase,Phase_autocorr_avg)
#plt.plot(k_phase[7:],Phase_autocorr_avg[7:])
plt.plot(np.zeros(len(Amp_autocorr_avg)))
#plt.plot(poly_phase, label= '$k^5$')
plt.plot(gauss_function(k_phase[0:22],pp[0],pp[1]),'*',label = '$\Phi_0$ $exp(-k^2/2\sigma_{k}^2)$')
plt.plot(new_k_phase,phase_ffit,'*',label= '14th-polynomail')


plt.ylabel('$\Phi$(k)')
plt.xlabel('k (MHz)' )
plt.legend(loc='best')
plt.tick_params(axis ='both')
plt.show()
#plt.title('Average over all antennas')
plt.plot(Amp_autocorr_avg)
#plt.plot(k_amp[8:],Amp_autocorr_avg[8:])

plt.plot(gauss_function(k_amp[0:22],pp2[0],pp2[1]),'*', label = '$G_0$ $exp(-k^2/2\sigma_{k}^2)$')
plt.plot(np.zeros(len(Amp_autocorr_avg)))
plt.plot(new_k_amp,amp_ffit,'*',label= '12th-polynomail')
plt.ylabel('G(k)')
plt.xlabel('k (MHz)')
plt.tick_params(axis ='both')
plt.legend(loc='best')
plt.show()


