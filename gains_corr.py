import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy.polynomial.polynomial as poly

#redundant data
gains, Gain_Sol = np.load('GAIN_SOL_beamvar_0.05_100.npy')[0], np.load('lincal_results_20_data_redundant88.npy')

Gain_Sol_20_bv = np.load('lincal_results_20_data_beamvar88.npy')
Gain_Sol_20_red = np.load('lincal_results_20_data_redundant88.npy')
gains_bv= np.load('GAIN_SOL_20_data_beamvar.npy')
gains_red = np.load('GAIN_SOL_20_data_redundant.npy')



#selecting gains
xdim =8
ydim =8
nu = np.arange(170.0,230,1.0)
#Gain_Sol = Gain_Sol_20_red
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
	ant_ii = tmp # np.column_stack(tmp)
	#print ant_ii.shape, ii

	GAIN_SOL[ii]=ant_ii


N_autocorr_amp={}
N_autocorr_phase={}
for ant_i in range(len(GAIN_SOL)):
	Autocorr_amp=[]
	Autocorr_phase =[]
	
	for ni in range(len(GAIN_SOL[ant_i])):
		GAIN_SOL_AVG = GAIN_SOL[ant_i][ni]/np.mean(GAIN_SOL[ant_i][ni])

		phase_autocorr ={}
		amp_autocorr ={}
		n_data = len(GAIN_SOL_AVG)
		for k in range(n_data-1):
			#phase = np.arctan2(np.real(GAIN_SOL_AVG[ant]),np.imag(GAIN_SOL_AVG[ant]))
			phase = np.angle(GAIN_SOL_AVG) 
			phase = phase - np.mean(phase)
			#amp = np.sqrt(np.real(GAIN_SOL_AVG[ant])**2 + np.imag(GAIN_SOL_AVG[ant])**2)
			amp = np.absolute(GAIN_SOL_AVG)
			amp = amp - np.mean(amp)
			tmp0=[]
			tmp1 =[]
			for n in range(n_data-k):
				tmp0.append(phase[n]*phase[n+k])
				tmp1.append(amp[n]*amp[n+k])
		
				
			phase_autocorr[k]=np.sum(tmp0)/(n_data-k)
			amp_autocorr[k]=np.sum(tmp1)/(n_data-k)	
		Autocorr_amp.append(amp_autocorr)
		Autocorr_phase.append(phase_autocorr)
	N_autocorr_amp[ant_i]= Autocorr_amp
	N_autocorr_phase[ant_i]= Autocorr_phase	


N_amp_autocorr =[]
N_phase_autocorr =[]

for i in range(len(N_autocorr_amp)):
	tmp1=[]
	tmp2=[]
	for k in range(len(N_autocorr_amp[i])):
		tmp1.append(N_autocorr_amp[i][k].values())
		tmp2.append(N_autocorr_phase[i][k].values())

	N_amp_autocorr.append(tmp1)
	N_phase_autocorr.append(tmp2)	

	
Amp_autocorr_avg=[]
Phase_autocorr_avg=[]	
for f in range(len(N_amp_autocorr)):
	g_0 = np.column_stack(N_amp_autocorr[f])
	phase_0 = np.column_stack(N_phase_autocorr[f])	
	N_amp_autocorr_avg=[]
	N_phase_autocorr_avg=[]
	for t in range(len(g_0)):
		N_amp_autocorr_avg.append(np.mean(g_0[t]))
		N_phase_autocorr_avg.append(np.mean(phase_0[t]))

	Amp_autocorr_avg.append(N_amp_autocorr_avg)
	Phase_autocorr_avg.append(N_phase_autocorr_avg)
Amp_autocorr_avg= np.column_stack(Amp_autocorr_avg)
Phase_autocorr_avg= np.column_stack(Phase_autocorr_avg)
AMP_AUTOCORR_AVG =[]
PHASE_AUTOCORR_AVG =[]
for p in range(n_data-1):
	AMP_AUTOCORR_AVG.append(np.mean(Amp_autocorr_avg[p]))
	PHASE_AUTOCORR_AVG.append(np.mean(Phase_autocorr_avg[p]))


def gauss_function(x, a, sigma):
    return a*np.exp(-(x)**2/(2*sigma**2))
k_phase = np.arange(len(PHASE_AUTOCORR_AVG))
k_amp = np.arange(len(AMP_AUTOCORR_AVG))
pp, cov = curve_fit(gauss_function,k_phase[0:10],PHASE_AUTOCORR_AVG[0:10])

pp2, cov2 = curve_fit(gauss_function,k_amp[0:10],AMP_AUTOCORR_AVG[0:10])

plt.figure()
plt.plot(AMP_AUTOCORR_AVG)
#plt.plot(gauss_function(k_amp[0:10],pp2[0],pp2[1]),'*', label = '$G_0$ $exp(-k^2/2\sigma_{k}^2)$')

plt.ylabel('Gain Amplitude Autocorrelation')
plt.xlabel('k (MHz)')
plt.tick_params(axis ='both')
plt.legend(loc='best')
plt.grid(True)
plt.show()
plt.figure()
plt.plot(PHASE_AUTOCORR_AVG)
plt.ylabel('Gain Phase Autocorrelation')
#plt.plot(gauss_function(k_phase[0:10],pp[0],pp[1]),'*',label = '$\Phi_0$ $exp(-k^2/2\sigma_{k}^2)$')
plt.xlabel('k (MHz)' )
plt.legend(loc='best')
plt.tick_params(axis ='both')
plt.grid(True)
plt.show()

"""
for ant_j in range(len(GAIN_SOL)):
	tmp0 = np.column_stack(GAIN_SOL[ant_j])
	
	tmp=[]
	for n_j in range(60):
		tmp.append(np.mean(tmp0[n_j]))
	


	plt.plot(nu,np.arctan2(np.real(tmp/np.mean(tmp)), np.imag(tmp/np.mean(tmp)))- np.arctan2(np.real(gains_red[ant_j]/np.mean(gains_red[ant_j])), np.imag(gains_red[ant_j]/np.mean(gains_red[ant_j]))),'o')
	plt.xlabel('Frequency (MHz)')
	plt.ylabel('Gain Phase Residual (radians)')
		
"""	

	
