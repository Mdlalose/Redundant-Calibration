import numpy as np
import matplotlib.pyplot as plt

gaain_true, Gain_Sol = np.load('GAIN_SOL_beamvar_0.05_50.npy')[0], np.load('GAIN_SOL_beamvar_0.05_50.npy')[1]

Gain=[]

for v in range(len(Gain_Sol)):
	Gain.append(Gain_Sol[v]/np.mean(Gain_Sol[v]))


Gain_Sol = Gain
	

n_data =  np.array(Gain_Sol[0]).size
phase_autcorr = np.zeros(n_data-1)
amp_autcorr = np.zeros(n_data-1)
Phase_autocorr = []
Amp_autocorr =[]
for ant in range(len(Gain_Sol)):

	phase_autocorr ={}
	amp_autocorr ={}
	for k in range(n_data-1):
		phase = np.arctan2(np.real(Gain_Sol[ant]),np.imag(Gain_Sol[ant])) 
		phase = phase - np.mean(phase)
		amp = np.sqrt(np.real(Gain_Sol[ant])**2 + np.imag(Gain_Sol[ant])**2)
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
	
	


plt.title('Average over all antennas')	
plt.plot(Phase_autocorr_avg)
plt.ylabel('Phase Autocorr')
plt.show()
plt.title('Average over all antennas')
plt.plot(Amp_autocorr_avg)
plt.ylabel('Amplitude Autocorr')
plt.show()


