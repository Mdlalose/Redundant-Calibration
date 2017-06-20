import numpy as np, math
import matplotlib.pyplot as plt
from scipy.optimize import minimize as min_newton_cg
import get_chi2_func
import lincal
import sys
import ml_lincal as ml

xdim = int(sys.argv[1])
ydim =  int(sys.argv[2])

def chi_func(x,*data):

     gains_p = x[0:gains.size]
     sky_p = x[gains.size:x.size]
     return get_chi2_func.get_chisqd(np.array(data),gains_p,sky_p,ant1,ant2,vis_map)


def chi_func_grad(x,*data):
     get_chi2_func.get_grad_chisqd(np.array(data),gains_p,sky_p,ant1,ant2,vis_map)
     gains_p = x[0:gains.size]
     sky_p = x[gains.size:x.size]
     return
lambda_0=0.01
l_down =9.0
lambda_0 = np.maximum(lambda_0/l_down, 10**-7)
def chi_func_curv(x,*data):

     gains_p = x[0:gains.size]
     sky_p = x[gains.size:x.size]
     return get_chi2_func.get_curv_chisqd(np.array(data),gains_p,sky_p,ant1,ant2,vis_map) + lambda_0*np.diag(get_chi2_func.get_curv_chisqd(np.array(data),gains_p,sky_p,ant1,ant2,vis_map))
Nfeval = 1
def callbackF(x,*data):
    #chi2=[]
    global Nfeval
    #print '{0:4d} {1:3.6f}'.format(Nfeval,chi_func(x,np.array(data)))
    #chi2.append(chi_func(x)
    Nfeval +=1

gains = np.zeros(ydim**2,dtype='complex')
eta_0, phi_0 = np.zeros(ydim**2), np.zeros(ydim**2)
for g in range(ydim**2):
     eta= np.random.normal(0.0,1.0)
     eta_0[g]=eta

     amp = np.random.uniform(0.1,1.0)
     phase = np.random.uniform(0.0,math.pi)
     phi_0[g] = phase
     gains[g]=amp*(np.cos(phase)+ 1j*np.sin(phase))

data_vis= np.load('obstrue_vis_8x8.npy')
ant1,ant2, vis_map = data_vis[1], data_vis[2], data_vis[3]
data_vis_0_05= np.load('obstrue_vis_8x820_data_redundant.npy')


noise_frac_gains = float(sys.argv[3])
noise_frac_sky = float(sys.argv[4])
#data = np.array(np.conj(gains[ant1])*gains[ant2]*data_vis_0_05[0][0])
nu = np.arange(170.0,230,1.0)
N_sol =[]
for n in range(1):

	data =  []  #[]#np.array(np.conj(gains[ant1])*gains[ant2]*sky[0][vis_map])#np.zeros(len(data_vis[0]))		
	for v_nu in range(1):
   		data.append(np.array(np.conj(gains[ant1])*gains[ant2]*data_vis_0_05[0][n][v_nu]))

	ncal_sol =[]
	n_iter = 20

	for i in range(1):
		g_0=gains + noise_frac_gains*np.random.randn(len(gains))
		s_0 = data_vis[0][i] + noise_frac_sky*np.random.randn(len(data_vis[0][i]))
		fitp= np.concatenate((g_0,s_0))
		tmp = []
		for ii in range(n_iter):
			chi2 = get_chi2_func.get_chisqd(data[i],fitp[0:g_0.size],fitp[g_0.size:],ant1,ant2,vis_map)
        		tmp.append(chi2)
			print chi2
			if chi2< 1e-14:
				break
			grad = get_chi2_func.get_grad_chisqd(data[i],fitp[0:g_0.size],fitp[g_0.size:],ant1,ant2,vis_map)
			curv = get_chi2_func.get_curv_chisqd(data[i],fitp[0:g_0.size],fitp[g_0.size:],ant1,ant2,vis_map)
		
			if ii<5:
				fac =0.25
				fitp=fitp - fac*grad/np.diag(curv)

			else:

				fac =0.8
				fitp = fitp - fac*grad/np.diag(curv)

		
			#print chi2
        	ncal_sol.append([fitp[0:g_0.size]/np.mean(fitp[0:g_0.size]),fitp[g_0.size:],chi2])
	N_sol.append(ncal_sol)		

#np.save('lincal_results_20_data_redundant' + repr(xdim) + repr(ydim) + '.npy',N_sol)



sky = data_vis[0][0]
best_fit_sky = ncal_sol[0][1]
best_fit_gains = ncal_sol[0][0]




gains = gains/np.average(gains)
#best_fit_gains = best_fit_gains/np.average(best_fit_gains)

plt.figure()
plt.quiver(np.real(sky), np.imag(sky), np.real(best_fit_sky-sky), np.imag(best_fit_sky-sky),angles='xy', scale_units='xy', scale=1)
plt.title('1000% offset at 170MHz ')
plt.plot(np.real(sky), np.imag(sky),'r.',label='True Visibilities')
plt.plot(np.real(best_fit_sky), np.imag(best_fit_sky),'*',label='Visibility Solution' )
plt.legend(loc='best')
plt.xlabel('Real Part')
plt.ylabel('Imaginary Part')
plt.grid(True)


plt.figure()
plt.title('1000% offset at 170MHz ')
plt.plot(np.real(gains),np.real(best_fit_gains),'r.')
#plt.plot(np.real(best_fit_sky), np.imag(best_fit_sky),'*',label='Visibility Solution' )
plt.legend(loc='best')
plt.xlabel('input Amplitude')
plt.ylabel('Output Amplitude')
plt.grid(True)

plt.figure()
plt.title('1000% offset at 170MHz ')
plt.plot(np.angle(gains),np.angle(best_fit_gains),'r.')
#plt.plot(np.real(best_fit_sky), np.imag(best_fit_sky),'*',label='Visibility Solution' )
plt.legend(loc='best')
plt.xlabel('input Phase (radians)')
plt.ylabel('Output Phase (radians)')
plt.grid(True)

"""
print 'finish calculating lincal'
Gains_Sol = []
SKY_sol =[]
#for q_i in range(len(N_sample_sol)):

Gain_Sol ={}
Sky_sol ={}
ant = np.arange(xdim*ydim)
n_unique = np.arange(len(data_vis[0][0]))
for ii in range(len(ant)):
		tmp=[]
		for nu_i in range(len(nu)):
			for i in range(len(ncal_sol[nu_i][0])):
				if ii == i:
					tmp.append(ncal_sol[nu_i][0][i])
		Gain_Sol[ii]= tmp





np.save('GAIN_SOL_20_redu.npy',[gains,Gain_Sol])
#np.save('SKY_SOL_20_redu.npy',SKY_sol)
"""			
"""
plt.figure()
plt.title('ANt0 50% offset')
plt.plot(nu, np.sqrt(np.real(Gain_Sol[0])**2+ np.imag(Gain_Sol[0])**2)-np.sqrt(np.real(gains[0])**2+ np.imag(gains[0])**2),'o')	
plt.ylabel('Amplitude Residual')
plt.xlabel('Frequency [MHz]')


plt.figure()
plt.title('ANt0 50% offset')
plt.plot(nu, np.arctan2(np.real(Gain_Sol[0]),np.imag(Gain_Sol[0]))-np.arctan2(np.real(gains[0]),np.imag(gains[0])),'o')	
plt.ylabel('Phase Residual [radians]')
plt.xlabel('Frequency [MHz]')
"""


"""
ax2.plot(nu, np.sqrt(np.real(Gain_Sol[7])**2+ np.imag(Gain_Sol[7])**2)-np.sqrt(np.real(gains[7])**2+ np.imag(gains[7])**2),'o')	

	


ax3.plot(nu, np.sqrt(np.real(Gain_Sol[56])**2+ np.imag(Gain_Sol[56])**2)-np.sqrt(np.real(gains[56])**2+ np.imag(gains[56])**2),'o')	



ax4.plot(nu, np.sqrt(np.real(Gain_Sol[63])**2+ np.imag(Gain_Sol[63])**2)-np.sqrt(np.real(gains[63])**2+ np.imag(gains[63])**2),'o')	

ax.set_xlabel('Amplitude Residual')
ax.set_ylabel('Frequency [MHz]')
ax1.set_title('#Ant0')
ax2.set_title('#Ant7')
ax3.set_title('#Ant56')
ax4.set_title('#Ant63')

fig =plt.figure()
ax = fig.add_subplot(111)    # The big subplot
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)




ax1.plot(nu, np.arctan2(np.real(Gain_Sol[0]),np.imag(Gain_Sol[0]))-np.arctan2(np.real(gains[0]),np.imag(gains[0])),'o')	



ax2.plot(nu, np.arctan2(np.real(Gain_Sol[7]),np.imag(Gain_Sol[7]))-np.arctan2(np.real(gains[7]),np.imag(gains[7])),'o')	



ax3.plot(nu, np.arctan2(np.real(Gain_Sol[56]),np.imag(Gain_Sol[56]))-np.arctan2(np.real(gains[56]),np.imag(gains[56])),'o')	




ax4.plot(nu, np.arctan2(np.real(Gain_Sol[63]),np.imag(Gain_Sol[63]))-np.arctan2(np.real(gains[63]),np.imag(gains[63])),'o')	
ax.set_xlabel('Phase Residual [radians]')
ax.set_ylabel('Frequency [MHz]')
ax1.set_title('#Ant0')
ax2.set_title('#Ant7')
ax3.set_title('#Ant56')
ax4.set_title('#Ant63')
"""


