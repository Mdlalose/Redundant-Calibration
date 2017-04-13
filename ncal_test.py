import numpy as np, math
import matplotlib.pyplot as plt
from scipy.optimize import minimize as min_newton_cg
import get_chi2_func
import lincal
import sys


xdim = int(sys.argv[1])
ydim =  int(sys.argv[2])
    
def chi_func(x,*data):
     
     gains_p = x[0:gains.size]
     sky_p = x[gains.size:x.size]
     return get_chi2_func.get_chisqd(np.array(data),gains_p,sky_p,ant1,ant2,vis_map)
    

def chi_func_grad(x,*data):
     
     gains_p = x[0:gains.size]
     sky_p = x[gains.size:x.size]
     return get_chi2_func.get_grad_chisqd(np.array(data),gains_p,sky_p,ant1,ant2,vis_map)
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
data = [] #np.zeros(len(data_vis[0]))
for v_nu in range(len(data_vis[0])):
	data.append(np.array(np.conj(gains[data_vis[1]])*gains[data_vis[2]]*data_vis[0][v_nu][data_vis[3]])) 

noise_frac_gains = float(sys.argv[3])
noise_frac_sky = float(sys.argv[4])

nu = np.arange(170.0,231,2.0)
ncal_sol =[]
for i in range(len(data)):
	x0= np.concatenate((gains + noise_frac_gains*np.ones(gains.size),data_vis[0][i] + noise_frac_sky*np.ones(data_vis[0][i].size)))
	x_fits = min_newton_cg(chi_func,x0,args=tuple(data[i]),callback=callbackF,method='Newton-CG',jac=chi_func_grad,hess=chi_func_curv,tol = 1e-6,options={'maxiter':100,'disp': True})
	best_fit_gains = x_fits.x[0:gains.size]
        best_fit_sky = x_fits.x[gains.size:x_fits.x.size]
	chi2 = x_fits.fun
        ncal_sol.append([best_fit_gains,best_fit_sky,chi2])
	plt.title('100% Noise')
	plt.plot(nu[i],chi2,'*')
	plt.xlabel('Frequency')
	plt.ylabel('$\chi^2$')


sky = data_vis[0][0]
best_fit_sky = ncal_sol[0][1]
best_fit_gains = ncal_sol[0][0]

gains = gains/np.average(gains)
best_fit_gains = best_fit_gains/np.average(best_fit_gains)

plt.figure()
plt.title('100% Noise Level at 170 MHz')
plt.quiver(np.real(sky), np.imag(sky), np.real(best_fit_sky-sky), np.imag(best_fit_sky-sky),angles='xy', scale_units='xy', scale=1)

plt.plot(np.real(sky), np.imag(sky),'r.',label='True Visibilities')
plt.plot(np.real(best_fit_sky), np.imag(best_fit_sky),'b.',label='Visibility Solution' )
plt.xlabel('Real Part')
plt.ylabel('Imaginary Part')
plt.legend(loc='best')	

plt.figure()
plt.title('100% Noise Level at 170 MHz')
plt.quiver(np.real(gains), np.imag(gains), np.real(best_fit_gains-gains), np.imag(best_fit_gains-gains),angles='xy', scale_units='xy', scale=1)
#plt.quiver(np.real(sky), np.imag(sky), np.real(best_fit_sky_err-sky), np.imag(best_fit_sky_err-sky),angles='xy', scale_units='xy', scale=1)
plt.plot(np.real(gains), np.imag(gains),'r.',label='True Gains')
plt.plot(np.real(best_fit_gains), np.imag(best_fit_gains),'b.',label='Gain Solution' )
plt.xlabel('Real Part')
plt.ylabel('Imaginary Part')
plt.legend(loc='best')	

plt.figure()
plt.title('100% Noise Level 170 MHz')
#plt.plot(np.abs(trueVis), np.abs(logcalVisSols))
plt.plot(np.abs(gains),np.abs(best_fit_gains),'.')
#plt.plot(np.abs(gains),np.abs(best_fit_gains_err),'.',label='$N(0.0,05cm)$')
plt.xlabel('Amp(True Gains)'); plt.ylabel('Amp(Gain Solution)')
plt.plot([0, 1.1*np.max(np.abs(gains))], [0, 1.1*np.max(np.abs(gains))], '--k')
plt.legend(loc='best')

plt.figure()
plt.title('100% Noise Level at 170 MHz')
#plt.plot(np.abs(trueVis), np.abs(logcalVisSols))
plt.plot(np.abs(gains)-np.abs(best_fit_gains),'.')
#plt.plot(np.abs(gains),np.abs(best_fit_gains_err),'.',label='$N(0.0,05cm)$')
plt.ylabel('Amp(True Gains)-Amp(Gain Solution)')
#plt.plot([0, 1.1*np.max(np.abs(gains))], [0, 1.1*np.max(np.abs(gains))], '--k')
plt.legend(loc='best')	

