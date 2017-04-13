import numpy as np, math
import matplotlib.pyplot as plt
from scipy.optimize import minimize as min_newton_cg
import get_chi2_func
import lincal
import sys

def chi_func(x):
     
     gains_p = x[0:gains.size]
     sky_p = x[gains.size:x.size]
     return get_chi2_func.get_chisqd(data,gains_p,sky_p,ant1,ant2,vis_map)
    

def chi_func_grad(x):
     
     gains_p = x[0:gains.size]
     sky_p = x[gains.size:x.size]
     return get_chi2_func.get_grad_chisqd(data,gains_p,sky_p,ant1,ant2,vis_map)

def chi_func_curv(x):
     
     gains_p = x[0:gains.size]
     sky_p = x[gains.size:x.size]
     return get_chi2_func.get_curv_chisqd(data,gains_p,sky_p,ant1,ant2,vis_map)


### ##################################### data test#############################################


xdim = int(sys.argv[1])
ydim =  int(sys.argv[2])
    



gains = np.zeros(ydim**2,dtype='complex')
eta_0, phi_0 = np.zeros(ydim**2), np.zeros(ydim**2)
for g in range(ydim**2):
     eta= np.random.normal(0.0,1.0)
     eta_0[g]=eta
     
     amp = np.random.uniform(0.1,1.0)
     phase = np.random.uniform(0.0,math.pi)
     phi_0[g] = phase
     gains[g]=amp*(np.cos(phase)+ 1j*np.sin(phase))






""" 
sky = np.zeros(n_unique,dtype='complex')

for s in range(n_unique):
     amp= np.random.normal(0.0,0.5)
     phase = np.random.uniform(0.0,math.pi)
     sky[s]=amp*(np.cos(phase)+ 1j*np.sin(phase))
     

"""
#data =np.conj(gains[ant1])*gains[ant2]*sky[vis_map] #+ 0.1*np.random.randn(q.size)

#gains_0 =gains + 0.1*np.random.randn(gains.size)
#sky_0 = sky + 0.1*np.random.randn(sky.size)

trueVis,ant1,ant2,vis_map= np.load('obstrue_vis_5x5.npy')[0],np.load('obstrue_vis_5x5.npy')[1],np.load('obstrue_vis_5x5.npy')[2],np.load('obstrue_vis_5x5.npy')[3]
#data=  np.conj(gains[ant1])*gains[ant2]*trueVis[vis_map]
obsVis = np.load('obstrue_vis_5x5_1.0m.npy')[0]
data = np.conj(gains[ant1])*gains[ant2]*obsVis

	

#data = np.load('obsvis_10x10_100cm.npy')
#gains=  np.load('gains_10x10_100cm.npy')
noise_frac_gains = float(sys.argv[3])
noise_frac_sky = float(sys.argv[4])
sky=trueVis
 #np.load('data_10x10.npy')



gains_0 = gains  + noise_frac_gains*np.ones(gains.size) # gain_0 
#g_0 = g_0/np.mean(g_0) 
g_0= gains
sky_0 = sky + noise_frac_sky*np.ones(sky.size)
x0= np.concatenate((gains_0,sky_0))
Nfeval =1

def callbackF(x):
    #chi2=[]
    global Nfeval
    print '{0:4d} {1:3.6f}'.format(Nfeval,chi_func(x))
    #chi2.append(chi_func(x)
    Nfeval +=1
    
    

print 'full chi square'
x_fits = min_newton_cg(chi_func,x0,callback=callbackF,method='Newton-CG',jac=chi_func_grad,hess=chi_func_curv,options={'xtol':1e-6,'disp': True})
print 'linearize chi sqaured'
#np.save('x_fits_err.npy',x_fits)
#x_lin_fits = min_newton_cg(get_lin_chi,x0,callback=callbackF,method='Newton-CG',jac=get_lin_grad,hess=get_lin_curv,options={'xtol':1e-10,'disp':True})
#print get_lin_chi(x0)
#plt.ion()

#x_fits = np.load('x_fits.npy')
#x_fits_err = np.load('x_fits_err.npy')
best_fit_gains = x_fits.x[0:gains.size]
best_fit_sky = x_fits.x[gains.size:x_fits.x.size]
#np.save('best_fit_err_5cm.npy',[best_fit_gains,best_fit_sky])
best_fits_data = np.conj(best_fit_gains[ant1])*best_fit_gains[ant2]*best_fit_sky[vis_map]

scatter_gains = np.array([np.std(np.load('best_fit_0.05m.npy')[0]),np.std(np.load('best_fit_0.15m.npy')[0]),np.std(np.load('best_fit_0.25m.npy')[0]),np.std(np.load('best_fit_0.35m.npy')[0]),np.std(np.load('best_fit_0.45m.npy')[0]),np.std(np.load('best_fit_0.55m.npy')[0]),np.std(np.load('best_fit_1.0m.npy')[0])])
scatter_sky = np.array([np.std(np.load('best_fit_0.05m.npy')[1]),np.std(np.load('best_fit_0.15m.npy')[1]),np.std(np.load('best_fit_0.25m.npy')[1]),np.std(np.load('best_fit_0.35m.npy')[1]),np.std(np.load('best_fit_0.45m.npy')[1]),np.std(np.load('best_fit_0.55m.npy')[1]),np.std(np.load('best_fit_1.0m.npy')[1])])

dp =np.array([0.05,0.15,0.25,0.35,0.45,0.55,1.0])

plt.figure()
plt.plot(dp,scatter_gains,'*')
plt.xlabel('$\delta \mathbf{b}$ [m]')
plt.ylabel('$\sigma (gain solution)$')
plt.figure()
plt.plot(dp,scatter_sky,'*')
plt.xlabel('$\delta \mathbf{b}$ [m]')
plt.ylabel('$\sigma (true visibility solution)$')


"""

p_err = np.sqrt(np.diag(np.linalg.pinv(get_chi2_func.get_curv_chisqd(data,best_fit_gains,best_fit_sky,ant1,ant2,vis_map))))
print 'variance', (np.std(best_fits_data-data)/np.std(best_fits_data))**2
np.save('best_fit_1.0m.npy',[best_fit_gains,best_fit_sky,p_err])
plt.figure()
plt.title('$\sigma_{antloc} = 5cm$')
plt.quiver(np.real(sky), np.imag(sky), np.real(best_fit_sky-sky), np.imag(best_fit_sky-sky),angles='xy', scale_units='xy', scale=1)

plt.plot(np.real(sky), np.imag(sky),'r.',label='True Visibilities')
plt.plot(np.real(best_fit_sky), np.imag(best_fit_sky),'b.',label='Visibility Solution' )
plt.xlabel('Real Part')
plt.ylabel('Imaginary Part')
plt.legend(loc='best')						
		

plt.figure()
plt.title(' $\sigma_{antloc} = 5cm$ ')
plt.quiver(np.real(gains), np.imag(gains), np.real(best_fit_gains-gains), np.imag(best_fit_gains-gains),angles='xy', scale_units='xy', scale=1)
#plt.quiver(np.real(sky), np.imag(sky), np.real(best_fit_sky_err-sky), np.imag(best_fit_sky_err-sky),angles='xy', scale_units='xy', scale=1)
plt.plot(np.real(gains), np.imag(gains),'r.',label='True Gains')
plt.plot(np.real(best_fit_gains), np.imag(best_fit_gains),'b.',label='Gain Solution' )

#plt.plot(np.real(best_fit_sky_err), np.imag(best_fit_sky_err),'k.',label='Calibrated Visibilities $N(0.0,5cm)$')
#plt.plot(np.real(visSols2), np.imag(visSols2),'g.',label='Calibrated Visibilities 2')
plt.xlabel('Real Part')
plt.ylabel('Imaginary Part')
plt.legend(loc='best')						
	
plt.figure()
#plt.plot(np.abs(trueVis), np.abs(logcalVisSols))
#plt.plot(np.angle(sky), np.angle(best_fit_sky),'.')
plt.plot(np.angle(best_fit_sky)-np.angle(sky),'.',label =' $\sigma_{antloc}(5cm)$')
#plt.plot(np.angle(sky),np.angle(best_fit_sky_err),'.',label='$N(0.0,05cm)$')
#plt.plot([-np.pi, np.pi], [-np.pi, np.pi], '--k')
plt.ylabel('Phase(Visibility Solution)-Phase(True Visibility)')
plt.legend(loc='best')

gains = gains/np.average(gains)
best_fit_gains = best_fit_gains/np.average(best_fit_gains)
plt.figure()
plt.title('10% Noise Level')
#plt.plot(np.abs(trueVis), np.abs(logcalVisSols))
plt.plot(np.abs(gains),np.abs(best_fit_gains),'.',label =' $\sigma_{antloc}(5cm)$')
#plt.plot(np.abs(gains),np.abs(best_fit_gains_err),'.',label='$N(0.0,05cm)$')
plt.xlabel('Amp(True Gains)'); plt.ylabel('Amp(Gain Solution)')
plt.plot([0, 1.1*np.max(np.abs(gains))], [0, 1.1*np.max(np.abs(gains))], '--k')
plt.legend(loc='best')

plt.figure()
#plt.plot(np.abs(trueVis), np.abs(logcalVisSols))
plt.plot(np.abs(gains)-np.abs(best_fit_gains),'.',label =' $\sigma_{antloc}(5cm)$')
#plt.plot(np.abs(gains),np.abs(best_fit_gains_err),'.',label='$N(0.0,05cm)$')
plt.ylabel('Amp(True Gains)-Amp(Gain Solution)')
#plt.plot([0, 1.1*np.max(np.abs(gains))], [0, 1.1*np.max(np.abs(gains))], '--k')
plt.legend(loc='best')	
"""
			
