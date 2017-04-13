import numpy as np, math 
from numpy.linalg import inv
from matplotlib import pyplot as plt
from scipy import linalg as lin
import get_chi2_func
import sys


def get_config_matrix(case,n_ants,ant1,ant2,vis_map,gain_param,sky_param,q,n_unique):
    # real vis, gains and sky'R==1'
    # configaration matrix for real and complex cases
    if case is 1:
        # real data
        A = np.zeros((q.size, n_ants + n_unique))
        for vis_ in range(q.size):
            A[vis_][n_unique + ant1[vis_]]   = gain_param[ant2[vis_]]*sky_param[vis_map[vis_]]
            A[vis_][n_unique + ant2[vis_]]  = gain_param[ant1[vis_]]*sky_param[vis_map[vis_]]
            A[vis_][vis_map[vis_]]          = gain_param[ant1[vis_]]*gain_param[ant2[vis_]]
        return A
        
    else:
        #complex vis, gains and sky
        A = np.zeros((q.size, n_ants + n_unique),dtype ="complex")
        for vis_ in range(q.size):
            A[vis_][ n_unique + ant1[vis_]]   = gain_param[ant2[vis_]]*sky_param[vis_map[vis_]]
            A[vis_][n_unique + ant2[vis_]]  = np.conj(gain_param[ant1[vis_]])*sky_param[vis_map[vis_]]
            A[vis_][vis_map[vis_]]          =      np.conj(gain_param[ant1[vis_]])*gain_param[ant2[vis_]]
        return A


def get_ml_func(data,g,s,ant1,ant2,vis_map,n_steps,lambda_0=0.01,eps=0.1,l_up= 11.0,l_down =9.0):
	S=[]
        G=[]
        Chi2=[]
	for step in range(n_steps):
		grad = get_chi2_func.get_grad_chisqd(np.array(data),g,s,ant1,ant2,vis_map)
                Curv = get_chi2_func.get_curv_chisqd(np.array(data),g,s,ant1,ant2,vis_map)
                lambda_0 = lambda_0*np.max(np.diag(Curv))
		#print lambda_0,np.max(np.diag(Curv))
                
                mod_curv = Curv + lambda_0*np.diag(Curv)
           
                h_lm = np.linalg.pinv(mod_curv).dot(grad.T)
       
                g_1 = g + h_lm[s.size:len(h_lm)]
		s_1 = s + h_lm[0:s.size]
		#print g_1, s_1, h_lm
                #print h
                
		delta_chi2 = (get_chi2_func.get_chisqd(data,g,s,ant1,ant2,vis_map)-get_chi2_func.get_chisqd(data,g_1,s_1,ant1,ant2,vis_map))
		delta_h_lm = h_lm.T.dot(lambda_0*np.diag(Curv).dot(h_lm) + grad.T)
		#print delta_chi2
		x0 = np.concatenate((g,s))
		x1 = np.concatenate((g_1,s_1))
		alpha = 0.1
                if np.linalg.norm(x1-x0)<= 1e-4:
			g= g + alpha*h_lm[s.size:len(h_lm)]
			s = s + alpha*h_lm[0:s.size]
			#print 'yes'
                        #S.append(s)
                        #G.append(g)
                        #Chi2.append(get_chi2_func.get_chisqd(data,g,s,ant1,ant2,vis_map))

			#ambda_new = np.maximum(lambda_0/(1.0 + alpha),10**-7)
			lambda_0 = np.maximum(lambda_0/l_down, 10**-7)
                        
		else:
			g_1= g + alpha*h_lm[s.size:len(h_lm)]
			s_1 = s + alpha*h_lm[0:s.size] 
			#print 'no'

			#ambda_0 = lambda_0 + np.linalg.norm(get_chi2_func.get_chisqd(data,g_1,s_1,ant1,ant2,vis_map) -get_chi2_func.get_chisqd(data,g,s,ant1,ant2,vis_map))/(2.0*alpha)
                        lambda_0 = np.minimum(lambda_0*l_up,10**7)

                
		
                     
               

                     

      

        return [g_1,s_1,get_chi2_func.get_chisqd(data,g_1,s_1,ant1,ant2,vis_map)] 

def  model_vis(g,sky,vis_map):
  
    return np.conj(g[ant1])*g[ant2]*sky[vis_map]
###########################testing ml ##############################################################
xdim = int(sys.argv[1])
ydim = int(sys.argv[2])
gains = np.zeros(ydim**2,dtype='complex')
eta_0, phi_0 = np.zeros(ydim**2), np.zeros(ydim**2)
for g in range(ydim**2):
     eta= np.random.normal(0.0,1.0)
     eta_0[g]=eta
     
     amp = np.random.uniform(0.1,1.0)
     phase = np.random.uniform(0.0,math.pi)
     phi_0[g] = phase
     gains[g]=amp*(np.cos(phase)+ 1j*np.sin(phase))

noise_frac_gains =float(sys.argv[3])
noise_frac_sky = float(sys.argv[4])


n_iter = int(sys.argv[5])
#fit_p = get_ml_func(Vis_data,gains,sky_true,ant1,ant2,vis_map,n_iter)

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

nu = np.arange(170.0,231.0,2.0)
ncal_sol =[]
for i in range(len(data)):
	g_0 = gains + noise_frac_gains*np.ones(gains.size)
	sky_0 = data_vis[0][i] + noise_frac_sky*np.ones(data_vis[0][i].size)
	fit_p = get_ml_func(np.array(data[i]),g_0,sky_0,ant1,ant2,vis_map,n_iter)
        ncal_sol.append([fit_p[0],fit_p[1],fit_p[2]])
	plt.title('100% Noise')
	plt.plot(nu[i],fit_p[2],'*')
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


