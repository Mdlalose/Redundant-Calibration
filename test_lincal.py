import numpy as np, math 
from numpy.linalg import inv
from matplotlib import pyplot as plt
from scipy import linalg as lin
import lincal as omnical
#import simulated_vis as vis


ant1,ant2,vis_map = np.load('ants_data_test_.npy')[0],np.load('ants_data_test_.npy')[1],np.load('ants_data_test_.npy')[2]
Vis_uniq = np.load('data_vis_uniq_test.npy')
Vis_data,g_true = np.load('data_vis_test.npy'), np.load('data_vis_test_data.npy')[1]
n_ants = np.arange(g_true.size)



n_unique = Vis_uniq.size
q = Vis_data


def  model_vis(g,sky,vis_map):
    #g = np.exp(eta_0)*(np.cos(phi_0)+1j*np.sin(phi_0))
    return np.conj(g[ant1])*g[ant2]*sky[vis_map]



def linncal_func(Vis,g_0,s_0,N_steps):
    Parameter ={}
    #chi_lin =np.zeros(N_steps)
    chi=[]
    for j in range(N_steps):
        P = omnical.lincal(Vis,model_vis(g_0,s_0,vis_map),g_0,s_0,omnical.config_matrix(2,n_ants.size,ant1,ant2,vis_map,g_0,s_0,q,n_unique),ant1,ant2,vis_map)
        Parameter[j] =P
        #print P
        
        s_0 = P[1]
        g_0 = P[0] #/np.mean(P[0])
        chi.append(P[3])
     
    return np.array([Parameter,chi])



def linncal2_func(Vis,eta_0,phi_0,s_0,N_steps):
    Parameter ={}
    
    for j in range(N_steps):
        P = omnical.lincal2(Vis,eta_0,phi_0,s_0,vis_map,ant1,ant2,n_ants)
        Parameter[j] =P
        #print P
        
        eta_0 = P[0]
        phi_0 = P[1]
        s_0 = P[2] #/np.mean(P[0])

    return Parameter

gain_0 = np.zeros(g_true.size,dtype ="complex")

for ant in range(g_true.size):
    
            eta= np.random.normal(0.0,1.0)
            amp = np.exp(eta)
            phase = np.random.uniform(0.0,2.0*math.pi)
            gain_0[ant]= amp*(np.cos(phase) +1j*np.sin(phase))
            
    
    
gain_0 = gain_0/np.mean(gain_0)      


g_0 = g_true  + 0.0001*np.random.randn(g_true.size) # gain_0 
g_0 = g_0/np.mean(g_0) 

s_0 = Vis_uniq  + 0.00001*np.random.randn(Vis_uniq.size)


def noise_dB(output,input):
    results = np.zeros(input.size)
    for j in range(input.size):
        results[j] = 10*math.log(abs(output[j]/input[j]))
    return results


n_iter = 100
data_inter = linncal_func(Vis_data,g_0,s_0,n_iter)[0]#nncal2_func(Vis_data,eta_0,phi_0,s_0,n_iter) 
chisqd  = np.array(linncal_func(Vis_data,g_0,s_0,n_iter)[1])

plt.plot(model_vis(data_inter[1][0],data_inter[1][1],vis_map).real/Vis_data.real,'.',label='$First Iteration $')


plt.plot(np.ones(len(model_vis(data_inter[n_iter-1][0],data_inter[n_iter-1][1],vis_map))),'r',label='CRf=1.0')

plt.ylabel(r'Correlation Residual Fraction')


plt.legend(loc = 'best')
plt.show()



plt.plot(model_vis(data_inter[n_iter-1][0],data_inter[n_iter-1][1],vis_map).real/Vis_data.real,'.',label='100th Iteration ')


plt.plot(np.ones(len(model_vis(data_inter[n_iter-1][0],data_inter[n_iter-1][1],vis_map))),'r',label='CRF=1.0')

plt.ylabel(r'Visibility Residual Fraction ')
#plt.axis([0,120,-20,20])
plt.title('Real Part')
plt.legend(loc = 'best')
plt.show()

plt.plot(model_vis(data_inter[n_iter-1][0],data_inter[n_iter-1][1],vis_map).imag/Vis_data.imag,'.',label='100th Iteration')


plt.plot(np.ones(len(model_vis(data_inter[n_iter-1][0],data_inter[n_iter-1][1],vis_map))),'r',label='CRF=1.0')

plt.ylabel(r'Visibility Residual Fraction')
plt.title('Imaginary Part')
#plt.axis([0,120,-20,20])
plt.legend(loc = 'best')
plt.show()

#plt.plot(model_vis(data_inter[n_iter-1][0],data_inter[n_iter-1][1],vis_map)[16:len(model_vis(data_inter[n_iter-1][0],data_inter[n_iter-1][1],vis_map))].real/Vis_data[16:len(model_vis(data_inter[n_iter-1][0],data_inter[n_iter-1][1],vis_map))].real,'.',label='$n_{iter}=3$')

#plt.plot(np.ones(len(model_vis(data_inter[n_iter-1][0],data_inter[n_iter-1][1],vis_map)[16:len(model_vis(data_inter[n_iter-1][0],data_inter[n_iter-1][1],vis_map))])),label='line_y = 1')

#plt.title('Change ant0 D=6.0 to D=6.5 $\lambda$ = 0.5 to  $\lambda$ = 0.3 ')
#plt.title('D=6.0m to (D_0,D_03,D_12,D_15)=(6.10,6.10,6.10,6.10)m')
#plt.title('Perfect Redundancy')
#plt.title('D=6.0 to D_0=6.10')
#plt.title(r'$V_{i,j}, i,j for (1,...,15)$')
#plt.ylabel(r'Re($V_{ouptput}/V_{sim}$)')
#plt.xlabel('Visibility Points')

#plt.legend(loc = 'best')
#plt.show()

plt.plot(Vis_data.real,Vis_data.imag,'<r',label='Simulated Visibilities')
plt.plot(model_vis(data_inter[1][0],data_inter[n_iter-1][1],vis_map).real,model_vis(data_inter[1][0],data_inter[1][1],vis_map).imag,'*',label='Best Fit Visibities$ ')
#plt.errorbar(model_vis(data_inter[1][0],data_inter[n_iter-1][1],vis_map).real,model_vis(data_inter[1][0],data_inter[1][1],vis_map).imag, yerr= diag, xerr=diag)
plt.ylabel(r'Imaginary Part ')
plt.xlabel('Real Part')
plt.title('First Iteration')

plt.legend(loc = 'best')
plt.grid(True)
plt.axis('equal')
plt.show()


plt.plot(Vis_data.real,Vis_data.imag,'>r',label='Simulated Visibilities')
plt.plot(model_vis(data_inter[n_iter-1][0],data_inter[n_iter-1][1],vis_map).real,model_vis(data_inter[n_iter-1][0],data_inter[n_iter-1][1],vis_map).imag,'.k',label='Best Fit Visibilities 100th Iterations')
#plt.errorbar(model_vis(data_inter[n_iter-1][0],data_inter[n_iter-1][1],vis_map).real,model_vis(data_inter[n_iter-1][0],data_inter[n_iter-1][1],vis_map).imag, yerr= diag ,xerr=diag)

plt.ylabel(r'Imaginary Part ')
plt.xlabel('Real Part')
#plt.title('3rd Iteration')

plt.legend(loc = 'best')
plt.grid(True)
plt.axis('equal')
plt.show()

#plt.plot(Vis_data.real[0:16],Vis_data[0:16].imag,'.',label='Simulated Visibility')

#plt.plot(model_vis(data_inter[n_iter-1][0],data_inter[n_iter-1][1],vis_map)[16:len(model_vis(data_inter[n_iter-1][0],data_inter[n_iter-1][1],vis_map))].real,model_vis(data_inter[n_iter-1][0],data_inter[n_iter-1][1],vis_map)[16:len(model_vis(data_inter[n_iter-1][0],data_inter[n_iter-1][1],vis_map))].imag,'.',label='$n_{iter}=3$ $V_{i,j}=i,j(1,15)$ ')

#plt.ylabel(r'$Im(V)$')
#plt.xlabel('Re(V)')
#plt.title('Perfect Redundancy')

#plt.legend(loc = 'best')
#plt.show()


plt.plot(chisqd.real)
plt.xlabel('N_iterations')
plt.ylabel(r'$\chi^2$')
#plt.title('Perfect Redundancy')
#plt.title('D=6.0m to (D_0,D_03,D_12,D_15)=(6.10,6.10,6.10,6.10)m')
#plt.title('D=6.0 to D_0=6.10')
#plt.title('Change ant0 D=6.0 to D=6.5 & $\lambda$ = 0.5 to  $\lambda$ = 0.3 ')
plt.show()

plt.plot(g_true.real, g_true.imag,'<', label ='Simulated Input gains factors')
plt.plot(data_inter[n_iter-1][0].real,data_inter[n_iter-1][0].imag,'.', label ='Best Fit gain factors')

plt.xlabel('Real Part')
plt.ylabel('Imaginary Part')
plt.legend(loc="best")
plt.show()

plt.plot(Vis_uniq.real,Vis_uniq.imag,'*',label='Simulated True Correlations')
plt.plot(data_inter[n_iter-1][1].real,data_inter[n_iter-1][1].imag,'v',label='Best Fit True Correlations')
plt.xlabel('Real Part')
plt.ylabel('Imaginary Part')
plt.legend(loc='best')
plt.show()

def rec_sim(r_,i_):
    
     n_r=[]
     n_i=[]
     n_ri=[]
     for i in range(r_.size):
            if 0.5<= r_[i] <=1.5:
                k=1
                n_r.append(k)
                
            else:
                pass
                
                
     for j in range(i_.size):
          if 0.5<= i_[i] <= 2.5:
                k=1
                n_i.append(k)
          else:
                pass
     for c in range(r_.size):
        if 0.5< r_[i] <2.5 and 0.5<= i_[i] <= 2.5:
                k=1
                n_ri.append(k)
                
        else:
            pass
             
            
    
     return  [float(len(n_r))/float(len(r_)),float(len(n_i))/float(len(i_)),float(len(n_ri))/float(len(r_))]  
res_corr_i= model_vis(data_inter[n_iter-1][0],data_inter[n_iter-1][1],vis_map).imag/Vis_data.imag
res_corr_r= model_vis(data_inter[n_iter-1][0],data_inter[n_iter-1][1],vis_map).real/Vis_data.real

print rec_sim(res_corr_r,res_corr_i)



