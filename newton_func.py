#Newton multavaraite method
import numpy as np, math
import matplotlib.pyplot as plt
import get_chi2_func


def get_fits_param(n_steps,tol,data,gains,sky,ant1,ant2,vis_map):
    """ Netown's Multivaraite Method"""
    p=[]
    for step in range(n_steps):
              param_old= np.concatenate((gains,sky))
              #print param_old[0]
              
              grad = get_chi2_func.get_grad_chisqd(data,param_old[0:gains.size],param_old[gains.size:len(param_old)],ant1,ant2,vis_map)
              curv = get_chi2_func.get_curv_chisqd(data,param_old[0:gains.size],param_old[gains.size:len(param_old)],ant1,ant2,vis_map)
              param_new = param_old - np.linalg.pinv(curv).dot(grad)
	      print param_new
              if np.linalg.norm(param_new-param_old) <= tol:
                  p.append(param_new)
                  print param_new[0], 'good'
                  break
              else:
                  param_old = param_new
                  #print param_new[0], 'bad'
   
          
    return p   
 

def get_nwr_param(data,g,s,ant1,ant2,vis_map,n_steps,lambda_0=0.01,eps=0.1):
         
        
	S=[]
        G=[]
        Chi2=[]
	for step in range(n_steps):
		grad_chi2 = np.diag(get_chi2_func.get_grad_chisqd(data,g,s,ant1,ant2,vis_map))
                curv_chi2 = get_chi2_func.get_curv_chisqd(data,g,s,ant1,ant2,vis_map)
                lambda_0 = lambda_0*np.max(np.diag(curv_chi2))
		
                mod_curv = curv_chi2 + lambda_0*np.diag(curv_chi2)
           
                h_lm = np.linalg.pinv(mod_curv).dot(grad_chi2)
          
                g_1 = g + h_lm[s.size:len(h_lm)]
		s_1 = s + h_lm[0:s.size]
                #print h_lm
                
		delta_chi2 = (get_chi2_func.get_chisqd(data,g,s,ant1,ant2,vis_map)-get_chi2_func.get_chisqd(data,g_1,s_1,ant1,ant2,vis_map))
		delta_h_lm = h_lm.T.dot(lambda_0*np.diag(curv_chi2).dot(h_lm) + grad_chi2.T.dot(data-np.conj(g[ant1])*g[ant2]*s[vis_map]))
		
		rho_h_lm = delta_chi2/delta_h_lm
                #print rho_h_lm
		
                
		alpha_up = (grad_chi2.T.dot(data-np.conj(g[ant1])*g[ant2]*s[vis_map])).T.dot(h_lm)
                
                
		alpha_d = (get_chi2_func.get_chisqd(data,g_1,s_1,ant1,ant2,vis_map) -get_chi2_func.get_chisqd(data,g,s,ant1,ant2,vis_map))/2.0 + 2.0*(grad_chi2.T.dot(data-np.conj(g[ant1])*g[ant2]*s[vis_map])).T.dot(h_lm)
                 
                
		alpha =  alpha_up/alpha_d
                #print alpha
		
                if rho_h_lm> eps:
			g= g + alpha*h_lm[s.size:len(h_lm)]
			s = s + alpha*h_lm[0:s.size]
                        S.append(s)
                        G.append(g)
                        Chi2.append(get_chi2_func.get_chisqd(data,g,s,ant1,ant2,vis_map))

			lambda_new = np.maximum(lambda_0/(1.0 + alpha),10**-7)
			#lambda_0 = np.maximum(lambda_0/l_down, 10**-7)
                        
		else:
			g_1= g + alpha*h_lm[s.size:len(h_lm)]
			s_1 = s + alpha*h_lm[0:s.size]

			lambda_0 = lambda_0 + np.linalg.norm(get_chi2_func.get_chisqd(data,g_1,s_1,ant1,ant2,vis_map) -get_chi2_func.get_chisqd(data,g,s,ant1,ant2,vis_map))/(2.0*alpha)
                        #lambda_0 = np.minimum(lambda_0*l_up,10**7)

                
		
                     
               

                     

      

        return [G,S,Chi2] 

