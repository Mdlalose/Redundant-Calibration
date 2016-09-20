#Classical omnical

import numpy as np, math 
from numpy.linalg import inv
from matplotlib import pyplot as plt
from scipy import linalg as lin
#import simulated_vis as sim_data


def config_matrix(case,n_ants,ant1,ant2,vis_map,gain_param,sky_param,q,n_unique):
    # real vis, gains and sky'R==1'
    if case is 1:
        A = np.zeros((q.size, n_ants + n_unique))
        for vis_ in range(q.size):
            A[vis_][n_unique +ant1[vis_]]   = gain_param[ant2[vis_]]*sky_param[vis_map[vis_]]
            A[vis_][n_unique + ant2[vis_]]  = gain_param[ant1[vis_]]*sky_param[vis_map[vis_]]
            A[vis_][vis_map[vis_]]          = gain_param[ant1[vis_]]*gain_param[ant2[vis_]]
        return A
        
    else:
        #complex vis, gains and sky
        A = np.zeros((q.size, n_ants + n_unique),dtype ="complex")
        for vis_ in range(q.size):
            A[vis_][ant1[vis_]]   = gain_param[ant2[vis_]]*sky_param[vis_map[vis_]]
            A[vis_][ant2[vis_]]  = gain_param[ant1[vis_]]*sky_param[vis_map[vis_]]
            A[vis_][n_ants + vis_map[vis_]]          =      np.conj(gain_param[ant1[vis_]])*gain_param[ant2[vis_]]
        return A

def chi_squared(data,model_sky):
    return (data -model_sky).T.dot(data-model_sky)


def chi_squared_lin(data,model_sky,A,X):
    return (data -(model_sky + A.dot(X))).T.dot(data-(model_sky + A.dot(X)))

def Pinv(Matrix,eigen_threshold=10**-6):
    #Matrix = Matrix.conjugate()
    U,s_value,V= lin.svd(Matrix)
    S_inv= np.zeros((len(s_value),len(s_value)))
    
    for ss in range(len(s_value)):
        #print s_value[ss]
        if s_value[ss] <= eigen_threshold :
            #print s_value[ss]
            S_inv[ss][ss] = 0.0

        else:
            s_inv =1.0/s_value[ss]
            #print s_inv
            S_inv[ss][ss]= s_inv
            

    


    return V.T.dot(S_inv).dot(U.T)

def  model_vis(g,sky,ant1,ant2,vis_map):
    return np.conj(g[ant1])*g[ant2]*sky[vis_map]



def lincal(data,model_sky,g_0,s_0,A,ant1,ant2,vis_u_map):
     
     Curvature_matrix = A.T.dot(A)
     #delta_X = Pinv(Curvature_matrix).dot(A.T).dot(data-model_sky)
     delta_X = lin.pinv(Curvature_matrix).dot(A.T).dot(data-model_sky)
     g_1 = g_0 + delta_X[0:g_0.size]
     s_1 = s_0 +  delta_X[g_0.size:len(delta_X)]
     chi= get_chisqd(g_1,data,s_1,ant1,ant2,vis_u_map)
     return np.array([g_1,s_1,delta_X,chi])




def get_chisqd(gain_0,data,v_true,ant1,ant2,vis_u_map,fac=1.0):
      gain_0 =gain_0/fac
      chi = 0.0
      for i in range(data.size):
         
         chi += np.linalg.norm((data[i] - gain_0[ant1[i]]*gain_0[ant2[i]]*v_true[vis_u_map[i]]))**2
            
      return chi
            

        
def get_grad(gain_0,data,v_true,ant1,ant2,fac=1.0):
    gain_0 =gain_0/fac
    grad =np.zeros(np.array(gain_0).size,dtype = "complex")
    for i in range(data.size):
        grad[ant1[i]] +=  -2.0*(v_true[i]-gain_0[ant1[i]]*np.conj(gain_0[ant2[i]])*data[i])*np.conj(gain_0[ant2[i]])*data[i]
        grad[ant2[i]] +=  -2.0*(v_true[i]-gain_0[ant1[i]]*np.conj(gain_0[ant2[i]])*data[i])*np.conj(gain_0[ant1[i]])*data[i]
        
    return grad
        
     
def newton_cal(n_ants,ant1,ant2,vis_map,gain_param,sky_param,q,n_unique,tol,n_iter):
    
        x_0 = np.concatenate(gain_param,sky_param)
        for k in range(n_iter):
            grad_matrix     = config_matrix(n_ants,ant1,ant2,vis_map,x_0[0:n_ants.size],x_0[n_ants.size:len(x_0)],q,n_unique)
            Curvature_matrix = config_matrix(n_ants,ant1,ant2,x_0[0:n_ants.size],x_0[n_ants.size:len(x_0)],q,n_unique).T.dot(config_matrix(n_ants,ant1,ant2,vis_map,gain_param,sky_param,q,n_unique))
            x_n= x_0 - lin.pinv(Curvature_matrix ).dot(grad_)
           
            
            if lin.norm(x_n- x_0)<= tol:
                return x_n
            else:
                x_0=x_n
            

            
            
             


    
def B_matrix(n_ants,ant1,ant2,vis_map,gain_param,sky_param,q,n_unique): 
        A = np.zeros((q.size, 2.0*n_ants.size + n_unique),dtype ="complex")
        for vis_ in range(q.size):
            A[vis_][vis_map[vis_]]          =  np.conj(gain_param[ant1[vis_]])*gain_param[ant2[vis_]]
            A[vis_][n_unique +ant1[vis_]]   = np.conj(gain_param[ant1[vis_]])*gain_param[ant2[vis_]]*sky_param[vis_map[vis_]]
            
            A[vis_][n_unique + ant2[vis_]]  = np.conj(gain_param[ant1[vis_]])*gain_param[ant2[vis_]]*sky_param[vis_map[vis_]]
            A[vis_][n_unique + n_ants +ant1[vis_]]   = -1j*np.conj(gain_param[ant1[vis_]])*gain_param[ant2[vis_]]*sky_param[vis_map[vis_]]
            A[vis_][n_unique + n_ants +ant2[vis_]]   = 1j*np.conj(gain_param[ant1[vis_]])*gain_param[ant2[vis_]]*sky_param[vis_map[vis_]]
            
            
           
        return A
    
    
    
def lincal2(data,eta_0,phi_0,s_0,vis_u_map,ant1,ant2,n_ants):
    n_unique = s_0.size
    gain_0 = np.exp(eta_0)*(np.cos(phi_0)+1j*np.sin(phi_0))
    #gain_0 = gain_0/np.mean(gain_0)
    model =  np.conj(gain_0[ant1])*gain_0[ant2]*s_0[vis_u_map]
    
    B = B_matrix(n_ants,ant1,ant2,vis_u_map,gain_0,s_0,data,n_unique)
    Curvature = B.T.dot(B)
    
    delta_x = np.linalg.pinv(Curvature).dot(B.T).dot(data-model)
    s_1 = s_0 + delta_x[0:s_0.size]
    eta_1 = eta_0 + delta_x[s_0.size: s_0.size+ n_ants.size]
    phi_1 = phi_0 + delta_x[s_0.size +n_ants.size:len(delta_x)]
    
    return np.array([eta_1,phi_1,s_1,delta_x])
    
      