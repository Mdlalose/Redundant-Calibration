#Classical omnical

import numpy as np, math 
from numpy.linalg import inv
from matplotlib import pyplot as plt
from scipy import linalg as lin
import simulated_vis as sim_data


def config_matrix(n_ants,ant1,ant2,vis_map,gain_param,sky_param,q,n_unique):
    
    A = np.zeros((q.size, n_ants + n_unique))
    for vis_ in range(q.size):
    
        # contructing A
        A[vis_][n_unique +ant1[vis_]]  = gain_param[ant2[vis_]]*sky_param[vis_map[vis_]]
        
        A[vis_][n_unique + ant2[vis_]]  = gain_param[ant1[vis_]]*sky_param[vis_map[vis_]]
        A[vis_][vis_map[vis_]]     = gain_param[ant1[vis_]]*gain_param[ant2[vis_]]

    return A

def chi_squared(data,model_sky):
    return (data -model_sky).T.dot(data-model_sky)


def chi_squared_lin(data,model_sky,A,X):
    return (data -(model_sky + A.dot(X))).T.dot(data-(model_sky + A.dot(X)))

def Pinv(Matrix,eigen_threshold=10**-5):
    #Matrix = Matrix.conjugate()
    U,s_value,V= lin.svd(Matrix)
    S_inv= np.zeros((len(s_value),len(s_value)))
    
    for ss in range(len(s_value)):
        if s_value[ss] <= eigen_threshold :
            #print s_value[ss]
            S_inv[ss][ss] = 0.0

        else:
            s_inv =1.0/s_value[ss]
            #print s_inv
            S_inv[ss][ss]= s_inv
            

    


    return V.T.dot(S_inv).dot(U.T)

def  model_vis(g,sky,ant1,ant2,vis_map):
    return g[ant1]*g[ant2]*sky[vis_map]



def lincal(data,model_sky,g_0,s_0,A):
     Curvature_matrix = A.T.dot(A)
     delta_X = Pinv(Curvature_matrix).dot(A.T).dot(data-model_sky)
     g_0 = g_0 + delta_X[s_0.size:len(delta_X)]
     s_0 = s_0 +  delta_X[0:s_0.size]
     
     
     return [g_0,s_0,delta_X,chi_squared(data,model_sky),chi_squared_lin(data,model_sky,A,delta_X)]






    
    
