#Classical omnical

import numpy as np
from numpy.linalg import inv
from matplotlib import pyplot as plt
from scipy import linalg as lin
import simulated_vis as sim_data
xdim =3
ydim = 3
n_ants = xdim*ydim
sim_redu_baseline = sim_data.sim_redundant_baseline(xdim,ydim)

ant1 = sim_redu_baseline[0]
ant2 = sim_redu_baseline[1]
q= sim_redu_baseline[3]
vis_map = sim_redu_baseline[2]
q_unique = sim_redu_baseline[5]
n_unique = sim_redu_baseline[-1]
sim_vis_data = sim_data.sim_vis(n_ants,q,ant1,ant2,q_unique,vis_map)

gains_data = sim_vis_data[3]
sky_data =sim_vis_data[2]
vis_data = sim_vis_data[0]
cov_matrix =sim_vis_data[1]


def config_matrix(n_ants,ant1,ant2,vis_map,gain_param,sky_param,q,n_unique):
    
    A = np.zeros((q.size, n_ants + n_unique))
    for vis_ in range(q.size):
    
        # contructing A
        A[vis_][ant1[vis_]]  = gain_param[ant2[vis_]]*sky_param[vis_map[vis_]]
        
        A[vis_][ant2[vis_]]  = gain_param[ant1[vis_]]*sky_param[vis_map[vis_]]
        A[vis_][n_ants + vis_map[vis_]]     = gain_param[ant1[vis_]]*gain_param[ant2[vis_]]

    return A

def chi_squared(data,model_sky):
    return (data -model_sky).T.dot(data-model_sky)

sky_0= np.random.normal(0.0,0.5)*np.random.randn(q_unique.size)
gains_0= np.random.uniform(0.8,1.2)*np.random.randn(n_ants)

model_vis_data = gains_0[ant1]*gains_0[ant2]*sky_0[vis_map]

def chi_squared_lin(data,model_sky,A,X):
    return (data -(model_sky + A.dot(X))).T.dot(data-(model_sky + A.dot(X)))

def least_square(data,model_sky,A):
     Curvature_matrix = A.T.dot(A)
     delta_X = lin.pinv(Curvature_matrix).dot(A.T).dot(data-model_sky)
     
     return [delta_X,chi_squared(data,model_sky),chi_squared_lin(data,model_sky,A,delta_X)]


def  model_vis(g,sky,):
    return g[ant1]*g[ant2]*sky[vis_map]



    
    
def omnical_method(data,g_0,s_0,N_steps):

    for n_i in range(N_steps):
        delta_x = least_square(data,model_vis(g_0,s_0),config_matrix(n_ants,ant1,ant2,vis_map,g_0,s_0,q,n_unique))[0]
        g_update = g_0 + delta_x[0:n_ants]
        s_update = s_0 + delta_x[n_ants:len(delta_x)]
        av_g = np.mean(g_update)
        g_0 =  g_update/av_g
        s_0 =  s_update
    return [g_0,s_0,least_square(data,model_vis(g_0,s_0),config_matrix(n_ants,ant1,ant2,vis_map,g_0,s_0,q,n_unique))[0], least_square(data,model_vis(g_0,s_0),config_matrix(n_ants,ant1,ant2,vis_map,g_0,s_0,q,n_unique))[1], least_square(data,model_vis(g_0,s_0),config_matrix(n_ants,ant1,ant2,vis_map,g_0,s_0,q,n_unique))[2]]
