#simulated Visibilitity
import numpy as np
from numpy.linalg import inv
from matplotlib import pyplot as plt





def sim_redundant_baseline(xdim,ydim,plot_qu=False,n_unique_baseline =False):
    """" This function simulated redundant baseline for a given N= xdim*ydim Anttenna Array
    This function return an a column array ant1,ant2 index, vis_map an array which contains indes which labe visibility according to their redundant unique group set  e.g first redundant set are has an index 0. It also returns baline q and u, and unique baseline q_unique and, u_nique"""

    
    xx=np.arange(xdim)
    yy=np.arange(ydim)
    x,y=np.meshgrid(xx,yy) # put antennnas position into a xdim by ydim grid
    x=np.reshape(x,[x.size,1])
    y=np.reshape(y,[y.size,1])
    ant=x*xdim + y     # making a 1d list of antenna indexes from 0 to 48
    x1,x2=np.meshgrid(x,x) # putting  position x into a N by N grid
    ant1,ant2=np.meshgrid(ant,ant) # putting antennas into a grid
    y1,y2=np.meshgrid(y,y) 
    # computing baseline in x and y direction

    q=y2-y1
    u=x2-x1
    q=np.reshape(q,[q.size,1])
    u=np.reshape(u,[u.size,1])
    ant1=np.reshape(ant1,[ant1.size,1]) # reshaping ant1 into 1d of len 4
    ant2=np.reshape(ant2,[ant2.size,1]) #reshaping ant2 into 1d of len 4
    # removing repeated antennas
    isgood=ant1>ant2
    #selecting q,u are not repeated and corresponding ant1 and ant2
    q=q[isgood]
    u=u[isgood]
    ant1=ant1[isgood]
    ant2=ant2[isgood]
    if plot_qu is True:
        
        plt.ion()
        plt.plot(q,u,'*')
        plt.xlabel('q')
        plt.ylabel('u')
       #plt.axis([-10,10,-10,10])
    else:
        pass
        


    qu=q+np.complex(0,1)*u
    qu=np.unique(qu) # selecting unique baselines

    n_unique=(xdim)**2-1+(xdim-1)**2
    if n_unique_baseline is True:
        print "unique baseline" , repr(n_unique)
    else:
        pass

    
    q_unique=qu.real
    u_unique=qu.imag
    # this is map of visibilities which group them  according to their redundant unique group set"
    vis_map=0*q
    for ind in range(q.size):
          #print ind, q_unique, q[ind]
          unique_ind=np.where( (q[ind]==q_unique) & (u[ind]==u_unique))
          #print myind, ind
          vis_map[[ind]]=unique_ind

    

    return [ant1,ant2,vis_map,q,u,q_unique,q_unique,n_unique]




# simulated visibilities form random sky
def sim_vis(n_ants,q,ant1,ant2,q_unique,vis_map,g_low = 0.8,g_high = 1.2, sigma_sky= 0.5, noise_fac = 0.01):

    """ This function simulate visibility using true sky as Guassian random field with sigma =0.5, gains [0.8,1.2] and with add a random noise of noise_fac = 0.01"""
    """ This function returns an array of real  simulated visibility, diagonal covaraince matrix,simulated true sky visibity and gains """
   
    sky= np.random.normal(0.0,sigma_sky)*np.random.randn(q_unique.size)
    g_param= np.random.uniform(g_low,g_high)*np.random.randn(n_ants)
    g_param =np.array(g_param)
    sky =np.array(sky)
    ant1 = np.array(ant1)
    ant2 = np.array(ant2)
    vis_map =np.array(vis_map)
    Vis =[g_param[ant1]*g_param[ant2]*sky[vis_map]+ noise_fac*np.random.randn(q.size), noise_fac*np.ones(q.size)]
    cov_matrix = np.zeros([q.size,q.size])
    for k in range(q.size):
        cov_matrix[k][k] =Vis[1][k]

    return [Vis[0],cov_matrix,sky,g_param]
   

