import numpy as np, math 
from numpy.linalg import inv
from matplotlib import pyplot as plt
from scipy import linalg as lin
import classical_omnical as omnical

xdim =3
ydim =3
n_ants =xdim*ydim

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

        


qu=q+np.complex(0,1)*u
qu=np.unique(qu) # selecting unique baselines

n_unique=(xdim)**2-1+(xdim-1)**2
   


    
q_unique=qu.real
u_unique=qu.imag
# this is map of visibilities which group them  according to their redundant unique group set"
vis_map=0*q
for ind in range(q.size):
          #print ind, q_unique, q[ind]
          unique_ind=np.where( (q[ind]==q_unique) & (u[ind]==u_unique))
          #print myind, ind
          vis_map[[ind]]=unique_ind


sigma_sky = 0.5
Tsys = 50 #K
CHIME_B = 10**6 #  1MHz
tau = 21 # integration time in seconds

def noise_sky(Tsys,tau,bandwidth):
    return Tsys/(math.sqrt(tau*bandwidth))


noise_fac = noise_sky(Tsys,CHIME_B,tau)
sky= np.random.normal(0.0,sigma_sky)*np.random.randn(q_unique.size)
amp= 0.3
s_n = ((0.3**2)*0.5)/noise_fac

print 'S/N', s_n


g_param= amp*np.random.randn(n_ants)


g_param =np.array(g_param)/np.mean(g_param)
sky =np.array(sky)
ant1 = np.array(ant1)
ant2 = np.array(ant2)
vis_map =np.array(vis_map)
Vis =[g_param[ant1]*g_param[ant2]*sky[vis_map]+ noise_fac*np.random.randn(q.size), noise_fac*np.ones(q.size)]

B= omnical.config_matrix(n_ants,ant1,ant2,vis_map,g_param,sky,q,n_unique)
def  model_vis(g,sky):
    return g[ant1]*g[ant2]*sky[vis_map]

s_0 =  sky + 0.1*np.random.randn(sky.size)
g_0 = g_param + 0.1*np.random.randn(g_param.size)

def linncal_func(Vis,g_0,s_0,N_steps):
    Parameter ={}
    for j in range(N_steps):
        P = omnical.lincal(Vis[0],model_vis(g_0,s_0),g_0,s_0,omnical.config_matrix(n_ants,ant1,ant2,vis_map,g_0,s_0,q,n_unique))
        Parameter[j] =P
        s_0 = P[1]
        g_0 = P[0]

    return Parameter




plt.plot(g_param,linncal_func(Vis,g_0,s_0,3)[2][0],'.',label=' 1st iteration')
plt.xlabel('input gains')
plt.ylabel('output gains')
plt.legend(loc = 'upper left')
plt.show()



plt.plot(sky,linncal_func(Vis,g_0,s_0,3)[2][1],'.',label=' 1st iteration')
plt.xlabel('input true sky ')
plt.ylabel('output true sky')

plt.legend(loc = 'upper left')
plt.show()


"""
"""
plt.plot(Vis[0],model_vis(linncal_func(Vis,g_0,s_0,3)[2][0],linncal_func(Vis,g_0,s_0,2)[1][1]),'.',label=' 1 iteration')
plt.xlabel('simulated vis ')
plt.ylabel('output vis')

plt.legend(loc = 'upper left')
plt.show()

def noise_dB(output,input):
    results = np.zeros(input.size)
    for j in range(input.size):
        results[j] = 10*math.log(abs(output[j]/input[j]))
    return results

        

plt.plot(model_vis(linncal_func(Vis,g_0,s_0,3)[2][0],linncal_func(Vis,g_0,s_0,2)[1][1])/Vis[0],'.',label=' 1 iteration')

plt.ylabel(r'$vis_{output}/vis_{input}$')

plt.legend(loc = 'upper left')
plt.show()


plt.plot(noise_dB(model_vis(linncal_func(Vis,g_0,s_0,3)[2][0],linncal_func(Vis,g_0,s_0,2)[1][1]),Vis[0]),'.',label=' 1 iteration')



plt.ylabel(r'$dB$')

plt.legend(loc = 'upper left')
plt.show()

def  sigma_ration( output,input):
    return  np.std(output -input)/np.std(input)



