#simulated Visibilitity
import numpy as np, math
from numpy.linalg import norm
from matplotlib import pyplot as plt
import random as rand
from scipy.optimize import fmin_cg
import lincal as lin



xdim =4
ydim = 4
    
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

#print   repr(n_unique)
    
    
q_unique=qu.real
u_unique=qu.imag
 # this is map of visibilities which group them  according to their redundant unique group set"
vis_map=0*q
for ind in range(q.size):
          #print ind, q_unique, q[ind]
          unique_ind=np.where( (q[ind]==q_unique) & (u[ind]==u_unique))
          
          #print myind, ind
          vis_map[[ind]]=unique_ind

    
    

def R_vec(uniq_block):
    R_={}
    for (k1,vis_real_imag) in uniq_block.iteritems():

         vis_real_imag[1] = np.array(vis_real_imag[1])
         # R vector for real visibilities    
         R_r=np.zeros(vis_real_imag[1].size)
         R_r[::2]=1
         vis_uniq_r= np.mean(vis_real_imag[1][::2])
        
         # R vector for imaginary visibilities
         R_i=np.zeros(vis_real_imag[1].size)
         R_i[1::2]=1
         vis_uniq_imag= np.mean(vis_real_imag[1][1::2])           
         R_[k1]=[vis_uniq_r,R_r,vis_uniq_imag,R_i]
         
    return R_  



print " done calculating u,q, vis_map"
   


def src(n_src,l,m):
    """ This function generate n_src  dictionary"""
    src_vec ={}
    for n_i in range(n_src):
        #amp = np.random.normal(0.0,0.5)
        amp = 0.5
        
        #phase= np.random.uniform(-math.pi,math.pi)
        phase = math.pi/10
        #src_vec[n_i] = np.array([amp,phase,rand.choice(l),rand.choice(m)])
        src_vec[n_i] = np.array([amp,phase,l,m])
        
    
    return src_vec

 
def V_src(n_src,l,m):
    """ This function generate n_src  dictionary"""
    src_vec ={}
    for n_i in range(n_src):
        amp = np.random.normal(0.0,0.5)
        #amp = 1.0
        
        phase= np.random.uniform(-math.pi,math.pi)
        #phase = math.pi/10
        src_vec[n_i] = np.array([amp,phase,rand.choice(l),rand.choice(m)])
        #src_vec[n_i] = np.array([amp,phase,l,m])
        
    
    return src_vec    


def I_src(l,m,src_dict):
    " This function generate intensity map given src vector"
    I=np.zeros((l.size,m.size),dtype ="complex")
    for key in src_dict:
        i = 956.0*(l[-1] +src_dict[key][2])
        j = 956.0*(l[-1] + src_dict[key][3])
        #print src_dict[key][0],i,j
        
        
        I[i][j] =src_dict[key][0]*complex(np.cos(src_dict[key][1]),math.sin(src_dict[key][1]))
        
    return I
def I_src_ref(l,m,src_dict):
    " This function generate intensity map given src vector"
    I=np.zeros((l.size,m.size),dtype ="complex")
    for key in src_dict:
        i = 956.0*(l[-1] +src_dict[key][2])
        j = 956.0*(l[-1] + src_dict[key][3])
        #print src_dict[key][0],i,j
        
        
        I[i][j] =src_dict[key][0]**2
    return I

def model_visib(I,A,l,m,u,q,n_ants):
      
        V_uq =[]  #np.zeros(u.size,dtype ="complex")
        for u_i,q_i,n_ant in zip(u,q,n_ants):
            exp_phase = np.cos(2.0*np.math.pi*(u_i*l+q_i*m)) -1j* np.sin(2.0*np.math.pi*(u_i*l+q_i*m)) 
            tmp= I.dot(A[n_ant])
            tmp = tmp*exp_phase
            tmp = tmp.sum(axis=1,dtype="complex")
            V_uq.append(sum(tmp))
            
        return V_uq    
    

        
        
    


    
    
def sigma(l_wave,dia_m):
    return (0.5*(1.22*l_wave))/dia_m

def gaussian_beam(sigma,l,m):
    " This is a 2d Gaussian circular aperture field pattern"

    norm = 1.0/(2.0*math.pi*sigma**2)
    sigma_ =sigma**2
    Beam_p =np.zeros((l.size,m.size))
    for i in range(l.size):
        for j in range(m.size):
        
           Beam_p[i][j] =np.exp(-0.5*(l[i]**2+m[j]**2)/sigma_)
            
     
    return Beam_p

def antennas_beam(n_ants,sigma_vec,l,m):
        beams ={}
        for i in range(n_ants):
            beams[i]= gaussian_beam(sigma_vec[i],l,m)
            
        return beams
    
    
def get_rows(ydim):
    row =[]
    row_i=0
    row_vec = np.zeros(ydim-1,'int')
    while row_i< ydim-1:
        row.append(row_vec + row_i)
        row_i +=1
        
    return np.array(row)
        
    
def get_cols(xdim):
    col =[]
    col_i=1
    while col_i<xdim:
        col.append(np.arange(1,xdim))
        col_i +=1
    return np.array(col)
        
    
         
        
            
         

def noise(T_sys,t,B):
    return T_sys/np.sqrt(t*B)

def get_ants_pos_dev(x_vec,y_vec,df):
    x_=np.zeros(x_vec.size)
    y_=np.zeros(x_vec.size)
    
    for i in  range(x_vec.size):
        if i !=0:
            x_[i]= x_vec[i][0] +df*np.random.randn()
            y_[i]= y_vec[i][0] + df*np.random.randn()
        else:
            pass
    return [x_,y_]


        
def get_sim_vis(I_lm,A_lm,l,m,x_ants,y_ants,n_ants,ant1,ant2,xdim,ydim):
    # this function simulate visibity power output in   given sky  with E =A*exp(itheta) in physical coordinate l,m
    Voltage =[]#np.zeros(n_ants, dtype ="complex")
    Voltage_true =[]
    Gains =[]
    Init_phi=[]
    Init_eta=[]
    noise_frac = noise(50,600,np.power(400,6))
    for x_i, y_i, n_ant in zip(x_ants,y_ants,n_ants):
          
          #phese delay
          delta_phi = x_i*np.sin(l)+ y_i*np.sin(m)
          #delta_phi_ant = x_i*np.sin(l[l.size/2 +3])+ y_i*np.sin(m[m.size/2 +3])
          
            
          phase_shift = np.cos(delta_phi) + 1j*np.sin(delta_phi)
          #phase_shift_ant = np.cos(delta_phi_ant) + 1j*np.sin(delta_phi_ant)
          #print delta_phi_ant
          #print phase_shift_ant
          
          if  norm(delta_phi) == 0.0:
          #if x_i == x_ants[0] & y_i ==y_ants[0]:
            #reference ant
            #print n_ant,A_lm[n_ant]
            
            V_ant=I_lm.dot(A_lm[n_ant])
    
            tmp =V_ant.sum(axis=1,dtype ="complex")
        
            eta= np.random.normal(0.0,1.0)
            amp = np.exp(eta)
            phase = np.random.uniform(0.0,2.0*math.pi)
            gain= amp*(np.cos(phase) +1j*np.sin(phase))
            Init_phi.append(phase)
            Init_eta.append(eta) 
            Gains.append(gain)
            Voltage.append(sum(tmp)*gain)
            Voltage_true.append(sum(tmp))
           
            
            
                
                
          else:
            I_lm_shifted =I_lm*phase_shift
            V_ant= I_lm_shifted.dot(A_lm[n_ant])
            tmp = V_ant.sum(axis=1,dtype ="complex")    
            #Voltage[n_ant] = sum(tmp)
            eta= np.random.normal(0.0,1.0)
            amp = np.exp(eta)
            phase = np.random.uniform(0.0,2.0*math.pi)
            gain= amp*(np.cos(phase) +1j*np.sin(phase))
            Gains.append(gain)
            Init_phi.append(phase)
            Init_eta.append(eta)
            Voltage.append(sum(tmp)*gain) 
            Voltage_true.append(sum(tmp))
            
            
    # cross correlation voltages
    Voltage= np.array(Voltage)
    Voltage_true= np.array(Voltage_true)
    Gains =np.array(Gains)
    #power_output = np.conj(Gains[ant1]*Voltage[ant1])*Gains[ant2]*Voltage[ant2]
    power_output =  np.conj(Voltage_true[ant1])*Voltage_true[ant2] 
    #diag =noise_frac*np.random.randn(power_output.size)
    
                
    
    
              
    
    return  [power_output,Gains] #[np.conj(Voltage[0])*Voltage[1]]
    



print " calculating visibilities"
m,l= np.arange(-math.pi/10.0,math.pi/10.0,0.001), np.arange(-math.pi/10.0,math.pi/10.0,0.001)
n_ants =np.arange(xdim*ydim)
n_ants_uniq =np.arange(q_unique.size)
sigma_vec = np.array([sigma(0.5,6.0),sigma(0.5,6.0),sigma(0.5,6.0),sigma(0.5,6.0),sigma(0.5,6.0),sigma(0.5,6.0),sigma(0.5,6.0), sigma(0.5,6.0),sigma(0.5,6.0),sigma(0.5,6.0),sigma(0.5,6.0),sigma(0.5,6.0),sigma(0.6,6.0),sigma(0.5,6.0),sigma(0.5,6.0), sigma(0.5,6.0),sigma(0.5,6.0),sigma(0.5,6.0),sigma(0.5,6.0),sigma(0.5,6.0),sigma(0.5,6.0),sigma(0.5,6.0),sigma(0.5,6.0), sigma(0.5,6.0),sigma(0.5,6.0),sigma(0.5,6.0),sigma(0.5,6.0),sigma(0.5,6.0),sigma(0.5,6.0),sigma(0.5,6.0),sigma(0.5,6.0), sigma(0.5,6.0)])
Vis_data = np.array([get_sim_vis(I_src(l,m,V_src(1,l,m)),antennas_beam(n_ants.size,sigma_vec,l,m),l,m,get_ants_pos_dev(x,y,0.0)[0]*200, get_ants_pos_dev(x,y,0.0)[1]*200,n_ants,ant1,ant2,xdim,ydim)])[0]

#.01, .02, .05, .1,
np.save('data_vis_test_data.npy',Vis_data)
def unique_blocks(data,u,q,u_nique,q_unique):
    "This fuction sort visibilities in redundablocks, each vector in vecs_unique contains all visibililies with same baseline "

    
    vecs_unique =[]
    for u_uniq, q_uniq in zip(u_nique,q_unique):
        tmp = []
        
        
        for i in range(u.size):
            if u[i] ==u_uniq and q[i] == q_uniq:
                tmp.append(data[i])
                #tmp.append(data[i].real) # real visibility
                #tmp.append(data[i].imag) # imaginary visibility
                
        vecs_unique.append(np.array(tmp))

    return np.array(vecs_unique) 

uniq_vis= np.zeros(unique_blocks(Vis_data[0],u,q,u_unique,q_unique).size,dtype="complex")

for i in range(unique_blocks(Vis_data[0],u,q,u_unique,q_unique).size):
    uniq_vis[i]=unique_blocks(Vis_data[0],u,q,u_unique,q_unique)[i][0]
    #print unique_blocks(Vis_data[1],u,q,u_unique,q_unique)[i][0]
    
    

#Vis_uniq =np.array(model_visib(I_src_ref(l,m,V_src(1,l,m)),antennas_beam(n_ants_uniq.size,sigma_vec,l,m),l,m,absu_unique,q_unique,n_ants_uniq))
np.save('data_vis_uniq_test.npy',uniq_vis)




#vis_u, ind = np.unique(Vis_data[1],return_index ="True")
#vis_u= vis_u[np.argsort(ind)]
#vis_u = vis_u[0:n_unique]
#np.save('data_vis_uniq_test.npy',vis_u)
#vis_u_map = np.zeros(Vis_data[1].size,dtype="int")
#for i in range(Vis_data[1].size):
    #index = np.where(Vis_data[1][i] == vis_u)

    #vis_u_map[i]=index[0]


#print "done calculating visibilities"  

ant = np.array([ant1,ant2,vis_map])
np.save('ants_data_test_.npy',ant)
    
def  model_vis(g,sky,vis_map):   
    return np.conj(g[ant1])*g[ant2]*sky[vis_map] 
#g_0 = Vis_data[2]/np.mean(Vis_data)
#s_0 = vis_u

data_gains = Vis_data[1]
true_sky_vis= uniq_vis

sim_data = np.conj(data_gains[ant1])*data_gains[ant2]*true_sky_vis[vis_map]

np.save('data_vis_test.npy',sim_data)

plt.plot(sim_data.real-model_vis(np.array(Vis_data[1]),uniq_vis,vis_map).real,'.')
plt.xlabel('Visibility Points')
plt.ylabel('Data -Model')
plt.show()



#fac=1.0;

#asdf=fmin_cg(lin.get_chi_phase,gain_phase*fac,lin.get_grad_phase,(v_data,v_model,ant1,ant2,fac))

#fit_gains=asdf/fac


#Beam = antennas_beam(n_ants_uniq.size,sigma_vec,l,m)


#plt.imshow(Beam[1]-Beam[0],cmap="Blues")
#plt.title(r"Beam Residual for  D=6.00m &  D= 6.01 m")
#plt.ylabel(r'$956.0(\pi/10 + m)$')
#plt.xlabel(r'$956.0(\pi/10 + l)$')
#plt.colorbar()
#plt.show()




        
        


cf =(200*0.5)/(2.0*math.pi),
plt.plot(x*200,y*200,'*',label="Perfect Redundant")
plt.plot( get_ants_pos_dev(x,y,0.3)[0]*200, get_ants_pos_dev(x,y,0.3)[1]*200,'.',label="Quasi-Redundant")
plt.ylabel('North-South Antenna Position (m)')
plt.xlabel('East-West Antenna Position (m)')

plt.legend(loc='best')
plt.grid(True)
plt.show()
"""
V =[]
x_ant,y_ant =np.array([0.0,0.0]),np.array([200.0,0.0])
n_ants= np.arange(2)
Beam = antennas_beam(n_ants.size,sigma_vec,l,m)
for j in range(m.size):      
      
      #print m_src[j],j
      V.append(get_sim_vis(I_src(l,m,src(1,l[0],m[j])),Beam,l,m,x_ant,y_ant,n_ants,ant1,ant2))
     
      #plt.ion()
      #plt.imshow(I_src(l,m,src(1,l[0],m[j])).real)
      #print  I_src(l,m,src(1,l_src,m_src)), l_src,m_src,i
      #plt.show()
    
    
 





plt.imshow(Beam[0],cmap="Greys")


  
    
m_src= m
plt.imshow(Beam[0],cmap="Blues")
plt.title(r"Gaussian Beam for $\sigma(\lambda =0.5 m, D= 6 m)$")
plt.ylabel(r'$956.0(\pi/10 + m)$')
plt.xlabel(r'$956.0(\pi/10 + l)$')
plt.colorbar()
plt.show()

plt.plot(m_src,np.array(V).real)
plt.ylabel('$Re(\mathbb{V})$ $[{Voltage}^2]$')
plt.xlabel('Radains')

plt.show()

plt.plot(m_src,np.array(V).imag)

plt.ylabel('$Im(\mathbb{V})$ $[{Voltage}^2]$')
plt.xlabel('Radians')
plt.show()

plt.plot(m_src,np.sqrt((np.array(V).real)**2 + (np.array(V).imag)**2))
plt.xlabel('Radians')
plt.ylabel('$\mid \mathbb{V} \mid$ [Voltage]')
plt.show()

plt.plot(m_src,np.arctan2(np.array(V).imag,np.array(V).real))
plt.ylabel(r'$\phi_{\mathbb{V}}$')
plt.xlabel('Radians')
plt.show()
"""

#if __name__ == __main__:
    #pass
