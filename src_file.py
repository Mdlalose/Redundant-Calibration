from astropy.io import fits
from astropy import units as u
import numpy as np
import math

"""
xdim =4
ydim = 4
#zdim = 3    
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

    

print "calculating u v w"
#get uv plane
space_d =10.0
h0=22
dec=-30
#phasecenter
mm = np.array([[np.sin(h0), np.cos(h0), 0],
      [-math.sin(dec)*np.cos(h0), math.sin(dec)*np.sin(h0), math.cos(dec)],
      [math.cos(dec)*np.cos(h0), -math.cos(dec)*np.sin(h0), math.sin(dec)]]) 
w = np.zeros(len(q_unique)) + abs(0.0000001*np.random.randn(len(q_unique)))
blz = np.zeros(len(u))
u_bl = np.outer(mm[0,0],(u_unique*space_d)/0.5) + np.outer(mm[0,1],(q_unique*space_d)/0.5) + np.outer(mm[0,2],w/0.5)
v_bl = np.outer(mm[1,0],(u_unique*space_d)/0.5) + np.outer(mm[1,1],(q_unique*space_d)/0.5) + np.outer(mm[1,2],w/0.5)

w_bl = np.outer(mm[2,0],(u_unique*space_d)/0.5) + np.outer(mm[2,1],(q_unique*space_d)/0.5) + np.outer(mm[2,2],w/0.5)        
 



u_bl, v_bl, w_bl= [ x.flatten() for x in (u_bl, v_bl, w_bl) ]
#w = np.zeros(len(v_bl))

# unique baselines
uniq_blz = np.zeros(len(u_unique))        
u_bl_uniq = np.outer(mm[0,0],u_unique) + np.outer(mm[0,1],q_unique) + np.outer(mm[0,2],uniq_blz)
v_bl_uniq = np.outer(mm[1,0],u_unique) + np.outer(mm[1,1],q_unique) + np.outer(mm[1,2],uniq_blz)

"""

hdulist = fits.open('GLEAM_EGC.fits')
scidata = hdulist[1].data

RA_src =[]
Dec_src =[]
src_flux =[] #jy
src_flux_unc =[]
for i in scidata:
	        #if (10*u.degree)<=i[12]<=(30*u.degree):
		RA_src.append(i[5])
		Dec_src.append(i[7])
	        src_flux.append(i[11])
		src_flux_unc.append(i[12])
	

np.save('src_pos.npy',np.array([RA_src,Dec_src]))
np.save('src_flux.npy',np.array([src_flux,src_flux_unc]))


"""
#select scr from +30,+10

#phasecenter 22,30
src_0 = [22,-30]

def get_lmn(s_0,src_pos):
          M_s0 = np.zeros((3,3))
          M_s0[0][0], M_s0[0][1] = np.sin(np.radians(s_0[0])), np.cos(np.radians(s_0[1]))
	  M_s0[1][0], M_s0[1][1], M_s0[1][2] = -np.sin(np.radians(s_0[1]))*np.cos(np.radians(s_0[0])), np.sin(np.radians(s_0[1]))*np.sin(np.radians(s_0[0])), np.cos(np.radians(s_0[1]))
	  M_s0[2][0], M_s0[2][1], M_s0[2][2] = np.cos(np.radians(s_0[1]))*np.cos(np.radians(s_0[0])),-np.cos(np.radians(s_0[1]))*np.sin(np.radians(s_0[0])), np.sin(np.radians(s_0[1]))
          l=[]
          m =[]
	  for pos in src_pos:
                  tmp = np.array([np.cos(np.radians(pos[1]))*np.sin(np.radians(pos[0])),-np.cos(np.radians(pos[1]))*np.sin(np.radians(pos[0])), np.sin(np.radians(pos[1]))])
                  l.append(M_s0.dot(tmp)[0])
		  m.append(M_s0.dot(tmp)[1])
		
             
	  
          return np.array([l,m])

src_pos = np.load('src_pos.npy')
src_flux = np.load('src_flux.npy')
print "calculating src l,m"

l_src, m_src =get_lmn(src_0,src_pos)

#l_max, l_min = np.max(l_src), np.min(l_src)
#m_max, m_min = np.max(m_src), np.min(m_src)

#l = np.arange(l_min,l_max)
#m = np.arange(m_min, m_max)

def I_src(l,m,l_src,m_src,src_flux):
        index_l = np.arange
	I = np.zeros(l.size,m.size)
	for i in range(len(l_src)):
             
		l_i = 10*(l[-1] + l_src[i])
		m_i = 10*(m[-1] + m_src[i])
		I[l_i][m_i] = src_flux[i]



def ant_pos(x,y,h,d):
        cosdcosh = np.cos(np.radians(d))*np.cos(np.radians(h))
	x_pos = x*cosdcosh
        y_pos = y*(-np.cos(np.radians(d))*np.sin(np.radians(h)))
        return [x_pox,y_pos]





#l_src, m_src, src_flux, src_pos = np.array(l_src)[0:10],np.array(m_src)[0:10],np.array(src_flux)[0:10], np.array(src_pos)[0:10]



def sigma(l_wave,dia_m):
    return (0.5*(1.22*l_wave))/dia_m
def get_V_uvw(u_bl,v_bl,w,l_src,m_src,src_flux,src_pos,l_wave,dia_m):
	V_uvw = []
        gain_ =[]
        n=0	
	for (u,v,w_) in zip(u_bl, v_bl, w):
	       tmp = 0.0
               for (l,m,flux,src_p) in zip(l_src,m_src,src_flux,src_pos):

			phase_shift = 2.0*np.pi*np.inner([u,v,w_],[l,m,np.sqrt(1-(l**2+ m**2))])
			exp_phase_shift = np.cos(phase_shift) -1j*np.sin(phase_shift)
               		#print exp_phase_shift
			sigma_ = sigma(l_wave,dia_m)
			if np.cos(np.radians(src_p[0])) >= 0.0 or np.cos(np.radians(src_p[1])) >= 0.0:  
				Beam = np.exp(-0.5*(np.cos(np.radians(src_p[0]))**2 + np.cos(np.radians(src_p[1]))**2)/sigma_)
				#print  Beam*flux[0]*exp_phase_shift
				tmp += Beam*flux[0]*exp_phase_shift
			else:
			   	pass
	       gain_i=  np.random.uniform(0.1,1.0)*(np.cos(np.random.uniform(0.0,2.0*np.pi)) +1j*np.sin(np.random.uniform(0.0,2.0*np.pi)))
	       gain_j = np.random.uniform(0.1,1.0)*(np.cos(np.random.uniform(0.0,2.0*np.pi)) +1j*np.sin(np.random.uniform(0.0,2.0*np.pi)))
	       obs_vis = tmp
	       n +=1
               print n
	       V_uvw.append(obs_vis)



	return np.array(V_uvw)

print "calculating observed visibility "
              
np.save('obs_vis_4x4.npy',get_V_uvw(u_bl,v_bl,w_bl,l_src,m_src,src_flux,src_pos,0.5,6.0))

print " successfully calculated visibilities"
"""

