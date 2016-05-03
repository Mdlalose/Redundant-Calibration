#omnical using linear algebra
import numpy as np
from numpy.linalg import inv
from matplotlib import pyplot as plt
from scipy import linalg as lin
import simulated_vis as sim_vis
#contruct configation matrix A


xdim= sim_vis.xdim
ydim= sim_vis.ydim
vis_map=sim_vis.vis_map
ant1= sim_vis.ant1
ant2 = sim_vis.ant2

def B_matrix(Vis,gain_param,sky_param):
    xdim= sim_vis.xdim
    ydim= sim_vis.ydim
    vis_map=sim_vis.vis_map
    ant1= sim_vis.ant1
    ant2 = sim_vis.ant2
    A = np.zeros((sim_vis.q.size,xdim*ydim + sim_vis.n_unique))
    for vis_ in range(len(Vis)):
        #print 'Vis',Vis[0][vis_], 'ant1 j',ant1[vis_] ,'ant2 j',ant2[vis_], 'vis',vis_map[vis_]
    
    
        # contructing A
        A[vis_][ant1[vis_]]  = gain_param[ant2[vis_]]*sky_param[vis_map[vis_]]
        
        A[vis_][ant2[vis_]]  = gain_param[ant1[vis_]]*sky_param[vis_map[vis_]]
        A[vis_][xdim*ydim + vis_map[vis_]]     = gain_param[ant1[vis_]]*gain_param[ant2[vis_]]
        #print ant2[vis_],vis_map[vis_],gain_param[ant2[vis_]]*sky_param[vis_map[vis_]],ant1[vis_],vis_map[vis_], gain_param[ant1[vis_]]*sky_param[vis_map[vis_]],ant1[vis_],ant2[vis_],  gain_param[ant1[vis_]]*gain_param[ant2[vis_]],xdim*ydim+ vis_map[vis_]

        #print vis_,ant1[vis_]
    return A











def S_inv(s_value):
    S= np.zeros((len(s_value),len(s_value)))
    for ss in range(len(s_value)):
        if s_value[ss] <= 10**-5:
            #print s_value[ss]
            S[ss][ss] = 0.0

        else:
            s_inv =1.0/s_value[ss]
            #print s_inv
            S[ss][ss]= s_inv


    return S
            
    














def least_sqaured_param(data,model,conf_matr,gain_param,sky_param):
    delta_vis_model = data[0] - model
     
    B= conf_matr.transpose().dot(conf_matr)
    curv_inv =lin.pinv(B)
    #U,s,V =lin.svd(B)
    #curv_inv =   (V.dot(S_inv(s))).dot(U)
    A_T_data =  conf_matr.T.dot(delta_vis_model)
    delta_param =  curv_inv.dot(A_T_data)
    return delta_param


# find removinh the singularity in the curv matrix by  adding
#a Curv(Chi^2)= Curv(Chi^2_old) + 2lambda

def curv_pern(Cur_matr,l):
    U,s,V = lin.svd(Cur_matr)
    
    delta_lambda = np.zeros((len(s),len(s)))
    for s_value in range(len(s)):
           if s[s_value] <= 10**-10:
                 delta_lambda[ s_value][s_value]= l


           else:
                pass

    Lambda = lin.diagsvd(s,len(s),len(s))

    

    
    # curv pen
    new_curv =func.matrixmult(func.matrixmult (U,Lambda),V) + func.matrixmult(func.matrixmult(U,delta_lambda),V)

    return [delta_lambda,Lambda,new_curv]        
# SIMULATED DATA
Vis_data= sim_vis.Vis
gains_data = sim_vis.g_param
sky_data =sim_vis.sky

#intianial guess of gains and true sky vis and visibility model
g_low = 0.8
g_high = 1.2
sigma_sky= 0.5
sky_guess= np.random.normal(0.0,sigma_sky)*np.random.randn(sim_vis.q_unique.size)


g_param_guess= np.random.uniform(g_low,g_high)*np.random.randn(sim_vis.xdim*sim_vis.ydim)

#Vis_model =g_param_guess[sim_vis.ant1]*g_param_guess[sim_vis.ant2]*sky_guess[sim_vis.vis_map]

def chi2_fun( Data, model,A,deltax):
     B=  A.dot(deltax)
     C = model + B
     chi2 = 0.0
     for d in range(len(model)):
         c = (Data[0][d] - C[d])**2
         chi2 +=c
         

     
    
     return chi2

N_steps = 3
l= 0.0001
for k in range(N_steps):
   Vis_model =g_param_guess[sim_vis.ant1]*g_param_guess[sim_vis.ant2]*sky_guess[sim_vis.vis_map]
   Delta_X = least_sqaured_param(Vis_data,Vis_model,B_matrix(Vis_model,g_param_guess,sky_guess),g_param_guess,sky_guess)
   #Updating GUESSES
   
   g_param_guess =g_param_guess + Delta_X[0:sim_vis.xdim*sim_vis.ydim]
   sky_guess = sky_guess + Delta_X[sim_vis.xdim*sim_vis.ydim:len(Delta_X)]
   avg_g = np.mean(g_param_guess)
   #g_parma_guess = (g_param_guess)/avg_g
   
   
   print 'chi2' , chi2_fun(Vis_data,Vis_model,B_matrix(Vis_model,g_param_guess,sky_guess),Delta_X)
   
   print 'delta_x', Delta_X, g_param_guess

"""
Vis_model = g_param_guess[sim_vis.ant1]*g_param_guess[sim_vis.ant2]*sky_guess[sim_vis.vis_map]
delta_v =   B_matrix(Vis_model,g_param_guess,sky_guess).dot(least_sqaured_param(Vis_data,Vis_model,B_matrix(Vis_model,g_param_guess,sky_guess),g_param_guess,sky_guess))

Delta_V = Vis_data[0]- Vis_model
plt.ion()
plt.plot(delta_v,Delta_V,'.')
plt.xlabel(r' A$\Delta X$')
plt.ylabel(r'$ V_{ij} - g_ig_jy_{ij}$') 
                                                         
"""   
"""

CHI_PRED =[]
CHI_CALC =[]
CHI_PRED_2nd=[]
grad =[]
DELTA=[]
Grad_pre=[]
n=1000
H_matrix= 2.0*(B_matrix(Vis_data,gains_data,sky_data).T.dot( B_matrix(Vis_data,gains_data,sky_data)))
#grad_chi2_old = -2.0*(B_matrix(Vis,g_param,sky).T.dot(Vis[0]-B_matrix(Vis,g_param,sky).dot(x_old)))
for k in range(n):
    step_size= 0.5
    sky_new= sky_data #+ step_size*np.random.randn(n_unique)
    rand =np.zeros(xdim*ydim)
    rand[0]=step_size*np.random.randn(xdim*ydim)[0] # varying on variable
    g_new=  gains_data + rand
    x_old = np.concatenate((gains_data,sky_data)) 
    x_new= np.concatenate((g_new,sky_new))
    
    # chi2 and grad_chi2 at the centre
    chi2_old = (Vis_data[0]- B_matrix(Vis_data,gains_data,sky_data).dot(x_old)).T.dot(Vis_data[0]-B_matrix(Vis_data,gains_data,sky_data).dot(x_old))
    grad_chi2_old = -2.0*(B_matrix(Vis_data,gains_data,sky_data).T.dot(Vis_data[0]-B_matrix(Vis_data,gains_data,sky_data).dot(x_old)))

    # chi2_new calcualted at new point 
    
    chi2_calcu = (Vis_data[0]- B_matrix(Vis_data,g_new,sky_new).dot(x_new)).T.dot(Vis_data[0]-B_matrix(Vis_data,g_new,sky_new).dot(x_new))

    #chi2 predicted from chi2_old 
   
    chi2_pred= chi2_old + (x_new-x_old).T.dot(grad_chi2_old)
    chi2_pred_2nd =  chi2_old + (x_new-x_old).T.dot(grad_chi2_old) + 0.5*((x_new-x_old).T.dot(H_matrix)).dot(x_new-x_old) 
    CHI_PRED.append(chi2_pred)
    CHI_PRED_2nd.append(chi2_pred_2nd)
    CHI_CALC.append(chi2_calcu)
    

    ## grad vs delta
    Delta_X = x_new-x_old
    
    DELTA.append(Delta_X[0])
    
    
    


    


    #print chi2_old, chi2_calcu, chi2_pred, grad_chi2_old
   

    
1.003856618200047

plt.plot(DELTA,CHI_CALC,'.',label ='$\chi^2_{calc}$')
plt.plot(DELTA,CHI_PRED,'.',label ='$\chi^2_{1st}$')
plt.plot(DELTA,CHI_PRED_2nd,'.',label ='$\chi^2_{2nd}$')
plt.plot(DELTA,np.array(CHI_CALC)-np.array(CHI_PRED),'.',label ='$\chi^2_{resid1st}$')
plt.plot(DELTA,np.array(CHI_CALC)-np.array(CHI_PRED_2nd),'.',label ='$\chi^2_{resid_2nd}$')
plt.xlabel(r'$\delta$')
plt.ylabel(r'$\chi^2$')
plt.legend(loc ='best')
#plt.axis([min(DELTA), max(DELTA), -4*min(CHI_CALC), 4*max(CHI_CALC)])
plt.show()  
"""
