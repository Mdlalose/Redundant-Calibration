# Chi squared function, Gradient & Curvature
import numpy as np
from matplotlib import pyplot as plt


# chi squared function

def get_chisqd(data,gains,sky,ant1,ant2,vis_map):
    chi=0.0
    for vis in range(data.size):
            error= data[vis] - np.conj(gains[ant1[vis]])*gains[ant2[vis]]*sky[vis_map[vis]] 
            chi += (np.conj(error)*error)
    
    return chi


def get_grad_chisqd(data,gains,sky,ant1,ant2,vis_map):
     
     grad= np.zeros(gains.size + sky.size,dtype="complex")
     for vis in range(data.size):
            error= data[vis] - np.conj(gains[ant1[vis]])*gains[ant2[vis]]*sky[vis_map[vis]]
            # grad chi2 wrt gains
            # ant1
            grad[ant1[vis]] += -np.conj(gains[ant2[vis]])*np.conj(sky[vis_map[vis]])*error-np.conj(error)*gains[ant2[vis]]*sky[vis_map[vis]]
            
            #ant2
            grad[ant2[vis]] += -gains[ant1[vis]]*np.conj(sky[vis_map[vis]])*error-np.conj(error)*np.conj(gains[ant1[vis]])*sky[vis_map[vis]]
            # grad chi2 wrt sky
            grad[gains.size + vis_map[vis]] += -gains[ant1[vis]]*np.conj(gains[ant2[vis]])*error-np.conj(error)*np.conj(gains[ant1[vis]])*gains[ant2[vis]]
           
     
     return grad

def get_curv_chisqd(data,gains,sky,ant1,ant2,vis_map):
     curv = np.zeros((gains.size + sky.size,gains.size + sky.size),dtype="complex")
     for vis in range(data.size):
         error= data[vis] - np.conj(gains[ant1[vis]])*gains[ant2[vis]]*sky[vis_map[vis]]
         curv[ant1[vis]][ant1[vis]]+= 2.0*np.conj(gains[ant2[vis]])*np.conj(sky[vis_map[vis]])*sky[vis_map[vis]]*gains[ant2[vis]]
         curv[ant2[vis]][ant1[vis]] += -np.conj(sky[vis_map[vis]])*error + np.conj(gains[ant2[vis]])*np.conj(sky[vis_map[vis]])*np.conj(gains[ant1[vis]])*sky[vis_map[vis]] + gains[ant1[vis]]*np.conj(sky[vis_map[vis]])*gains[ant2[vis]]*sky[vis_map[vis]]-np.conj(error)*sky[vis_map[vis]]
         curv[ant1[vis]][ant2[vis]] += -np.conj(sky[vis_map[vis]])*error + np.conj(gains[ant2[vis]])*np.conj(sky[vis_map[vis]])*np.conj(gains[ant1[vis]])*sky[vis_map[vis]] + gains[ant1[vis]]*np.conj(sky[vis_map[vis]])*gains[ant2[vis]]*sky[vis_map[vis]] - np.conj(error)*sky[vis_map[vis]]
         curv[gains.size + vis_map[vis]][ant1[vis]] += -np.conj(gains[ant2[vis]])*error + np.conj(gains[ant2[vis]])*np.conj(sky[vis_map[vis]])*np.conj(gains[ant1[vis]])*gains[ant2[vis]] + gains[ant1[vis]]*sky[vis_map[vis]]*np.conj(gains[ant2[vis]])*sky[vis_map[vis]] - np.conj(error)*gains[ant2[vis]]
         curv[ant1[vis]][gains.size + vis_map[vis]] += -np.conj(gains[ant2[vis]])*error + np.conj(gains[ant2[vis]])*np.conj(sky[vis_map[vis]])*np.conj(gains[ant1[vis]])*gains[ant2[vis]] + gains[ant1[vis]]*sky[vis_map[vis]]*np.conj(gains[ant2[vis]])*sky[vis_map[vis]] - np.conj(error)*gains[ant2[vis]]
         curv[ant2[vis]][ant2[vis]] += 2.0*gains[ant1[vis]]*np.conj(sky[vis_map[vis]])*np.conj(gains[ant1[vis]])*sky[vis_map[vis]]
         curv[gains.size + vis_map[vis]][ant2[vis]] += -gains[ant1[vis]]*error + gains[ant1[vis]]*np.conj(sky[vis_map[vis]])*np.conj(gains[ant1[vis]])*gains[ant2[vis]]  -error*np.conj(gains[ant1[vis]])+  gains[ant1[vis]]*np.conj(gains[ant2[vis]])*np.conj(gains[ant2[vis]])*sky[vis_map[vis]]
         curv[ant2[vis]][gains.size + vis_map[vis]] += -gains[ant1[vis]]*error + gains[ant1[vis]]*np.conj(sky[vis_map[vis]])*np.conj(gains[ant1[vis]])*gains[ant2[vis]]  -error*np.conj(gains[ant1[vis]]) +  gains[ant1[vis]]*np.conj(gains[ant2[vis]])*np.conj(gains[ant2[vis]])*sky[vis_map[vis]]
         curv[gains.size + vis_map[vis]][gains.size + vis_map[vis]] += 2.0*gains[ant1[vis]]*np.conj(gains[ant2[vis]])*np.conj(gains[ant1[vis]])*gains[ant2[vis]]
         
     return curv     

     
           
            
            
            




