# omnical @ J DILSON
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
import anta_pos as ant
import sys

xdim = int(sys.argv[1])
ydim = int(sys.argv[2])
long_0= 21.4278 # deg
lalt_0=  -30.7224 # deg
ant_positions= ant.get_antenna_pos(xdim,ydim,long_0,lalt_0)
"""
#def get_ant_pos(ant_separ,precFrac,hexNum):
ant_separ= 14.0
precFac = 100000
hexNum = 4
pos=[]
for row in range(hexNum-1,-(hexNum),-1):
	for col in range(2*hexNum-abs(row)-1):
		#print col, 'col', row, 'row'
		xPos = ((-(2*hexNum-abs(row))+2)/2.0 + col)*ant_separ;
        	yPos = row*ant_separ*3**.5/2;
        	pos.append([xPos, yPos, 0])
positions = np.array(pos)


#Calculating Baselines    
nAntennas = len(positions)

blVectors = []
blPairs = []
for ant1 in range(nAntennas):
    for ant2 in range(ant1+1,nAntennas):
        delta = np.array([int(np.round(precFac*(positions[ant1][i]-positions[ant2][i]))) for i in range(3)])
	print delta[0], delta[1]
        if delta[1] > 0 or (delta[1] == 0 and delta[0] > 0): 
            blVectors.append(tuple(delta))
            blPairs.append((ant1, ant2))
        else: 
            blVectors.append(tuple(-delta))
            blPairs.append((ant2, ant1))
        
#Calculating Unique Baselines
ublDict = {}
for b in range(len(blVectors)):
    if ublDict.has_key(blVectors[b]): ublDict[blVectors[b]].append(blPairs[b])
    else: ublDict[blVectors[b]] = [blPairs[b]]
ublIndexDict = {antPair: i for i,antPairs in enumerate(ublDict.values()) for antPair in antPairs }
ublVectors = np.array([positions[antList[0][0]]-positions[antList[0][1]] for antList in ublDict.values()])
print "With", len(positions), "antennas there are", len(ublDict.items()), "unique baselines."

nant, nbl, nubl = len(positions), len(blPairs), len(ublVectors)

#Simulating visibilities

xdim = 3 #int(sys.argv[1])
ydim =  3 #int(sys.argv[2])
    
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

"""
gains = np.zeros(ydim**2,dtype='complex')
eta_0, phi_0 = np.zeros(ydim**2), np.zeros(ydim**2)
for g in range(ydim**2):
     eta= np.random.normal(0.0,1.0)
     eta_0[g]=eta
     
     amp = np.random.uniform(0.1,1.0)
     phase = np.random.uniform(0.0,2.0*math.pi)
     phi_0[g] = phase
     gains[g]=amp*(np.cos(phase)+ 1j*np.sin(phase))









noise_frac_gains = float(sys.argv[3])
noise_frac_sky = float(sys.argv[4])
data_vis= np.load('obstrue_vis_8x8.npy')
ant1,ant2, vis_map = data_vis[1], data_vis[2], data_vis[3]
data_vis_0_05= np.load('obstrue_vis_8x820_data_redundant.npy')
trueVis = data_vis[0][0]
gains = gains + noise_frac_gains*np.random.randn(len(gains))
trueVis = trueVis + noise_frac_sky*np.random.randn(len(trueVis))
obsVis = gains[ant1] * np.conj(gains[ant2]) * trueVis[vis_map]

noise_frac_gains = float(sys.argv[3])
noise_frac_sky = float(sys.argv[4])
#data = np.array(np.conj(gains[ant1])*gains[ant2]*data_vis_0_05[0][0])
nu = np.arange(170.0,230,1.0)
#obsVis = [ for (ant1,ant2) in blPairs]
#noiselessVis = np.array(obsVis)
#obsVis = np.array(obsVis) + noise

blVectors, blPairs = ant.get_blvectors(ant_positions)
	#np.save('redund_blVectors.npy',[blVectors,blPairs])
	
nant= xdim*ydim
nubl= len(trueVis)
ublDict = {}
for b in range(len(blVectors)):
     	# grouping togather all blPairs with same blVector
    	if ublDict.has_key(blVectors[b]): ublDict[blVectors[b]].append(blPairs[b])
    	else: ublDict[blVectors[b]] = [blPairs[b]]

ublIndexDict = {antPair: i for i,antPairs in enumerate(ublDict.values()) for antPair in antPairs }
ublVectors = np.array([ant_positions[antList[0][0]]-ant_positions[antList[0][1]] for antList in ublDict.values()])

def DegeneracyProjectionMatrices(AtA):
    #Find orthonormal vectors that span the null-space
    evals, evecs = np.linalg.eigh(AtA)
    zeroEVecs = (evecs[:,evals < 1e-12*np.max(evals)]).T
    
#     for i in range(len(zeroEVecs)):
#         for j in range(i+1,len(zeroEVecs)):
#             zeroEVecs[j,:] -= zeroEVecs[i,:].dot(zeroEVecs[j,:]) * zeroEVecs[i,:]
#             zeroEVecs[j,:] /= np.linalg.norm(zeroEVecs[j,:])
    
    Proj = np.eye(len(zeroEVecs[0,:])) - zeroEVecs.T.dot(zeroEVecs)
    IMinusProj = np.eye(len(zeroEVecs[0,:])) - Proj
    return Proj, IMinusProj
nbl = len(obsVis)
# logcal
Acoeffs, Bcoeffs, rowIndices, colIndices = [np.zeros(nbl*3) for i in range(4)]
for n,(ant1,ant2) in enumerate(blPairs):
    rowIndices[3*n:3*n+3] = n
    colIndices[3*n:3*n+3] = [ant1, ant2, nant + ublIndexDict[(ant1,ant2)]]
    Acoeffs[3*n:3*n+3] = [1.0, 1.0, 1.0]
    Bcoeffs[3*n:3*n+3] = [1.0, -1.0, 1.0]

logcalA = csr_matrix((Acoeffs,(rowIndices,colIndices)), shape=(nbl, nant + nubl))
logcalB = csr_matrix((Bcoeffs,(rowIndices,colIndices)), shape=(nbl, nant + nubl))
AtA = (logcalA.conj().T.dot(logcalA)).toarray()
BtB = (logcalB.conj().T.dot(logcalB)).toarray()



#degeneracies projection amtrixies
AProj, IMinusAProj = DegeneracyProjectionMatrices(AtA)
BProj, IMinusBProj = DegeneracyProjectionMatrices(BtB)

# delta x for real and imaginary part
xReal = (np.linalg.pinv(AtA)).dot(logcalA.conj().T.dot(np.real(np.log(obsVis))))
xImag = (np.linalg.pinv(BtB)).dot(logcalB.conj().T.dot(np.imag(np.log(obsVis))))

# removing the degeneracies
xRealTrue = np.real(np.log(np.append(gains, trueVis)))
xImagTrue = np.imag(np.log(np.append(gains, trueVis)))
xReal = AProj.dot(xReal) + IMinusAProj.dot(xRealTrue)
xImag = BProj.dot(xImag) + IMinusBProj.dot(xImagTrue)

# hat x
xHat = np.exp(xReal + 1.0j*xImag)
logcalGainSols, logcalVisSols = xHat[0:nant], xHat[nant:]

calVis = [obs / (logcalGainSols[ant1]*np.conj(logcalGainSols[ant2])) for (ant1,ant2),obs in zip(blPairs,obsVis)]

print 'Zero Eigenvalues: ', [np.sum(np.abs(np.linalg.eigvals(XtX)<1e-10)) for XtX in [AtA, BtB]]


# plots the logcal results
plt.figure()
plt.quiver(np.real(trueVis), np.imag(trueVis), np.real(logcalVisSols-trueVis), np.imag(logcalVisSols-trueVis),
           angles='xy', scale_units='xy', scale=1)
plt.plot(np.real(trueVis), np.imag(trueVis),'r.',label='True Visibilities')
plt.plot(np.real(logcalVisSols), np.imag(logcalVisSols),'b.',label='Log-Calibrated Visibilities')
plt.legend(loc='best')

#plt.figure()
#plt.plot(np.abs(trueVis), np.abs(logcalVisSols))
#plt.plot(np.angle(trueVis), np.angle(logcalVisSols),'.')
#plt.xlabel('Phase(True Visibility)'); plt.ylabel('Phase(Visibility Solution)')
#lincal

gainSols = logcalGainSols.copy()
visSols = logcalVisSols.copy()

allProj = []

for iteration in range(100):

    Acoeffs, rowIndices, colIndices = [np.zeros(nbl*12) for i in range(3)]
    for n,(ant1,ant2) in enumerate(blPairs):
        rowIndices[12*n:12*n+6] = 2*n
        rowIndices[12*n+6:12*n+12] = 2*n+1
        ublIndex = ublIndexDict[(ant1,ant2)]
        colIndices[12*n:12*n+6] = [2*ant1, 2*ant1+1, 2*ant2, 2*ant2+1, 2*nant+2*ublIndex, 2*nant+2*ublIndex+1]
        colIndices[12*n+6:12*n+12] = [2*ant1, 2*ant1+1, 2*ant2, 2*ant2+1, 2*nant+2*ublIndex, 2*nant+2*ublIndex+1]
        #Compute coefficients
        gi0V0 = gainSols[ant1] * visSols[ublIndex]
        gj0starV0 = np.conj(gainSols[ant2]) * visSols[ublIndex]
        gi0gj0star = gainSols[ant1] * np.conj(gainSols[ant2])
        Acoeffs[12*n:12*n+6] = [np.real(gj0starV0), -np.imag(gj0starV0), np.real(gi0V0), 
                                np.imag(gi0V0), np.real(gi0gj0star), -np.imag(gi0gj0star)]
        Acoeffs[12*n+6:12*n+12] = [np.imag(gj0starV0), np.real(gj0starV0), np.imag(gi0V0), 
                                   -np.real(gi0V0), np.imag(gi0gj0star), np.real(gi0gj0star)]

    A = csr_matrix((Acoeffs,(rowIndices,colIndices)), shape=(2*nbl, 2*nant + 2*nubl))
    AtA = (A.conj().T.dot(A)).toarray()

    modelObs = np.array([gainSols[ant1] * np.conj(gainSols[ant2]) * visSols[ublIndexDict[(ant1,ant2)]]  
                         for (ant1,ant2),obs in zip(blPairs,obsVis)])
    error = obsVis - modelObs
    y = np.dstack((np.real(error),np.imag(error))).flatten()
    xHat = np.linalg.pinv(AtA).dot(A.conj().T.dot(y))
    
    Proj, IMinusProj = DegeneracyProjectionMatrices(AtA)
    allProj.append(Proj) 
    xCurrent = np.dstack((np.real(np.append(gainSols, visSols)),np.imag(np.append(gainSols, visSols)))).flatten()
    xTrue = np.dstack((np.real(np.append(gains, trueVis)),np.imag(np.append(gains, trueVis)))).flatten()    

    xHat = Proj.dot(xHat) + IMinusProj.dot(xTrue - xCurrent) 
    
    updates = xHat[0::2] + 1.0j*xHat[1::2]
    alpha = 1
    gainSols += alpha*updates[0:nant]
    visSols += alpha*updates[nant:]

    modelObs = np.array([gainSols[ant1] * np.conj(gainSols[ant2]) * obs  for (ant1,ant2),obs in zip(blPairs,obsVis)])
    
    print 'Relative change in solutions:', np.linalg.norm(updates) / np.linalg.norm(np.append(gainSols,visSols))
    if np.linalg.norm(updates) / np.linalg.norm(np.append(gainSols,visSols)) < 1e-14: break



############################################################
# phase/amplitude
gainSols2 = logcalGainSols.copy()
#gainSols2 = gainSols2/np.mean(gainSols2)
visSols2 = logcalVisSols.copy()

allProj2 = []

for iteration in range(100):

    Acoeffs, rowIndices, colIndices = [np.zeros(nbl*12) for i in range(3)]
    for n,(ant1,ant2) in enumerate(blPairs):
        rowIndices[12*n:12*n+6] = 2*n
        rowIndices[12*n+6:12*n+12] = 2*n+1
        ublIndex = ublIndexDict[(ant1,ant2)]
        colIndices[12*n:12*n+6] = [2*ant1, 2*ant1+1, 2*ant2, 2*ant2+1, 2*nant+2*ublIndex, 2*nant+2*ublIndex+1]
        colIndices[12*n+6:12*n+12] = [2*ant1, 2*ant1+1, 2*ant2, 2*ant2+1, 2*nant+2*ublIndex, 2*nant+2*ublIndex+1]
        #Compute coefficients
        gi0gj0star = gainSols2[ant1] * np.conj(gainSols2[ant2])
        gi0gj0starVij0 = gi0gj0star * visSols[ublIndex]
        Acoeffs[12*n:12*n+6] = [np.real(gi0gj0starVij0), -np.imag(gi0gj0starVij0), np.real(gi0gj0starVij0), 
                                np.imag(gi0gj0starVij0), np.real(gi0gj0star), -np.imag(gi0gj0star)]
        Acoeffs[12*n+6:12*n+12] = [np.imag(gi0gj0starVij0), np.real(gi0gj0starVij0), np.imag(gi0gj0starVij0), 
                                   -np.real(gi0gj0starVij0), np.imag(gi0gj0star), np.real(gi0gj0star)]

    A = csr_matrix((Acoeffs,(rowIndices,colIndices)), shape=(2*nbl, 2*nant + 2*nubl))
    AtA = (A.conj().T.dot(A)).toarray()

    modelObs = np.array([gainSols2[ant1] * np.conj(gainSols2[ant2]) * visSols2[ublIndexDict[(ant1,ant2)]]  
                         for (ant1,ant2),obs in zip(blPairs,obsVis)])
    error = obsVis - modelObs
    y = np.dstack((np.real(error),np.imag(error))).flatten()
    xHat = np.linalg.pinv(AtA).dot(A.conj().T.dot(y))
    
    Proj, IMinusProj = DegeneracyProjectionMatrices(AtA)
    allProj2.append(Proj)
    deltaGains = gains/gainSols2
    deltaVis = trueVis - visSols2
    deltax = np.dstack((np.append(np.abs(np.log(deltaGains)),np.real(deltaVis)), 
                        np.append(np.angle(deltaGains), np.imag(deltaVis)))).flatten()
    xHat = Proj.dot(xHat) + IMinusProj.dot(deltax)                  
    updates = xHat[0::2] + 1.0j*xHat[1::2]
    
    gainSols2 *= (1.0+updates[0:nant])
    visSols2 += updates[nant:]
    modelObs = np.array([gainSols[ant1] * np.conj(gainSols[ant2]) * obs  for (ant1,ant2),obs in zip(blPairs,obsVis)])
    
    print 'Relative change in solutions:', np.linalg.norm(updates) / np.linalg.norm(np.append(gainSols2,visSols2))
    if np.linalg.norm(updates) / np.linalg.norm(np.append(gainSols2,visSols2)) < 1e-14: break


plt.figure()
plt.title('100% offset at 170MHz ')
plt.quiver(np.real(trueVis), np.imag(trueVis), np.real(visSols2-trueVis), np.imag(visSols2-trueVis),
           angles='xy', scale_units='xy', scale=1)
plt.plot(np.real(trueVis), np.imag(trueVis),'r.',label='True Visibilities')
#plt.plot(np.real(visSols), np.imag(visSols),'b.',label='Calibrated Visibilities')
plt.plot(np.real(visSols2), np.imag(visSols2),'g.',label='Calibrated Visibilities 2')
plt.legend(loc='best')
plt.grid(True)

gainSols2 = gainSols2/np.mean(gainSols2)
gains = gains/np.mean(gains)


plt.figure()
plt.title('100% offset at 170MHz ')
#plt.plot(np.abs(trueVis), np.abs(logcalVisSols))
plt.plot(np.abs(gains), np.abs(logcalGainSols),'.', label='Logcal')
plt.plot(np.abs(gains), np.abs(gainSols2),'.', label='Lincal')
plt.xlabel('Amplitude(True Gains)'); plt.ylabel('Amplitude(Gain Solution)')
#plt.plot([0, 1.1*np.max(np.abs(gains))], [0, 1.1*np.max(np.abs(gains))], '--k')
plt.legend(loc='best')
plt.grid(True)

plt.figure()
plt.title('100% offset at 170MHz ')
#plt.plot(np.abs(trueVis), np.abs(logcalVisSols))
plt.plot(np.angle(gains), np.angle(logcalGainSols),'.', label='Logcal')
plt.plot(np.angle(gains), np.angle(gainSols2),'.', label='Lincal')
plt.xlabel('Phase(True Gains)'); plt.ylabel('Phase(Gain Solution)')
#plt.plot([0, 1.1*np.max(np.abs(gains))], [0, 1.1*np.max(np.abs(gains))], '--k')
plt.legend(loc='best')
plt.grid(True)


