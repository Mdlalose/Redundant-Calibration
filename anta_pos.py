import numpy as np, math
from numpy.linalg import norm
from matplotlib import pyplot as plt



def get_antenna_pos(xdim,ydim,long_,lalt,space =15, R_E=6e6,a=6378140,h=1000,f=1.0/298.257):
	# geocentric coordinate of ref ant at lat & long
	C=((np.cos(lalt))**2+(1-f)**2*(np.sin(lalt))**2)**-0.5
        S=(1-f)**2*C
	x_g = (a*C+h)*np.cos(np.radians(lalt))*np.cos(np.radians(long_))
	y_g = (a*C+h)*np.cos(np.radians(lalt))*np.sin(np.radians(long_))
	z_g = (a*S+h)*np.sin(np.radians(lalt))
        XYZ_g = np.array([x_g,y_g,z_g])



        # Rotation matrix R
        mm = np.array([[-np.sin(np.radians(long_)), np.cos(np.radians(long_)), 0],
      [-math.sin(np.radians(lalt))*np.cos(np.radians(long_)), -math.sin(np.radians(lalt))*np.sin(np.radians(long_)), math.cos(np.radians(lalt))],
      [math.cos(np.radians(lalt))*np.cos(np.radians(long_)), math.cos(np.radians(lalt))*np.sin(np.radians(long_)), math.sin(np.radians(lalt))]])
	# transforming from geodetic-geocentric coordinate to topocentric coord system
	XYZ_t = XYZ_g

	xx, yy, zz =[],[],[]
	j=0
	while j <=xdim-1:
		xx.append( XYZ_t[0]+ space)
		zz.append( XYZ_t[1] + space)
		yy.append( XYZ_t[2])
		space +=15 # meters
                j +=1
	x,z = np.meshgrid(xx,zz)
        x = np.reshape(x,[x.size,1])
	z = np.reshape(z,[z.size,1])
	y = np.ones(len(z))*XYZ_t[2]

        ant_p =[]
        for i in range(len(x)):
		ant_p.append([int(x[i]),int(y[i]),int(z[i])])



        return np.array(ant_p)





def get_blvectors(Antenna_positions):
	ant_p = np.array(Antenna_positions)
	nAntennas = len(ant_p)

	blVectors = []
	blPairs = []
	for ant1 in range(nAntennas):
   		 for ant2 in range(ant1+1,nAntennas):
        	    if ant2 >ant1:
        		delta = np.array([float(np.round((ant_p[ant1][i]-ant_p[ant2][i]),5)) for i in range(3)])
			
			#print delta, ant1, ant2
			blVectors.append(tuple(delta))
           		blPairs.append((ant1, ant2))
	return [blVectors,blPairs]

def get_ublvectors(Antenna_positions,blVectors,blPairs):
	ublDict = {}
	for b in range(len(blVectors)):
     		# grouping togather all blPairs with same blVector
    		if ublDict.has_key(blVectors[b]): ublDict[blVectors[b]].append(blPairs[b])
    		else: ublDict[blVectors[b]] = [blPairs[b]]

	ublIndexDict = {antPair: i for i,antPairs in enumerate(ublDict.values()) for antPair in antPairs }
	ublVectors = np.array([ant_p[antList[0][0]]-ant_p[antList[0][2]] for antList in ublDict.values()])
	print "With", len(ant_p), "antennas there are", len(ublDict.items()), "unique baselines."
	return [ublVectors,ublIndexDict]
