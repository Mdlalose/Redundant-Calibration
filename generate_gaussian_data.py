import numpy
from matplotlib import pyplot as plt
xdim = 8
ydim = 8
dnu=1
nfreq=100
corr_len=10/2.35

freq=numpy.arange(-nfreq,nfreq)*dnu;
mycorr=numpy.exp(-0.5*freq**2/corr_len**2)
#plt.ion()
#plt.clf()


mycorr_ft=numpy.abs(numpy.fft.rfft(mycorr))


mysig=numpy.sqrt(mycorr_ft)

Ant_eps = []
nn_sample = 100
for n_i in range(nn_sample):
	ant_eps = []
	for i in range(64):
		datft=numpy.random.randn(len(mysig))+numpy.complex(0,1)*numpy.random.randn(len(mysig))
		datft=datft*mysig
		datft[0]=numpy.real(datft[0])
		datft[-1]=numpy.real(datft[-1])

		dat=numpy.fft.irfft(datft)

		dat = dat[len(dat)/2:]
		dat = (dat/numpy.std(dat))*0.05
		ant_eps.append(dat)
		plt.plot(dat)
	


 	ant_eps = numpy.column_stack(ant_eps)
	Ant_eps.append(ant_eps)


numpy.save('Ant_epsolon0.05.npy',Ant_eps)

"""
RA_src, Dec_src = numpy.load('src_pos.npy')[0], numpy.load('src_pos.npy')[1]
src_flux = numpy.load('src_flux.npy')[1]
x, y, z = numpy.cos(numpy.radians(Dec_src))*numpy.cos(numpy.radians(RA_src)), numpy.cos(numpy.radians(Dec_src))*numpy.sin(numpy.radians(RA_src)), numpy.sin(numpy.radians(Dec_src))
zen_vec =  [0.0,1.0,0.0]
#xyz = np.concatenate((x,y,z))
mydot = x*zen_vec[0]+ y*zen_vec[1] + z*zen_vec[2]
r=numpy.sqrt(x**2 + y**2 + z**2)
theta = numpy.arccos(mydot)
sigma = 0.5*(1.22*(0.5/6.0))
beam = numpy.exp(-0.5*(theta/sigma)**2)
ii =  beam > 0.1

RA_src= RA_src[ii]
Dec_src= Dec_src[ii]
src_flux = src_flux[ii]
x, y, z = numpy.cos(numpy.radians(Dec_src))*numpy.cos(numpy.radians(RA_src)), numpy.cos(numpy.radians(Dec_src))*numpy.sin(numpy.radians(RA_src)), numpy.sin(numpy.radians(Dec_src))
r=numpy.sqrt(x**2 + y**2 + z**2)
mydot = (x/r)*zen_vec[0] + (y/r)*zen_vec[1] + (z/r)*zen_vec[2]
theta = numpy.arccos(mydot)
for k in range(len(dat)):

	beam1 = numpy.exp(-0.5*(theta/(sigma))**2)
	beam12 = numpy.exp((-0.5*(theta)**2)*(1.0/(sigma + numpy.random.choice(dat)) + 1.0/(sigma + numpy.random.choice(dat)))**2)
	plt.title('Primary Beam Variations')
	plt.plot(theta,beam12, '.')
	plt.plot(theta,beam1,'*')

"""


