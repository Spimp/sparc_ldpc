import numpy as np
from random import random 
import py.ldpc as ldpc
from pylab import *

sim_param = [
    ("802.16","1/2",100,"A"),
    ("802.16","2/3",100,"A"),
    ("802.16","2/3",100,"B"),
    ("802.16","3/4",3,"A"),
    ("802.16","3/4",3,"B"),
    ("802.16","5/6",3,"A"),
    ("802.16","1/2",27,"A"),
    ("802.16","2/3",27,"A"),
    ("802.16","2/3",27,"B"),
    ("802.16","3/4",27,"A"),
    ("802.16","3/4",27,"B"),
    ("802.16","5/6",27,"A"),
    ("802.16","1/2",54,"A"),
    ("802.16","2/3",54,"A"),
    ("802.16","2/3",54,"B"),
    ("802.16","3/4",54,"A"),
    ("802.16","3/4",54,"B"),
    ("802.16","5/6",54,"A"),
    ("802.16","1/2",81,"A"),
    ("802.16","2/3",81,"A"),
    ("802.16","2/3",81,"B"),
    ("802.16","3/4",81,"A"),
    ("802.16","3/4",81,"B"),
    ("802.16","5/6",81,"A"),
    ("802.11n","1/2",100,"A"),
    ("802.11n","2/3",100,"A"),
    ("802.11n","3/4",27,"A"),
    ("802.11n","5/6",27,"A"),
    ("802.11n","1/2",54,"A"),
    ("802.11n","2/3",54,"A"),
    ("802.11n","3/4",54,"A"),
    ("802.11n","5/6",54,"A"),
    ("802.11n","1/2",81,"A"),
    ("802.11n","2/3",81,"A"),
    ("802.11n","3/4",81,"A"),
    ("802.11n","5/6",81,"A"),
]
# x is passed through binary sysmetric channel and returned 
def bsc(z, p):
	for i in range(len(z)):
		# flip roughly p of the bits in x
		if random() < p:
			z[i] = 1 - z[i]
	return z

# ch y to LLR
def ch2llr(ch, p):
	LLR = (1-ch)*log((1-p)/p) - ch*log((1-p)/p)

	return LLR

# Do I also need the bpsk (binary phase shift keying) function? Or is it ok to just transmit bits as 0 and 1.
# shouldn't need it if it's not required for LLR equations.

# p is the bsc cross over probability. 
# start with p that gives capicity equal to rate and then reduce p to increase the capacity. 
def sim(standard, rate, z, ptype='A'):

	if rate == "1/2":
		p=0.11
	elif rate == "2/3":
		p=0.06
	else:
		raise NameError("Rate unsupported")

	mycode = ldpc.code(standard, rate, z, ptype)
	N = mycode.N # assuming I can get N in the same way I can get K
	K = mycode.K

	repeats = 200

	print(rate)
	print(N)
	print(mycode.K/N)

	res = []
	parray = np.linspace(p, 0.01, int(p*100))
	print(parray)
	for q in parray:
		#print("q is ", q)
		errors = []
		print("q is ", q)
		for i in range(repeats):

			x = np.zeros(N)
			#u = np.zeros(N)

			# note that the value of x has changed. 
			y = bsc(x, q)

			yl = ch2llr(y, q)
			(app, it) = mycode.decode(yl)
			xh = (app<0.0)
			nbiterror = np.sum(xh)
			errors.append(nbiterror/N)
		#print("errors ", np.sum(errors))
		res.append(np.sum(errors)/repeats)

	return res, parray



if __name__ == "__main__":

	res, parray = sim(*sim_param[1])
	print("The size of parray is ", parray.size)
	print("The size of res is", res)
	plot(parray, res)
	title(sim_param[0])
	xlabel('p')
	ylabel('Bit errors')
	savefig("sim_param2.png")
	show()




