import numpy as np
import py.ldpc as ldpc


# x is the input 
# xp is x predicted, i.e. the prediction of the input to the channel based on the channel output.
# calculate the errors between the two.
def error_count(x, xp):
	errors = [(xp!=x)]
	errors = np.reshape(errors, -1)
	return sum(errors)

# standard is the ldpc standard to use. Options are '802.11n' or '802.16'
# rate can be '1/2', '2/3', '3/4', '5/6'
# z is the size of each matrix that replaces an entry in a protograph matrix
# ptype is 'A' or 'B' (only needed for 802.16, rate 2/3 or 3/4)

standard = '802.16'
rate = '1/2'
z = 3
ptype = 'A'
#standard deviation of gaussian noise.
sigma = 0.8

# initialise an ldpc object c
c = ldpc.code()#standard, rate, z, ptype)

# u is the information vector to be encoded
u = np.random.randint(0,2,c.K)

# Encode u using the ldpc object
x = c.encode(u)

#y = 10*(.5-x)
#generate gaussian noise.
n = np.random.normal(loc=0.0, scale=sigma, size=len(x))
# Add noise to the codeword x to get the channel output.
y = x + n 

# note the log likelihoods should be what gets passed into decoder. L(y_j) = 2(y_j)/sigma^2
# log likelihood of y
Ly = 2.0*y/(sigma**2)
app, it = c.decode(Ly)

# prediction of x. 
xp = [app<0]
xp = np.reshape(xp, -1)

#print(x.shape)
#print(xp.shape)
errors = error_count(x, xp)
print("Number of errors: ", errors)
print("Error rate: ", (errors/len(x)))
print(it)