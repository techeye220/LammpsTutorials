import sys
import numpy as np

# Read all lines from the input and add a last one for lambda=1
lines = open(sys.argv[1],"r").readlines()
lines.append("0.0 1.0 0.0")
# Convert it to a numpy array
data = np.array([line.strip().split() for line in lines[2:]],dtype=float)
lam = data[:,1]
grad = data[:,2]
# Linear interpolation to lambda=1 from the last two lambda values simulated
grad[-1] = np.polyfit(lam[-3:-1],grad[-3:-1],1).sum()
# Integrates the gradient using trapezoid
print "Free energy = %.3f kcal/mol"%np.trapz(grad,x=lam)