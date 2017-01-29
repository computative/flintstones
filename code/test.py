from numpy import *
n = 3
state = zeros((n*(n+1),4))

for i in range(n): # for hvert skall
    for j in range(2*(i+1)): # for hver tilstand i skall
        state[i*(i+1) + j][0] = (j - j%2)/2.          # set nx
        state[i*(i+1) + j][1] = i - (j - j%2)/2.  # set ny
        state[i*(i+1) + j][2] = j % 2      # set spinprojeksjon
        state[i*(i+1) + j][3] = i+1      # set energiniva

print state
        
i = 3 # E -1. 
nx = 3
sigma = 1

print i*(i+1) + 2*nx + sigma+1
