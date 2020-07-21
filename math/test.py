import math_API
import time
import numpy as np

l = 12 

m = -2 

theta = np.cos(0.5) 

phi = 0.5 

cos_sin_phi = np.array([np.cos(phi), np.sin(phi)]) 

#math_API.compute_associated_legendre_poly(l, m, theta) 

start = time.time() 

for i in range(10000): 

    #gg = math_API.compute_spherical_harmonics(l, theta, cos_sin_phi)

    dd = math_API.compute_optimized_Y12(np.sin(0.5), theta, cos_sin_phi[0], cos_sin_phi[1])

end = time.time() 

#print ("normal any l: ", end -start)

print ("expand l: ", end -start) 
