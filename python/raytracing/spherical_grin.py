import numpy as np

def raytrace_spherical(params):
    # Perform raytracing calculations for a spherical GRIN lens

    #   Simulation parameters
    N = 2^9     # Number of grid points
    NN = 2^9    # Number of propagation steps
    L = 10^(-2) # Step size for ray 

    #   Phisical parameters
    R = 1   # Radius of sphere
    h = 1   # x coordinate of sphere's center 
    k = 1   # y coordinate of sphere's center
    n1 = 1  # Refractive index of outside medium
    lens = 'Luneburg' # Lens type

    if lens == 'Luneburg':
        n = lambda r : np.sqrt(2 - r^2)

    elif lens  == 'Gutman':
        f = 0.75 
        n = lambda r : np.sqrt((1+f^2)-r^2) / f 

    elif lens == 'Non Normalized Gutman':
        a = 0.74
        f = 0.75
        n = lambda r : np.sqrt((R+f)^2 - a*r^2) / (f)

    elif lens == 'Maxwell fisheye':
        n = lambda r : 2 / (1+r^2)
    
    elif lens == 'Concentrator':
        n = lambda r : 1 / r
    
    else: 
        print('Invalid lens type')
        return
    
    #   Make 2D Sphere (matrices)
    [x,y]=np.meshgrid(np.linspace(-R,R,N),np.linspace(-R,R,N)) 
    r = np.sqrt(x^2 + y^2) # Radius of each point in space
    n = n(r) # Refractive index of each point in space
    lens_boundary = lambda y : - np.sqrt(R^2 - (y-k)^2) + h # Boundary of the lens (here y is relative to the center of sphere [-R,R])

    #   Initialize ray
    x0 = -10 + h ; y0 = k # Origin of ray (y constrained to k for e0 calculation)
    H0 = 0.5*R + y0 # Incident height of ray

    x = np.zeros(NN); y = np.zeros(NN)  # x, y coordinates of ray
    r = np.zeros(NN); theta = np.zeros(NN) # r, theta coordinates of ray
    gamma = np.zeros(NN); phi = np.zeros(NN) # gamma, phi propagation angles of ray


    y[0] = H0
    x[0] = lens_boundary(y[0])
    r[0] = np.sqrt(x[0]^2 + y[0]^2)
    theta[0] = np.arcsin(y[0]-y0 / r[0])
    
    #   Raytracing 
    epsilon0 = abs(np.arctan((y[0]-y0)/(x[0]-x0)))
    epsilon1 = epsilon0 + theta[0]
    epsilon2 = np.arcsin(n1/n(r[0]) * np.sin(epsilon1)) # Snell's law for refracted incident ray

    if H0 - y0 > 0:
        phi[0] = np.pi - epsilon2
    elif H0 - y0 < 0:
        phi[0] = epsilon2
        theta[0] = 2*np.pi-theta[0]

    K = n(r[0]) * r[0] * np.sin(phi[0]) # Constant for raytracing

    gamma[0] = theta[0] + phi[0] - np.pi

    for i in range(1,NN):
        x[i] = x[i-1] + np.sign(H0-y0)* L * np.cos(gamma[i-1])
        y[i] = y[i-1] - np.sign(H0-y0)* L * np.sin(gamma[i-1])
        r[i] = np.sqrt((x[i]-h)^2 + (y[i]-k)^2)

        if x[i]-h < 0 and y[i]-k > 0: # Cuadrant 2 
            theta[i] = np.arctan((y[i]-k)/(x[i]-h)) 
        elif x[i]-h > 0 and y[i]-k > 0: # Cuadrant 1
            theta[i] = np.pi - np.arctan((y[i]-k)/(x[i]-h))

        elif x[i]-h < 0 and y[i]-k < 0: # Cuadrant 3
            theta[i] = 2* np.pi - np.arctan((y[i]-k)/(x[i]-h))
        elif x[i]-h > 0 and y[i]-k < 0: # Cuadrant 4
            theta[i] = np.pi+ np.arctan((y[i]-k)/(x[i]-h))
        else:
            print('Error in theta calculation')
            return
        
        if r[i] > R:
            print('Ray has left the lens')
            n = lambda r: n1
        
        if phi[i-1] > np.pi/2 - 0.01: # + 90 degrees
            phi[i] = np.pi - np.arcsin(K/(r[i]*n(r[i])))
        elif phi[i-1] < np.pi/2 + 0.01: # - 90 degrees
            phi[i] = np.arcsin(K/(r[i]*n(r[i])))
        else:
            phi[i] = phi[i-1]
            counter = counter + 1
            print(f'Error in phi calculation {counter} times')

        gamma[i] = theta[i] + phi[i] - np.pi


    return x,y,r,theta,phi,gamma