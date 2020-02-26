from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt

loads = np.genfromtxt("load_A380.dat", delimiter =",")
testloads = np.ones([81, 41])

zset = []
xset = []
ca = 0.547
la = 2.771
for i in range(1,82):
    thetaz = (i-1)*np.pi/81
    thetaz1 = i*np.pi/81
    z = -0.5*(0.5*ca*(1-np.cos(thetaz)) + 0.5*ca*(1-np.cos(thetaz1)))
    zset.append(z)

for i in range(1,42):
    thetaz = (i-1)*np.pi/41
    thetaz1 = i*np.pi/41
    x = 0.5*(0.5*la*(1-np.cos(thetaz)) + 0.5*la*(1-np.cos(thetaz1)))
    xset.append(x)

zstuff = np.array(zset)
xstuff = np.array(xset)

##x, z = np.meshgrid(xstuff, zstuff)

##fig = plt.figure()
##ax = plt.axes(projection='3d')
##ax.plot_wireframe(z, x, testloads, color='black')
##ax.plot_surface(z, x, loads, rstride=1, cstride=1,
                #cmap='viridis', edgecolor='none')

##plt.show()

'---------------------- Cubic spline interpolation ----------------------------'
##def jacobi(A, b, x0, tol, n_iterations=300):
##    """
##    Performs Jacobi iterations to solve the line system of
##    equations, Ax=b, starting from an initial guess, ``x0``.
##    
##    Returns:
##    x, the estimated solution
##    """
##    
##    n = A.shape[0]
##    x = x0.copy()
##    x_prev = x0.copy()
##    counter = 0
##    x_diff = tol+1
##    
##    while (x_diff > tol) and (counter < n_iterations): #iteration level
##        for i in range(0, n): #element wise level for x
##            s = 0
##            for j in range(0,n): #summation for i !=j
##                if i != j:
##                    s += A[i,j] * x_prev[j] 
##            
##            x[i] = (b[i] - s) / A[i,i]
##        #update values
##        counter += 1
##        x_diff = (np.sum((x-x_prev)**2))**0.5 
##        x_prev = x.copy() #use new x for next iteration
##        
##    
##    print("Number of Iterations: ", counter)
##    print("Norm of Difference: ", x_diff)
##    return x
##
##
##def cubic_spline(x, y, tol = 1e-100):
##    """
##    Interpolate using natural cubic splines.
##    
##    Generates a strictly diagonal dominant matrix then applies Jacobi's method.
##    
##    Returns coefficients:
##    b, coefficient of x of degree 1
##    c, coefficient of x of degree 2
##    d, coefficient of x of degree 3
##    """ 
##    x = np.array(x)
##    y = np.array(y)
##    ### check if sorted
##    if np.any(np.diff(x) < 0):
##        idx = np.argsort(x)
##        x = x[idx]
##        y = y[idx]
##
##    size = len(x)
##    delta_x = np.diff(x)
##    delta_y = np.diff(y)
##    
##    ### Get matrix A
##    A = np.zeros(shape = (size,size))
##    b = np.zeros(shape=(size,1))
##    A[0,0] = 1
##    A[-1,-1] = 1
##    
##    for i in range(1,size-1):
##        A[i, i-1] = delta_x[i-1]
##        A[i, i+1] = delta_x[i]
##        A[i,i] = 2*(delta_x[i-1]+delta_x[i])
##    ### Get matrix b
##        b[i,0] = 3*(delta_y[i]/delta_x[i] - delta_y[i-1]/delta_x[i-1])
##        
##    ### Solves for c in Ac = b
##    print('Jacobi Method Output:')
##    c = jacobi(A, b, np.zeros(len(A)), tol = tol, n_iterations=1000)
##    
##    ### Solves for d and b
##    d = np.zeros(shape = (size-1,1))
##    b = np.zeros(shape = (size-1,1))
##    for i in range(0,len(d)):
##        d[i] = (c[i+1] - c[i]) / (3*delta_x[i])
##        b[i] = (delta_y[i]/delta_x[i]) - (delta_x[i]/3)*(2*c[i] + c[i+1])    
##    
##    return b.squeeze(), c.squeeze(), d.squeeze()
##
##
##coeffs = cubic_spline(zstuff, loads.T[0])
##
##z = np.linspace(zstuff[0], zstuff[-1], 10000000)  #np.linspace(0, -ca, -0.01)
##loadiny =0
##loadsety = []
##j = 0
##realz = []
##
##for i in z:
##    #print(f"going for {i}")
##    #print(zstuff[j], zstuff[j+1], zstuff[j] - zstuff[j+1], zstuff.shape)
##    if i <= zstuff[j] and i > zstuff[j+1]:
##            #print(i)
##            loadiny = loads.T[0][j] + coeffs[0][j]*(i-zstuff[j])+coeffs[1][j]*(i-zstuff[j])**2+coeffs[2][j]*(i-zstuff[j])**3
##            loadsety.append(loadiny)
##            realz.append(i)
##            
##    else: j = j+1
##
##
##plt.plot(realz, loadsety)
##plt.scatter(zstuff, loads.T[0])
##plt.show()


"---------- Linear spline interpolation -----------"

##z=[]
##for i in range(len(zstuff)):
##    z.append(zstuff[len(zstuff)-1-i])
h_z = np.diff(zstuff)
h_aload = np.diff(loads.T[0])

coeffs = np.zeros((80,2))

coeffs.T[0] = loads.T[0][:-1]
coeffs.T[1] = h_aload/h_z

##z = np.linspace(zstuff[0], zstuff[-1], 10000000)  #np.linspace(0, -ca, 10000)
##loadiny =0
##loadsety = []
##j = 0
##realz = []
##
##for i in z:
##    #print(f"going for {i}")
##    #print(zstuff[j], zstuff[j+1], zstuff[j] - zstuff[j+1], zstuff.shape)
##    if i <= zstuff[j] and i > zstuff[j+1]:
##            #print(i)
##            loadiny = coeffs.T[0][j] +coeffs.T[0][j]*(i-zstuff[j])
##            loadsety.append(loadiny)
##            realz.append(i)
##            
##    else: j = j+1
zvalue = zstuff[0]
for i in range(len(zstuff) -1):
    loadsety = []
    realz = []
    while zvalue > zstuff[i+1]:
        loadiny = coeffs.T[0][i] + coeffs.T[1][i]*(zvalue-zstuff[i])
        #print(loadiny)
        loadsety.append(loadiny)
        realz.append(zvalue)
        zvalue -= 0.00005
    plt.plot(realz,loadsety)

    
plt.scatter(zstuff, loads.T[0])
plt.show()
