import numpy as np

def SOR(A, b, x0, omega):
    x = np.zeros((3, 1))
    R = np.zeros((3, 1))
    D = np.diag(A)
    n = len(A)
      
    # for loop for 3 times as to calculate x, y , z
    for j in range(0, n):
        R[j] = b[j]
        
        # to calculate respective xi, yi, zi
        for i in range(0, n):
            R[j] -= A[j, i]*x[i]
            
        # updating the value of our solution
        x[j] = x[j] + (omega / D[j]) * R[j]
        
    return x

omega = 1.2
epsilon = 1e-3

x0 = np.array([0,0,0]).reshape(3, 1)


A = np.array([[3, 1, 1], 
             [3, 3, 1], 
             [1, 1, 2]])

b = np.array([4, 0, 5])


x = np.zeros((3, 1))
R = np.zeros((3, 1))
D = np.diag(A)
n = len(A)


iter1 = 0

for k in range(100):
    
    iter1 = iter1 + 1
    x_old = x.copy()
    for j in range(0, n):
        R[j] = b[j]
        
        # to calculate respective xi, yi, zi
        for i in range(0, n):
            R[j] -= A[j, i]*x[i]
            
        # updating the value of our solution
        x[j] = x[j] + (omega / D[j]) * R[j]
        
        # Maxixumum relative variation 
        xold_norm = np.linalg.norm(x_old.ravel(), np.inf)
        
        residual = (np.linalg.norm((x.ravel() - x_old.ravel()), np.inf) / 
                    xold_norm)
        
            
            
    if residual < epsilon:   
        break
    else:
        print(k, x.ravel(), residual)
        
        
    
