import sympy as sp
import numpy as np



P = [[1,2],[2,9],[3,4],[7,9]]

def PolynomeInterpolateur(D):
    x = sp.var("x")
    P = 0
    L = 1
    for i in range(len(D)):
        for j in range(len(D)):
            if j!=i:
                L = L*((x-D[j][0])/(D[i][0] - D[j][0]))
        P = P + D[i][1]*L
        L = 1
    return P.expand()

def Spline_Cubique_Naturel(D,x):
    points = sorted(D, key=lambda point: point[0])
    H = []
    F = []
    for i in range(len(D)-1):
        H.append(points[i+1][0] - points[i][0])
    for i in range(len(D)):
        if (i == 0) or (i == (len(D) - 1)):
            F.append(0)
        else:
            F.append((D[i+1][1] - D[i][1])/H[i] - (D[i][1] - D[i-1][1])/H[i-1])
    R = np.zeros([len(D),len(D)])
    R[0][0] = 1
    R[-1][-1] = 1
    for i in range(len(D)):
        if (i != 0) and (i != (len(D) - 1)):
            R[i][i] = (H[i-1] + H[i])/3
            R[i][i+1] = H[i]/6
            R[i][i-1] = H[i-1]/6
    M = np.dot(np.linalg.inv(R),F)
    C = []
    Cp = []
    for i in range(len(D)-1):
        C.append((D[i+1][1] - D[i][1])/H[i] - (H[i]/6)*(M[i+1] - M[i])) 
        Cp.append(D[i][1] - M[i]*H[i]**2/6) 
    S = 0
    for i in range(len(D)):
        if (points[i][0] <= x) and (x <= points[i+1][0]):
            S = i
    return M[S]*(D[S+1][0] - x)**3/(6*H[S]) - M[S+1]*(D[S][0] - x)**3/(6*H[S]) + C[S]*(x - D[S][0]) + Cp[S]

