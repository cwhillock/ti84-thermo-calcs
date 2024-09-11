from math import sqrt,pow,cos,acos,pi,log,exp
print(__name__)
s = 1 + sqrt(2)
e = 1 - sqrt(2)
omega = 0.0778
psi = 0.45724

def calcB0(Tr):
    return 0.083 - 0.422/pow(Tr,1.6)

def calcdB0(Tr):
    return 0.675 / pow(Tr,2.6)

def calcB1(Tr):
    return 0.139 - 0.172/pow(Tr,4.2)

def calcdB1(Tr):
    return 0.722 / pow(Tr,5.2)

def calcBhat(Tr,W):
    return calcB0(Tr) + W * calcB1(Tr)

def calcB(T,Tc,W,Pc):
    Tr = T/Tc
    return calcBhat(Tr,W) * 83.14 * Tc / Pc

def calcPRa(T,Tc,Pc,W):
    Tr = T/Tc
    return psi * alpha(Tr,W) * 83.14**2 * Tc**2 / Pc

def calcPRb(Tc,Pc):
    return omega * 83.14 * Tc / Pc

def coef(w):
    return 0.37464 + 1.54226 * w - 0.26992 * (w**2)

def alpha(Tr,w):
    return (1 + coef(w) * (1-sqrt(Tr)))**2

def beta(Tr,Pr):
    return omega * Pr/Tr

def queue(alpha,Tr):
    return (psi * alpha) / (omega * Tr)

def calcZ_PR(beta,Q,phase):
    A = (s + e - 1) * beta - 1
    B = (s * e - s - e) * beta**2 + (Q - s - e) * beta
    C =  -(s * e * (1 + beta) + Q) * beta**2

    K = (A**2 - B * 3) / 9
    L = (2 * A**3 - 9 * A * B + 27 * C) / 54
    M = L**2 - K**3

    if M > 0:
        'math.pow(abs(x),float(1)/3) * (1,-1)[x<0]'
        term1 = sqrt(M) - L
        term2 = -L-sqrt(M)
        return pow(abs(term1),float(1)/3) * (1,-1)[term1<0] + pow(abs(term2),float(1)/3) * (1,-1)[term2<0] - A/3
    
    theta = acos(L/pow(K,1.5))
    solutions = [-2*sqrt(K)*cos(theta/3)-A/3,-2*sqrt(K)*cos((theta+2*pi)/3)-A/3,-2*sqrt(K)*cos((theta-2*pi)/3)-A/3]
    if phase == "l":
        return min(solutions)
    else:
        return max(solutions)
    
def calcI(Z,beta):
    return (1/(s-e)) * log((Z+s*beta)/(Z+e*beta))

def calclnterm(alpha,Tr,coeff):
    return -coeff*sqrt(Tr/alpha)

def GroverRT(Z,beta,Q,I):
    return Z - 1 - log(Z-beta) - Q*I

def HroverRT(Z,lnterm,Q,I):
    return Z - 1 + (lnterm - 1) * Q * I

def SroverR(Z,beta,lnterm,Q,I):
    return log(Z-beta) + lnterm * Q * I