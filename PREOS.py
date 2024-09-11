from funcs import coef,alpha,beta,queue,calcZ_PR,calcI,calclnterm,GroverRT,HroverRT,SroverR
from math import sqrt,pow,cos,acos,pi,log,exp
print(__name__)

Tc = float(input("Tc:"))
T = float(input("temperature:"))
Pc = float(input("Pc:"))
P = float(input("pressure:"))
W = float(input("w:"))
phase = input('phase: "v" or "l":')
Tr = T/Tc
Pr = P/Pc
C = coef(W)
alpha = alpha(Tr,W)
beta = beta(Tr,Pr)
Q = queue(alpha,Tr)
Z = calcZ_PR(beta,Q,phase)
print("C=",C)
print("alpha=",alpha)
print("beta=",beta)
print("Q=",Q)
print("Z=",Z)
input("Enter to continue")
I = calcI(Z,beta)
lnterm = calclnterm(alpha,Tr,C)
Groverstuff = GroverRT(Z,beta,Q,I)
Hroverstuff = GroverRT(Z,lnterm,Q,I)
Sroverstuff = SroverR(Z,beta,lnterm,Q,I)
print("I=",I)
print("d(ln(alpha))/d(Tr)=",lnterm)
print("Gres/RT=",Groverstuff)
print("Hr/RT=",Hroverstuff)
print("Sr/R=",Sroverstuff)
input("Enter to continue")
print("fugacity coefficient=",exp(Groverstuff))
print("fugacity=",exp(Groverstuff) * P)
print("Hres=",Hroverstuff * 8.314 * T)
print("Sres=",Sroverstuff * 8.314)
