from funcs import calcPRa,calcPRb,calcZ_PR,calcI
from math import sqrt,pow,exp,log
print(__name__)

num = int(input("Num components:"))
T = float(input("Temperature:"))
P = float(input("Pressure:"))
Tc = []
Pc = []
W = []
X = []
for i in range(1,num+1):
    Tc.append(float(input("component {} Tc:".format(i))))
    Pc.append(float(input("component {} Pc:".format(i))))
    W.append(float(input("component {} w:".format(i))))
    X.append(float(input("component {} mol frac:".format(i))))

PRa = []
PRb = []
overallPRb = 0
for i in range(num):
    PRa.append(calcPRa(T,Tc[i],Pc[i],W[i]))
    PRb.append(calcPRb(Tc[i],Pc[i]))
    overallPRb += X[i] * PRb[i]

mixesPRa = {"name":"mixture a's"}
for i in range(1,num):
    for j in range(1,num+1):
        if i >= j:
            continue 
        stringindex = (str(i),str(j))
        mixesPRa[stringindex] = sqrt(PRa[i-1]*PRa[j-1])

overallPRa = 0
for i in range(1,num+1):
    for j in range(1,num + 1):
        if i > j:
            continue
        if i == j:
            overallPRa += X[i-1]**2 * PRa[i-1]
        else:
            overallPRa += 2 * X[i-1] * X[j-1] * mixesPRa[str(i),str(j)]

overallq = overallPRa / (overallPRb * 83.14 * T)
overallbeta = overallPRb * P / (83.14 * T)
abars = []

for k in range(1,num+1):
    tempsum = 2 * X[k-1] * PRa[k-1]
    for i in range(1,num+1):
        if i == k:
            continue
        if i > k:
            stringindex = (str(k),str(i))
        else:
            stringindex = (str(i),str(k))
        tempsum += 2 * X[i-1] * mixesPRa[stringindex]
    abars.append(tempsum - overallPRa)

qbars = []
Z = calcZ_PR(overallbeta,overallq,"v")
I = calcI(Z,overallbeta)
lnFugacityCoefficients = []
fugacityCoefficients = []

for i in range(num):
    qbars.append(overallq * (1 + abars[i]/overallPRa - PRb[i]/overallPRb))
    fugacityCoefficients.append(exp(PRb[i]/overallPRb * (Z -1) - log(Z-overallbeta)-qbars[i]*I))

print("Component a's: [1,2,3,etc]")
for x in PRa: print(x)
input("enter to continue")
print("Component b's: [1,2,3,etc]")
for x in PRb: print(x)
input("enter to continue")
print(mixesPRa["name"])
mixesPRa.pop("name")
for x in mixesPRa: print(x,"=",mixesPRa[x])
input("enter to continue")
print("a:",overallPRa)
print("b:",overallPRb)
print("q:",overallq)
print("beta:",overallbeta)
print("Z:",Z)
print("I:",I)
input("enter to continue")
print("abars [1,2,3,etc]:")
for x in abars: print(x)
input("enter to continue")
print("qbars:")
for x in qbars: print(x)
input("enter to continue")
print("component fugacity coefficients: [1,2,3,etc]")
for x in fugacityCoefficients: print(x)