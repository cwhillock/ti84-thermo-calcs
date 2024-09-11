from funcs import calcBhat,calcB
from math import sqrt,pow,exp
print(__name__)

num = int(input("Num components:"))
T = float(input("Temperature:"))
P = float(input("Pressure:"))
Tc = []
Pc = []
W = []
Vc = []
Zc = []
Y = []
for i in range(1,num+1):
    Tc.append(float(input("component {} Tc:".format(i))))
    Pc.append(float(input("component {} Pc:".format(i))))
    W.append(float(input("component {} w:".format(i))))
    Zc.append(float(input("component {} Zc:".format(i))))
    Vc.append(float(input("component {} Vc:".format(i))))
    Y.append(float(input("component {} Y:".format(i))))

B = []
for i in range(num):
    B.append(calcB(T,Tc[i],W[i],Pc[i]))

mixesW = {"name":"mixture w's"}
mixesTc = {"name":"mixture Tc's"}
mixesZc = {"name":"mixture Zc's"}
mixesVc = {"name":"mixture Vc's"}
mixesPc = {"name":"mixture Pc's"}
mixesbhat = {"name":"mixture bhats"}
mixesB = {"name":"mixture B's"}
deltas = {"name":"mixture deltas"}

for i in range(1,num):
    for j in range(1,num+1):
        if i >= j:
            continue 
        stringindex = (str(i),str(j))
        mixesW[stringindex]= ((W[i-1]+W[j-1])/2)
        mixesTc[stringindex]= ((sqrt(Tc[i-1]*Tc[j-1])))
        mixesZc[stringindex]= ((Zc[i-1]+Zc[j-1])/2)
        mixesVc[stringindex]= ( pow( (pow(Vc[i-1],1/3) + pow(Vc[j-1],1/3)) /2 ,3) )
        mixesPc[stringindex] = mixesZc[stringindex] * 83.14 * mixesTc[stringindex] / mixesVc[stringindex]
        mixesbhat[stringindex] = calcBhat(T/mixesTc[stringindex],mixesW[stringindex])
        mixesB[stringindex] = mixesbhat[stringindex] * 83.14 * mixesTc[stringindex] / mixesPc[stringindex]
        deltas[stringindex] = 2 * mixesB[stringindex] - B[i-1] - B[j-1]

lnFugacityCoefficients = []
fugacityCoefficients = []
for k in range(1,num+1):
    tempsum = 0
    for i in range(1,num+1):
        if i == k:
            continue
        for j in range(1,num+1):
            if j == k or i>j:
                continue
            if i > k:
                stringindex = (str(k),str(i))
            else:
                stringindex = (str(i),str(k))
            if i == j:
                tempsum += Y[i-1]**2 * (2*deltas[stringindex])
            else:
                tempsum += Y[i-1] * Y[j-1] * (2*deltas[stringindex]-deltas[str(i),str(j)])
    fugacityCoefficients.append(exp(P/83.14/T * (B[k-1]+0.5*tempsum)))

print("Component B's: [1,2,3,etc]")
for x in B: print(x)
input("enter to continue")
for i in [mixesW,mixesTc,mixesZc,mixesVc,mixesPc,mixesbhat,mixesB,deltas]:
    print(i["name"])
    i.pop("name")
    for x in i:
        print(x,"=",i[x])
    input("enter to continue")
print("component fugacity coefficients: [1,2,3,etc]")
for x in fugacityCoefficients: print(x)
