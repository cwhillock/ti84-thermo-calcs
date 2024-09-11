# SINGLE COMPONENT TWO TERM VIRIAL CALCULATION

from math import pow,exp
print(__name__)
Tc = float(input("Tc:"))
T = float(input("temperature:"))
Pc = float(input("Pc:"))
P = float(input("pressure:"))
W = float(input("w:"))
Tr = T/Tc
Pr = P/Pc
B0 = 0.083 - 0.422/pow(Tr,1.6)
B1 = 0.139 - 0.172/pow(Tr,4.2)
Bhat = B0 + W * B1
B = Bhat * 83.14 * Tc / Pc
Z = 1 + B * P / (83.14*T)

print("B0=",B0)
print("B1=",B1)
print("Bhat=",Bhat)
print("B=",B)
print("Z=",Z)
input("Enter to continue")
dB0 = 0.675 / pow(Tr,2.6)
dB1 = 0.722 / pow(Tr,5.2)
fugc = exp(Z-1)
Hres = Pr * (B0 - Tr * dB0 + W * (B1 - Tr * dB1)) * 8.314 * Tc
Sres = -Pr * (dB0 + W * dB1) * 8.314
print("dB0=",dB0)
print("dB1=",dB1)
print("Hr=",Hres)
print("Sr=",Sres)
print("fugacity coefficient=",fugc)
print("fugacity=",fugc*P)
