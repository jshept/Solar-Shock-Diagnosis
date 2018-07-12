########################################################################
# Shock diagnosing tool
# Developed by Wenzhi Ruan, Jiansen He, Limei Yan
# Email: wruan@pku.edu.cn; jshept@pku.edu.cn; lmyan@mail.iggcas.ac.cn
# This code derives the parameters of the observed shocks based on
# the method presented by Ruan, Yan, He et al., (2018)
########################################################################

import sys
import numpy as np
import matplotlib.pyplot as plt
import idlsave
import math


#Constants
gamma=5.0/3
kB=1.38e-23     
mp=1.6726e-27   #mass of proton


#---------------------------------------------------------------------#
#-----------------------Parameters input begin------------------------#
#---------------------------------------------------------------------#
#Input
I1=2321         #upstream line intensity
D1=12.2e3       #[m/s],upstream Doppler velocity
I2=3481         #downstream line intensity
D2=-3.4e3       #[m/s], downstream Doppler velocity
vbeta=2.1e4     #[m/s], component of propagation speed perpendicular 
                #to the LOS direction


#Assumptions
mpr=1.1*mp/2.0   #average mass of charged particles (ions + electrons) 
                 #local acoustic speed is assumped to be 
                 #sqrt(gamma*kB*T/mpr)

#Emission line
Line='SiIV'

#Temperature of the emission line
T0=80000        #[K]

#Step for iteration. A small step lead to accurate results but a long 
#running time
step=500.0     #[m/s]

#---------------------------------------------------------------------#
#-----------------------Parameters input end--------------------------#
#---------------------------------------------------------------------#




#---------------------------------------------------------------------#
#----------------------------Data base begin--------------------------#
#---------------------------------------------------------------------#

#Contribution functions of Si IV emission from CHIANTI Atomic Database
T_SiIV=[4.,         4.1,        4.2,        4.3,       
        4.40000001, 4.50000001, 4.60000001, 4.70000001, 
        4.80000001, 4.90000001, 5.00000001, 5.10000002,
        5.20000002, 5.30000002, 5.40000002, 5.50000002, 
        5.60000002, 5.70000003, 5.80000003, 5.90000003, 
        6.00000003, 6.10000003, 6.20000003, 6.30000003,
        6.40000004, 6.50000004, 6.60000004, 6.70000004, 
        6.80000004, 6.90000004, 7.00000004, 7.10000005, 
        7.20000005, 7.30000005, 7.40000005, 7.50000005,
        7.60000005, 7.70000006, 7.80000006, 7.90000006, 
        8.00000006]

GT_SiIV=[0.00000000e+00, 1.93533655e-40, 7.15981215e-36, 2.04686128e-32,
         5.72640405e-30, 3.23092896e-28, 6.85244696e-27, 7.55682283e-26,
         4.67101473e-25, 9.55142072e-25, 4.33380741e-25, 1.45387404e-25,
         5.81942295e-26, 3.07693806e-26, 2.01894530e-26, 1.18998796e-26,
         4.49471491e-27, 9.76143303e-28, 1.29855510e-28, 1.11544509e-29,
         6.23585371e-31, 2.22202348e-32, 4.31947134e-34, 3.48537142e-36,
         0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
         0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
         0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
         0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
         0.00000000e+00]

#---------------------------------------------------------------------#
#----------------------------Data base end----------------------------#
#---------------------------------------------------------------------#





#---------------------------------------------------------------------#
#-------------------------Define functions begin----------------------#
#---------------------------------------------------------------------#
#Contribution function GT(T)
def GT(Te):
  global Line
  global Tl, GTl

  logTe=math.log10(Te)
  if (logTe<min(Tl) or logTe>max(Tl)):
    print 'ERROR! T out of range of G(T)! T={}'.format(Te)

  for iT in range(0,len(Tl)):
    if (Tl[iT]<=logTe and Tl[iT+1]>logTe):
      iTe=iT
  GTe=(GTl[iTe+1]-GTl[iTe])*(logTe-Tl[iTe])/(Tl[iTe+1]-Tl[iTe]) + GTl[iTe]
  return GTe*1e24



#Local acoustic speed Cs(T)
def Cs(Te):
  global gamma, kB, mpr

  CsT=np.sqrt(gamma*kB*Te/mpr)
  return CsT



#Square upstream Mach number M1^2(va)
def M1s(va):
  global D1, D2, gamma

  temp1=2*(D1+va)
  temp2=(gamma+1)*(D2+va)-(gamma-1)*(D1+va)
  if (temp2==0):
    print 'ERROR! temp2==0 in M1s. va={}'.format(va)
    sys.exit()
  M1sVa=temp1/temp2
  return M1sVa



#Upstream temperature T1(va)
def T1(va):
  global D1, vbeta, gamma, kB, mpr

  temp1=(va**2+vbeta**2)*(D1+va)**2
  temp2=mpr/kB
  temp3=gamma*M1s(va)*va**2
  T1Va=temp1*temp2/temp3
  return T1Va



#Downstream temperature T2(va)
def T2(va):
  global gamma

  temp1=2*(gamma-1)/((gamma+1)**2)
  temp2=(gamma*M1s(va)+1)/M1s(va)
  temp3=M1s(va)-1
  T2Va=(1+temp1*temp2*temp3)*T1(va)
  return T2Va



#Ratio of density R(va)=rho2/rho1=sqrt(I2/GT(T2))/sqrt(I1/GT(T1))
def R(va):
  global I1, I2

  rho1=np.sqrt(I1/GT(T1(va)))
  rho2=np.sqrt(I2/GT(T2(va)))
  RVa=rho2/rho1
  return RVa



#Function for iteration va=f(va)
def f(va):
  global I1, I2, D1, D2

  rho1=np.sqrt(I1/GT(T1(va)))
  rho2=np.sqrt(I2/GT(T2(va)))
  fVa=rho2/rho1 - (D1+va)/(D2+va)
  return fVa



#Gradient of f Gradf(va,dva)=f'(va)
def Gradf(va,dva):
  dfdva=(f(va+dva)-f(va-dva))/(2*dva)
  return dfdva
  


#---------------------------------------------------------------------#
#-------------------------Define functions end------------------------#
#---------------------------------------------------------------------#




#---------------------------------------------------------------------#
#-----------------------------------Main------------------------------#
#---------------------------------------------------------------------#

#Check database
if (Line=='SiIV'):
  Tl=T_SiIV
  GTl=GT_SiIV
else:
  print 'ERROR! No data for the emission line {}'.format(Line)
  sys.exit()


#Check the domain of GT(va)
T_min=Tl[0]
T_max=Tl[len(Tl)-1]
for iT in range(0,len(Tl)-1):
  if (GTl[iT]==0 and GTl[iT+1]>0):
    T_min=Tl[iT+1]
  if (GTl[iT]>0 and GTl[iT+1]==0):
    T_max=Tl[iT]
T_min=10**T_min
T_max=10**T_max

vaP=((gamma-1)*D1-(gamma+1)*D2)/2.0   #Pole of M1(va)
vaZ=-D1   #Zero of M1(va)
dva=10.0
va=max(vaP,vaZ)+dva
while (va<10*Cs(T0)):
  va_next=va+dva
  if (T1(va)<T_min or T2(va)<T_min):
    if (T1(va_next)>=T_min and T2(va_next)>=T_min):
      va_min=va_next
  if (T1(va)<=T_max and T2(va)<=T_max):
    if (T1(va_next)>T_max or T2(va_next)>T_max):
      va_max=va
  va=va_next


#Looking for the solution via Iteration
va=Cs(T0)
if (va<va_min or va>va_max):
  va=(va_max+va_min)/2.0
dva=step
f_current=f(va)
while (f_current>1.0e-4):
  va=va-dva*Gradf(va,dva)
  f_current=f(va)


#------------------------ Results ------------------------------------#
#Two significant digits for v_alpha, Tu and Td
v_alpha=round(va/1.0e3)          #km/s
Tu=round(T1(va)/1.0e3)/10.0         #10^4 K
Td=round(T2(va)/1.0e3)/10.0         #10^4 K

#Two significant digits for vp
v_beta=round(vbeta/1.0e3)        #km/s
vp=round(math.sqrt(1.0*v_alpha**2+1.0*v_beta**2))      #km/s

#Two significant digits for theta
theta=math.atan(1.0*v_beta/v_alpha)
theta=round(theta*180/np.pi)    #degree

#One decimal place for vu and vd
vu=-D1/math.cos(theta/180.0*np.pi)
vu=round(vu/100)/10.0   #km/s
vd=-D2/math.cos(theta/180.0*np.pi)
vd=round(vd/100)/10.0   #km/s

#Two significant digits for M1 and M2
M1=abs(vu*1.0e3-vp*1.0e3)/Cs(Tu*1.0e4)
M2=abs(vd*1.0e3-vp*1.0e3)/Cs(Td*1.0e4)
M1=round(M1*10)/10.0
M2=round(M2*10)/10.0


print '#----------------------------------------------------------------#'
print '#---------------------------  INPUT  ----------------------------#'
print '#----------------------------------------------------------------#'
print 'Emission line:'
print '{}'.format(Line)
print 'Upstream emission intensity:'
print 'I1 = {}'.format(I1)
print 'Upstream Doppler shift:'
print 'D1 = {:.1f} km/s'.format(D1/1000.)
print 'Downstream emission intensity:'
print 'I1 = {}'.format(I2)
print 'Downstream Doppler shift:'
print 'D1 = {:.1f} km/s'.format(D2/1000.)
print 'Component of propagation speed perpendicular to the LOS direction:'
print 'v_beta = {:.0f} km/s'.format(v_beta)
print ' '
print ' '

print '#----------------------------------------------------------------#'
print '#---------------------------  OUTPUT  ---------------------------#'
print '#----------------------------------------------------------------#'
print 'Component of propagation speed parallel to the LOS direction:'
print 'v_alpha = {:.0f} km/s'.format(v_alpha)
print 'Upstream temperature:'
print 'T1 = {:.1f}x10^4 K'.format(Tu)
print 'Downstream temperature:'
print 'T2 = {:.1f}x10^4 K'.format(Td)
print ' '

print 'Propagation speed:'
print 'vp = {:.0f} km/s'.format(vp)
print 'Angle between the LOS direction and the propogation direction:'
print 'theta = {:.0f} degree'.format(theta)
print ' '

print 'Upstream velocity:'
print 'v1 = {:.1f} km/s'.format(vu)
print 'Downstream velocity:'
print 'v2 = {:.1f} km/s'.format(vd)
print ' '

print 'Upstream Mach number:'
print 'M1 = {:.1f}'.format(M1)
print 'Downstream Mach number:'
print 'M2 = {:.1f}'.format(M2)
