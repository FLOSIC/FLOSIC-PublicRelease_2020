import os
import math
import subprocess
import numpy as np
import matplotlib
import tkinter
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


"""
Description: This python program generates the linear plot
for related-density functions:
   ELF    Electron localization function
   D      Tau -1/8*grad(n)^2/n (as used in the ELF)
   Dv2    Tau/2 -1/8*grad(n)^2/n (for testing)
   DENS   Electron density
   TAU    Positive kinetic energy density
   VWKE   Von Weiszacker kinetic energy density  

Copyright by Juan E. Peralta and Duyen B. Nguyen, 2021
Please contact peral1j@cmich.edu (Dr.Peralta) or nguye8t@cmich.edu
for any questions, updates, modifications, and bug reports.
"""

#**********  get mesh rid ********** 
start   = np.array([0, 0, -2])
end     = np.array([0, 0,  2])
npoints = 500
delta   = np.add(end,-start)/npoints


xyz = []
for i in range(npoints):
  xyz.append(np.add(start,i*delta))


print("generating points done.")
xyz = np.asarray(xyz)



#********** read atomic geometry from CLUSTER file ************* 
file_name = open("CLUSTER", "r")
lines     = file_name.readlines()
data      = []
n      = lines[2].split()
natoms = int(n[0].strip())
s      = lines[3+natoms].split()
spin   = float(s[1].strip())



for i in range(3, 3+natoms):
    st = lines[i].split()
    x  = float(st[0].strip())
    y  = float(st[1].strip())
    z  = float(st[2].strip())
    data.append([x,y,z])
    reference= np.asarray(data)



#********** calling the cube_exe with points generated as arguments ********** 
cube_string = ' cat XYZ.dat   | ~/git_8/PublicRelease_2020/flosic/cube_exe '

xyz_txt = '%s \n' % npoints 
for j in range(len( xyz )):
       xyz_txt += "%s  %s   %s  \n" % ( xyz[j,0], xyz[j,1], xyz[j,2] )

open('XYZ.dat', "wb").write(xyz_txt.encode())

proc = subprocess.Popen(cube_string, shell=True, stdout=subprocess.PIPE)
(output,err) = proc.communicate()

print("calling cube_exe done.")

output =  output.decode().split('\n')


#********** get all qualities ********** 
rho_up = []
rho_dx_up   = []
rho_dy_up   = []
rho_dz_up   = []
rho_dx2_up  = []
rho_dy2_up  = []
rho_dz2_up  = []
rho_dxdy_up = []
rho_dxdz_up = []
rho_dydz_up = []
tau_up = []

rho_down = []
rho_dx_down   = []
rho_dy_down   = []
rho_dz_down   = []
rho_dx2_down  = []
rho_dy2_down  = []
rho_dz2_down  = []
rho_dxdy_down = []
rho_dxdz_down = []
rho_dydz_down = []
tau_down = []


if spin == 0:
  for line in output:
    if "dens=" in line:
        rho_up.append(float(line.split()[4]))
        rho_dx_up.append(float(line.split()[5]))
        rho_dy_up.append(float(line.split()[6]))
        rho_dz_up.append(float(line.split()[7]))
        rho_dx2_up.append(float(line.split()[8]))
        rho_dy2_up.append(float(line.split()[9]))
        rho_dz2_up.append(float(line.split()[10]))
        rho_dxdy_up.append(float(line.split()[11]))
        rho_dxdz_up.append(float(line.split()[12]))
        rho_dydz_up.append(float(line.split()[13]))
        tau_up.append(float(line.split()[14]))   
else:
  for line in output:
    if "dens=" in line:
        rho_up.append(float(line.split()[4]))
        rho_dx_up.append(float(line.split()[5]))
        rho_dy_up.append(float(line.split()[6]))
        rho_dz_up.append(float(line.split()[7]))
        rho_dx2_up.append(float(line.split()[8]))
        rho_dy2_up.append(float(line.split()[9]))
        rho_dz2_up.append(float(line.split()[10]))
        rho_dxdy_up.append(float(line.split()[11]))
        rho_dxdz_up.append(float(line.split()[12]))
        rho_dydz_up.append(float(line.split()[13]))
        tau_up.append(float(line.split()[14]))
        rho_down.append(float(line.split()[15]))
        rho_dx_down.append(float(line.split()[16]))
        rho_dy_down.append(float(line.split()[17]))
        rho_dz_down.append(float(line.split()[18]))
        rho_dx2_down.append(float(line.split()[19]))
        rho_dy2_down.append(float(line.split()[20]))
        rho_dz2_down.append(float(line.split()[21]))
        rho_dxdy_down.append(float(line.split()[22]))
        rho_dxdz_down.append(float(line.split()[23]))
        rho_dydz_down.append(float(line.split()[24]))
        tau_down.append(float(line.split()[25]))


#********** alpha density ********** 
rho_up      = np.array(rho_up)
rho_dx_up   = np.array(rho_dx_up)
rho_dy_up   = np.array(rho_dy_up)
rho_dz_up   = np.array(rho_dz_up)
rho_dx2_up  = np.array(rho_dx2_up)
rho_dy2_up  = np.array(rho_dy2_up)
rho_dz2_up  = np.array(rho_dz2_up)
rho_dxdy_up = np.array(rho_dxdy_up)
rho_dxdz_up = np.array(rho_dxdz_up)
rho_dydz_up = np.array(rho_dydz_up)
tau_up      = np.array(tau_up)

#********** beta density ********** 
rho_down      = np.array(rho_down)
rho_dx_down   = np.array(rho_dx_down)
rho_dy_down   = np.array(rho_dy_down)
rho_dz_down   = np.array(rho_dz_down)
rho_dx2_down  = np.array(rho_dx2_down)
rho_dy2_down  = np.array(rho_dy2_down)
rho_dz2_down  = np.array(rho_dz2_down)
rho_dxdy_down = np.array(rho_dxdy_down)
rho_dxdz_down = np.array(rho_dxdz_down)
rho_dydz_down = np.array(rho_dydz_down)
tau_down      = np.array(tau_down)

#********** Functions **********
#Von-Weizsacker Kinetic Energy density 
wke_up   = (rho_dx_up**2  +rho_dy_up**2. +rho_dz_up**2)  /(8*rho_up)
wke_down = (rho_dx_down**2+rho_dy_down**2+rho_dz_down**2)/(8*rho_down)

#The Mobility Function
mofunc_up   = tau_up - wke_up
mofunc_down = tau_down - wke_down
Dn = (3/10)*(math.pi**2)**(5/3)*rho_up**(5/3)

#Electron Localization Function
ELF_up   = 1/(1+(mofunc_up/Dn)**2)
Fhole_up = 8*tau_up/rho_up - (rho_dx_up**2+ rho_dy_up**2+ rho_dz_up**2)/rho_up**2


#print the maximumn value
#ELF_max1=max(rho_up)
#max_x1  = xyz[rho_up.argmax()]

#print(ELF_max1)
#print(max_x1)


#********** plotting ********** 
fig, ax= plt.subplots()
ax.plot( xyz.T[2], rho_up ,color='red', label='density')
plt.legend()
ax.set(xlabel='d (Bohr)')

print("generating plot done.")
plt.show()
