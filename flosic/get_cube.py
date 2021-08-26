#!/usr/bin/env python
from __future__ import print_function
from os import listdir, environ
from os.path import isfile, join
import math
import sys
import subprocess
import numpy as np
import argparse
"""
Description: This python program generates cube files as a post-processing tool to the 
FLOSIC code. Cubes can be higly customized. It can generate pre-coded and arbitrary 
functions with minimal changes.  Use python 'get_cube.py --help' for help.
This code must be used in conjunction with the utility cube_exe in the FLOSIC code.

Existing functions so far:
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



get_cube_help = """get_cube is a python script to generate diverse cube files from the FLOSIC code"""




# Default options
options = {"FOLDER"         : [None, "FLOSIC calculation folder"], 
           "CUBE_EXE"       : [None, "cube_exe executable folder"], 
           "FUNCTION"       : [None,"Function to cube"],
           "SPIN"           : [None,"Spin ID"],   # 1=up; 2=down; 3=up+down; 4=up-down
           "THRESHOLD"      : [None,"Maximum value"],  # not used yet
           "CUBENAME"       : [None,"cube file name"]}



# Default parameters (these are the cube parameters only)
parameters = {"XMIN"           : [None,"x min"],
              "XMAX"           : [None,"x max"],
              "YMIN"           : [None,"y min"],
              "YMAX"           : [None,"y max"],
              "ZMIN"           : [None,"z min"],
              "ZMAX"           : [None,"z max"],
              "NXPOINTS"       : [None,"Nx points"],
              "NYPOINTS"       : [None,"Ny points"],
              "NZPOINTS"       : [None,"Nz points"],
              "STYPE"          : [None,"Spin type"]}

# atoms and coordinates holders
atoms  = []
coords = []

#############################################################################################
#
# cube_params
#
#############################################################################################

def cube_params(options,parameters,atoms,coords):
    """
    Read CLUSTER file to get molecular geometry 
    We get back atoms (atomic number list as integers) and coords (x,y,z in Bohr)
    This can probably be done in a more elegant way. All units are Bohr
    """ 
    file_name = join( options["FOLDER"][0] , 'CLUSTER')
    f = open(file_name, 'r') 
    lines =  f.readlines()
    n=lines[2].split()
    natoms = int(n[0].strip())
#    s = lines[3+natoms].split()
#    spin = int( float((s[1].strip())  )
    for i in range(3, 3+natoms):
        string = lines[i].split()
        x = float(string[0].strip())
        y = float(string[1].strip())
        z = float(string[2].strip())
        at = int(string[3].strip())
        coords.append([x,y,z])
        atoms.append(at)
    coords = np.array(coords)
    atoms  = np.array(atoms)    

    # Now that we have the atomic coordinates, possibly change the cube parameters
    change = False

    if (parameters["XMIN"][0] == -10.0):
        parameters["XMIN"][0]  = min(coords.T[0]) - 4.0
        change = True


    if (parameters["XMAX"][0] == +10.0):
        parameters["XMAX"][0]  = max(coords.T[0]) + 4.0
        change = True


    if (parameters["YMIN"][0] == -10.0):
        parameters["YMIN"][0]  = min(coords.T[1]) - 4.0
        change = True


    if (parameters["YMAX"][0] == +10.0):
        parameters["YMAX"][0]  = max(coords.T[1]) + 4.0
        change = True


    if (parameters["ZMIN"][0] == -10.0):
        parameters["ZMIN"][0]  = min(coords.T[2]) - 4.0
        change = True


    if (parameters["ZMAX"][0] == +10.0):
        parameters["ZMAX"][0]  = max(coords.T[2]) + 4.0
        change = True



    print("Cube coordinates (Bohr, min max points)")
    print("X:     %11.7f    %11.7f     %5d" % ( parameters["XMIN"][0], parameters["XMAX"][0], parameters["NXPOINTS"][0]  )  )
    print("Y:     %11.7f    %11.7f     %5d" % ( parameters["YMIN"][0], parameters["YMAX"][0], parameters["NYPOINTS"][0] )  )
    print("Z:     %11.7f    %11.7f     %5d" % ( parameters["ZMIN"][0], parameters["ZMAX"][0], parameters["NZPOINTS"][0] )  )









#############################################################################################
#
# read_data 
#
#############################################################################################



def read_data(options,parameters):
    """
    Reads the command line parameters using the parser
    Here we read the parameters and place them in either options for the cube options or 
    in parameters for the grid parameters, just to keep them separate. Every argument has a 
    default, a type, and a help line.
    """
 
    parser = argparse.ArgumentParser(description=get_cube_help)
    parser.add_argument('folder', metavar='<calculation folder>', type=str, nargs='?',default=".",
                   help='directory containing the FLOSIC calculation files')

    parser.add_argument('cube_exe', metavar='<cube_exe executable folder>', type=str, nargs='?',default="~/PHY790/code/modify/0.2a",
                   help='directory containing the FLOSIC calculation files')
 
    parser.add_argument('--function', metavar='<function>', type=str, nargs='?',default="ELF",
                   help='function to be cubed (string, default = ELF)')

    parser.add_argument('--spin', metavar='<integer>', type=int, nargs='?',default=1,
                   help='spin index (integer, default = 1) alpha=1/beta=2/total=3/spin=4')
                   
    parser.add_argument('--threshold', metavar='<threshold>', type=float, nargs='?',default=1.0E20,
                   help='cutoff above which the funcion is converted to zero (float, default = 1.E20)')

    parser.add_argument('--cube', metavar='<file name>', type=str, nargs='?',default="cube",
                   help='the file name of the output cube (string, default = cube)')
                   
    parser.add_argument('--xmin', metavar='<x min>', type=float, nargs='?',default=-10.0,
                   help='minimum x (float, default = auto)')

    parser.add_argument('--xmax', metavar='<x max>', type=float, nargs='?',default=+10.0,
                   help='maximum x (float, default = auto)')

    parser.add_argument('--ymin', metavar='<y min>', type=float, nargs='?',default=-10.0,
                   help='minimum y (float, default = auto)')

    parser.add_argument('--ymax', metavar='<y max>', type=float, nargs='?',default=+10.0,
                   help='maximum y (float, default = auto)')

    parser.add_argument('--zmin', metavar='<z min>', type=float, nargs='?',default=-10.0,
                   help='minimum z (float, default = auto)')

    parser.add_argument('--zmax', metavar='<z max>', type=float, nargs='?',default=+10.0,
                   help='maximum z (float, default = auto)')

    parser.add_argument('--xpts', metavar='<x points>', type=int, nargs='?',default=40,
                   help='points in x (integer, default = 40)')

    parser.add_argument('--ypts', metavar='<y points>', type=int, nargs='?',default=40,
                   help='points in y (integer, default = 40)')

    parser.add_argument('--zpts', metavar='<z points>', type=int, nargs='?',default=40,
                   help='points in z (integer, default = 40)')

    args = parser.parse_args()
    options["FOLDER"][0] = str(args.folder)
    options["CUBE_EXE"][0] = str(args.cube_exe)
    options["FUNCTION"][0] = str(args.function)
    options["SPIN"][0] = int(args.spin)
    options["THRESHOLD"][0] = float(args.threshold)
    options["CUBENAME"][0] = str(args.cube) 

    parameters["XMIN"][0] = float(args.xmin)
    parameters["XMAX"][0] = float(args.xmax)
    parameters["YMIN"][0] = float(args.ymin)
    parameters["YMAX"][0] = float(args.ymax)
    parameters["ZMIN"][0] = float(args.zmin)
    parameters["ZMAX"][0] = float(args.zmax)
    parameters["NXPOINTS"][0] = int(args.xpts)
    parameters["NYPOINTS"][0] = int(args.ypts)
    parameters["NZPOINTS"][0] = int(args.zpts)


    print("Options:")
    sorted_options = sorted(options.keys())
    for k in sorted_options:
        print(" %-20s  %20s" % (options[k][1],options[k][0]))







#############################################################################################
#
# gen_cube_data 
#
#############################################################################################


def gen_cube_data(options,parameters,atoms,coords):
    """
    This is the main worker. It does:
    1) Sets up a 3D grid from the deta in parameters.
    2) Calls cube_exe to evaluate the 11 quatities on that grid
    3) Process these quantities to obtain the wanted function to put un the cube
    4) Writes the cube file 
    """


    #
    # set up the grid
    #   
    folder = options["CUBE_EXE"][0]
    cube_string = ' cat XYZ.dat   | ' + join( folder,  'cube_exe  ') 
    xmin = parameters["XMIN"][0]
    xmax = parameters["XMAX"][0]
    ymin = parameters["YMIN"][0]
    ymax = parameters["YMAX"][0]
    zmin = parameters["ZMIN"][0]
    zmax = parameters["ZMAX"][0]
    nx   = parameters["NXPOINTS"][0] 
    ny   = parameters["NYPOINTS"][0] 
    nz   = parameters["NZPOINTS"][0] 


    xpts = np.linspace(xmin, xmax, nx, endpoint=True  ) 
    ypts = np.linspace(ymin, ymax, ny, endpoint=True  ) 
    zpts = np.linspace(zmin, zmax, nz, endpoint=True  )  
    xyz = []
    for i in xpts:
         for j in ypts:
             for k in zpts:
                xyz.append([i,j,k])


    xyz_txt = '%s \n' % len(xyz)
    for j in range(len(xyz)):
          xyz_txt += "%s  %s   %s  \n" % ( xyz[j][0], xyz[j][1], xyz[j][2] )

    open('XYZ.dat', "wb").write(xyz_txt.encode())
#
#   call cube_exe
#
    print("Calling cube_exe")
    proc = subprocess.Popen(cube_string, shell=True, stdout=subprocess.PIPE)
    (output,err) = proc.communicate()
    output =  output.decode().split('\n')
#    print(output)


#
#   start processing quantities
#
    n    = [ [],[] ]
    dx   = [ [],[] ]
    dy   = [ [],[] ]
    dz   = [ [],[] ]
    dx2  = [ [],[] ]
    dy2  = [ [],[] ]
    dz2  = [ [],[] ]
    dxdy = [ [],[] ]
    dxdz = [ [],[] ]
    dydz = [ [],[] ]
    tau  = [ [],[] ]


    for line in output:
        if "dens=" in line:
          for isp in (0,1):
            n[isp].append(float(line.split()[4+11*isp]))
            dx[isp].append(float(line.split()[5+11*isp]))
            dy[isp].append(float(line.split()[6+11*isp]))
            dz[isp].append(float(line.split()[7+11*isp]))
            dx2[isp].append(float(line.split()[8+11*isp]))
            dy2[isp].append(float(line.split()[9+11*isp]))
            dz2[isp].append(float(line.split()[10+11*isp]))
            dxdy[isp].append(float(line.split()[11+11*isp]))
            dxdz[isp].append(float(line.split()[12+11*isp]))
            dydz[isp].append(float(line.split()[13+11*isp]))
            tau[isp].append(float(line.split()[14+11*isp]))

#   decide what is found for spin  0=colosed-shell; 1=open-shell
#   do some processing to decode the spin options

    spin = 1
    if (max(n[1]) == 0.0 and min(n[1]) == 0.0): spin = 0
    parameters["STYPE"][0]  = spin




    if (spin == 0) : 
       print("Found a closed-shell calculation")
    if (spin == 1) : 
       print("Found a spin-unrestricted calculation")




    for isp in (0,1):
      n[isp] = np.array(n[isp])
      dx[isp] = np.array(dx[isp])
      dy[isp] = np.array(dy[isp])
      dz[isp] = np.array(dz[isp])
      dx2[isp] = np.array(dx2[isp])
      dy2[isp] = np.array(dy2[isp])
      dz2[isp] = np.array(dz2[isp])
      dxdy[isp] = np.array(dxdy[isp])
      dxdz[isp] = np.array(dxdz[isp])
      dydz[isp] = np.array(dydz[isp])
      tau[isp] = np.array(tau[isp])


#
#   this is that we want in the cube file
    what_spin = options["SPIN"][0]
#
#   now decode some options
#   s_loop is a tuple that tells what spin indeces need to be used
#
    factor =[0.0,0.0]
    if (spin == 0 and what_spin == 1): 
        factor[0] = 1.0
        factor[1] = 0.0
        s_loop = (0,) 
    elif (spin == 0 and what_spin == 2):
        sys.exit('Spin down not available in this calculation') 
    elif (spin == 0 and what_spin == 3):
        factor[0] = 2.0
        factor[1] = 0.0
        s_loop = (0,) 
    elif (spin == 0 and what_spin == 4):
        sys.exit('Spin down not available in this calculation') 
    elif (spin == 1 and what_spin == 1):
        factor[0] = 1.0
        factor[1] = 0.0
        s_loop = (0,) 
    elif (spin == 1 and what_spin == 2):
        factor[0] = 0.0
        factor[1] = 1.0
        s_loop = (1,) 
    elif (spin == 1 and what_spin == 3):
        factor[0] = 1.0
        factor[1] = 1.0
        s_loop = (0, 1) 
    elif (spin == 1 and what_spin == 4):
        factor[0] = 1.0
        factor[1] = -1.0
        s_loop = (0, 1) 
    else:
        sys.exit('Spin option cannot be understood') 
     
#
#   here we decode the functions. More can be added
#   be careful to use [i] only for quantities from cube_exe
#   more functions should be added
#
    what_func = options["FUNCTION"][0]
    function = np.zeros(nx*ny*nz)
    if (what_func == "ELF"): 
        for i in s_loop:
           grad2 = dx[i]**2 + dy[i]**2 + dz[i]**2
           vwke  = grad2/8.0/n[i]
           Dh     = 3.0/10.0*(3*math.pi**2)**(5./3.)*n[i]**(5./3.)
           D = tau[i] - vwke
           function +=   factor[i] * 1.0/( 1.0 + (D/Dh)**2 )
    elif (what_func == "D"):
        for i in s_loop:
           grad2 = dx[i]**2 + dy[i]**2 + dz[i]**2
           vwke  = grad2/8.0/n[i]
           D = tau[i] - vwke
           function +=   factor[i] *D 
    elif (what_func == "Dv2"): 
        for i in s_loop:
           grad2 = dx[i]**2 + dy[i]**2 + dz[i]**2
           vwke  = grad2/4.0/n[i]
           D = tau[i] - vwke
           function +=   factor[i] *D 
    elif (what_func == "DENS"): 
        for i in s_loop:
           function +=   factor[i] *n[i]
    elif (what_func == "TAU"): 
        for i in s_loop:
           function +=   factor[i] *tau[i]
    elif (what_func == "VWKE"): 
        for i in s_loop:
           grad2 = dx[i]**2 + dy[i]**2 + dz[i]**2
           vwke  = grad2/8.0/n[i]
           function +=   factor[i] *vwke

#
#    uncomment these lines to add a custom function with a custom name
#    new arrays do not need index
#    factor contains the spin prefactor information, it depends of 
#    the type of calculation and the type of spin integration requested, 
#    and should not be changed
#    
#    elif (what_func == "MYFUNC"): 
#        for i in s_loop:
#           grad2 = dx[i]**2 + dy[i]**2 + dz[i]**2
#           vwke  = grad2/8.0/n[i]
#           function +=   factor[i] *vwke


    else:
        sys.exit('Function ' + what_func +' cannot be understood') 




    # now we print the cube file 


    cname = options["CUBENAME"][0] + "."  + what_func +".cube"
    function = function.reshape((nx,ny,nz))
   


    xvec =  xpts[1] - xpts[0] 
    yvec =  ypts[1] - ypts[0] 
    zvec =  zpts[1] - zpts[0]
    print("Writting cube file "+ cname)
    with open(cname, "w") as cube:

        # first two lines are comments

        cube.write("Cubefile created by get_cube.py\nCentral Michigan University\n")
        natm = len(atoms)
        cube.write(_cubeline(natm, xmin,ymin,zmin  )) 
        cube.write(_cubeline(nx,   xvec, 0.0 ,0.0 ) )
        cube.write(_cubeline(ny,   0.0, yvec, 0.0 ) )
        cube.write(_cubeline(nz,   0.0 ,0.0, zvec)  )

        # and  finaly write the data   



        for i in range(natm):
            cube.write(_cubeline( atoms[i], atoms[i] , coords[i][0] , coords[i][1], coords[i][2]  )) 
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if (i or j or k) and k%6==0:
                        cube.write("\n")
                    cube.write(" {0: .5E}".format(function[i,j,k]))





#############################################################################################
#
# cubeline
#
#############################################################################################

def _cubeline(*args):
    """
    Generates a line to be written to a cube file 
      *args: first arg is  int and all remaining floats
       returns a formatted string to be written to a cube file
    """
    q = "{0:^ 8d}".format(args[0])
    q += "".join("{0:< 12.6f}".format(arg) for arg in args[1:])
    return q + "\n"


#############################################################################################
#
# main routine
#
#############################################################################################


def main(argv):
    read_data(options,parameters)
    cube_params(options,parameters,atoms, coords)
    gen_cube_data(options,parameters,atoms, coords)
    print("Cube done.")







if __name__ == '__main__':
    main(sys.argv)

