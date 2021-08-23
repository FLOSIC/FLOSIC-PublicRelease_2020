
**FLOSIC_CUBE2020**
------------
Cube utility generates cube files of arbitrary functions as a post-processing tool for FLOSIC calculations. 
Please use ./get_cube.py --help for more information.   

**Brief instructions**
------------

To use cube utility:
- Compile the code (hpcc):  
 Make -f Makefile.hpcc  
 Make -f Makefile.hpcc cube  
 (other Makefiles need to be updated)

- Mark the get_cube.py file as executable:  
 chmod +x get_cube.py  

- Load necessary modules (hpcc):  
ml -* iccifort/2019.5.281 impi/2018.5.288 SciPy-bundle/2019.10-Python-3.7.4  
module load matplotlib  
export OMP_NUM_THREADS=8  (to use OMP parallelization)

**Support** 
------------
For any question regarding the cube utility,  
please contact peral1j@cmich.edu (Dr. Peralta) or nguye8t@cmich.edu
