#!/usr/bin/python

import os

logfile='gitlog.txt'
os.system('git log | head -3 > {0}'.format(logfile))
outfile=open('gitversion.f90','w')
outfile.write('SUBROUTINE GITVERSION\n\n')
if os.path.isfile('{0}'.format(logfile)):
    infile=open('{0}'.format(logfile),'r')
    line=infile.readline()
    line=line.strip()
    outfile.write('PRINT \'(A)\',\'{0}\'\n'.format(line))
    line=infile.readline()
    line=line.strip()
    outfile.write('PRINT \'(A)\',\'{0}\'\n'.format(line))
    line=infile.readline()
    line=line.strip()
    outfile.write('PRINT \'(A)\',\'{0}\'\n'.format(line))
    infile.close()
    os.system('rm {0}'.format(logfile))
else:
    outfile.write('PRINT \'(A)\',\'GIT REPOSITORY INFORMATION UNAVAILABLE\'\n')
outfile.write('\nEND SUBROUTINE GITVERSION\n')
outfile.close()
