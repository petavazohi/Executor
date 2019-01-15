#!/usr/bin/env python3

import os
import subprocess
import argparse
import time
import re
import numpy as np 
import pychemia


def create_kpoints(length,work_dir):
    length = str(length)+'\n'
    comment = 'Automatic mesh'+'\n'
    nkp     = str(0)+'\n'
    method  = 'Auto'+'\n'
    wf = open(work_dir+os.sep+"KPOINTS",'w')
    wf.write(comment)
    wf.write(nkp    )
    wf.write(method )
    wf.write(length )
    wf.close()
    return

def create_incar(mode, sys_name,work_dir):
    if mode == "kpoint_convergence":
        incar = pychemia.code.vasp.VaspInput()
        incar['ENCUT']  = 400
        incar['SYSTEM'] = sys_name
        incar['EDIFF' ] = 1e-5
        incar['NWRITE'] = 2
        incar['PREC'  ] = 'Accurate'
        incar['NCORE']  = 4
        incar.write(work_dir+os.sep+"INCAR")
    elif mode == "workfunction":
        incar = pychemia.code.vasp.VaspInput()
        incar.set_encut(1.4,POTCAR='POTCAR')
        incar['SYSTEM'] = sys_name
        incar['IBRION'] = -1
        incar['ISTART'] = 0
        incar['EDIFF'] = 1e-6
        incar['NELMIN'] = 6
        incar['NCORE'] = 4
        incar['LREAL'] = 'Auto'
        incar['LVTOT'] =  True
        incar.write(work_dir+os.sep+"INCAR")
    return 

def kpoint_convergence(e_threshold,step,executable,nparal):
    address_kpoint = os.getcwd() 
    toten = []
    kmesh = []
    wf = open(address_kpoint+os.sep+"convergence",'w')
    create_incar(mode='kpoint_convergence', sys_name=address_kpoint.split('/')[-1],work_dir=address_kpoint)
    # create potcar here 
    for klength in np.arange(1,11)*step:
        create_kpoints(klength,work_dir=address_kpoint)
        runtime = execute(nparal,address_kpoint,executable)
        rf = open(address_kpoint+os.sep+"OUTCAR",'r')
        data = rf.read()
        rf.close()
        toten.append(float(re.findall("TOTEN\s*=\s*([-+0-9.]*)\s*eV",data)[-1]))
        kmesh.append(re.findall("generate k-points for.*",data)[0])
        wf.write("kpoint length = %i ,kmesh = %s, TOTEN =%f \n" %(klength,kmesh[-1],toten[-1]))
        if klength/10 > 1 : # this is to check if it's not doing the calculation for the first time
            change = abs(toten[-2]-toten[-1])
            if change < e_threshold :
                print("VASP calculations converged with k points length %i and kmesh %s " % (klength,kmesh[-1]))
                wf.close()
                break
            else : 
                os.remove(address_kpoint+os.sep+"CHG")
                os.remove(address_kpoint+os.sep+"CHGCAR")
                os.remove(address_kpoint+os.sep+"WAVECAR")
                
        else : 
            os.remove(address_kpoint+os.sep+"CHG")
            os.remove(address_kpoint+os.sep+"CHGCAR")
            os.remove(address_kpoint+os.sep+"WAVECAR")
    return 


def execute(nparal,work_dir,vasp_exe):

     wf=open(work_dir+os.sep+'RUNNING','w')
     wf.write(time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime()))
     wf.close()
     start_time=time.time()
     status = subprocess.call("mpirun -np {} {}".format(nparal,vasp_exe),cwd=work_dir ,shell=True)
     status = 0
     end_time=time.time()
     runtime=end_time-start_time
     if status== 0:
         print("VASP execution completed with returcode: %d runtime: %d secs" % (status, runtime))
     else:
         print("VASP execution failed with returcode: %d runtime: %d secs" % (status, runtime))
     os.remove(work_dir+os.sep+'RUNNING')
     return runtime

if __name__ == "__main__" : 
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--structure",dest="structure",type=str,help='poscar structure that you want to run',default = 'POSCAR')
    parser.add_argument("-np" ,dest="np",type=int ,action="store", help="Number of MPI processes for the code")
    parser.add_argument("--mode",dest="mode",type=str,action="store",help="type of calculation")
    parser.add_argument("--executable",dest="executable",type=str,action="store",help="vasp executable",default="vasp_std")
    args = parser.parse_args()
    if args.mode == "kpoint_convergence":
        kpoint_convergence(e_threshold=10e-3,step=10,executable=args.executable,nparal=args.np)

