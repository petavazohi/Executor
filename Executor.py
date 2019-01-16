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

def kpoint_convergence(e_threshold,start,end,step,executable,nparal):
    kpoint_convergence(e_threshold=args.Ethreshold,start=args.Kstart,end=args.Kend,step=args.Kstep,executable=args.executable,nparal=args.np)
    address_kpoint = os.getcwd() 
    toten = []
    kmesh = []
    wf = open(address_kpoint+os.sep+"kpoint_convergence",'w')
    incar = pychemia.code.vasp.VaspInput()
    incar.set_encut(1.4,POTCAR='POTCAR')
    incar['EDIFF' ] = 1e-4
    incar['NWRITE'] = 2
    incar['PREC'  ] = 'Accurate'
    incar['NCORE' ] = 4
    incar['SYSTEM'] = address_kpoint.replace("/",'-')
    incar.write(work_dir+os.sep+"INCAR")

    # create potcar here 
    for klength in np.arange(1,11)*step:
        create_kpoints(klength,work_dir=address_kpoint)
        runtime = execute(args.np,address_kpoint,executable)
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
    return 

def encut_convergence(e_threshold,step,executable):
    address_encut = os.getcwd()
    toten = []
    kmesh = []
    wf = open(address_encut+os.sep+"encut_convergence",'w')
    if not os.path.exists(address_encut+os.sep+'KPOINTS'):
        create_kpoints(10,workdir=address_encut)
    rf = open(address_encut+os.sep+'POTCAR')
    potcar = rf.read()
    rf.close()
    encut_init = round(max([float(x) for x in re.findall('ENMAX\s*=\s*([0-9.]*);',potcar)])*1.3)
    encuts = np.arange(encut_init,1000,step)
    for iencut in encuts:
        incar = pychemia.code.vasp.VaspInput()
        incar.set_encut(ENCUT=iencut)
        incar['SYSTEM'] = sys_name
        incar['EDIFF' ] = 1e-5
        incar['NWRITE'] = 2
        incar['PREC'  ] = 'Accurate'
        incar['NCORE' ] = 4
        incar.write(work_dir+os.sep+"INCAR")
        runtime = execute(args.np,address_kpoint,executable)
        rf = open(address_encut+os.sep+"OUTCAR",'r')
        data = rf.read()
        rf.close()
        toten.append(float(re.findall("TOTEN\s*=\s*([-+0-9.]*)\s*eV",data)[-1]))
        kmesh.append(re.findall("generate k-points for.*",data)[0])
        wf.write("kpoint length = %i ,kmesh = %s, TOTEN =%f \n" %(klength,kmesh[-1],toten[-1]))
        if iencut != encut_init : # this is to check if it's not doing the calculation for the first time                                                                                                  
            change = abs(toten[-2]-toten[-1])
            if change < e_threshold :
                print("VASP calculations converged with encut " % (iencut))
                wf.close()
                break
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
    parser.add_argument("-np" ,dest="np",type=int ,action="store", help="Number of MPI processes for the code",default = '1')
    subparsers = parser.add_subparsers()
    parser_kpnt = subparsers.add_parser('kpoint_convergence')
    parser_kpnt.add_argument('--Kstart',default=10)
    parser_kpnt.add_argument('--Kend',default=100)
    parser_kpnt.add_argument('--Kstep',default=10)
    parser_kpnt.add_argument('--Ethreshold',default=1e-3)
    parser_encut = subparsers.add_parser('encut_convergence')
    parser_encut.add_argument('--Estart',default=400)
    parser_encut.add_argument('--Eend',default=1200)
    parser_encut.add_argument('--Estep',default=50)
    parser_encut.add_argument('--Ethreshold',default=1e-3)
    parser.add_argument("--executable",dest="executable",type=str,action="store",help="vasp executable",default="vasp_std")
    args = parser.parse_args()
    if 'Kstart' in args :
        kpoint_convergence(e_threshold=args.Ethreshold,start=args.Kstart,end=args.Kend,step=args.Kstep,executable=args.executable,nparal=args.np)
    elif 'Estart' in args :
        encut_convergence(e_threshold=args.Ethreshold,start=args.Estart,end=args.Eend,step=args.Estep,executable=args.executable,nparal=args.np)

