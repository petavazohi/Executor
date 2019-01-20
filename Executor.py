#!/usr/bin/env python3

import os
import subprocess
import argparse
import time
import re
import numpy as np 
import pychemia

def create_kpoints(length):
    length = str(length)+'\n'
    comment = 'Automatic mesh'+'\n'
    nkp     = str(0)+'\n'
    method  = 'Auto'+'\n'
    wf = open("KPOINTS",'w')
    wf.write(comment)
    wf.write(nkp    )
    wf.write(method )
    wf.write(length )
    wf.close()
    return

def kpoint_convergence(e_threshold,start,end,step,executable,nparal):
    address = os.getcwd() 
    toten = []
    kmesh = []
    if not os.path.exists('INCAR'):
        incar = pychemia.code.vasp.VaspInput()
        incar.set_encut(1.4,POTCAR='POTCAR')
        incar['EDIFF' ] = 1e-4
        incar['NWRITE'] = 2
        incar['PREC'  ] = 'Accurate'
        incar['NCORE' ] = 4
        incar['SYSTEM'] = '-'.join(address.split('/')[-2:])
        incar.write("INCAR")
    # create potcar here 
    klengths = np.arange(start,end,step)
    for klength in klengths:
        create_kpoints(klength)
        print("===================================================================================")
        print("===================================================================================")
        print("===================================================================================")
        runtime = execute(nparal,address,executable)
        rf = open("OUTCAR",'r')
        data = rf.read()
        rf.close()
        toten.append(float(re.findall("TOTEN\s*=\s*([-+0-9.]*)\s*eV",data)[-1]))
        kmesh.append(re.findall("generate k-points for.*",data)[0])
        wf = open("kpoint_convergence",'a')        
        wf.write("kpoint length = %i ,kmesh = %s, TOTEN =%f \n" %(klength,kmesh[-1],toten[-1]))
        wf.close()
        if klength != start : # this is to check if it's not doing the calculation for the first time
            change = abs(toten[-2]-toten[-1])
            if change < e_threshold : # need to add more comparision here
                print("VASP calculations converged with k points length %i and kmesh %s " % (klength,kmesh[-1]))
                break
    print("VASP calculations converged with k points length %i and kmesh %s " % (klengths[conv_idx],kmesh[conv_idx]))
    return 

def encut_convergence(e_threshold,start,end,step,executable,nparal):
    address = os.getcwd()
    toten = []
    if not os.path.exists('KPOINTS'):
        create_kpoints(10)
    rf = open('POTCAR')
    potcar = rf.read()
    rf.close()
    encut_init = round(max([float(x) for x in re.findall('ENMAX\s*=\s*([0-9.]*);',potcar)])*1.3)
    if start < encut_init :
        print('Initial value provided for ENUCT is less than 1.3*ENMAX pseudo potential, replacing Estart with %f'% encut_init)
        start = round(encut_init,-2)
    encuts = np.arange(start,end,step)
    for iencut in encuts:
        incar = pychemia.code.vasp.VaspInput()
        incar['ENCUT']  = iencut
        incar['SYSTEM'] = '-'.join(address.split('/')[-2:])
        incar['EDIFF' ] = 1e-5
        incar['NWRITE'] = 2
        incar['PREC'  ] = 'Accurate'
        incar['NCORE' ] = 4
        incar.write("INCAR")
        print("===================================================================================")
        print("===================================================================================")
        print("===================================================================================")
        print('Running vasp with ENCUT = {}'.format(iencut))
        runtime = execute(args.np,address,executable)
        rf = open("OUTCAR",'r')
        data = rf.read()
        rf.close()
        toten.append(float(re.findall("TOTEN\s*=\s*([-+0-9.]*)\s*eV",data)[-1]))
        wf = open("encut_convergence",'a')
        wf.write("encut = %i , TOTEN =%f \n" %(iencut,toten[-1]))
        wf.close()
    toten = np.array(toten)
    encut_idx = toten - toten[-1] < e_threshold
    best_encut = encuts[encut_idx][0]
    wf = open('best_encut','w')
    wf.write('best_cut = {}'.format(best_encut))
    wf.close()
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
    parser_kpnt.add_argument('--Kstart',type=float,default=10)
    parser_kpnt.add_argument('--Kend',type=float,default=100)
    parser_kpnt.add_argument('--Kstep',type=float,default=10)
    parser_kpnt.add_argument('--Ethreshold',type=float,default=1e-3)
    parser_encut = subparsers.add_parser('encut_convergence')
    parser_encut.add_argument('--Estart',type=float,default=100)
    parser_encut.add_argument('--Eend',type=float,default=900)
    parser_encut.add_argument('--Estep',type=float,default=50)
    parser_encut.add_argument('--Ethreshold',type=float,default=1e-3)
    parser.add_argument("--executable",dest="executable",type=str,action="store",help="vasp executable",default="vasp_std")
    args = parser.parse_args()
    if 'Kstart' in args :
        kpoint_convergence(e_threshold=args.Ethreshold,start=args.Kstart,end=args.Kend,step=args.Kstep,executable=args.executable,nparal=args.np)
    elif 'Estart' in args :
        encut_convergence(e_threshold=args.Ethreshold,start=args.Estart,end=args.Eend,step=args.Estep,executable=args.executable,nparal=args.np)

