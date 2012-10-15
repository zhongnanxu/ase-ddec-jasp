#!/usr/bin/env python
import sys
import os
import shutil
from os.path import isfile
from myvasp import check_calc_complete
from subprocess import Popen, PIPE

def calculate(time='short', ppn=1, name=os.getcwd()):
    assert isfile('running') == False
    assert check_calc_complete() == True
    assert check_ddec_complete() == False
    CWD = os.getcwd()
    shutil.copy('/home/zhongnanxu/lib/python2.6/site-packages/myddec/chargemol_job.m', CWD)
    script = '''#!/bin/bash
cd $PBS_O_WORKDIR
matlab -nodesktop -r chargemol_job > chargemol_output.txt
'''
    if time == 'short':
        hours = 24
    else:
        hours = 168
    resources = '-l walltime=%i:00:00,nodes=1:ppn=%i' % (hours,ppn)
    p = Popen(['qsub', '-joe', '-N', '%s' % name, resources],
              stdin=PIPE, stdout=PIPE, stderr=PIPE)
    out, err = p.communicate(script)
    return

def check_ddec_complete():
    '''This function checks whether a calculation needs to be done
    
    The relevant files in a DDEC charge calculation is the chargemol_output.txt.
    If the calculatio is complete, it will have a line that says 
    'Exiting Chargemol'. We will open this file and search for this line'''
    if isfile('chargemol_output.txt') == False:
        return False
    for line in open('chargemol_output.txt', 'r'):
        if line.startswith('Exiting'):
            return True
    return False
