#!/usr/bin/env python
from jasp import *
from ddec_exceptions import *
import sys
import os
import shutil
from os.path import isfile, isdir
from myvasp import check_calc_complete
from subprocess import Popen, PIPE

CWD = os.getcwd()

def check_installation():
    """Checks if the module and DDEC are installed correctly"""
    try:
        ddec_path = os.environ['DDEC_PATH']
    except KeyError:
        raise DDECInitialize('Make sure your DDEC_PATH is defined correctly')
    if (not isfile(ddec_path + '/chargemol_job.m') 
        and not isdir(ddec_path + 'chargemol_10_10_2012')):
        raise DDECInitialize('Please make sure the DDEC program is installed')
    return True

def get_charges(self):
    check_installation()

    # If the job has never been started, we need to submit a calculation
    if not isfile('chargemol_output.txt'):
        ddec_run(self)

    # If the job has been started but not running, in the queue, or completed, we need to
    # clean up the files and restart the calculation
    elif (not check_ddec_complete()
          and not ddec_job_in_queue()):
        ddec_cleanup()
        ddec_run(self)
    
    # If the job is done, we should just read it and remove the running file
    elif check_ddec_complete():
        try:
            os.remove('ddec-running')
        except OSError:
            pass
        return read_output(self.atoms)
    
    # If the job is running, we should pass
    elif (not check_ddec_complete() 
          and ddec_job_in_queue()
          and os.path.exists('ddec-running')):
        raise DDECRunning()

    else:
        raise DDECUnknownState, 'I do not recognize the state of this directory'

Vasp.get_charges = get_charges

def ddec_run(calc):
    check_vasp_complete(calc)
    ddec_path = os.environ['DDEC_PATH']
    shutil.copy(ddec_path + 'chargemol_job.m', CWD)
    script = '''#!/bin/bash
cd $PBS_O_WORKDIR
matlab -nodesktop -r chargemol_job > chargemol_output.txt
'''
    resources = '-l walltime=168:00:00'
    p = Popen(['qsub', '-joe', '-N', calc.vaspdir + '-DDEC', resources],
              stdin=PIPE, stdout=PIPE, stderr=PIPE)
    out, err = p.communicate(script)
    f = open('ddec-running', 'w')
    f.close
    f = open('ddec-jobid', 'w')
    f.write(out)
    f.close
    raise DDECSubmitted(out)
    return
    
def read_output(atoms):
    '''The output we are interested in is reading the xyz file with all of the charges'''
    f = open('DDEC_net_atomic_charges.xyz', 'r')
    lines = f.readlines()
    f.close()
    n = 0
    for line in lines:
        if line.startswith('The following XYZ'):
            charges = []
            for i in range(len(atoms)):
                charges.append(float(lines[n+2+i].split()[4]))
        n += 1
    return np.array(charges)[atoms.calc.resort]

def ddec_job_in_queue():
    ''' return True or False if the directory has a job in the queue'''
    if not os.path.exists('ddec-jobid'):
        return False
    else:
        # get the jobid
        jobid = open('ddec-jobid').readline().strip()

        # see if jobid is in queue
        jobids_in_queue = commands.getoutput('qselect').split('\n')
        if jobid in jobids_in_queue:
            # get details on specific jobid
            status, output = commands.getstatusoutput('qstat %s' % jobid)
            if status == 0:
                lines = output.split('\n')
                fields = lines[2].split()
                job_status = fields[4]
                if job_status == 'C':
                    return False
                else:
                    return True
        else:
            return False

def check_ddec_complete():
    '''This function checks whether a calculation needs to be done
    
    One relevant files in a DDEC charge calculation is the chargemol_output.txt.
    If the calculation is complete, it will have a line that says 
    'Exiting Chargemol'. We will open this file and search for this line'''
    if isfile('chargemol_output.txt') == False:
        return False
    for line in open('chargemol_output.txt', 'r'):
        if line.startswith('Exiting'):
            return True
    return False

def check_vasp_complete(calc):
    '''This function checks whether a prior VASP calculation has been performed correctly'''
    # First check if the initial settings of the VASP calculation were
    # appropriate for a DDEC calculation
    if not ((calc.int_params['nsw'] == 0 or calc.int_params['nsw'] == None)
            and calc.bool_params['lcharg'] == True
            and calc.bool_params['laechg'] == True
            and calc.string_params['prec'] == 'Accurate'):
        raise DDECInitialize("""Please assure VASP calculation is correctly initalized with
nsw = 0
lcharg = True
laechg = True
prec = 'Accurate'""")
    # Now check to see if the calculation done and the right files are there
    if not (check_calc_complete()
            and isfile('POTCAR')
            and isfile('AECCAR0')
            and isfile('AECCAR2')
            and isfile('CHG')):
        raise DDECInitialize("""Please assure VASP calculation has completed all the files are present
POTCAR
AECCAR0
AECCAR2
CHG""")
    return
            
def ddec_cleanup():
    '''This function cleans the directory of the files needed to run the job'''
    files = ('add_missing_core_density.m',
             'atomic_number_to_symbol.m',
             'atomic_symbol_to_number.m',
             'calculate_inverse_Xi.m',
             'calculate_theta_scalar.m',
             'calculate_theta_vector.m',
             'calculate_Xi.m',
             'calculate_Xi_derivative.m',
             'charge_center_positions.m',
             'chargemol.m',
             'check_grid_spacing.m',
             'check_noncollinear_XC_functional.m',
             'cloud_penetration.m',
             'collinear_spin_moments_iterator.m',
             'compute_center_of_mass.m',
             'compute_dominant_atom_volumes.m',
             'constants.m',
             'core_iterator.m',
             'fast_calculate_inverse_Xi.m',
             'fast_calculate_Xi.m',
             'format_gaussian_cube_densities.m',
             'format_nongaussian_cube_densities.m',
             'format_total_cube_density.m',
             'format_vasp_densities.m',
             'format_xsf_densities.m',
             'gaussian_integration.m',
             'gaussian_overlap_integral.m',
             'gaussian_value.m',
             'generate_density_grids_from_gaussian_basis_set_coefficients.m',
             'generate_spin_lookup_tables.m',
             'generate_spin_magnetic_moment_file.m',
             'generate_xyzfile.m',
             'initialize_atomic_densities.m',
             'local_multipole_moment_analysis.m',
             'new_update_atomic_densities.m',
             'new_valence_iterator.m',
             'noble_gas_core.m',
             'noncollinear_spin_moments_iterator.m',
             'oxidation_density.m',
             'parallelpiped.m',
             'print_atomic_densities.m',
             'quote.m',
             'read_wfx_file.m',
             'run_valence_core_densities.m',
             'standard_atomic_weights.m',
             'total_multipole_moment_analysis.m',
             'file_cleanup.m')
    for f in files:
        try:
            os.remove(f)
        except OSError:
            pass
