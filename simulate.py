# -*- coding: utf-8 -*-
"""
To be used with something like this:
./nC.sh -python /home/ucbtepi/code/network/src/simulate.py SimpleParameterSpacePoint(150,6,2.90,4,28.74,0.5,0,50,10,10,10,100,50) 9 matlem
"""
import random
import time
import sys
import tempfile
import shutil
import subprocess
import os
import os.path

from java.lang import System, Float
from java.io import File
from java.util import Vector, ArrayList

from ucl.physiol.neuroconstruct.project import ProjectManager
from ucl.physiol.neuroconstruct.nmodleditor.processes import ProcessManager
from ucl.physiol.neuroconstruct.neuron import NeuronFileManager

from utils.pure import SimpleParameterSpacePoint
from utils.network import generate_nC_network, generate_nC_saves, generate_nC_stimuli

point = eval(sys.argv[1].replace('+', ','))
stim_pattern_number = int(sys.argv[2])
system = sys.argv[3]

scripts_path = '/home/ucbtepi/code/network/src/scripts/'
project_path = '/home/ucbtepi/nC_projects/if_gl/'
project_filename = 'if_gl.ncx' # neuroConstruct project file name
if system == 'matlem'
    # we're on matlem
    sim_path = '/scratch0/ucbtepi/'+os.environ['JOB_ID']+'.'+os.environ['SGE_TASK_ID']
    temp_dir = sim_path
    os.makedirs(temp_dir)
elif system == 'legion':
    # we're on legion
    sim_path = os.environ['TMPDIR']
    temp_dir = sim_path
else:
    # this normally means we're on crugiat
    sim_path = '/home/ucbtepi/nC_projects/if_gl/simulations'
    temp_dir = tempfile.mkdtemp(dir=sim_path)
sim_config_name = 'Default Simulation Configuration'
nC_seed = 1234

# copy nC project in temp folder
shutil.copy2(project_path+project_filename, temp_dir+"/"+project_filename)
shutil.copytree(project_path+"morphologies", temp_dir+"/morphologies")
shutil.copytree(project_path+"cellMechanisms", temp_dir+"/cellMechanisms")

# load project and initialise
project_file = File(temp_dir + "/" + project_filename)
print ('Loading project from file: ' + project_file.getAbsolutePath() + ", exists: " + str(project_file.exists()))
pm = ProjectManager(None, None)
project = pm.loadProject(project_file)
print 'Loaded project: ' + project.getProjectName()

# load existing simulation configurations and set sim duration and a
# couple of neuron options
sim_config = project.simConfigInfo.getSimConfig(sim_config_name)
sim_config.setSimDuration(point.sim_duration)
project.neuronSettings.setNoConsole()
project.neuronSettings.setCopySimFiles(1)

# generate network
generate_nC_network(point, pm, project, sim_config)

# generate saves
generate_nC_saves(point, project)

# load stimulation patterns
print(point.stim_pattern_filename)
spf = open(point.stim_pattern_filename, "r")
stim_patterns = [[int(mf) for mf in line.split(' ')[0:-1]] for line in spf.readlines()]
spf.close()

# select the pattern I have been assigned to simulate
stim_pattern = stim_patterns[stim_pattern_number]

refs_list = [] # used to keep track of the last simulation that is run

# main loop
for trial in range(point.n_trials):
    simulation_trial_is_successful = False
    while not simulation_trial_is_successful:
        simulator_seed = random.getrandbits(32)
        # innermost loop: determine the simulation reference name
        sim_ref = point.get_simulation_reference(stim_pattern_number, trial)
        refs_list.append(sim_ref)
        project.simulationParameters.setReference(sim_ref)

        # generate stimuli in nC
        generate_nC_stimuli(point, project, sim_config, stim_pattern)

        # generate and compile neuron files
        print "Generating NEURON scripts..."
        project.neuronFileManager.generateTheNeuronFiles(sim_config, None, NeuronFileManager.RUN_HOC,simulator_seed)
        compile_process = ProcessManager(project.neuronFileManager.getMainHocFile())
        compile_success = compile_process.compileFileWithNeuron(0,0)
        # simulate
        if compile_success:
            print "...success."
            print "Simulating: simulation reference %s" % str(sim_ref)
            pm.doRunNeuron(sim_config)
            timefile_path = temp_dir+"/simulations/"+sim_ref+"/time.dat"
            print "NEURON should have started the simulation. Waiting for timefile to appear at "+timefile_path
            while os.path.exists(timefile_path)==0:
                time.sleep(2)

        hdf5_file_name = temp_dir+'/simulations/'+sim_ref+'/'+sim_ref+'_.h5'

        # wait for NEURON to finish writing the results on disk
        while not os.path.exists(hdf5_file_name):
            print ("Waiting for NEURON to finish writing to disk...")
            time.sleep(1)

        # test whether the output data archive has been created
        # properly (sometimes it's just an empty archive, with no
        # groups or datasets!)
        simulation_trial_is_successful = not bool(subprocess.call("python "+scripts_path+'test_trial_output_data.py '+hdf5_file_name, shell=True))

        if simulation_trial_is_successful:
            # copy results to main data folder
            print "Moving "+ hdf5_file_name + " to " + point.data_folder_path
            shutil.move(hdf5_file_name, point.data_folder_path)
            # perform an additional check to see if the
            # simulations results file has been successfully
            # copied
            simulation_trial_is_successful = not bool(subprocess.call("python "+scripts_path+'test_trial_output_data.py '+point.data_folder_path+'/'+sim_ref+'_.h5', shell=True))
            if not simulation_trial_is_successful:
                print("WARNING: archive" + hdf5_file_name + "was not successfuly copied to " + point.data_folder_path + "! restarting simulation.")
        else:
            print ("WARNING: archive" + hdf5_file_name + "was not found to contain simulation results! restarting simulation.")

        # delete useless files left over by neuroConstruct
        try:
            shutil.rmtree(temp_dir+'/simulations/'+sim_ref)
        except OSError:
            print('Unable to delete %s' % temp_dir+'/simulations/'+sim_ref)

# remove temp directory and exit
if system is not 'legion':
    print('removing job-specific temporary directory %s' % temp_dir)
    delete_attempts = 0
    while delete_attempts < 10:
        try:
            shutil.rmtree(temp_dir)
            print('temporary directory removed')
            break
        except OSError:
            delete_attempts += 1
            print('waiting and retrying to remove temporary directory')
            time.sleep(20)

print("done simulating pattern number " + str(stim_pattern_number) + ". closing job.")
System.exit(0)
