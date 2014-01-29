# -*- coding: utf-8 -*-
"""
To be used with something like this:
./nC.sh -python /home/ucbtepi/code/network/src/simulate.py SimpleParameterSpacePoint(4+0+0+0.5+1+0+0+80+0+10+0+128+50+200) 0
"""
import random
import time
import sys
import tarfile
import shutil
import subprocess
import os
import os.path

import java.lang
import java.io

from ucl.physiol.neuroconstruct.project import ProjectManager
from ucl.physiol.neuroconstruct.nmodleditor.processes import ProcessManager
from ucl.physiol.neuroconstruct.neuron import NeuronFileManager

from utils.cluster_system import ClusterSystem
from utils.pure import SimpleParameterSpacePoint
from utils.network import generate_nC_network, generate_nC_saves, generate_nC_stimuli, set_tonic_GABA

point = eval(sys.argv[1].replace('+', ','))
stim_pattern_number = int(sys.argv[2])

# first of all, check if a temporary tar archive is already available
# for the results of the simulation of this stimulus pattern. If this
# is the case there is no need for simulating everything again, so
# just print a message and exit.
if os.path.isfile(point.get_tar_simulation_archive_path(stim_pattern_number)):
    print("WARNING: tar spikes archive found on disk at path " + point.get_tar_simulation_archive_path(stim_pattern_number) + ".\nNot running simulations for stim pattern number " + str(stim_pattern_number) + ". closing job.")
    java.lang.System.exit(0)

print("Simulating point:\n" + str(point) + "\nSimulation results will be stored in " + point.get_tar_simulation_archive_path(stim_pattern_number))

scripts_path = '/home/ucbtepi/code/network/src/scripts/'
project_path = '/home/ucbtepi/nC_projects/if_gl/'
project_filename = 'if_gl.ncx' # neuroConstruct project file name
sim_config_name = 'Default Simulation Configuration'
nC_seed = 1234

with ClusterSystem() as system:
    work_dir = system.work_dir
    # copy nC project in temp folder
    shutil.copy2(project_path+project_filename, work_dir+"/"+project_filename)
    shutil.copytree(project_path+"morphologies", work_dir+"/morphologies")
    shutil.copytree(project_path+"cellMechanisms", work_dir+"/cellMechanisms")

    # prepare temporary tar archive for simulation data
    tar_archive_path = work_dir+'/temporary_archive.tar'
    tar_archive = tarfile.open(tar_archive_path, 'w')

    # set level of inhibition by modifying GrC model. This is somewhat
    # involved, to follow what was done with the binary model. If
    # DTA=0, then the total amount of inhibitory conductance is
    # independent of th number of dendrites, and it's given by the
    # default value of gGABA_base multiplied by the gaba_scale
    # factor. If DTA is not zero, then the total amount of conductance
    # is gaba_scale*gGABA_base*(d/4)(1+DTA*p(MF)). Note that this
    # means that when DTA is not zero the total GABA conductance
    # scales with the number of dendrites, and conversely it is not
    # possible - strictly speaking - to have the GABA conductance
    # depend on the number of dendrites but not on p(MF).
    gGABA_base = 438.
    if not point.dta:
        total_GABA_conductance_in_pS = gGABA_base * point.gaba_scale
    else:
        total_GABA_conductance_in_pS = gGABA_base * point.gaba_scale * (float(point.n_grc_dend)/4) * (1 + (point.active_mf_fraction * point.dta))
    total_GABA_conductance_in_nS = total_GABA_conductance_in_pS/1000.
    set_tonic_GABA(work_dir, total_GABA_conductance_in_nS)

    # load project and initialise
    project_file = java.io.File(work_dir + "/" + project_filename)
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
        simulation_attempt = 0
        simulation_trial_is_successful = False
        while simulation_attempt < 3 and not simulation_trial_is_successful:
            simulation_attempt += 1
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
                timefile_path = work_dir+"/simulations/"+sim_ref+"/time.dat"
                print "NEURON should have started the simulation. Waiting for timefile to appear at "+timefile_path
                while os.path.exists(timefile_path)==0:
                    time.sleep(2)

            hdf5_file_name = work_dir+'/simulations/'+sim_ref+'/'+sim_ref+'_.h5'

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
                print "Archiving "+ hdf5_file_name + " to temporary tar file"
                tar_archive.add(hdf5_file_name, arcname=sim_ref+'_.h5')
            else:
                read_results_attempt = 0 
                while read_results_attempt < 10 and not simulation_trial_is_successful:
                    read_results_attempt += 1
                    print("WARNING: simulation results not found. Waiting and retrying: attempt "+str(read_results_attempt))
                    time.sleep(60)
                    simulation_trial_is_successful = not bool(subprocess.call("python "+scripts_path+'test_trial_output_data.py '+hdf5_file_name, shell=True))
                if not simulation_trial_is_successful:
                    print ("WARNING: archive " + hdf5_file_name + " was not found to contain simulation results! restarting simulation.")

            # delete useless files left over by neuroConstruct
            try:
                shutil.rmtree(work_dir+'/simulations/'+sim_ref)
            except OSError:
                print('Unable to delete %s' % work_dir+'/simulations/'+sim_ref)

        if not simulation_trial_is_successful:
            raise Exception("Simulation of trial {} was attempted several times, but it always failed. Something must be wrong; giving up.".format(trial))

    # move tar archive to main data directory
    tar_archive.close()
    shutil.move(tar_archive_path,
                point.get_tar_simulation_archive_path(stim_pattern_number))

print("done simulating pattern number " + str(stim_pattern_number) + ". closing job.")
java.lang.System.exit(0)
