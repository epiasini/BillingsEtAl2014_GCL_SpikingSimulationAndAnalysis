# -*- coding: utf-8 -*-
"""
To be used with something like this:
./nC.sh -python /home/ucbtepi/code/network/src/simulate.py SimpleParameterSpacePoint(4,0,0,0.5,1,0.3,0,1,0,80,0,10,0,128,50,200) 0
"""
import sys
import subprocess
import shutil

import java.lang
import java.io

from ucl.physiol.neuroconstruct.project import ProjectManager
from ucl.physiol.neuroconstruct.nmodleditor.processes import ProcessManager

sys.path.append('/home/ucbtepi/code/network/src')
sys.path.append('/home/ucbtepi/code/networkx')

from utils.cluster_system import ClusterSystem
from utils.pure import SimpleParameterSpacePoint
from utils.network import generate_nC_network, generate_nC_saves, generate_nC_stimuli, set_tonic_GABA, scale_excitatory_conductances, export_to_neuroml

point = SimpleParameterSpacePoint(4,0,0,0.5,1,0.3,0,1,0,80,0,10,0,128,50,200)
stim_pattern_number = 0
trial = 0

print("Exporting static NeuroML for parameter space point:\n" + str(point))

scripts_path = '/home/ucbtepi/code/network/src/scripts/'
project_path = '/home/ucbtepi/nC_projects/if_gl/'
project_filename = 'if_gl.ncx' # neuroConstruct project file name
sim_config_name = 'Default Simulation Configuration'
nC_seed = 1234

with ClusterSystem() as system:
    work_dir = system.work_dir
    out_dir = "/tmp/generatedNeuroML2"
    # copy nC project in temp folder
    shutil.copy2(project_path+project_filename, work_dir+"/"+project_filename)
    shutil.copytree(project_path+"morphologies", work_dir+"/morphologies")
    shutil.copytree(project_path+"cellMechanisms", work_dir+"/cellMechanisms")

    # set level of inhibition by modifying the GrC model. If
    # inh_cond_scaling=0, then the total amount of inhibitory
    # conductance is independent of the number of dendrites. If it is
    # not zero, it's proportional to the number of dendrites divided
    # by 4, with inh_cond_scaling as the proportionality constant. The
    # DTA parameter sets how the total conductance scales with p(MF)
    gGABA_base = 438.
    dta_factor = 1 + (point.active_mf_fraction * point.dta)
    if point.inh_cond_scaling:
        inh_scaling_factor = point.inh_cond_scaling * (float(point.n_grc_dend)/4)
    else:
        inh_scaling_factor = 1
    total_GABA_conductance_in_pS = gGABA_base * point.gaba_scale * dta_factor * inh_scaling_factor
    total_GABA_conductance_in_nS = total_GABA_conductance_in_pS/1000.
    set_tonic_GABA(work_dir, total_GABA_conductance_in_nS)

    # scale excitatory conductance by modifying the AMPA and NMDA
    # synaptic models. In general, the amplitude of the excitatory
    # conductances is inversely proportional to the number of
    # dendrites divided by 4, with the exc_cond_scaling parameter as
    # the proportionality constant. If this is set to zero, the
    # conductances are not scaled.
    if point.exc_cond_scaling:
        exc_scaling_factor = point.exc_cond_scaling / (float(point.n_grc_dend)/4)
        scale_excitatory_conductances(work_dir, exc_scaling_factor)

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

    # select the pattern I have been assigned to export
    stim_pattern = stim_patterns[stim_pattern_number]

    # determine the simulation reference name
    sim_ref = point.get_simulation_reference(stim_pattern_number, trial)
    project.simulationParameters.setReference(sim_ref)

    # generate stimuli in nC
    generate_nC_stimuli(point, project, sim_config, stim_pattern)

    # generate NeuroML2 files
    export_to_neuroml(project, sim_config, out_dir)

print("done exporting network")
java.lang.System.exit(0)
