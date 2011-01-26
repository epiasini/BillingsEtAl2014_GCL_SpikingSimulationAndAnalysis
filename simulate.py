# -*- coding: utf-8 -*-
"""
To be used with something like this:
./nC.sh -python ~/data/eugenio/network/trunk/simulate.py base_name
"""
import random
import time
import sys
import os.path

from java.lang import System
from java.io import File

from ucl.physiol.neuroconstruct.project import ProjectManager
from ucl.physiol.neuroconstruct.utils import NumberGenerator
from ucl.physiol.neuroconstruct.nmodleditor.processes import ProcessManager
from ucl.physiol.neuroconstruct.neuron import NeuronFileManager

# simulation control parameters
base_name = sys.argv[1]
ntrials = 5
sim_path = '/home/ucgbgbi/data/eugenio/nC_projects/if_gl/'
project_filename = 'if_gl.ncx'
sim_config_name = 'Default Simulation Configuration'
existing_stim = 'MF_stim'
nC_seed = 1234
active_mf_fraction = 0.3
sim_duration = 1000.0
n_stim_patterns = 1
n_conn_patterns = 1
bias_values = [-0.005]

def simulate(base_name):
    # load project and initialise
    project_file = File(sim_path + project_filename)
    print ('Loading project from file: ' + project_file.getAbsolutePath() + ", exists: " + str(project_file.exists()))
    pm = ProjectManager(None, None)
    project = pm.loadProject(project_file)
    print 'Loaded project: ' + project.getProjectName()

    # generate network
    pm.doGenerate(sim_config_name, nC_seed)
    while pm.isGenerating():
        print 'Waiting for the project to be generated...'
        time.sleep(2)

    # load existing simulation configurations and set sim duration and a couple of neuron options
    sim_config = project.simConfigInfo.getSimConfig(sim_config_name)
    sim_config.setSimDuration(sim_duration)
    project.neuronSettings.setNoConsole()
    project.neuronSettings.setCopySimFiles(1)

    # count generated cells
    n_mf=project.generatedCellPositions.getNumberInCellGroup('MFs')
    n_gr=project.generatedCellPositions.getNumberInCellGroup('GrCs')
    n_generated_cells = project.generatedCellPositions.getNumberInAllCellGroups()
    print "Number of cells generated: " + str(n_generated_cells)

    # load connection patterns
    conn_patterns = [[random.sample(range(n_mf),4) for each in range(n_gr)] for each in range(n_conn_patterns)]
    #inv_connection_patterns = [[gr for gr in range(n_gr) if mf in connection_patterns[gr]] for mf in range(n_mf)]

    # generate random stimulation pattern (at fixed sparsity)
    active_mf_number = int(round(n_mf*active_mf_fraction))
    stim_patterns = [random.sample(range(n_mf), active_mf_number) for each in range (n_stim_patterns)]

    # simulate and record spiketimes
    if n_generated_cells > 0:
        for cpn, cp in enumerate(conn_patterns):
            # delete all existing connections
            project.generatedNetworkConnections.reset()

            for spn, sp in enumerate(stim_patterns):
                # delete all existing stimuli
                project.generatedElecInputs.reset()

                for bias in bias_values:
                    for trial in range(ntrials):
                        simulator_seed = random.getrandbits(32)
                        # innermost loop: determine the simulation reference name
                        sim_ref = base_name + "_cp" + str(cpn) + "_sp" + str(spn) + "_b" + str(bias) + "_t" + str(trials)
                        project.simulationParameters.setReference(sim_ref)

                        #Set the thresholding current
                        bias_input = project.elecInputInfo.getStim("bias")
                        bias_input.setAmp(NumberGenerator(bias))
                        bias_input.setDur(NumberGenerator(sim_duration))
                        #project.elecInputInfo.updateStim(bias_input)

                        # regenerate network
                        pm.doGenerate(sim_config_name, nC_seed)
                        while pm.isGenerating():
                            print 'Waiting for the project to be regenerated...'
                            time.sleep(1)

                        conn_pattern_file=open(sim_path+"simulations/"+sim_ref+'_conn.txt',"w")
                        for gr in range(n_gr):
                            for mf in cp[gr]:
                                # create connections, following the current connection pattern
                                project.generatedNetworkConnections.addSynapticConnection('NetConn_MFs_GrCs', mf, gr)
                                conn_pattern_file.write(str(mf) + " ")
                            conn_pattern_file.write("\n")
                        conn_pattern_file.close()

                        stim_pattern_file=open(sim_path+"simulations/"+sim_ref+'_stim.txt',"w")
                        for mf in sp:
                            # create random spiking stimuli on active MFs, following the current stimulus pattern
                            project.generatedElecInputs.addSingleInput(existing_stim,'RandomSpikeTrain','MFs',mf,0,0,None)
                            stim_pattern_file.write(str(mf) + " ")
                        stim_pattern_file.close()

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
                            propsfile_path=sim_path+"simulations/"+sim_ref+"/simulation.props"
                            while os.path.exists(propsfile_path)==0:
                                time.sleep(2)
                            print("")
    else:
        print "No cells generated."

    System.exit(0)

# the "main"
simulate(base_name)
