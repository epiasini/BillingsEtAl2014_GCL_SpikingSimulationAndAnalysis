# -*- coding: utf-8 -*-
import random
import time
import os.path

from java.lang import System
from java.io import File

from ucl.physiol.neuroconstruct.project import ProjectManager
from ucl.physiol.neuroconstruct.utils import NumberGenerator
from ucl.physiol.neuroconstruct.nmodleditor.processes import ProcessManager
from ucl.physiol.neuroconstruct.neuron import NeuronFileManager

# simulation control parameters
sim_path = '/home/eugenio/phd/nC_projects/if_network/'
project_filename = 'if_network.ncx'
sim_config_name = 'Default Simulation Configuration'
existing_stim = 'MF_stim'
nC_seed = 1234
simulator_seed = random.getrandbits(32)
active_mf_fraction = 0.3
sim_duration = 100.0
n_stim_patterns = 1
n_conn_patterns = 1
bias_values = [-0.005]
# timestamp
now = time.localtime(time.time())
tstr = time.strftime("%Y-%m-%d.%H-%M-%S", now)

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
                # innermost loop: determine the simulation reference name
                sim_ref = tstr + "_cp" + str(cpn) + "_sp" + str(spn) + "_b" + str(bias)
                project.simulationParameters.setReference(sim_ref)
                pattern_file=open(sim_path+"simulations/"+sim_ref+'_patterns.txt',"w")
                
                #Set the thresholding current
                bias_input = project.elecInputInfo.getStim("bias")
                bias_input.setAmp(NumberGenerator(bias))
                bias_input.setDur(NumberGenerator(sim_duration))
                #project.elecInputInfo.updateStim(bias_input)

                pm.doGenerate(sim_config_name, nC_seed)
                while pm.isGenerating():
                    print 'Waiting for the project to be regenerated...'
                    time.sleep(1)
                
                pattern_file.write(str(cp) + "\n")
                for gr in range(n_gr):
                    for mf in cp[gr]:
                        # create connections, following the current connection pattern
                        project.generatedNetworkConnections.addSynapticConnection('NetConn_MFs_GrCs', mf, gr)

                pattern_file.write(str(sp))
                for mf in sp:
                    # create random spiking stimuli on active MFs, following the current stimulus pattern
                    project.generatedElecInputs.addSingleInput(existing_stim,'RandomSpikeTrain','MFs',mf,0,0,None)

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
                pattern_file.close()
else:
    print "No cells generated."

System.exit(0)
