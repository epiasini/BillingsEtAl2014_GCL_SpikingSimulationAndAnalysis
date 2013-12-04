# -*- coding: utf-8 -*-
"""
To be used with something like this:
./nC.sh -python ~/home/ucbtepi/code/network/trunk/simulate.py SimpleParameterSpacePoint(300,6,2.00,4,5.00,0.5,-20,120,30,30,10,20,200) 9
"""
import random
import time
import sys
import tempfile
import shutil
import subprocess
import os.path

from java.lang import System, Float
from java.io import File
from java.util import Vector, ArrayList

from ucl.physiol.neuroconstruct.project import ProjectManager, CellGroupsInfo, SingleElectricalInput, SearchPattern, MaxMinLength, SynapticProperties, ConnectivityConditions
from ucl.physiol.neuroconstruct.project.stimulation import RandomSpikeTrainInstanceProps
from ucl.physiol.neuroconstruct.project.cellchoice import FixedNumberCells, AllCells
from ucl.physiol.neuroconstruct.utils import NumberGenerator
from ucl.physiol.neuroconstruct.simulation import RandomSpikeTrainSettings, IClampSettings
from ucl.physiol.neuroconstruct.nmodleditor.processes import ProcessManager
from ucl.physiol.neuroconstruct.neuron import NeuronFileManager

from utils.pure import SimpleParameterSpacePoint, plast_correction_factor

point = eval(sys.argv[1])
rank = int(sys.argv[2])

scripts_path = '/home/ucbtepi/code/network/src/scripts/'
project_path = '/home/ucbtepi/nC_projects/if_gl/'
project_filename = 'if_gl.ncx' # neuroConstruct project file name
sim_path = '/home/ucbtepi/nC_projects/if_gl/simulations'
sim_config_name = 'Default Simulation Configuration'
nC_seed = 1234

temp_dir = tempfile.mkdtemp(dir=sim_path)
shutil.copy2(project_path+project_filename, temp_dir+"/"+project_filename)
shutil.copytree(project_path+"morphologies", temp_dir+"/morphologies")
shutil.copytree(project_path+"cellMechanisms", temp_dir+"/cellMechanisms")

# load project and initialise
project_file = File(temp_dir + "/" + project_filename)
print ('Loading project from file: ' + project_file.getAbsolutePath() + ", exists: " + str(project_file.exists()))
pm = ProjectManager(None, None)
project = pm.loadProject(project_file)
print 'Loaded project: ' + project.getProjectName()

# set network size
group_info = project.cellGroupsInfo
mf_pack_adapter = group_info.getCellPackingAdapter('MFs')
gr_pack_adapter = group_info.getCellPackingAdapter('GrCs')
mf_pack_adapter.setMaxNumberCells(point.n_mf)
gr_pack_adapter.setMaxNumberCells(point.n_grc)

# generate network
i = 0
while True:
    try:
        pm.doGenerate(sim_config_name, nC_seed)
        while pm.isGenerating():        
            print 'Waiting for the project to be generated...'
            time.sleep(2)
        break
    except java.util.ConcurrentModificationException:
        i = i + 1
        print "ConcurrentModificationException raised by nC; retrying network generation."
        if i>5:
            raise java.util.ConcurrentModificationException

# load existing simulation configurations and set sim duration and a couple of neuron options
sim_config = project.simConfigInfo.getSimConfig(sim_config_name)
sim_config.setSimDuration(point.sim_duration)
project.neuronSettings.setNoConsole()
project.neuronSettings.setCopySimFiles(1)

# count generated cells
n_mf=project.generatedCellPositions.getNumberInCellGroup('MFs')
n_gr=project.generatedCellPositions.getNumberInCellGroup('GrCs')
n_generated_cells = project.generatedCellPositions.getNumberInAllCellGroups()
print "Number of cells generated: " + str(n_generated_cells)

# load connection pattern
cpf = open(point.conn_pattern_filename, "r")
print(cpf)
conn_pattern = [[int(mf) for mf in line.split(' ')[0:-1]] for line in cpf.readlines()]
cpf.close()

# load stimulation patterns
print(point.stim_pattern_filename)
spf = open(point.stim_pattern_filename, "r")
stim_patterns = [[int(mf) for mf in line.split(' ')[0:-1]] for line in spf.readlines()]
spf.close()

# calculate which patterns are mine to simulate
patterns_per_chunk = point.n_stim_patterns/point.SIZE_PER_SIMULATION
my_stim_lower_bound = rank*patterns_per_chunk
if rank != point.SIZE_PER_SIMULATION-1:
    my_stim_upper_bound = (rank + 1) * patterns_per_chunk
else:
    my_stim_upper_bound = point.n_stim_patterns

print (point.n_stim_patterns, patterns_per_chunk)
print (rank, my_stim_lower_bound, my_stim_upper_bound)

refs_list = [] # used to keep track of the last simulation that is run

# main loop
for spn, sp in list(enumerate(stim_patterns))[my_stim_lower_bound: my_stim_upper_bound]:
    for trial in range(point.n_trials):
        simulation_trial_is_successful = False
        while not simulation_trial_is_successful:
            # delete all existing connections
            project.generatedNetworkConnections.reset()
            project.morphNetworkConnectionsInfo.deleteAllNetConns()
            sim_config.setNetConns(ArrayList())        
            # delete all existing stimuli
            project.generatedElecInputs.reset()
            project.elecInputInfo.deleteAllStims()
            sim_config.setInputs(ArrayList())

            simulator_seed = random.getrandbits(32)
            # innermost loop: determine the simulation reference name
            sim_ref = point.get_simulation_reference(spn, trial)
            refs_list.append(sim_ref)
            project.simulationParameters.setReference(sim_ref)

            #Set the thresholding current
            bias_in_nA = 0.001 * point.bias
            bias_input = IClampSettings('my_bias', 'GrCs', AllCells(), 0, 0, point.sim_duration, bias_in_nA, False)
            project.elecInputInfo.addStim(bias_input)
            sim_config.addInput(bias_input.getReference())

            # regenerate network
            i = 0
            while True:
                try:
                    pm.doGenerate(sim_config_name, nC_seed)
                    while pm.isGenerating():        
                        print 'Waiting for the project to be regenerated...'
                        time.sleep(1)
                    break
                except java.util.ConcurrentModificationException:
                    i = i + 1
                    print "ConcurrentModificationException raised by nC; retrying network generation."
                    if i>5:
                        raise java.util.ConcurrentModificationException

            for mf in range(point.n_mf):
                if mf in sp:
                    # create random spiking stimuli on active MFs, following the current stimulus pattern
                    rate = max(0.1, random.gauss(point.stim_rate_mu, point.stim_rate_sigma))
                else:
                    rate = max(0.1, random.gauss(point.noise_rate_mu, point.noise_rate_sigma))
                rate_in_khz = rate/1000.
                stim = RandomSpikeTrainSettings('MF_stim_'+str(mf), 'MFs', FixedNumberCells(0), 0, NumberGenerator(rate_in_khz), 'FastSynInput')
                project.elecInputInfo.addStim(stim)
                sim_config.addInput(stim.getReference())
                project.generatedElecInputs.addSingleInput(stim.getReference(), 'RandomSpikeTrain', 'MFs', mf, 0, 0, None)
                # netconn set-up (remember that all the connections originating from the same MF are identical)
                for syn_type in ['RothmanMFToGrCAMPA', 'RothmanMFToGrCNMDA']:
                    syn_props_list = Vector([SynapticProperties(syn_type)])
                    syn_props_list[0].setFixedDelay(0)
                    syn_props_list[0].setThreshold(0)
                    syn_props_list[0].setWeightsGenerator(NumberGenerator(1))
                    conn_conditions = ConnectivityConditions()
                    conn_conditions.setNumConnsInitiatingCellGroup(NumberGenerator(0))
                    project.morphNetworkConnectionsInfo.addRow('conn_'+str(mf)+'_'+syn_type, 'MFs', 'GrCs', syn_props_list, SearchPattern.getRandomSearchPattern(), MaxMinLength(Float.MAX_VALUE, 0, 'r', 100), conn_conditions, Float.MAX_VALUE)

            for gr in range(n_gr):
                # set up thresholding current
                project.generatedElecInputs.addSingleInput('bias', 'IClamp', 'GrCs', gr, 0, 0, None)
                for mf in conn_pattern[gr]:
                    # create connections, following the current connection pattern
                    for syn_type in ['RothmanMFToGrCAMPA', 'RothmanMFToGrCNMDA']:
                        conn_name = 'conn_'+str(mf)+'_'+syn_type
                        sim_config.addNetConn(conn_name)
                        project.generatedNetworkConnections.addSynapticConnection(conn_name, mf, gr)


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

            # delete useless files left over by neuroConstruct
            try:
                shutil.rmtree(temp_dir+'/simulations/'+sim_ref)
            except OSError:
                print('Unable to delete %s' % temp_dir+'/simulations/'+sim_ref)

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
System.exit(0)
