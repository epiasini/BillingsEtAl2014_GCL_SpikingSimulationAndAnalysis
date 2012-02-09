# -*- coding: utf-8 -*-
"""
To be used with something like this:
./nC.sh -python ~/home/ucbtepi/code/network/trunk/simulate.py min_mf_number grc_mf_ratio n_grc_dend network_scale active_mf_fraction bias stim_rate_mu stim_rate_sigma noise_rate_mu noise_rate_sigma n_stim_patterns n_trials size rank
"""
import random
import time
import sys
import glob
import tempfile
import shutil
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

from utils.simulation import SimulationPoint, plast_correction_factor

min_mf_number = int(sys.argv[1])
grc_mf_ratio = float(sys.argv[2])
n_grc_dend = int(sys.argv[3])
network_scale = float(sys.argv[4])
active_mf_fraction = float(sys.argv[5])
bias = float(sys.argv[6])
stim_rate_mu = float(sys.argv[7])
stim_rate_sigma = float(sys.argv[8])
noise_rate_mu = float(sys.argv[9])
noise_rate_sigma = float(sys.argv[10])
n_stim_patterns = int(sys.argv[11])
n_trials = int(sys.argv[12])
size = int(sys.argv[13])
rank = int(sys.argv[14])

project_path = '/home/ucbtepi/nC_projects/if_gl/'
project_filename = 'if_gl.ncx' # neuroConstruct project file name
sim_path = '/home/ucbtepi/nC_projects/if_gl/simulations'
sim_config_name = 'Default Simulation Configuration'
nC_seed = 1234

stim_rate_mu = 120
stim_rate_sigma = 30
noise_rate_mu = 0
noise_rate_sigma = 10

sim_duration = 300.0

point = SimulationPoint(sim_duration,
                 min_mf_number,
                 grc_mf_ratio,
                 n_grc_dend,
                 network_scale,
                 active_mf_fraction,
                 bias,
                 stim_rate_mu,
                 stim_rate_sigma,
                 noise_rate_mu,
                 noise_rate_sigma,
                 n_stim_patterns,
                 n_trials)

mf_number = int(round(min_mf_number * network_scale))
gr_number = int(round(mf_number * grc_mf_ratio))

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
mf_pack_adapter.setMaxNumberCells(mf_number)
gr_pack_adapter.setMaxNumberCells(gr_number)

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
sim_config.setSimDuration(sim_duration)
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
spf = open(point.stim_pattern_filename, "r")
stim_patterns = [[int(mf) for mf in line.split(' ')[0:-1]] for line in spf.readlines()]
spf.close()

# calculate which patterns are mine to simulate
patterns_per_chunk = n_stim_patterns/size
my_stim_lower_bound = rank*patterns_per_chunk
if rank != size-1:
    my_stim_upper_bound = (rank + 1) * patterns_per_chunk
else:
    my_stim_upper_bound = n_stim_patterns

print (n_stim_patterns, patterns_per_chunk)
print (rank, my_stim_lower_bound, my_stim_upper_bound)

refs_list = [] # used to keep track of the last simulation that is run

# main loop
for spn, sp in list(enumerate(stim_patterns))[my_stim_lower_bound: my_stim_upper_bound]:
    for trial in range(n_trials):
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
        sim_ref = point.get_ref(spn, trial)
        refs_list.append(sim_ref)
        project.simulationParameters.setReference(sim_ref)

        #Set the thresholding current
        bias_in_nA = 0.001 * bias
        bias_input = IClampSettings('my_bias', 'GrCs', AllCells(), 0, 0, sim_duration, bias_in_nA, False)
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
        
        for mf in range(mf_number):
            if mf in sp:
                # create random spiking stimuli on active MFs, following the current stimulus pattern
                rate = max(0.1, random.gauss(stim_rate_mu, stim_rate_sigma))
            else:
                rate = max(0.1, random.gauss(noise_rate_mu, noise_rate_sigma))
            rate_in_khz = rate/1000.
            stim = RandomSpikeTrainSettings('MF_stim_'+str(mf), 'MFs', FixedNumberCells(0), 0, NumberGenerator(rate_in_khz), 'FastSynInput')
            project.elecInputInfo.addStim(stim)
            sim_config.addInput(stim.getReference())
            project.generatedElecInputs.addSingleInput(stim.getReference(), 'RandomSpikeTrain', 'MFs', mf, 0, 0, None)
            # netconn set-up (remember that all the connections originating from the same MF are identical)
            syn_props_list = Vector([SynapticProperties('NMDA'), SynapticProperties('AMPA_direct'), SynapticProperties('AMPA_spillover')])
            for synp in syn_props_list:
                synp.setFixedDelay(0)
                synp.setThreshold(0)
            syn_props_list[0].setWeightsGenerator(NumberGenerator(plast_correction_factor(rate, 'NMDA')))
            syn_props_list[1].setWeightsGenerator(NumberGenerator(plast_correction_factor(rate, 'AMPA')))
            syn_props_list[2].setWeightsGenerator(NumberGenerator(plast_correction_factor(rate, 'AMPA')))
            conn_conditions = ConnectivityConditions()
            conn_conditions.setNumConnsInitiatingCellGroup(NumberGenerator(0))
            project.morphNetworkConnectionsInfo.addRow('conn_'+str(mf), 'MFs', 'GrCs', syn_props_list, SearchPattern.getRandomSearchPattern(), MaxMinLength(Float.MAX_VALUE, 0, 'r', 100), conn_conditions, Float.MAX_VALUE)

        for gr in range(n_gr):
            # set up thresholding current
            project.generatedElecInputs.addSingleInput('bias', 'IClamp', 'GrCs', gr, 0, 0, None)
            for mf in conn_pattern[gr]:
                # create connections, following the current connection pattern
                conn_name = 'conn_'+str(mf)
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
            timefile_path=temp_dir+"/simulations/"+sim_ref+"/time.dat"
            while os.path.exists(timefile_path)==0:
                time.sleep(2)

        # wait for NEURON to finish writing the results on disk
        while len(glob.glob(temp_dir+'/simulations/'+sim_ref+'/*.h5')) < 2:
            print ("Waiting for NEURON to finish writing to disk...")
            time.sleep(1)
        
        # copy results to main data folder and delete useless nC files
        destination = point.data_folder_path + '/' + sim_ref
        if os.path.isdir(destination):
            shutil.rmtree(destination)
        os.mkdir(destination)
        print "Moving "+ temp_dir+"/simulations/"+sim_ref + " to " + destination
        for source in glob.glob(temp_dir+'/simulations/'+sim_ref+'/*.h5'):
            shutil.move(source, destination)
        shutil.rmtree(temp_dir+'/simulations/'+sim_ref)

shutil.rmtree(temp_dir)
System.exit(0)
