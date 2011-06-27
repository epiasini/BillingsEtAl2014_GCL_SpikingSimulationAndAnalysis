# -*- coding: utf-8 -*-
"""
To be used with something like this:
./nC.sh -python ~/data/eugenio/network/trunk/simulate.py min_mf_number grc_mf_ratio n_grc_dend network_scale active_mf_fraction bias n_stim_patterns n_trials size rank
"""
import random
import time
import sys
import glob
import tempfile
import shutil
import os.path

from java.lang import System
from java.io import File

from ucl.physiol.neuroconstruct.project import ProjectManager, CellGroupsInfo
from ucl.physiol.neuroconstruct.utils import NumberGenerator
from ucl.physiol.neuroconstruct.nmodleditor.processes import ProcessManager
from ucl.physiol.neuroconstruct.neuron import NeuronFileManager

from utils import ref_ctor, conn_pattern_filename, stim_pattern_filename, data_folder_path_ctor

min_mf_number = int(sys.argv[1])
grc_mf_ratio = float(sys.argv[2])
n_grc_dend = int(sys.argv[3])
network_scale = float(sys.argv[4])
active_mf_fraction = float(sys.argv[5])
bias = float(sys.argv[6])
n_stim_patterns = int(sys.argv[7])
n_trials = int(sys.argv[8])
size = int(sys.argv[9])
rank = int(sys.argv[10])

project_path = '/home/ucbtepi/nC_projects/if_gl/'
project_filename = 'if_gl.ncx' # neuroConstruct project file name
sim_path = '/home/ucbtepi/nC_projects/if_gl/simulations'
sim_config_name = 'Default Simulation Configuration'
nC_seed = 1234

sim_duration = 300.0

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
cpf = open(conn_pattern_filename(grc_mf_ratio, n_grc_dend, network_scale), "r")
print(cpf)
conn_pattern = [[int(mf) for mf in line.split(' ')[0:-1]] for line in cpf.readlines()]
cpf.close()

# load stimulation patterns
spf = open(stim_pattern_filename(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, n_stim_patterns), "r")
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

# delete all existing connections
project.generatedNetworkConnections.reset()

# main loop
for spn, sp in list(enumerate(stim_patterns))[my_stim_lower_bound: my_stim_upper_bound]:
    # delete all existing stimuli
    project.generatedElecInputs.reset()

    for trial in range(n_trials):
        simulator_seed = random.getrandbits(32)
        # innermost loop: determine the simulation reference name
        sim_ref = ref_ctor(n_stim_patterns, n_trials, spn, trial)
        refs_list.append(sim_ref)
        project.simulationParameters.setReference(sim_ref)

        #Set the thresholding current
        bias_in_nA = 0.001 * bias
        bias_input = project.elecInputInfo.getStim("bias")
        bias_input.setAmp(NumberGenerator(bias_in_nA))
        bias_input.setDur(NumberGenerator(sim_duration))
        project.elecInputInfo.updateStim(bias_input)

        # regenerate network
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
                
        for gr in range(n_gr):
            # set up thresholding current
            project.generatedElecInputs.addSingleInput('bias', 'IClamp', 'GrCs', gr, 0, 0, None)
            for mf in conn_pattern[gr]:
                # create connections, following the current connection pattern
                project.generatedNetworkConnections.addSynapticConnection('NetConn_MFs_GrCs', mf, gr)

        for mf in sp:
            # create random spiking stimuli on active MFs, following the current stimulus pattern
            project.generatedElecInputs.addSingleInput("MF_stim",'RandomSpikeTrain','MFs',mf,0,0,None)

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

# the writing of the output files to disk can take a while after the simulations have finished, so this is to avoid the script exiting before all the results have been saved. This while loop puts the process to sleep until the results of the last simulation begin to be written to disk.
while len(glob.glob(temp_dir+"/simulations/"+refs_list[-1]+"/*.h5")) < 2:
    print "checking " + temp_dir+"/simulations/"+refs_list[-1]
    print ("Waiting for NEURON to finish writing to disk...")
    time.sleep(2)
old_file_number = 0
while len(glob.glob(temp_dir+"/simulations/"+refs_list[-1]+"/*")) > old_file_number:
    print ("Waiting for NEURON to finish writing the results of the last simulation...")
    # this is meant to prevent the case in which the script exits while the results in the last folder are being written.
    old_file_number = len(glob.glob(temp_dir+"/simulations/"+refs_list[-1]+"/*"))
    time.sleep(1)

for dn in refs_list:
    destination = data_folder_path_ctor(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias) + '/' + dn
    if os.path.isdir(destination):
        shutil.rmtree(destination)
    os.mkdir(destination)
    print "Moving "+ temp_dir+"/simulations/"+dn + " to " + destination
    for source in glob.glob(temp_dir+'/simulations/'+dn+'/*.h5'):
        shutil.move(source, destination)

shutil.rmtree(temp_dir)
System.exit(0)
