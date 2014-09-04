# -*- coding: utf-8 -*-
"""
To be used with something like this:
./nC.sh -python /home/ucbtepi/code/network/src/export_NeuroML_network_structure.py 4
"""
import sys
import subprocess
import shutil
from os.path import join

import java.lang
import java.io

from ucl.physiol.neuroconstruct.project import ProjectManager
from ucl.physiol.neuroconstruct.nmodleditor.processes import ProcessManager

sys.path.append('/home/ucbtepi/code/network/src')
sys.path.append('/home/eugenio/src/decorator-3.4.0/build/lib')
sys.path.append('/home/eugenio/src/networkx-1.9')

from utils.cluster_system import ClusterSystem
from utils.pure import SimpleParameterSpacePoint
from utils.network import generate_nC_network, generate_nC_saves, generate_nC_stimuli, set_tonic_GABA, scale_excitatory_conductances, export_to_neuroml

n_grc_dend = int(sys.argv[1])
connectivity_rule = 0 # 0: tissue model, 1: random bipartite graph

scripts_path = '/home/ucbtepi/code/network/src/scripts/'
project_path = '/home/ucbtepi/nC_projects/if_gl/'
project_filename = 'if_gl.ncx' # neuroConstruct project file name
sim_config_name = 'Default Simulation Configuration'
nC_seed = 1234

network_structures_dir = '/home/ucbtepi/code/network/data/network_structures'
out_dir = join(network_structures_dir,
               'gd{}'.format(n_grc_dend))

def main():
    with ClusterSystem() as system:
        point = SimpleParameterSpacePoint(n_grc_dend,connectivity_rule,0,0.5,1,0.3,0,1,0,80,0,10,0,1024,60,180)

        print("Exporting static NeuroML for network with " + str(n_grc_dend) + " dendrites.")

        work_dir = system.work_dir
        nml_gen_dir = join(work_dir, 'NeuroML')
        # copy nC project in temp folder
        shutil.copy2(project_path+project_filename, work_dir+"/"+project_filename)
        shutil.copytree(project_path+"morphologies", work_dir+"/morphologies")
        shutil.copytree(project_path+"cellMechanisms", work_dir+"/cellMechanisms")

        # load project and initialise
        project_file = java.io.File(work_dir + "/" + project_filename)
        print ('Loading project from file: ' + project_file.getAbsolutePath() + ", exists: " + str(project_file.exists()))
        pm = ProjectManager(None, None)
        project = pm.loadProject(project_file)
        print 'Loaded project: ' + project.getProjectName()

        # load existing simulation configuration
        sim_config = project.simConfigInfo.getSimConfig(sim_config_name)

        # generate network
        generate_nC_network(point, pm, project, sim_config)

        # generate NeuroML2 files
        export_to_neuroml(project, sim_config, nml_gen_dir)

        # move NeuroML2 files to the network structure directory. If the
        # cell morphology and mechanisms files are needed, just add them
        # here.
        shutil.move(join(nml_gen_dir, 'if_gl.net.nml'),
                    join(out_dir, 'gd{}_cr{}.net.nml'.format(n_grc_dend,
                                                             connectivity_rule)))
    print("done exporting network")
    java.lang.System.exit(0)

if __name__ == "__main__":
    main()
