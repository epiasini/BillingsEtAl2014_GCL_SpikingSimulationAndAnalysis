#
#
#   A file which opens a neuroConstruct project, adds some cells and network connections
#   and then saves a NetworkML file with the net structure
#
#   Author: Padraig Gleeson
#
#   This file has been developed as part of the neuroConstruct project
#   This work has been funded by the Medical Research Council and the
#   Wellcome Trust
#
#

try:
	from java.io import File
	from java.lang import System
except ImportError:
	print "Note: this file should be run using ..\\nC.bat -python XXX.py' or './nC.sh -python XXX.py'"
	print "See http://www.neuroconstruct.org/docs/python.html for more details"
	quit()

from ucl.physiol.neuroconstruct.project import ProjectManager

from math import *


# Load an existing neuroConstruct project

projFile = File("/home/eugenio/phd/nC_projects/if_network/if_network.ncx")
print "Loading project from file: " + projFile.getAbsolutePath()+", exists: "+ str(projFile.exists())

pm = ProjectManager()
myProject = pm.loadProject(projFile)
print "Loaded project: " + myProject.getProjectName() 


# Add a number of cells to the generatedCellPositions, connections to generatedNetworkConnections
# and electrical inputs to generatedElecInputs
numCells = 12

for i in range(0, numCells) :
    x = 100 * sin(i * 2 *pi / numCells)
    y = 100 * cos(i * 2 *pi / numCells)
    myProject.generatedCellPositions.addPosition("SampleCellGroup", i, x,y,0)
    
    if i != numCells-1 : 
        myProject.generatedNetworkConnections.addSynapticConnection("NC1", i, i+1)


# Print details
print myProject.generatedCellPositions.details()
print myProject.generatedNetworkConnections.details()



# Save to a NetworkML file
myNetworkMLFile = File("TestPython/savedNetworks/nmlt.nml")

simConfig = myProject.simConfigInfo.getDefaultSimConfig() 

pm.saveNetworkStructureXML(myProject, myNetworkMLFile, 0, 0, simConfig.getName(), "Physiological Units")

print "Network structure saved to file: "+ myNetworkMLFile.getAbsolutePath()


System.exit(0)


