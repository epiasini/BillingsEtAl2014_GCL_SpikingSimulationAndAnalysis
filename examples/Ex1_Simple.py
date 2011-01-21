#
#
#   A simple example of loading a neuroConstruct project using the Jython scripting interface
#
#   To execute this type of file, type '..\nC.bat -python Ex1_Simple.py' (Windows)
#   or '../nC.sh -python Ex1_Simple.py' (Linux/Mac). Note: you may have to update the
#   NC_HOME and NC_MAX_MEMORY variables in nC.bat/nC.sh
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

file = File("/home/eugenio/phd/nC_projects/if_network/if_network.ncx")
print 'Loading project file: ', file.getAbsolutePath()

pm = ProjectManager()
myProject = pm.loadProject(file)

print pm.status()

#  Remove this line to remain in interactive mode
System.exit(0)
