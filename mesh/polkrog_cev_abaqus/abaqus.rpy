# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2021.HF12 replay file
# Internal Version: 2022_03_29-21.44.38 168405
# Run by Aljaz 2 on Wed Apr 17 20:36:51 2024
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=77.3874969482422, 
    height=34.2577819824219)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
openMdb('modelska_cev.cae')
#: The model database "C:\Users\Aljaz 2\Desktop\modelska_cev.cae" has been opened.
session.viewports['Viewport: 1'].setValues(displayedObject=None)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
p = mdb.models['Model-1'].parts['cel']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
session.viewports['Viewport: 1'].view.setValues(nearPlane=4.05866, 
    farPlane=4.88561, width=2.29458, height=1.51805, viewOffsetX=0.0418109, 
    viewOffsetY=-0.0780669)
session.viewports['Viewport: 1'].partDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=ON)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=OFF)
p1 = mdb.models['Model-1'].parts['polkrog']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
p1 = mdb.models['Model-1'].parts['cel']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
p1 = mdb.models['Model-1'].parts['polkrog']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    adaptiveMeshConstraints=ON)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=ON, 
    constraints=ON, connectors=ON, engineeringFeatures=ON, 
    adaptiveMeshConstraints=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON, 
    interactions=OFF, constraints=OFF, connectors=OFF, engineeringFeatures=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
a = mdb.models['Model-1'].rootAssembly
partInstances =(a.instances['cel-1'], )
a.deleteMesh(regions=partInstances)
a = mdb.models['Model-1'].rootAssembly
partInstances =(a.instances['cel-1'], )
a.deleteSeeds(regions=partInstances)
a = mdb.models['Model-1'].rootAssembly
partInstances =(a.instances['cel-1'], )
a.seedPartInstance(regions=partInstances, size=0.01, deviationFactor=0.1, 
    minSizeFactor=0.1)
a = mdb.models['Model-1'].rootAssembly
partInstances =(a.instances['cel-1'], )
a.generateMesh(regions=partInstances)
session.viewports['Viewport: 1'].view.setValues(nearPlane=4.22117, 
    farPlane=4.7231, width=1.38994, height=0.923434, viewOffsetX=-0.757967, 
    viewOffsetY=0.19255)
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.68677, 
    farPlane=4.96211, width=0.884697, height=0.587764, cameraPosition=(2.92223, 
    -1.80814, 2.12659), cameraUpVector=(0.301357, 0.952031, 0.0531239), 
    cameraTarget=(-0.380067, -0.608436, -0.640281), viewOffsetX=-0.482444, 
    viewOffsetY=0.122558)
session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF, 
    optimizationTasks=ON, geometricRestrictions=ON, stopConditions=ON)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
mdb.jobs.changeKey(fromName='mesh', toName='mesh_01')
session.viewports['Viewport: 1'].partDisplay.setValues(mesh=OFF)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=OFF)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
p = mdb.models['Model-1'].parts['polkrog']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=ON, 
    engineeringFeatures=ON)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=OFF)
mdb.models['Model-1'].Material(name='Material-1')
mdb.models['Model-1'].materials['Material-1'].Conductivity(table=((1.0, ), ))
mdb.models['Model-1'].HomogeneousSolidSection(name='Section-1', 
    material='Material-1', thickness=None)
p = mdb.models['Model-1'].parts['polkrog']
f = p.faces
faces = f.getSequenceFromMask(mask=('[#3 ]', ), )
region = regionToolset.Region(faces=faces)
p = mdb.models['Model-1'].parts['polkrog']
p.SectionAssignment(region=region, sectionName='Section-1', offset=0.0, 
    offsetType=MIDDLE_SURFACE, offsetField='', 
    thicknessAssignment=FROM_SECTION)
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    adaptiveMeshConstraints=ON)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=ON, 
    constraints=ON, connectors=ON, engineeringFeatures=ON, 
    adaptiveMeshConstraints=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
    predefinedFields=ON, interactions=OFF, constraints=OFF, 
    engineeringFeatures=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, 
    predefinedFields=OFF, connectors=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    adaptiveMeshConstraints=ON)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    adaptiveMeshConstraints=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=ON, 
    constraints=ON, connectors=ON, engineeringFeatures=ON)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
    predefinedFields=ON, interactions=OFF, constraints=OFF, 
    engineeringFeatures=OFF)
session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=OFF, 
    engineeringFeatures=OFF)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
p = mdb.models['Model-1'].parts['polkrog']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=ON, 
    engineeringFeatures=ON)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=OFF)
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
del mdb.models['Model-1'].boundaryConditions['boundary']
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, 
    predefinedFields=OFF, connectors=OFF)
p = mdb.models['Model-1'].parts['polkrog']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    adaptiveMeshConstraints=ON)
mdb.models['Model-1'].HeatTransferStep(name='Step-1', previous='Initial', 
    response=STEADY_STATE, maxNumInc=1, timeIncrementationMethod=FIXED, 
    amplitude=RAMP)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-1')
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
    predefinedFields=ON, connectors=ON, adaptiveMeshConstraints=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['cel-1'].edges
edges1 = e1.getSequenceFromMask(mask=('[#3 ]', ), )
region = a.Set(edges=edges1, name='boundary_points')
mdb.models['Model-1'].TemperatureBC(name='boundary', createStepName='Step-1', 
    region=region, fixed=OFF, distributionType=UNIFORM, fieldName='', 
    magnitude=0.0, amplitude=UNSET)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON, loads=OFF, 
    bcs=OFF, predefinedFields=OFF, connectors=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=OFF)
mdb.jobs['mesh_01'].writeInput(consistencyChecking=OFF)
#: The job input file has been written to "mesh_01.inp".
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
    predefinedFields=ON, connectors=ON)
a = mdb.models['Model-1'].rootAssembly
region = a.sets['boundary']
mdb.models['Model-1'].boundaryConditions['boundary'].setValues(region=region)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, 
    predefinedFields=OFF, connectors=OFF)
mdb.jobs['mesh_01'].writeInput(consistencyChecking=OFF)
#: The job input file has been written to "mesh_01.inp".
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
    predefinedFields=ON, connectors=ON)
del mdb.models['Model-1'].boundaryConditions['boundary']
a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['cel-1'].edges
edges1 = e1.getSequenceFromMask(mask=('[#3 ]', ), )
region = regionToolset.Region(edges=edges1)
mdb.models['Model-1'].TemperatureBC(name='boundary', createStepName='Step-1', 
    region=region, fixed=OFF, distributionType=UNIFORM, fieldName='', 
    magnitude=0.0, amplitude=UNSET)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, 
    predefinedFields=OFF, connectors=OFF)
mdb.jobs['mesh_01'].writeInput(consistencyChecking=OFF)
#: The job input file has been written to "mesh_01.inp".
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
    predefinedFields=ON, connectors=ON)
del mdb.models['Model-1'].boundaryConditions['boundary']
session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Initial')
a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['cel-1'].edges
edges1 = e1.getSequenceFromMask(mask=('[#3 ]', ), )
region = a.Set(edges=edges1, name='set_points')
mdb.models['Model-1'].TemperatureBC(name='boundary', createStepName='Initial', 
    region=region, distributionType=UNIFORM, fieldName='', magnitude=0.0)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, 
    predefinedFields=OFF, connectors=OFF)
mdb.jobs['mesh_01'].writeInput(consistencyChecking=OFF)
#: The job input file has been written to "mesh_01.inp".
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
    predefinedFields=ON, connectors=ON)
del mdb.models['Model-1'].boundaryConditions['boundary']
a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['cel-1'].edges
edges1 = e1.getSequenceFromMask(mask=('[#3 ]', ), )
region = regionToolset.Region(edges=edges1)
mdb.models['Model-1'].TemperatureBC(name='boundary', createStepName='Initial', 
    region=region, distributionType=UNIFORM, fieldName='', magnitude=0.0)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-1')
session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Initial')
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, 
    predefinedFields=OFF, connectors=OFF)
mdb.jobs['mesh_01'].writeInput(consistencyChecking=OFF)
#: The job input file has been written to "mesh_01.inp".
mdb.jobs['mesh_01'].submit(consistencyChecking=OFF)
#: The job input file "mesh_01.inp" has been submitted for analysis.
#: Error in job mesh_01: 32384 elements have missing property definitions. The elements have been identified in element set ErrElemMissingSection.
#: Job mesh_01: Analysis Input File Processor aborted due to errors.
#: Error in job mesh_01: Analysis Input File Processor exited with an error - Please see the  mesh_01.dat file for possible error messages if the file exists.
#: Job mesh_01 aborted due to errors.
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=OFF)
p = mdb.models['Model-1'].parts['polkrog']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
p = mdb.models['Model-1'].parts['polkrog']
f = p.faces
faces = f.getSequenceFromMask(mask=('[#3 ]', ), )
region = regionToolset.Region(faces=faces)
p = mdb.models['Model-1'].parts['polkrog']
p.SectionAssignment(region=region, sectionName='Section-1', offset=0.0, 
    offsetType=MIDDLE_SURFACE, offsetField='', 
    thicknessAssignment=FROM_SECTION)
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
mdb.jobs['mesh_01'].submit(consistencyChecking=OFF)
#: The job input file "mesh_01.inp" has been submitted for analysis.
#: Error in job mesh_01: 32384 elements have missing property definitions. The elements have been identified in element set ErrElemMissingSection.
#: Job mesh_01: Analysis Input File Processor aborted due to errors.
#: Error in job mesh_01: Analysis Input File Processor exited with an error - Please see the  mesh_01.dat file for possible error messages if the file exists.
#: Job mesh_01 aborted due to errors.
p = mdb.models['Model-1'].parts['polkrog']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
    predefinedFields=ON, connectors=ON)
del mdb.models['Model-1'].rootAssembly.sets['set_points']
del mdb.models['Model-1'].rootAssembly.sets['boundary_points']
del mdb.models['Model-1'].boundaryConditions['boundary']
#: Warning: Cannot continue yet--complete the step or cancel the procedure.
a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['cel-1'].edges
edges1 = e1.getSequenceFromMask(mask=('[#3 ]', ), )
region = a.Set(edges=edges1, name='boundary_data')
mdb.models['Model-1'].TemperatureBC(name='boundary', createStepName='Initial', 
    region=region, distributionType=UNIFORM, fieldName='', magnitude=0.0)
del mdb.models['Model-1'].rootAssembly.sets['boundary']
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, 
    predefinedFields=OFF, connectors=OFF)
mdb.jobs['mesh_01'].writeInput(consistencyChecking=OFF)
#: The job input file has been written to "mesh_01.inp".
session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=ON, 
    constraints=ON, connectors=ON, engineeringFeatures=ON)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
    predefinedFields=ON, interactions=OFF, constraints=OFF, 
    engineeringFeatures=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
#: Warning: Cannot continue yet--complete the step or cancel the procedure.
session.viewports['Viewport: 1'].view.setValues(nearPlane=4.37196, 
    farPlane=4.57231, width=0.494676, height=0.327268, viewOffsetX=-0.385936, 
    viewOffsetY=-0.271168)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF)
session.viewports['Viewport: 1'].view.setValues(nearPlane=4.08283, 
    farPlane=4.86144, width=2.13646, height=1.41344, viewOffsetX=-0.176941, 
    viewOffsetY=0.0862497)
a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['cel-1'].edges
edges1 = e1.getSequenceFromMask(mask=('[#3 ]', ), )
a.Set(edges=edges1, name='boundary')
#: The set 'boundary' has been created (2 edges).
mdb.models['Model-1'].rootAssembly.sets.changeKey(fromName='boundary', 
    toName='boundary_geo')
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, 
    predefinedFields=OFF, connectors=OFF)
mdb.jobs['mesh_01'].writeInput(consistencyChecking=OFF)
#: The job input file has been written to "mesh_01.inp".
mdb.jobs['mesh_01'].submit(consistencyChecking=OFF)
#: The job input file "mesh_01.inp" has been submitted for analysis.
#: Error in job mesh_01: 32384 elements have missing property definitions. The elements have been identified in element set ErrElemMissingSection.
#: Job mesh_01: Analysis Input File Processor aborted due to errors.
#: Error in job mesh_01: Analysis Input File Processor exited with an error - Please see the  mesh_01.dat file for possible error messages if the file exists.
#: Job mesh_01 aborted due to errors.
session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=OFF, 
    engineeringFeatures=OFF)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
p1 = mdb.models['Model-1'].parts['cel']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
del mdb.models['Model-1'].parts['cel']
p = mdb.models['Model-1'].parts['polkrog']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
a1 = mdb.models['Model-1'].rootAssembly
a1.suppressFeatures(featureNames=('cel-1', ))
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
a1 = mdb.models['Model-1'].rootAssembly
p = mdb.models['Model-1'].parts['polkrog']
a1.Instance(name='polkrog-1', part=p, dependent=ON)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    adaptiveMeshConstraints=ON)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=ON, 
    constraints=ON, connectors=ON, engineeringFeatures=ON, 
    adaptiveMeshConstraints=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
    predefinedFields=ON, interactions=OFF, constraints=OFF, 
    engineeringFeatures=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['polkrog-1'].edges
edges1 = e1.getSequenceFromMask(mask=('[#1e ]', ), )
region = regionToolset.Region(edges=edges1)
mdb.models['Model-1'].boundaryConditions['boundary'].setValues(region=region)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON, loads=OFF, 
    bcs=OFF, predefinedFields=OFF, connectors=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
p1 = mdb.models['Model-1'].parts['polkrog']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
session.viewports['Viewport: 1'].partDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=ON)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=OFF)
p = mdb.models['Model-1'].parts['polkrog']
p.checkGeometry()
#: Part 'polkrog' contains valid geometry and topology.
#: Part 'polkrog' is a shell part(2 shell faces, 5 edges, 4 vertices).
p = mdb.models['Model-1'].parts['polkrog']
p.deleteMesh()
p = mdb.models['Model-1'].parts['polkrog']
p.seedPart(size=0.01, deviationFactor=0.1, minSizeFactor=0.1)
p = mdb.models['Model-1'].parts['polkrog']
p.generateMesh()
session.viewports['Viewport: 1'].view.setValues(nearPlane=4.11717, 
    farPlane=4.82711, width=1.75165, height=1.16374, viewOffsetX=-0.00592616, 
    viewOffsetY=0.0178007)
a1 = mdb.models['Model-1'].rootAssembly
a1.regenerate()
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=OFF)
mdb.jobs['mesh_01'].writeInput(consistencyChecking=OFF)
#: The job input file has been written to "mesh_01.inp".
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
    predefinedFields=ON, connectors=ON)
a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['polkrog-1'].edges
edges1 = e1.getSequenceFromMask(mask=('[#1e ]', ), )
a.Set(edges=edges1, name='boundary_geo')
#: The set 'boundary_geo' has been edited (4 edges).
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, 
    predefinedFields=OFF, connectors=OFF)
mdb.jobs['mesh_01'].writeInput(consistencyChecking=OFF)
#: The job input file has been written to "mesh_01.inp".
mdb.jobs['mesh_01'].submit(consistencyChecking=OFF)
#: The job input file "mesh_01.inp" has been submitted for analysis.
#: Error in job mesh_01: STRESS-DISPLACEMENT ELEMENTS OR OTHER ELEMENTS WITHOUT TEMPERATURE DEGREE OF FREEDOM ARE NOT ALLOWED IN A HEAT TRANSFER ANALYSIS
#: Job mesh_01: Analysis Input File Processor aborted due to errors.
#: Error in job mesh_01: Analysis Input File Processor exited with an error - Please see the  mesh_01.dat file for possible error messages if the file exists.
#: Job mesh_01 aborted due to errors.
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
p1 = mdb.models['Model-1'].parts['polkrog']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
p = mdb.models['Model-1'].parts['polkrog']
p.checkGeometry()
#: Part 'polkrog' contains valid geometry and topology.
#: Part 'polkrog' is a shell part(2 shell faces, 5 edges, 4 vertices).
elemType1 = mesh.ElemType(elemCode=DC2D4, elemLibrary=STANDARD)
elemType2 = mesh.ElemType(elemCode=DC2D3, elemLibrary=STANDARD)
p = mdb.models['Model-1'].parts['polkrog']
f = p.faces
faces = f.getSequenceFromMask(mask=('[#3 ]', ), )
pickedRegions =(faces, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))
a1 = mdb.models['Model-1'].rootAssembly
a1.regenerate()
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=OFF)
mdb.jobs['mesh_01'].submit(consistencyChecking=OFF)
#: The job input file "mesh_01.inp" has been submitted for analysis.
#: Job mesh_01: Analysis Input File Processor completed successfully.
#: Job mesh_01: Abaqus/Standard completed successfully.
#: Job mesh_01 completed successfully. 
mdb.Job(name='mesh_001', model='Model-1', description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1, 
    numGPUs=0)
o3 = session.openOdb(name='C:/Users/Aljaz 2/Desktop/mesh_01.odb')
#: Model: C:/Users/Aljaz 2/Desktop/mesh_01.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       2
#: Number of Node Sets:          2
#: Number of Steps:              1
session.viewports['Viewport: 1'].setValues(displayedObject=o3)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].view.setValues(nearPlane=4.89279, 
    farPlane=7.75632, width=1.74375, height=1.1026, viewOffsetX=0.00615823, 
    viewOffsetY=0.0116891)
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='HFL', outputPosition=INTEGRATION_POINT, refinement=(
    INVARIANT, 'Magnitude'), )
session.viewports['Viewport: 1'].odbDisplay.display.setValues(
    plotState=CONTOURS_ON_DEF)
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='NT11', outputPosition=NODAL, )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='RFL11', outputPosition=NODAL, )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='GRADT', outputPosition=INTEGRATION_POINT, refinement=(
    INVARIANT, 'Magnitude'), )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='HFL', outputPosition=INTEGRATION_POINT, refinement=(
    INVARIANT, 'Magnitude'), )
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
    predefinedFields=ON, connectors=ON)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-1')
a = mdb.models['Model-1'].rootAssembly
f1 = a.instances['polkrog-1'].faces
faces1 = f1.getSequenceFromMask(mask=('[#3 ]', ), )
region = regionToolset.Region(faces=faces1)
mdb.models['Model-1'].BodyHeatFlux(name='vol_gen', createStepName='Step-1', 
    region=region, magnitude=1.0)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, 
    predefinedFields=OFF, connectors=OFF)
mdb.jobs['mesh_01'].submit(consistencyChecking=OFF)
#: The job input file "mesh_01.inp" has been submitted for analysis.
#: Job mesh_01: Analysis Input File Processor completed successfully.
#: Job mesh_01: Abaqus/Standard completed successfully.
#: Job mesh_01 completed successfully. 
session.viewports['Viewport: 1'].setValues(
    displayedObject=session.odbs['C:/Users/Aljaz 2/Desktop/mesh_01.odb'])
#: Warning: The output database 'C:/Users/Aljaz 2/Desktop/mesh_01.odb' disk file has changed.
#: 
#: The current plot operation has been canceled, re-open the file to view the results
o3 = session.openOdb(name='C:/Users/Aljaz 2/Desktop/mesh_01.odb')
#: Model: C:/Users/Aljaz 2/Desktop/mesh_01.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       2
#: Number of Node Sets:          2
#: Number of Steps:              1
session.viewports['Viewport: 1'].setValues(displayedObject=o3)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='HFL', outputPosition=INTEGRATION_POINT, refinement=(
    INVARIANT, 'Magnitude'), )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='NT11', outputPosition=NODAL, )
session.viewports['Viewport: 1'].view.setValues(nearPlane=4.7633, 
    farPlane=7.88581, width=2.17432, height=1.37486, viewOffsetX=0.007209, 
    viewOffsetY=-0.0275226)
session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
session.viewports['Viewport: 1'].view.setValues(nearPlane=4.8249, 
    farPlane=7.82421, width=1.94608, height=1.23053, viewOffsetX=-0.114355, 
    viewOffsetY=0.0122037)
session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(
    visibleEdges=FEATURE)
session.viewports['Viewport: 1'].odbDisplay.contourOptions.setValues(
    contourType=LINE, maxValue=0.0976177, minValue=0)
session.viewports['Viewport: 1'].odbDisplay.contourOptions.setValues(
    contourType=QUILT)
#* VisError: Quilt plotting cannot be performed using results at nodes. Select 
#* an element result or change the contour type.
session.viewports['Viewport: 1'].odbDisplay.contourOptions.setValues(
    contourType=BANDED)
session.viewports['Viewport: 1'].odbDisplay.contourOptions.setValues(
    minAutoCompute=OFF)
session.viewports['Viewport: 1'].odbDisplay.contourOptions.setValues(
    minAutoCompute=ON)
session.viewports['Viewport: 1'].viewportAnnotationOptions.setValues(
    legendBox=OFF)
session.viewports['Viewport: 1'].viewportAnnotationOptions.setValues(
    legendMinMax=ON)
session.viewports['Viewport: 1'].viewportAnnotationOptions.setValues(
    legendMinMax=OFF)
session.viewports['Viewport: 1'].viewportAnnotationOptions.setValues(
    legendTitle=OFF)
session.viewports['Viewport: 1'].viewportAnnotationOptions.setValues(
    legendTitle=ON)
session.viewports['Viewport: 1'].viewportAnnotationOptions.setValues(
    legendNumberFormat=FIXED)
session.viewports['Viewport: 1'].viewportAnnotationOptions.setValues(
    legendNumberFormat=ENGINEERING)
session.viewports['Viewport: 1'].viewportAnnotationOptions.setValues(
    legendNumberFormat=FIXED)
session.viewports['Viewport: 1'].view.setValues(nearPlane=4.78788, 
    farPlane=7.86123, width=2.32505, height=1.47016, viewOffsetX=-0.133937, 
    viewOffsetY=0.00243549)
session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
session.viewports['Viewport: 1'].view.setValues(nearPlane=4.7633, 
    farPlane=7.88581, width=2.46075, height=1.55597, viewOffsetX=0.00594437, 
    viewOffsetY=-0.0594895)
session.viewports['Viewport: 1'].view.setValues(session.views['User-1'])
session.viewports['Viewport: 1'].view.setValues(nearPlane=4.73743, 
    farPlane=7.91168, width=2.30055, height=1.45467, viewOffsetX=0.0492925, 
    viewOffsetY=-0.054233)
session.viewports['Viewport: 1'].view.setValues(nearPlane=4.76514, 
    farPlane=7.88397, width=2.314, height=1.46318, viewOffsetX=0.00830268, 
    viewOffsetY=-0.0497132)
session.View(name='User-1', nearPlane=4.7651, farPlane=7.884, width=2.314, 
    height=1.4632, projection=PERSPECTIVE, cameraPosition=(0, -0.5, 6.3246), 
    cameraUpVector=(0, 1, 0), cameraTarget=(0, -0.5, 0), viewOffsetX=0.0083027, 
    viewOffsetY=-0.049713, autoFit=OFF)
session.viewports['Viewport: 1'].view.setValues(session.views['User-1'])
session.viewports['Viewport: 1'].view.setValues(session.views['User-2'])
session.viewports['Viewport: 1'].view.setValues(session.views['User-1'])
session.psOptions.setValues(logo=OFF, date=OFF, resolution=DPI_300)
session.printToFile(fileName='mesh_01', format=PS, canvasObjects=(
    session.viewports['Viewport: 1'], ))
p = mdb.models['Model-1'].parts['polkrog']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
p = mdb.models['Model-1'].parts['polkrog']
p.deleteMesh()
p = mdb.models['Model-1'].parts['polkrog']
p.seedPart(size=0.001, deviationFactor=0.1, minSizeFactor=0.1)
p = mdb.models['Model-1'].parts['polkrog']
p.generateMesh()
#*p.generateMesh()
#* Command Interrupted
p = mdb.models['Model-1'].parts['polkrog']
p.seedPart(size=0.005, deviationFactor=0.1, minSizeFactor=0.1)
session.viewports['Viewport: 1'].view.setValues(nearPlane=4.45353, 
    farPlane=4.49074, width=0.10341, height=0.0687024, viewOffsetX=-0.770829, 
    viewOffsetY=0.477628)
p = mdb.models['Model-1'].parts['polkrog']
p.generateMesh()
session.viewports['Viewport: 1'].view.setValues(nearPlane=4.25149, 
    farPlane=4.69278, width=1.22255, height=0.812222, viewOffsetX=-0.811992, 
    viewOffsetY=0.554147)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.18623, 
    farPlane=4.73749, width=0.916227, height=0.608711, cameraPosition=(1.68685, 
    -1.70046, 3.38117), cameraUpVector=(-0.0081292, 0.943476, 0.33134), 
    cameraTarget=(-0.0634248, -0.350254, -0.506434), viewOffsetX=-0.608538, 
    viewOffsetY=0.415299)
session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
session.viewports['Viewport: 1'].view.setValues(session.views['User-1'])
session.viewports['Viewport: 1'].view.setValues(nearPlane=5.96065, 
    farPlane=6.68855, width=1.79453, height=1.19223, viewOffsetX=-0.0213797, 
    viewOffsetY=0.084692)
a1 = mdb.models['Model-1'].rootAssembly
a1.regenerate()
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
a = mdb.models['Model-1'].rootAssembly
a.SetFromNodeLabels(name='Set-3', nodeLabels=(('polkrog-1', (2, 204, 205, 206, 
    207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 
    222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 
    237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 
    252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 
    267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 
    282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 
    297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 
    312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 
    327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 
    342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 
    357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 
    372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 
    387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 
    402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 
    417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 
    432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 
    447, 448, 449, 450, 451, 452, 453, 454, 455, 456, 457, 458, 459, 460, 461, 
    462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 
    477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 
    492, 493, 494, 495, 496, 497, 498, 499, 500, 501, 502, 503, 504, 505, 506, 
    507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 3)), ('polkrog-1', (3, 
    517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 
    532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 
    547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561, 
    562, 563, 564, 565, 566, 567, 568, 569, 570, 571, 572, 573, 574, 575, 576, 
    577, 578, 579, 580, 581, 582, 583, 584, 585, 586, 587, 588, 589, 590, 591, 
    592, 593, 594, 595, 596, 597, 598, 599, 600, 601, 602, 603, 604, 605, 606, 
    607, 608, 609, 610, 611, 612, 613, 614, 615, 616, 617, 618, 619, 620, 621, 
    622, 623, 624, 625, 626, 627, 628, 629, 630, 631, 632, 633, 634, 635, 636, 
    637, 638, 639, 640, 641, 642, 643, 644, 645, 646, 647, 648, 649, 650, 651, 
    652, 653, 654, 655, 656, 657, 658, 659, 660, 661, 662, 663, 664, 665, 666, 
    667, 668, 669, 670, 671, 672, 673, 674, 675, 676, 677, 678, 679, 680, 681, 
    682, 683, 684, 685, 686, 687, 688, 689, 690, 691, 692, 693, 694, 695, 696, 
    697, 698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 
    712, 713, 714, 715, 1)), ('polkrog-1', (1, 716, 717, 718, 719, 720, 721, 
    722, 723, 724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 
    737, 738, 739, 740, 741, 742, 743, 744, 745, 746, 747, 748, 749, 750, 751, 
    752, 753, 754, 755, 756, 757, 758, 759, 760, 761, 762, 763, 764, 765, 766, 
    767, 768, 769, 770, 771, 772, 773, 774, 775, 776, 777, 778, 779, 780, 781, 
    782, 783, 784, 785, 786, 787, 788, 789, 790, 791, 792, 793, 794, 795, 796, 
    797, 798, 799, 800, 801, 802, 803, 804, 805, 806, 807, 808, 809, 810, 811, 
    812, 813, 814, 815, 816, 817, 818, 819, 820, 821, 822, 823, 824, 825, 826, 
    827, 828, 829, 830, 831, 832, 833, 834, 835, 836, 837, 838, 839, 840, 841, 
    842, 843, 844, 845, 846, 847, 848, 849, 850, 851, 852, 853, 854, 855, 856, 
    857, 858, 859, 860, 861, 862, 863, 864, 865, 866, 867, 868, 869, 870, 871, 
    872, 873, 874, 875, 876, 877, 878, 879, 880, 881, 882, 883, 884, 885, 886, 
    887, 888, 889, 890, 891, 892, 893, 894, 895, 896, 897, 898, 899, 900, 901, 
    902, 903, 904, 905, 906, 907, 908, 909, 910, 911, 912, 913, 914, 4)), (
    'polkrog-1', (4, 915, 916, 917, 918, 919, 920, 921, 922, 923, 924, 925, 
    926, 927, 928, 929, 930, 931, 932, 933, 934, 935, 936, 937, 938, 939, 940, 
    941, 942, 943, 944, 945, 946, 947, 948, 949, 950, 951, 952, 953, 954, 955, 
    956, 957, 958, 959, 960, 961, 962, 963, 964, 965, 966, 967, 968, 969, 970, 
    971, 972, 973, 974, 975, 976, 977, 978, 979, 980, 981, 982, 983, 984, 985, 
    986, 987, 988, 989, 990, 991, 992, 993, 994, 995, 996, 997, 998, 999, 1000, 
    1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 1011, 1012, 
    1013, 1014, 1015, 1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023, 1024, 
    1025, 1026, 1027, 1028, 1029, 1030, 1031, 1032, 1033, 1034, 1035, 1036, 
    1037, 1038, 1039, 1040, 1041, 1042, 1043, 1044, 1045, 1046, 1047, 1048, 
    1049, 1050, 1051, 1052, 1053, 1054, 1055, 1056, 1057, 1058, 1059, 1060, 
    1061, 1062, 1063, 1064, 1065, 1066, 1067, 1068, 1069, 1070, 1071, 1072, 
    1073, 1074, 1075, 1076, 1077, 1078, 1079, 1080, 1081, 1082, 1083, 1084, 
    1085, 1086, 1087, 1088, 1089, 1090, 1091, 1092, 1093, 1094, 1095, 1096, 
    1097, 1098, 1099, 1100, 1101, 1102, 1103, 1104, 1105, 1106, 1107, 1108, 
    1109, 1110, 1111, 1112, 1113, 1114, 1115, 1116, 1117, 1118, 1119, 1120, 
    1121, 1122, 1123, 1124, 1125, 1126, 1127, 1128, 1129, 1130, 1131, 1132, 
    1133, 1134, 1135, 1136, 1137, 1138, 1139, 1140, 1141, 1142, 1143, 1144, 
    1145, 1146, 1147, 1148, 1149, 1150, 1151, 1152, 1153, 1154, 1155, 1156, 
    1157, 1158, 1159, 1160, 1161, 1162, 1163, 1164, 1165, 1166, 1167, 1168, 
    1169, 1170, 1171, 1172, 1173, 1174, 1175, 1176, 1177, 1178, 1179, 1180, 
    1181, 1182, 1183, 1184, 1185, 1186, 1187, 1188, 1189, 1190, 1191, 1192, 
    1193, 1194, 1195, 1196, 1197, 1198, 1199, 1200, 1201, 1202, 1203, 1204, 
    1205, 1206, 1207, 1208, 1209, 1210, 1211, 1212, 1213, 1214, 1215, 1216, 
    1217, 1218, 1219, 1220, 1221, 1222, 1223, 1224, 1225, 1226, 1227, 2))), 
    unsorted=True)
del mdb.models['Model-1'].rootAssembly.sets['boundary_data']
del mdb.models['Model-1'].rootAssembly.sets['boundary_geo']
mdb.models['Model-1'].rootAssembly.sets.changeKey(fromName='Set-3', 
    toName='boundary_nodes')
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF, loads=ON, 
    bcs=ON, predefinedFields=ON, connectors=ON)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Initial')
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, 
    predefinedFields=OFF, connectors=OFF)
mdb.jobs.changeKey(fromName='mesh_001', toName='mesh_005')
mdb.jobs['mesh_005'].writeInput(consistencyChecking=OFF)
#: The job input file has been written to "mesh_005.inp".
import os
os.chdir(r"C:\Users\Aljaz 2\Desktop\polkrog_cev")
mdb.jobs['mesh_005'].writeInput(consistencyChecking=OFF)
#: The job input file has been written to "mesh_005.inp".
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
p1 = mdb.models['Model-1'].parts['polkrog']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
p = mdb.models['Model-1'].parts['polkrog']
p.checkGeometry()
#: Part 'polkrog' contains valid geometry and topology.
#: Part 'polkrog' is a shell part(2 shell faces, 5 edges, 4 vertices).
p = mdb.models['Model-1'].parts['polkrog']
p.deleteMesh()
p = mdb.models['Model-1'].parts['polkrog']
p.seedPart(size=0.01, deviationFactor=0.1, minSizeFactor=0.1)
p = mdb.models['Model-1'].parts['polkrog']
p.generateMesh()
a1 = mdb.models['Model-1'].rootAssembly
a1.regenerate()
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=OFF)
mdb.jobs['mesh_01'].writeInput(consistencyChecking=OFF)
#: The job input file has been written to "mesh_01.inp".
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
p1 = mdb.models['Model-1'].parts['polkrog']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
p = mdb.models['Model-1'].parts['polkrog']
p.checkGeometry()
#: Part 'polkrog' contains valid geometry and topology.
#: Part 'polkrog' is a shell part(2 shell faces, 5 edges, 4 vertices).
p = mdb.models['Model-1'].parts['polkrog']
p.deleteMesh()
p = mdb.models['Model-1'].parts['polkrog']
p.seedPart(size=0.0025, deviationFactor=0.1, minSizeFactor=0.1)
p = mdb.models['Model-1'].parts['polkrog']
p.generateMesh()
#* Command Interrupted
p = mdb.models['Model-1'].parts['polkrog']
p.seedPart(size=0.01, deviationFactor=0.1, minSizeFactor=0.1)
p = mdb.models['Model-1'].parts['polkrog']
p.generateMesh()
#: Warning: Cannot continue yet--complete the step or cancel the procedure.
p = mdb.models['Model-1'].parts['polkrog']
p.SetFromNodeLabels(name='Set-1', nodeLabels=(2, 104, 105, 106, 107, 108, 109, 
    110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 
    125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 
    140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 
    155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 
    170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 
    185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 
    200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 
    215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 
    230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 
    245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 
    3, 3, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 
    274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 
    289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 
    304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 
    319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 
    334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 
    349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 1, 1, 359, 360, 361, 362, 
    363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 
    378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 
    393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 
    408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 
    423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 
    438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 
    453, 454, 455, 456, 457, 4, 4, 458, 459, 460, 461, 462, 463, 464, 465, 466, 
    467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480, 481, 
    482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 
    497, 498, 499, 500, 501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511, 
    512, 513, 514, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 
    527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 
    542, 543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 
    557, 558, 559, 560, 561, 562, 563, 564, 565, 566, 567, 568, 569, 570, 571, 
    572, 573, 574, 575, 576, 577, 578, 579, 580, 581, 582, 583, 584, 585, 586, 
    587, 588, 589, 590, 591, 592, 593, 594, 595, 596, 597, 598, 599, 600, 601, 
    602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613, 2), 
    unsorted=True)
mdb.models['Model-1'].parts['polkrog'].sets.changeKey(fromName='Set-1', 
    toName='boundary_nodes')
a1 = mdb.models['Model-1'].rootAssembly
a1.regenerate()
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=OFF)
mdb.jobs['mesh_01'].writeInput(consistencyChecking=OFF)
#: The job input file has been written to "mesh_01.inp".
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
p1 = mdb.models['Model-1'].parts['polkrog']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
p = mdb.models['Model-1'].parts['polkrog']
p.checkGeometry()
#: Part 'polkrog' contains valid geometry and topology.
#: Part 'polkrog' is a shell part(2 shell faces, 5 edges, 4 vertices).
p = mdb.models['Model-1'].parts['polkrog']
p.deleteMesh()
p = mdb.models['Model-1'].parts['polkrog']
p.seedPart(size=0.0025, deviationFactor=0.1, minSizeFactor=0.1)
p = mdb.models['Model-1'].parts['polkrog']
p.generateMesh()
session.viewports['Viewport: 1'].view.setValues(nearPlane=5.78765, 
    farPlane=6.86155, width=3.00704, height=1.99778, viewOffsetX=0.145707, 
    viewOffsetY=-0.494298)
session.viewports['Viewport: 1'].view.setValues(nearPlane=5.21518, 
    farPlane=6.95126, width=2.70961, height=1.80018, cameraPosition=(0.419414, 
    2.39074, 5.33657), cameraUpVector=(0.0327077, 0.883068, -0.468104), 
    cameraTarget=(0.0180735, -0.553835, -0.246352), viewOffsetX=0.131295, 
    viewOffsetY=-0.445406)
session.viewports['Viewport: 1'].view.setValues(session.views['User-1'])
a1 = mdb.models['Model-1'].rootAssembly
a1.regenerate()
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=OFF)
mdb.Job(name='mesh_0025', model='Model-1', description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1, 
    numGPUs=0)
mdb.jobs['mesh_0025'].setValues(numCpus=4, numDomains=4)
mdb.jobs['mesh_0025'].writeInput(consistencyChecking=OFF)
#* Academic Teaching license is restricted to 250000 nodes.
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
p1 = mdb.models['Model-1'].parts['polkrog']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
p = mdb.models['Model-1'].parts['polkrog']
p.checkGeometry()
#: Part 'polkrog' contains valid geometry and topology.
#: Part 'polkrog' is a shell part(2 shell faces, 5 edges, 4 vertices).
p = mdb.models['Model-1'].parts['polkrog']
p.deleteMesh()
p = mdb.models['Model-1'].parts['polkrog']
p.seedPart(size=0.003, deviationFactor=0.1, minSizeFactor=0.1)
p = mdb.models['Model-1'].parts['polkrog']
p.generateMesh()
p = mdb.models['Model-1'].parts['polkrog']
p.deleteMesh()
p = mdb.models['Model-1'].parts['polkrog']
p.seedPart(size=0.004, deviationFactor=0.1, minSizeFactor=0.1)
p = mdb.models['Model-1'].parts['polkrog']
p.generateMesh()
a1 = mdb.models['Model-1'].rootAssembly
a1.regenerate()
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
p1 = mdb.models['Model-1'].parts['polkrog']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
p = mdb.models['Model-1'].parts['polkrog']
p.checkGeometry()
#: Part 'polkrog' contains valid geometry and topology.
#: Part 'polkrog' is a shell part(2 shell faces, 5 edges, 4 vertices).
p = mdb.models['Model-1'].parts['polkrog']
p.deleteMesh()
p = mdb.models['Model-1'].parts['polkrog']
p.seedPart(size=0.1, deviationFactor=0.1, minSizeFactor=0.1)
p = mdb.models['Model-1'].parts['polkrog']
p.generateMesh()
a1 = mdb.models['Model-1'].rootAssembly
a1.regenerate()
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=OFF)
mdb.jobs.changeKey(fromName='mesh_0025', toName='mesh_1')
mdb.jobs['mesh_1'].writeInput(consistencyChecking=OFF)
#: The job input file has been written to "mesh_1.inp".
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF, loads=ON, 
    bcs=ON, predefinedFields=ON, connectors=ON)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
a = mdb.models['Model-1'].rootAssembly
a.SetFromNodeLabels(name='Set-2', nodeLabels=(('polkrog-1', (2, 14, 15, 16, 17, 
    18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 3)), ('polkrog-1', (3, 29, 30, 
    31, 32, 33, 34, 35, 36, 37, 1)), ('polkrog-1', (1, 38, 39, 40, 41, 42, 43, 
    44, 45, 46, 4)), ('polkrog-1', (4, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 
    57, 58, 59, 60, 61, 2))), unsorted=True)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF)
del mdb.models['Model-1'].rootAssembly.sets['boundary_nodes']
mdb.models['Model-1'].rootAssembly.sets.changeKey(fromName='Set-2', 
    toName='boundary_nodes')
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, 
    predefinedFields=OFF, connectors=OFF)
mdb.jobs['mesh_1'].writeInput(consistencyChecking=OFF)
#: The job input file has been written to "mesh_1.inp".
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
p1 = mdb.models['Model-1'].parts['polkrog']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
p = mdb.models['Model-1'].parts['polkrog']
p.checkGeometry()
#: Part 'polkrog' contains valid geometry and topology.
#: Part 'polkrog' is a shell part(2 shell faces, 5 edges, 4 vertices).
p = mdb.models['Model-1'].parts['polkrog']
p.deleteMesh()
p = mdb.models['Model-1'].parts['polkrog']
p.seedPart(size=0.05, deviationFactor=0.1, minSizeFactor=0.1)
p = mdb.models['Model-1'].parts['polkrog']
p.generateMesh()
del mdb.models['Model-1'].parts['polkrog'].sets['boundary_nodes']
p = mdb.models['Model-1'].parts['polkrog']
p.SetFromNodeLabels(name='boundary_nodes', nodeLabels=(2, 24, 25, 26, 27, 28, 
    29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 
    48, 49, 50, 51, 52, 53, 3, 3, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 
    65, 66, 67, 68, 69, 70, 71, 72, 1, 1, 73, 74, 75, 76, 77, 78, 79, 80, 81, 
    82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 4, 4, 92, 93, 94, 95, 96, 97, 98, 
    99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 
    114, 115, 116, 117, 118, 119, 120, 121, 2), unsorted=True)
a1 = mdb.models['Model-1'].rootAssembly
a1.regenerate()
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=OFF)
mdb.jobs['mesh_1'].writeInput(consistencyChecking=OFF)
#: The job input file has been written to "mesh_1.inp".
mdb.Job(name='mesh_05', model='Model-1', description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1, 
    numGPUs=0)
mdb.jobs['mesh_05'].writeInput(consistencyChecking=OFF)
#: The job input file has been written to "mesh_05.inp".
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
p1 = mdb.models['Model-1'].parts['polkrog']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
p = mdb.models['Model-1'].parts['polkrog']
p.checkGeometry()
#: Part 'polkrog' contains valid geometry and topology.
#: Part 'polkrog' is a shell part(2 shell faces, 5 edges, 4 vertices).
p = mdb.models['Model-1'].parts['polkrog']
p.deleteMesh()
p = mdb.models['Model-1'].parts['polkrog']
p.seedPart(size=0.1, deviationFactor=0.1, minSizeFactor=0.1)
p = mdb.models['Model-1'].parts['polkrog']
p.generateMesh()
del mdb.models['Model-1'].parts['polkrog'].sets['boundary_nodes']
p = mdb.models['Model-1'].parts['polkrog']
p.SetFromNodeLabels(name='boundary_nodes', nodeLabels=(2, 14, 15, 16, 17, 18, 
    19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 3, 3, 29, 30, 31, 32, 33, 34, 35, 
    36, 37, 1, 1, 38, 39, 40, 41, 42, 43, 44, 45, 46, 4, 4, 47, 48, 49, 50, 51, 
    52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 2), unsorted=True)
a1 = mdb.models['Model-1'].rootAssembly
a1.regenerate()
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=OFF)
mdb.jobs['mesh_1'].writeInput(consistencyChecking=OFF)
#: The job input file has been written to "mesh_1.inp".
mdb.jobs['mesh_005'].setValues(numCpus=4, numDomains=4)
mdb.jobs['mesh_005'].submit(consistencyChecking=OFF)
#: The job input file "mesh_005.inp" has been submitted for analysis.
mdb.jobs['mesh_005'].kill()
#: Job mesh_005: Analysis Input File Processor completed successfully.
#: Job mesh_005: Abaqus/Standard completed successfully.
#: Job mesh_005 completed successfully. 
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
p1 = mdb.models['Model-1'].parts['polkrog']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
p = mdb.models['Model-1'].parts['polkrog']
p.checkGeometry()
#: Part 'polkrog' contains valid geometry and topology.
#: Part 'polkrog' is a shell part(2 shell faces, 5 edges, 4 vertices).
p = mdb.models['Model-1'].parts['polkrog']
p.deleteMesh()
p = mdb.models['Model-1'].parts['polkrog']
p.seedPart(size=0.005, deviationFactor=0.1, minSizeFactor=0.1)
p = mdb.models['Model-1'].parts['polkrog']
p.generateMesh()
del mdb.models['Model-1'].parts['polkrog'].sets['boundary_nodes']
p = mdb.models['Model-1'].parts['polkrog']
p.SetFromNodeLabels(name='boundary_nodes', nodeLabels=(2, 204, 205, 206, 207, 
    208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 
    223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 
    238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 
    253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 
    268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 
    283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 
    298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 
    313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 
    328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 
    343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 
    358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 
    373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 
    388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 
    403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 
    418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 
    433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 
    448, 449, 450, 451, 452, 453, 454, 455, 456, 457, 458, 459, 460, 461, 462, 
    463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 
    478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 
    493, 494, 495, 496, 497, 498, 499, 500, 501, 502, 503, 504, 505, 506, 507, 
    508, 509, 510, 511, 512, 513, 514, 515, 516, 3, 3, 517, 518, 519, 520, 521, 
    522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 
    537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550, 551, 
    552, 553, 554, 555, 556, 557, 558, 559, 560, 561, 562, 563, 564, 565, 566, 
    567, 568, 569, 570, 571, 572, 573, 574, 575, 576, 577, 578, 579, 580, 581, 
    582, 583, 584, 585, 586, 587, 588, 589, 590, 591, 592, 593, 594, 595, 596, 
    597, 598, 599, 600, 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 
    612, 613, 614, 615, 616, 617, 618, 619, 620, 621, 622, 623, 624, 625, 626, 
    627, 628, 629, 630, 631, 632, 633, 634, 635, 636, 637, 638, 639, 640, 641, 
    642, 643, 644, 645, 646, 647, 648, 649, 650, 651, 652, 653, 654, 655, 656, 
    657, 658, 659, 660, 661, 662, 663, 664, 665, 666, 667, 668, 669, 670, 671, 
    672, 673, 674, 675, 676, 677, 678, 679, 680, 681, 682, 683, 684, 685, 686, 
    687, 688, 689, 690, 691, 692, 693, 694, 695, 696, 697, 698, 699, 700, 701, 
    702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 1, 1, 
    716, 717, 718, 719, 720, 721, 722, 723, 724, 725, 726, 727, 728, 729, 730, 
    731, 732, 733, 734, 735, 736, 737, 738, 739, 740, 741, 742, 743, 744, 745, 
    746, 747, 748, 749, 750, 751, 752, 753, 754, 755, 756, 757, 758, 759, 760, 
    761, 762, 763, 764, 765, 766, 767, 768, 769, 770, 771, 772, 773, 774, 775, 
    776, 777, 778, 779, 780, 781, 782, 783, 784, 785, 786, 787, 788, 789, 790, 
    791, 792, 793, 794, 795, 796, 797, 798, 799, 800, 801, 802, 803, 804, 805, 
    806, 807, 808, 809, 810, 811, 812, 813, 814, 815, 816, 817, 818, 819, 820, 
    821, 822, 823, 824, 825, 826, 827, 828, 829, 830, 831, 832, 833, 834, 835, 
    836, 837, 838, 839, 840, 841, 842, 843, 844, 845, 846, 847, 848, 849, 850, 
    851, 852, 853, 854, 855, 856, 857, 858, 859, 860, 861, 862, 863, 864, 865, 
    866, 867, 868, 869, 870, 871, 872, 873, 874, 875, 876, 877, 878, 879, 880, 
    881, 882, 883, 884, 885, 886, 887, 888, 889, 890, 891, 892, 893, 894, 895, 
    896, 897, 898, 899, 900, 901, 902, 903, 904, 905, 906, 907, 908, 909, 910, 
    911, 912, 913, 914, 4, 4, 915, 916, 917, 918, 919, 920, 921, 922, 923, 924, 
    925, 926, 927, 928, 929, 930, 931, 932, 933, 934, 935, 936, 937, 938, 939, 
    940, 941, 942, 943, 944, 945, 946, 947, 948, 949, 950, 951, 952, 953, 954, 
    955, 956, 957, 958, 959, 960, 961, 962, 963, 964, 965, 966, 967, 968, 969, 
    970, 971, 972, 973, 974, 975, 976, 977, 978, 979, 980, 981, 982, 983, 984, 
    985, 986, 987, 988, 989, 990, 991, 992, 993, 994, 995, 996, 997, 998, 999, 
    1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 1011, 
    1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023, 
    1024, 1025, 1026, 1027, 1028, 1029, 1030, 1031, 1032, 1033, 1034, 1035, 
    1036, 1037, 1038, 1039, 1040, 1041, 1042, 1043, 1044, 1045, 1046, 1047, 
    1048, 1049, 1050, 1051, 1052, 1053, 1054, 1055, 1056, 1057, 1058, 1059, 
    1060, 1061, 1062, 1063, 1064, 1065, 1066, 1067, 1068, 1069, 1070, 1071, 
    1072, 1073, 1074, 1075, 1076, 1077, 1078, 1079, 1080, 1081, 1082, 1083, 
    1084, 1085, 1086, 1087, 1088, 1089, 1090, 1091, 1092, 1093, 1094, 1095, 
    1096, 1097, 1098, 1099, 1100, 1101, 1102, 1103, 1104, 1105, 1106, 1107, 
    1108, 1109, 1110, 1111, 1112, 1113, 1114, 1115, 1116, 1117, 1118, 1119, 
    1120, 1121, 1122, 1123, 1124, 1125, 1126, 1127, 1128, 1129, 1130, 1131, 
    1132, 1133, 1134, 1135, 1136, 1137, 1138, 1139, 1140, 1141, 1142, 1143, 
    1144, 1145, 1146, 1147, 1148, 1149, 1150, 1151, 1152, 1153, 1154, 1155, 
    1156, 1157, 1158, 1159, 1160, 1161, 1162, 1163, 1164, 1165, 1166, 1167, 
    1168, 1169, 1170, 1171, 1172, 1173, 1174, 1175, 1176, 1177, 1178, 1179, 
    1180, 1181, 1182, 1183, 1184, 1185, 1186, 1187, 1188, 1189, 1190, 1191, 
    1192, 1193, 1194, 1195, 1196, 1197, 1198, 1199, 1200, 1201, 1202, 1203, 
    1204, 1205, 1206, 1207, 1208, 1209, 1210, 1211, 1212, 1213, 1214, 1215, 
    1216, 1217, 1218, 1219, 1220, 1221, 1222, 1223, 1224, 1225, 1226, 1227, 2), 
    unsorted=True)
a1 = mdb.models['Model-1'].rootAssembly
a1.regenerate()
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=OFF)
mdb.jobs['mesh_005'].writeInput(consistencyChecking=OFF)
#: The job input file has been written to "mesh_005.inp".
mdb.jobs['mesh_005'].submit(consistencyChecking=OFF)
#: The job input file "mesh_005.inp" has been submitted for analysis.
#: Job mesh_005: Analysis Input File Processor completed successfully.
#: Job mesh_005: Abaqus/Standard completed successfully.
#: Job mesh_005 completed successfully. 
session.viewports['Viewport: 1'].setValues(
    displayedObject=session.odbs['C:/Users/Aljaz 2/Desktop/mesh_01.odb'])
o3 = session.openOdb(name='C:/Users/Aljaz 2/Desktop/polkrog_cev/mesh_005.odb')
#: Model: C:/Users/Aljaz 2/Desktop/polkrog_cev/mesh_005.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       1
#: Number of Node Sets:          2
#: Number of Steps:              1
session.viewports['Viewport: 1'].setValues(displayedObject=o3)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='NT11', outputPosition=NODAL, )
session.viewports['Viewport: 1'].odbDisplay.display.setValues(
    plotState=CONTOURS_ON_DEF)
session.viewports['Viewport: 1'].view.setValues(session.views['User-1'])
session.printToFile(fileName='mesh_005', format=PS, canvasObjects=(
    session.viewports['Viewport: 1'], ))
mdb.save()
#: The model database has been saved to "C:\Users\Aljaz 2\Desktop\modelska_cev.cae".
mdb.save()
#: The model database has been saved to "C:\Users\Aljaz 2\Desktop\modelska_cev.cae".
mdb.save()
#: The model database has been saved to "C:\Users\Aljaz 2\Desktop\modelska_cev.cae".
