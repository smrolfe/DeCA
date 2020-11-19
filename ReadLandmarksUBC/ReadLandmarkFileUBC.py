import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
import numpy as np

#
# ReadLandmarkFileUBC
#

class ReadLandmarkFileUBC(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "ReadLandmarkFileUBC" # TODO make this more human readable by adding spaces
    self.parent.categories = ["DeCA Toolbox"]
    self.parent.dependencies = []
    self.parent.contributors = ["Sara Rolfe (UW)"] # replace with "Firstname Lastname (Organization)"
    self.parent.helpText = """
This module imports a directory of CSV files containing landmark points and exports as Slicer format FCSV files.
"""
    self.parent.helpText += self.getDefaultModuleDocumentationLink()
    self.parent.acknowledgementText = """

""" # replace with organization, grant and thanks.
       

#
# ReadLandmarkFileUBCWidget
#

class ReadLandmarkFileUBCWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setup(self):
    ScriptedLoadableModuleWidget.setup(self)

    # Instantiate and connect widgets ...

    #
    # Parameters Area
    #
    parametersCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersCollapsibleButton.text = "Parameters"
    self.layout.addWidget(parametersCollapsibleButton)

    # Layout within the dummy collapsible button
    parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)

    #
    # Select landmark folder to import
    #
    self.inputFileSelector = ctk.ctkPathLineEdit()
    self.inputFileSelector.filters = ctk.ctkPathLineEdit.Dirs
    self.inputFileSelector.setToolTip( "Select folder containing landmark names and coordinates for import." )
    parametersFormLayout.addRow("Select landmark folder for import: ", self.inputFileSelector)
    
    #
    # Select landmark folder for export
    #
    self.outputFileSelector = ctk.ctkPathLineEdit()
    self.outputFileSelector.filters = ctk.ctkPathLineEdit.Dirs
    self.outputFileSelector.setToolTip( "Select folder where the FCSV files will be written." )
    parametersFormLayout.addRow("Select destination folder for export: ", self.outputFileSelector)
    
    #
    # Select landmark numbers to export
    #
    self.startPoint = ctk.ctkDoubleSpinBox()
    self.startPoint.minimum = 0
    self.startPoint.value = 0
    self.startPoint.singleStep = 1
    self.startPoint.setDecimals(0)
    self.startPoint.setToolTip("First landmark number to import:")
    parametersFormLayout.addRow("Start landmark number: ", self.startPoint)
    
    self.stopPoint = ctk.ctkDoubleSpinBox()
    self.stopPoint.minimum = 0
    self.startPoint.value = 0
    self.stopPoint.singleStep = 1
    self.stopPoint.setDecimals(0)
    self.stopPoint.setToolTip("Last landmark number to import")
    parametersFormLayout.addRow("Stop landmark number: ", self.stopPoint)
    
    #
    # Apply Button
    #
    self.applyButton = qt.QPushButton("Run")
    self.applyButton.toolTip = "Convert the landmarks."
    self.applyButton.enabled = False
    parametersFormLayout.addRow(self.applyButton)

    # connections
    self.applyButton.connect('clicked(bool)', self.onApplyButton)
    self.inputFileSelector.connect('validInputChanged(bool)', self.onSelect)
    self.outputFileSelector.connect('validInputChanged(bool)', self.onSelect)

    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh Apply button state
    self.onSelect()

  def cleanup(self):
    pass

  def onSelect(self):
    self.applyButton.enabled = bool(self.inputFileSelector.currentPath) and bool(self.outputFileSelector.currentPath)

  def onApplyButton(self):
    logic = ReadLandmarkFileUBCLogic()
    logic.run(self.inputFileSelector.currentPath, self.outputFileSelector.currentPath, 
    int(self.startPoint.value), int(self.stopPoint.value))

#
# ReadLandmarkFileUBCLogic
#

class ReadLandmarkFileUBCLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def run(self, inputDirectory, outputDirectory, startValue, stopValue):
    """
    Run the actual algorithm
    """
    extensionInput = ".csv"
    extensionOutput = ".fcsv"
    for file in os.listdir(inputDirectory):
      if file.endswith(extensionInput):
        inputFilePath = os.path.join(inputDirectory, file)  
        (landmarkFileBase, ext) = os.path.splitext(file)
        array = np.genfromtxt(inputFilePath, delimiter=',')  
        
        # Create a markups node for imported points
        fiducialNode = slicer.vtkMRMLMarkupsFiducialNode()
        slicer.mrmlScene.AddNode(fiducialNode)
        for i in range(startValue-1, stopValue):
          point = array[:,i]
          pointName = landmarkFileBase + '_' + str(i+1)
          fiducialNode.AddFiducialFromArray([-point[0],-point[1],point[2]],pointName)
    
        # save output
        outputFilePath = os.path.join(outputDirectory, landmarkFileBase + extensionOutput)
        slicer.util.saveNode(fiducialNode, outputFilePath)
        logging.info('Processing completed for ' + landmarkFileBase)
        
        # clean up
        slicer.mrmlScene.RemoveNode(fiducialNode)

    return True

class ReadLandmarkFileUBCTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_ReadLandmarkFileUBC1()

  def test_ReadLandmarkFileUBC1(self):
    """ Ideally you should have several levels of tests.  At the lowest level
    tests should exercise the functionality of the logic with different inputs
    (both valid and invalid).  At higher levels your tests should emulate the
    way the user would interact with your code and confirm that it still works
    the way you intended.
    One of the most important features of the tests is that it should alert other
    developers when their changes will have an impact on the behavior of your
    module.  For example, if a developer removes a feature that you depend on,
    your test should break so they know that the feature is needed.
    """

    self.delayDisplay("Starting the test")
    #
    # first, get some data
    #
    import urllib
    downloads = (
        ('http://slicer.kitware.com/midas3/download?items=5767', 'FA.nrrd', slicer.util.loadVolume),
        )

    for url,name,loader in downloads:
      filePath = slicer.app.temporaryPath + '/' + name
      if not os.path.exists(filePath) or os.stat(filePath).st_size == 0:
        logging.info('Requesting download %s from %s...\n' % (name, url))
        urllib.urlretrieve(url, filePath)
      if loader:
        logging.info('Loading %s...' % (name,))
        loader(filePath)
    self.delayDisplay('Finished with download and loading')

    volumeNode = slicer.util.getNode(pattern="FA")
    logic = ReadLandmarkFileUBCLogic()
    self.assertIsNotNone( logic.hasImageData(volumeNode) )
    self.delayDisplay('Test passed!')
