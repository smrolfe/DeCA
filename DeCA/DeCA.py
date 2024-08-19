import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
import fnmatch
import  numpy as np
import random
import math
from pathlib import Path
import re
import csv
import vtk.util.numpy_support as vtk_np

#
# DeCA
#

class DeCA(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "DeCA" # TODO make this more human readable by adding spaces
    self.parent.categories = ["DeCA Toolbox"]
    self.parent.dependencies = []
    self.parent.contributors = ["Sara Rolfe (UW)"] # replace with "Firstname Lastname (Organization)"
    self.parent.helpText = """
      This module provides several flexible workflows for finding and analyzing dense correspondence points between models.
      """
    self.parent.helpText += self.getDefaultModuleDocumentationLink()
    self.parent.acknowledgementText = """

      """ # replace with organization, grant and thanks.

#
# DeCAWidget
#

class DeCAWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """
  #def onSelect(self):
  #  self.applyButton.enabled = bool (self.meshDirectory.currentPath and self.landmarkDirectory.currentPath and
  #    self.baseMeshSelector.currentNode() and self.baseLMSelector.currentNode() and self.baseSLMSelect.currentNode()
  #    and self.outputDirectory.currentPath)

  def setup(self):
    ScriptedLoadableModuleWidget.setup(self)

    # Set up tabs to split workflow
    tabsWidget = qt.QTabWidget()
    alignTab = qt.QWidget()
    alignTabLayout = qt.QFormLayout(alignTab)
    generateMeanTab = qt.QWidget()
    generateMeanTabLayout = qt.QFormLayout(generateMeanTab)
    DeCATab = qt.QWidget()
    DeCATabLayout = qt.QFormLayout(DeCATab)
    DeCALTab = qt.QWidget()
    DeCALTabLayout = qt.QFormLayout(DeCALTab)
    visualizeTab = qt.QWidget()
    visualizeTabLayout = qt.QFormLayout(visualizeTab)
    symmetryTab = qt.QWidget()
    symmetryTabLayout = qt.QFormLayout(symmetryTab)

    tabsWidget.addTab(alignTab, "Rigid Alignment")
    tabsWidget.addTab(generateMeanTab, "Generate Mean")
    tabsWidget.addTab(symmetryTab, "Mirror Data")
    tabsWidget.addTab(DeCATab, "DeCA")
    tabsWidget.addTab(DeCALTab, "DeCAL")
    tabsWidget.addTab(visualizeTab, "Visualize Results")

    self.layout.addWidget(tabsWidget)

    ################################### Align Tab ###################################
    # Layout within the tab
    alignWidget=ctk.ctkCollapsibleButton()
    alignWidgetLayout = qt.QFormLayout(alignWidget)
    alignWidget.text = "Align all samples to base sample"
    alignTabLayout.addRow(alignWidget)

    #
    # Select base mesh
    #
    self.baseMeshSelector = ctk.ctkPathLineEdit()
    self.baseMeshSelector.filters  = ctk.ctkPathLineEdit().Files
    self.baseMeshSelector.nameFilters=["Model (*.ply *.stl *.obj *.vtk *.vtp)"]
    alignWidgetLayout.addRow("Base model: ", self.baseMeshSelector)

    #
    # Select base landmarks
    #
    self.baseLMSelector = ctk.ctkPathLineEdit()
    self.baseLMSelector.filters  = ctk.ctkPathLineEdit().Files
    self.baseLMSelector.nameFilters=["Point set (*.fcsv *.json *.mrk.json)"]
    alignWidgetLayout.addRow("Base landmarks: ", self.baseLMSelector)

    #
    # Select meshes directory
    #
    self.meshDirectory=ctk.ctkPathLineEdit()
    self.meshDirectory.filters = ctk.ctkPathLineEdit.Dirs
    self.meshDirectory.setToolTip( "Select directory containing meshes" )
    alignWidgetLayout.addRow("Mesh directory: ", self.meshDirectory)

    #
    # Select landmarks directory
    #
    self.landmarkDirectory=ctk.ctkPathLineEdit()
    self.landmarkDirectory.filters = ctk.ctkPathLineEdit.Dirs
    self.landmarkDirectory.setToolTip( "Select directory containing landmarks" )
    alignWidgetLayout.addRow("Landmark directory: ", self.landmarkDirectory)

    #
    # Select aligned mesh directory
    #
    self.alignedMeshDirectory=ctk.ctkPathLineEdit()
    self.alignedMeshDirectory.filters = ctk.ctkPathLineEdit.Dirs
    self.alignedMeshDirectory.setToolTip( "Select directory for aligned meshes: " )
    alignWidgetLayout.addRow("Aligned mesh directory: ", self.alignedMeshDirectory)

    #
    # Select aligned landmark directory
    #
    self.alignedLMDirectory=ctk.ctkPathLineEdit()
    self.alignedLMDirectory.filters = ctk.ctkPathLineEdit.Dirs
    self.alignedLMDirectory.setToolTip( "Select directory for aligned landmarks: " )
    alignWidgetLayout.addRow("Aligned landmark directory: ", self.alignedLMDirectory)

    #
    # Remove scale option
    #
    self.removeScaleCheckBox = qt.QCheckBox()
    self.removeScaleCheckBox.checked = False
    self.removeScaleCheckBox.setToolTip("If checked, alignment will include isotropic scaling.")
    alignWidgetLayout.addRow("Remove scale: ", self.removeScaleCheckBox)

    #
    # Hidden semilandmark options
    self.semilandmarkCollapsibleButton = ctk.ctkCollapsibleButton()
    self.semilandmarkCollapsibleButton.text = "Semi-landmark Options"
    self.semilandmarkCollapsibleButton.collapsed = True
    self.semilandmarkCollapsibleButton.enabled = True
    alignWidgetLayout.addRow(self.semilandmarkCollapsibleButton)
    semilandmarkOptionLayout = qt.QFormLayout(self.semilandmarkCollapsibleButton)

    self.semilandmarkSelector = ctk.ctkPathLineEdit()
    self.semilandmarkSelector.filters  = ctk.ctkPathLineEdit().Dirs
    self.semilandmarkSelector.enabled = True
    self.semilandmarkSelector.setToolTip( "Select directory with semi-landmark files" )
    semilandmarkOptionLayout.addRow("Semi-landmark directory: ", self.semilandmarkSelector)

    self.alignedSemilandmarkSelector = ctk.ctkPathLineEdit()
    self.alignedSemilandmarkSelector.filters  = ctk.ctkPathLineEdit().Dirs
    self.alignedSemilandmarkSelector.enabled = True
    self.alignedSemilandmarkSelector.setToolTip( "Select directory to output aligned semilandmark files" )
    semilandmarkOptionLayout.addRow("Aligned semi-landmark directory: ", self.alignedSemilandmarkSelector)

    #
    # Apply Button
    #
    self.applyButton = qt.QPushButton("Run alignment")
    self.applyButton.toolTip = "Run alignment to base subject"
    self.applyButton.enabled = False
    alignWidgetLayout.addRow(self.applyButton)

    # connections
    self.baseMeshSelector.connect('validInputChanged(bool)', self.onSelectAlignPaths)
    self.baseLMSelector.connect('validInputChanged(bool)', self.onSelectAlignPaths)
    self.meshDirectory.connect('validInputChanged(bool)', self.onSelectAlignPaths)
    self.landmarkDirectory.connect('validInputChanged(bool)', self.onSelectAlignPaths)
    self.alignedMeshDirectory.connect('validInputChanged(bool)', self.onSelectAlignPaths)
    self.alignedLMDirectory.connect('validInputChanged(bool)', self.onSelectAlignPaths)
    self.applyButton.connect('clicked(bool)', self.onApplyAlignmentButton)

    ################################### Mean Tab ###################################
    # Layout within Mean tab
    MeanWidget=ctk.ctkCollapsibleButton()
    MeanWidgetLayout = qt.QFormLayout(MeanWidget)
    MeanWidget.text = "Generate mean mesh"
    generateMeanTabLayout.addRow(MeanWidget)

    #
    # Select aligned mesh directory
    #
    self.meanMeshDirectory=ctk.ctkPathLineEdit()
    self.meanMeshDirectory.filters = ctk.ctkPathLineEdit.Dirs
    self.meanMeshDirectory.setToolTip( "Select directory for aligned meshes: " )
    MeanWidgetLayout.addRow("Aligned mesh directory: ", self.meanMeshDirectory)

    #
    # Select aligned landmark directory
    #
    self.meanLMDirectory=ctk.ctkPathLineEdit()
    self.meanLMDirectory.filters = ctk.ctkPathLineEdit.Dirs
    self.meanLMDirectory.setToolTip( "Select directory for aligned landmarks: " )
    MeanWidgetLayout.addRow("Aligned landmark directory: ", self.meanLMDirectory)

    #
    # Select output directory
    #
    self.meanOutputDirectory=ctk.ctkPathLineEdit()
    self.meanOutputDirectory.filters = ctk.ctkPathLineEdit.Dirs
    self.meanOutputDirectory.setToolTip( "Select directory for mean output" )
    MeanWidgetLayout.addRow("Mean output directory: ", self.meanOutputDirectory)

    #
    # Generate mean button
    #
    self.generateMeanButton = qt.QPushButton("Generate mean")
    self.generateMeanButton.toolTip = "Generate mean template from aligned subjects"
    self.generateMeanButton.enabled = False
    MeanWidgetLayout.addRow(self.generateMeanButton)

    #
    # connections
    #
    self.meanMeshDirectory.connect('validInputChanged(bool)', self.onMeanSelect)
    self.meanLMDirectory.connect('validInputChanged(bool)', self.onMeanSelect)
    self.meanOutputDirectory.connect('validInputChanged(bool)', self.onMeanSelect)
    self.generateMeanButton.connect('clicked(bool)', self.onGenerateMean)

    ################################### DeCA Tab ###################################
    # Layout within the DeCA tab
    DeCAWidget=ctk.ctkCollapsibleButton()
    DeCAWidgetLayout = qt.QFormLayout(DeCAWidget)
    DeCAWidget.text = "Dense Correspondence I/0"
    DeCATabLayout.addRow(DeCAWidget)

    #
    # Select base mesh
    #
    self.DCBaseModelSelector = ctk.ctkPathLineEdit()
    self.DCBaseModelSelector.filters  = ctk.ctkPathLineEdit().Files
    self.DCBaseModelSelector.nameFilters=["Model (*.ply *.stl *.obj *.vtk *.vtp)"]
    DeCAWidgetLayout.addRow("Base mesh: ", self.DCBaseModelSelector)

    #
    # Select base landmarks
    #
    self.DCBaseLMSelector = ctk.ctkPathLineEdit()
    self.DCBaseLMSelector.filters  = ctk.ctkPathLineEdit().Files
    self.DCBaseLMSelector.nameFilters=["Point set (*.fcsv *.json *.mrk.json"]
    DeCAWidgetLayout.addRow("Base landmarks: ", self.DCBaseLMSelector)

    #
    # Select rigidly aligned meshes directory
    #
    self.DCMeshDirectory=ctk.ctkPathLineEdit()
    self.DCMeshDirectory.filters = ctk.ctkPathLineEdit.Dirs
    self.DCMeshDirectory.setToolTip( "Select directory containing rigidly aligned meshes" )
    DeCAWidgetLayout.addRow("Rigidly Aligned Mesh directory: ", self.DCMeshDirectory)

    #
    # Select rigidly aligned landmark directory
    #
    self.DCLandmarkDirectory=ctk.ctkPathLineEdit()
    self.DCLandmarkDirectory.filters = ctk.ctkPathLineEdit.Dirs
    self.DCLandmarkDirectory.setToolTip( "Select directory containing rigidly aligned landmarks" )
    DeCAWidgetLayout.addRow("Rigidly Aligned Landmark directory: ", self.DCLandmarkDirectory)

    #
    # Select DeCA output directory
    #
    self.DCOutputDirectory=ctk.ctkPathLineEdit()
    self.DCOutputDirectory.filters = ctk.ctkPathLineEdit.Dirs
    self.DCOutputDirectory.setToolTip( "Select directory for DeCA output" )
    DeCAWidgetLayout.addRow("DeCA output directory: ", self.DCOutputDirectory)

    #
    # Use CPD registration
    #
    self.CPDCheckBox = qt.QCheckBox()
    self.CPDCheckBox.checked = False
    self.CPDCheckBox.enabled = False
    self.CPDCheckBox.setToolTip("If checked, DeCA will use CPD for point correspondences.")
    #DeCAWidgetLayout.addRow("Use CPD registration: ", self.CPDCheckBox)

    #
    # Write directory for error checking
    #
    self.WriteErrorCheckBox = qt.QCheckBox()
    self.WriteErrorCheckBox.checked = False
    self.WriteErrorCheckBox.setToolTip("If checked, DeCA will create a directory of results for use in estimating point correspondence error.")
    DeCAWidgetLayout.addRow("Create output for error checking: ", self.WriteErrorCheckBox)

    #
    # Write directory for error checking
    #
    self.WriteCorrPointsCheckBox= qt.QCheckBox()
    self.WriteCorrPointsCheckBox.checked = False
    self.WriteCorrPointsCheckBox.setToolTip("If checked, DeCA will write out corresponding meshes.")
    #DeCAWidgetLayout.addRow("Write out point correspondences: ", self.WriteCorrPointsCheckBox)

    #
    # Select Analysis Type
    #
    self.analysisTypeShape=qt.QRadioButton()
    self.analysisTypeShape.setChecked(True)
    DeCAWidgetLayout.addRow("Shape analysis: ", self.analysisTypeShape)
    self.analysisTypeSymmetry=qt.QRadioButton()
    DeCAWidgetLayout.addRow("Symmetry analysis: ", self.analysisTypeSymmetry)

    #
    # Hidden symmetry options
    self.symmetryCollapsibleButton = ctk.ctkCollapsibleButton()
    self.symmetryCollapsibleButton.text = "Symmetry Options"
    self.symmetryCollapsibleButton.collapsed = True
    self.symmetryCollapsibleButton.enabled = False
    DeCAWidgetLayout.addRow(self.symmetryCollapsibleButton)
    symmetryOptionLayout = qt.QFormLayout(self.symmetryCollapsibleButton)

    self.mirrorMeshSelector = ctk.ctkPathLineEdit()
    self.mirrorMeshSelector.filters  = ctk.ctkPathLineEdit().Dirs
    self.mirrorMeshSelector.enabled = False
    self.mirrorMeshSelector.setToolTip( "Select directory with mirrored mesh data" )
    symmetryOptionLayout.addRow("Mirror mesh directory: ", self.mirrorMeshSelector)

    self.mirrorLMSelector = ctk.ctkPathLineEdit()
    self.mirrorLMSelector.filters  = ctk.ctkPathLineEdit().Dirs
    self.mirrorLMSelector.enabled = False
    self.mirrorLMSelector.setToolTip( "Select directory with mirrored landmark data" )
    symmetryOptionLayout.addRow("Mirror landmark directory: ", self.mirrorLMSelector)

    #
    # Run DeCA Button
    #
    self.DCApplyButton = qt.QPushButton("Run DeCA")
    self.DCApplyButton.toolTip = "Run non-rigid alignment"
    self.DCApplyButton.enabled = False
    DeCAWidgetLayout.addRow(self.DCApplyButton)

    # Connections
    self.DCBaseModelSelector.connect('validInputChanged(bool)', self.onDCSelect)
    self.DCBaseLMSelector.connect('validInputChanged(bool)', self.onDCSelect)
    self.DCMeshDirectory.connect('validInputChanged(bool)', self.onDCSelect)
    self.DCLandmarkDirectory.connect('validInputChanged(bool)', self.onDCSelect)
    self.DCOutputDirectory.connect('validInputChanged(bool)', self.onDCSelect)
    self.analysisTypeShape.connect('toggled(bool)', self.onToggleAnalysis)
    self.analysisTypeSymmetry.connect('toggled(bool)', self.onToggleAnalysis)
    self.mirrorMeshSelector.connect('validInputChanged(bool)', self.onDCSelect)
    self.mirrorLMSelector.connect('validInputChanged(bool)', self.onDCSelect)
    self.DCApplyButton.connect('clicked(bool)', self.onDCApplyButton)

    ################################### Visualize Tab ###################################
    # Layout within the tab
    visualizeWidget=ctk.ctkCollapsibleButton()
    visualizeWidgetLayout = qt.QFormLayout(visualizeWidget)
    visualizeWidget.text = "Visualize the output feature heat maps"
    visualizeTabLayout.addRow(visualizeWidget)

    #
    # Select base mesh
    #
    self.meshSelect = slicer.qMRMLNodeComboBox()
    self.meshSelect.nodeTypes = ( ("vtkMRMLModelNode"), "" )
    self.meshSelect.setToolTip( "Select model node with result arrays" )
    self.meshSelect.selectNodeUponCreation = False
    self.meshSelect.noneEnabled = True
    self.meshSelect.addEnabled = False
    self.meshSelect.removeEnabled = False
    self.meshSelect.showHidden = False
    self.meshSelect.setMRMLScene( slicer.mrmlScene )
    visualizeWidgetLayout.addRow("Result Model: ", self.meshSelect)

    #
    # Select Subject ID
    #
    self.subjectIDBox=qt.QComboBox()
    self.subjectIDBox.enabled = False
    visualizeWidgetLayout.addRow("Subject ID: ", self.subjectIDBox)

    # Connections
    self.meshSelect.connect("currentNodeChanged(vtkMRMLNode*)", self.onVisualizeMeshSelect)
    self.subjectIDBox.connect("currentIndexChanged(int)", self.onSubjectIDSelect)

    # Add vertical spacer
    self.layout.addStretch(1)

    ################################### Symmetry Tab ###################################
    # Layout within the tab
    symmetryWidget=ctk.ctkCollapsibleButton()
    symmetryWidgetLayout = qt.QFormLayout(symmetryWidget)
    symmetryWidget.text = "Generate mirrored data for point-wise symmetry analysis"
    symmetryTabLayout.addRow(symmetryWidget)

    #
    # Select meshes directory
    #
    self.symMeshDirectory=ctk.ctkPathLineEdit()
    self.symMeshDirectory.filters = ctk.ctkPathLineEdit.Dirs
    self.symMeshDirectory.setToolTip( "Select directory containing aligned meshes" )
    symmetryWidgetLayout.addRow("Aligned mesh directory: ", self.symMeshDirectory)

    #
    # Select landmarks directory
    #
    self.symLandmarkDirectory=ctk.ctkPathLineEdit()
    self.symLandmarkDirectory.filters = ctk.ctkPathLineEdit.Dirs
    self.symLandmarkDirectory.setToolTip( "Select directory containing aligned landmarks" )
    symmetryWidgetLayout.addRow("Aligned landmark directory: ", self.symLandmarkDirectory)

    #
    # Select aligned mesh directory
    #
    self.mirrorMeshDirectory=ctk.ctkPathLineEdit()
    self.mirrorMeshDirectory.filters = ctk.ctkPathLineEdit.Dirs
    self.mirrorMeshDirectory.setToolTip( "Select output directory for mirrored meshes: " )
    symmetryWidgetLayout.addRow("Ouput mirror mesh directory: ", self.mirrorMeshDirectory)

    #
    # Select aligned landmark directory
    #
    self.mirrorLMDirectory=ctk.ctkPathLineEdit()
    self.mirrorLMDirectory.filters = ctk.ctkPathLineEdit.Dirs
    self.mirrorLMDirectory.setToolTip( "Select output directory for mirrored landmarks: " )
    symmetryWidgetLayout.addRow("Output mirror landmark directory: ", self.mirrorLMDirectory)

    #
    # Axis selection
    #
    self.xAxis=qt.QRadioButton()
    self.xAxis.setChecked(True)
    symmetryWidgetLayout.addRow('X-Axis', self.xAxis)
    self.yAxis=qt.QRadioButton()
    self.yAxis.setChecked(False)
    symmetryWidgetLayout.addRow('Y-Axis', self.yAxis)
    self.zAxis=qt.QRadioButton()
    self.zAxis.setChecked(False)
    symmetryWidgetLayout.addRow('Z-Axis', self.zAxis)

    #
    # Mirror Landmark Index
    #
    self.landmarkIndexText=qt.QLineEdit()
    self.landmarkIndexText.setToolTip("No spaces. Seperate numbers by commas.  Example:  2,1,3,5,4")
    symmetryWidgetLayout.addRow('Mirror landmark index', self.landmarkIndexText)

    #
    # Hidden semilandmark options
    #
    self.semilandmarkMirrorCollapsibleButton = ctk.ctkCollapsibleButton()
    self.semilandmarkMirrorCollapsibleButton.text = "Semi-landmark Options"
    self.semilandmarkMirrorCollapsibleButton.collapsed = True
    self.semilandmarkMirrorCollapsibleButton.enabled = True
    symmetryWidgetLayout.addRow(self.semilandmarkMirrorCollapsibleButton)
    semilandmarkMirrorOptionLayout = qt.QFormLayout(self.semilandmarkMirrorCollapsibleButton)

    self.semilandmarkMirrorSelector = ctk.ctkPathLineEdit()
    self.semilandmarkMirrorSelector.filters  = ctk.ctkPathLineEdit().Dirs
    self.semilandmarkMirrorSelector.enabled = True
    self.semilandmarkMirrorSelector.setToolTip( "Select directory with semi-landmark files" )
    semilandmarkMirrorOptionLayout.addRow("Semi-landmark directory: ", self.semilandmarkMirrorSelector)

    self.alignedSemilandmarkMirrorSelector = ctk.ctkPathLineEdit()
    self.alignedSemilandmarkMirrorSelector.filters  = ctk.ctkPathLineEdit().Dirs
    self.alignedSemilandmarkMirrorSelector.enabled = True
    self.alignedSemilandmarkMirrorSelector.setToolTip( "Select directory to output aligned semilandmark files" )
    semilandmarkMirrorOptionLayout.addRow("Mirrored semi-landmark directory: ", self.alignedSemilandmarkMirrorSelector)

    self.semilandmarkIndexMirrorText=qt.QLineEdit()
    self.semilandmarkIndexMirrorText.setToolTip("No spaces. Seperate numbers by commas.  Example:  2,1,3,5,4")
    semilandmarkMirrorOptionLayout.addRow('Mirror semi-landmark index', self.semilandmarkIndexMirrorText)

    #
    # Apply Button
    #
    self.mirrorButton = qt.QPushButton("Mirror")
    self.mirrorButton.toolTip = "Generate mirrored and aligned copies of models and landmarks"
    self.mirrorButton.enabled = False
    symmetryWidgetLayout.addRow(self.mirrorButton)

    # connections
    self.symMeshDirectory.connect('validInputChanged(bool)', self.onMirrorSelect)
    self.symLandmarkDirectory.connect('validInputChanged(bool)', self.onMirrorSelect)
    self.mirrorMeshDirectory.connect('validInputChanged(bool)', self.onMirrorSelect)
    self.mirrorLMDirectory.connect('validInputChanged(bool)', self.onMirrorSelect)
    self.mirrorButton.connect('clicked(bool)', self.onMirrorButton)

    ################################### DeCAL Tab ###################################
    # Layout within the DeCA tab
    DeCALWidget=ctk.ctkCollapsibleButton()
    DeCALWidgetLayout = qt.QFormLayout(DeCALWidget)
    DeCALWidget.text = "Dense Correspondence Landmarking"
    DeCALTabLayout.addRow(DeCALWidget)

    #
    # Select base mesh
    #
    self.DCLBaseModelSelector = ctk.ctkPathLineEdit()
    self.DCLBaseModelSelector.filters  = ctk.ctkPathLineEdit().Files
    self.DCLBaseModelSelector.nameFilters=["Model (*.ply *.stl *.obj *.vtk *.vtp)"]
    DeCALWidgetLayout.addRow("Base mesh: ", self.DCLBaseModelSelector)

    #
    # Select base landmarks
    #
    self.DCLBaseLMSelector = ctk.ctkPathLineEdit()
    self.DCLBaseLMSelector.filters  = ctk.ctkPathLineEdit().Files
    self.DCLBaseLMSelector.nameFilters=["Point set (*.fcsv *.json *.mrk.json"]
    DeCALWidgetLayout.addRow("Base landmarks: ", self.DCLBaseLMSelector)

    #
    # Select rigidly aligned meshes directory
    #
    self.DCLMeshDirectory=ctk.ctkPathLineEdit()
    self.DCLMeshDirectory.filters = ctk.ctkPathLineEdit.Dirs
    self.DCLMeshDirectory.setToolTip( "Select directory containing rigidly aligned meshes" )
    DeCALWidgetLayout.addRow("Rigidly Aligned Mesh directory: ", self.DCLMeshDirectory)

    #
    # Select rigidly aligned landmark directory
    #
    self.DCLLandmarkDirectory=ctk.ctkPathLineEdit()
    self.DCLLandmarkDirectory.filters = ctk.ctkPathLineEdit.Dirs
    self.DCLLandmarkDirectory.setToolTip( "Select directory containing rigidly aligned landmarks" )
    DeCALWidgetLayout.addRow("Rigidly Aligned Landmark directory: ", self.DCLLandmarkDirectory)

    #
    # Select DeCA output directory
    #
    self.DCLOutputDirectory=ctk.ctkPathLineEdit()
    self.DCLOutputDirectory.filters = ctk.ctkPathLineEdit.Dirs
    self.DCLOutputDirectory.setToolTip( "Select directory for DeCA output" )
    DeCALWidgetLayout.addRow("DeCA output directory: ", self.DCLOutputDirectory)

    #
    # Set spacing tolerance
    #
    self.spacingTolerance = ctk.ctkSliderWidget()
    self.spacingTolerance.singleStep = .1
    self.spacingTolerance.minimum = 0
    self.spacingTolerance.maximum = 10
    self.spacingTolerance.value = 4
    self.spacingTolerance.setToolTip("Set tolerance of spacing as a percentage of the image diagonal")
    DeCALWidgetLayout.addRow("Spacing tolerance: ", self.spacingTolerance)

    #
    # Get Subsample Rate Button
    #
    self.getPointNumberButton = qt.QPushButton("Get subsample number")
    self.getPointNumberButton.toolTip = "Output the number of points sampled on the template shape."
    self.getPointNumberButton.enabled = False
    DeCALWidgetLayout.addRow(self.getPointNumberButton)

    #
    # Subsample Information
    #
    self.subsampleInfo = qt.QPlainTextEdit()
    self.subsampleInfo.setPlaceholderText("Subsampling information")
    self.subsampleInfo.setReadOnly(True)
    DeCALWidgetLayout.addRow(self.subsampleInfo)

    #
    # Apply Button
    #
    self.DCLApplyButton = qt.QPushButton("Run DeCAL")
    self.DCLApplyButton.toolTip = "Generate a set of corresponding landmarks"
    self.DCLApplyButton.enabled = False
    DeCALWidgetLayout.addRow(self.DCLApplyButton)

    # connections
    self.DCLBaseModelSelector.connect('validInputChanged(bool)', self.testConnect)
    self.DCLBaseLMSelector.connect('validInputChanged(bool)', self.onDCLSelect)
    self.DCLMeshDirectory.connect('validInputChanged(bool)', self.onDCLSelect)
    self.DCLLandmarkDirectory.connect('validInputChanged(bool)', self.onDCLSelect)
    self.DCLOutputDirectory.connect('validInputChanged(bool)', self.onDCLSelect)
    self.DCLApplyButton.connect('clicked(bool)', self.onDCLApplyButton)
    self.getPointNumberButton.connect('clicked(bool)', self.onGetPointNumberButton)

### GUI SUpport Functions
  def testConnect(self):
    print("Got it!!")

  def onToggleAnalysis(self):
    if self.analysisTypeSymmetry.checked == True:
      self.mirrorMeshSelector.enabled = True
      self.mirrorLMSelector.enabled = True
      self.symmetryCollapsibleButton.collapsed = False
      self.symmetryCollapsibleButton.enabled = True
    else:
      self.mirrorMeshSelector.enabled = False
      self.mirrorLMSelector.enabled = False
      self.symmetryCollapsibleButton.collapsed = True
      self.symmetryCollapsibleButton.enabled = False

  def onSubjectIDSelect(self):
    try:
      subjectID = self.subjectIDBox.currentText
      self.resultNode.GetDisplayNode().SetActiveScalarName(subjectID)
      self.resultNode.GetDisplayNode().SetAndObserveColorNodeID('vtkMRMLColorTableNodeFilePlasma.txt')
      print(subjectID)
    except:
      print("Error: No array found")

  def onVisualizeMeshSelect(self):
    if bool(self.meshSelect.currentNode()):
      self.resultNode = self.meshSelect.currentNode()
      self.resultNode.GetDisplayNode().SetVisibility(True)
      self.resultNode.GetDisplayNode().SetScalarVisibility(True)
      resultData = self.resultNode.GetPolyData().GetPointData()
      self.subjectIDBox.enabled = True
      arrayNumber = resultData.GetNumberOfArrays()
      if arrayNumber > 0:
        for i in range(resultData.GetNumberOfArrays()):
          arrayName = resultData.GetArrayName(i)
          self.subjectIDBox.addItem(arrayName)
      else:
        self.subjectIDBox.clear()
        self.subjectIDBox.enabled = False

  def cleanup(self):
    pass

  def onSelectAlignPaths(self):
    self.applyButton.enabled = bool (self.meshDirectory.currentPath and self.landmarkDirectory.currentPath and
      self.baseMeshSelector.currentPath and self.baseLMSelector.currentPath
      and self.alignedMeshDirectory.currentPath and self.alignedLMDirectory.currentPath)

    self.meanMeshDirectory.currentPath = self.alignedMeshDirectory.currentPath
    self.meanLMDirectory.currentPath = self.alignedLMDirectory.currentPath

    self.DCMeshDirectory.currentPath = self.alignedMeshDirectory.currentPath
    self.DCLandmarkDirectory.currentPath = self.alignedLMDirectory.currentPath

  def onMirrorSelect(self):
    self.mirrorButton.enabled = bool (self.symMeshDirectory.currentPath and self.symLandmarkDirectory.currentPath
    and self.mirrorMeshDirectory.currentPath and self.mirrorLMDirectory.currentPath)

  def onApplyAlignmentButton(self):
    logic = DeCALogic()
    logic.runAlign(self.baseMeshSelector.currentPath, self.baseLMSelector.currentPath, self.meshDirectory.currentPath,
      self.landmarkDirectory.currentPath, self.alignedMeshDirectory.currentPath, self.alignedLMDirectory.currentPath, self.removeScaleCheckBox.checked,
      self.semilandmarkSelector.currentPath, self.alignedSemilandmarkSelector.currentPath)

  def onMeanSelect(self):
    self.generateMeanButton.enabled = bool (self.meanMeshDirectory.currentPath and self.meanLMDirectory.currentPath and self.meanOutputDirectory)

  def onGenerateMean(self):
    logic = DeCALogic()
    base, modelExt = os.path.splitext(self.baseMeshSelector.currentPath)
    logic.runMean(self.meanLMDirectory.currentPath, self.meanMeshDirectory.currentPath, modelExt, self.meanOutputDirectory.currentPath)

  def onDCApplyButton(self):
    logic = DeCALogic()
    if self.analysisTypeShape.checked == True:
      logic.runDCAlign(self.DCBaseModelSelector.currentPath, self.DCBaseLMSelector.currentPath, self.DCMeshDirectory.currentPath,
      self.DCLandmarkDirectory.currentPath, self.DCOutputDirectory.currentPath, self.CPDCheckBox.checked, self.WriteErrorCheckBox.checked)
    else:
      logic.runDCAlignSymmetric(self.DCBaseModelSelector.currentPath, self.DCBaseLMSelector.currentPath, self.DCMeshDirectory.currentPath,
      self.DCLandmarkDirectory.currentPath, self.mirrorMeshSelector.currentPath, self.mirrorLMSelector.currentPath, self.DCOutputDirectory.currentPath,
      self.CPDCheckBox.checked, self.WriteErrorCheckBox.checked, self.WriteCorrPointsCheckBox.checked)

  def onDCSelect(self):
    if self.analysisTypeShape.checked == True:
      self.DCApplyButton.enabled = bool (self.DCMeshDirectory.currentPath and self.DCLandmarkDirectory.currentPath
      and self.DCOutputDirectory.currentPath and self.DCBaseModelSelector.currentPath and self.DCBaseLMSelector.currentPath)
    else:
      self.DCApplyButton.enabled = bool (self.DCMeshDirectory.currentPath and self.DCLandmarkDirectory.currentPath
      and self.DCOutputDirectory.currentPath and self.DCBaseModelSelector.currentPath and self.DCBaseLMSelector.currentPath
      and self.mirrorMeshSelector.currentPath and self.mirrorLMSelector.currentPath)

  def onDCLSelect(self):
    self.getPointNumberButton.enabled = bool (self.DCLMeshDirectory.currentPath and self.DCLLandmarkDirectory.currentPath
    and self.DCLOutputDirectory.currentPath and self.DCLBaseModelSelector.currentPath and self.DCLBaseLMSelector.currentPath)
    print(bool (self.DCLMeshDirectory.currentPath ))
    print(bool (self.DCLLandmarkDirectory.currentPath) )
    print(bool (self.DCLOutputDirectory.currentPath ) )
    print(bool ( self.DCLBaseLMSelector.currentPath) )

  def onGetPointNumberButton(self):
    logic = DeCALogic()
    self.baseNode, self.templateModel, pointNumber = logic.runCheckPoints(self.DCLBaseModelSelector.currentPath, self.spacingTolerance.value)
    self.subsampleInfo.insertPlainText(f'The subsampled template has a total of {pointNumber} points. \n')
    self.DCLApplyButton.enabled = True

  def onDCLApplyButton(self):
    logic = DeCALogic()
    logic.runDeCAL(self.baseNode, self.templateModel, self.DCLBaseLMSelector.currentPath, self.DCLMeshDirectory.currentPath,
    self.DCLLandmarkDirectory.currentPath, self.DCLOutputDirectory.currentPath, self.spacingTolerance.value)

  def onMirrorButton(self):
    logic = DeCALogic()
    axis = [1,1,1]
    if self.xAxis.isChecked():
      axis[0] = -1
    elif self.yAxis.isChecked():
      axis[1] = -1
    else:
      axis[2] = -1

    logic.runMirroring(self.symMeshDirectory.currentPath, self.symLandmarkDirectory.currentPath, self.mirrorMeshDirectory.currentPath,
      self.mirrorLMDirectory.currentPath, axis, self.landmarkIndexText.text, self.semilandmarkMirrorSelector.currentPath, self.alignedSemilandmarkMirrorSelector.currentPath, self.semilandmarkIndexMirrorText.text)

#
# DeCALogic
#

class DeCALogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
    computation done by your module.  The interface
    should be such that other python code can import
    this class and make use of the functionality without
    requiring an instance of the Widget.
    Uses ScriptedLoadableModuleLogic base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """
  def runCheckPoints(self, baseMeshPath, spacingTolerance):
    baseNode = slicer.util.loadModel(baseMeshPath)
    indexArrayName = "indexArray"
    self.addIndexArray(baseNode, indexArrayName)
    spacingPercentage = spacingTolerance/100
    templateModel = self.downsampleModel(baseNode, spacingPercentage)
    templatePointNumber = templateModel.GetNumberOfPoints()
    #slicer.mrmlScene.RemoveNode(baseNode)
    return baseNode, templateModel, templatePointNumber

  def runDeCAL(self, baseNode, templateModel, baseLMPath, meshDirectory, landmarkDirectory, outputDirectory, spacingPercentage):
    baseLandmarks=self.fiducialNodeToPolyData(baseLMPath).GetPoints()
    modelExt=['ply','stl','vtp']
    self.modelNames, models = self.importMeshes(meshDirectory, modelExt)
    landmarks = self.importLandmarks(landmarkDirectory)
    self.outputDirectory = outputDirectory
    denseCorrespondenceGroup = self.denseCorrespondenceBaseMesh(landmarks, models, baseNode.GetPolyData(), baseLandmarks)

    # get base mesh index
    indexArrayName = "indexArray"
    templateIndex = templateModel.GetPointData().GetArray(indexArrayName)

    # saving point correspondences
    if(templateIndex):
      sampleNumber = denseCorrespondenceGroup.GetNumberOfBlocks()
      for i in range(sampleNumber):
        alignedMesh = denseCorrespondenceGroup.GetBlock(i)
        alignedPointNode= slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode',"alignedPoints")
        for j in range(templateIndex.GetNumberOfValues()):
          baseIndex = templateIndex.GetValue(j)
          alignedPoint = alignedMesh.GetPoint(baseIndex)
          alignedPointNode.AddControlPoint(alignedPoint)
        outputLMPath = os.path.join(outputDirectory, self.modelNames[i]+".mrk.json")
        slicer.util.saveNode(alignedPointNode, outputLMPath)
        slicer.mrmlScene.RemoveNode(alignedPointNode)

      # save base node correspondences
      basePointNode= slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode',"basePoints")
      for j in range(templateIndex.GetNumberOfValues()):
        baseIndex = templateIndex.GetValue(j)
        basePoint = baseNode.GetPolyData().GetPoint(baseIndex)
        basePointNode.AddControlPoint(basePoint)
      baseLMPath = os.path.join(outputDirectory, "baseModel.mrk.json")
      slicer.util.saveNode(basePointNode, baseLMPath)
      slicer.mrmlScene.RemoveNode(basePointNode)

  def downsampleModel(self, model, spacingPercentage):
    points=model.GetPolyData()
    cleanFilter=vtk.vtkCleanPolyData()
    cleanFilter.SetToleranceIsAbsolute(False)
    cleanFilter.SetTolerance(spacingPercentage)
    cleanFilter.SetInputData(points)
    cleanFilter.Update()
    return cleanFilter.GetOutput()

  def addIndexArray(self, mesh, arrayName):
    # Array of original index values
    indexArray = vtk.vtkIntArray()
    indexArray.SetNumberOfComponents(1)
    indexArray.SetName(arrayName)
    for i in range(mesh.GetPolyData().GetNumberOfPoints()):
      indexArray.InsertNextValue(i)
    mesh.GetPolyData().GetPointData().AddArray(indexArray)

  def computeNormals(self, inputModel):
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputData(inputModel.GetPolyData())
    normals.SetAutoOrientNormals(True)
    normals.Update()
    inputModel.SetAndObservePolyData(normals.GetOutput())

  def runMirroring(self, meshDirectory, lmDirectory, mirrorMeshDirectory, mirrorLMDirectory, mirrorAxis, mirrorIndexText, slmDirectory, outputSLMDirectory, mirrorSLMIndexText):
    mirrorMatrix = vtk.vtkMatrix4x4()
    mirrorMatrix.SetElement(0, 0, mirrorAxis[0])
    mirrorMatrix.SetElement(1, 1, mirrorAxis[1])
    mirrorMatrix.SetElement(2, 2, mirrorAxis[2])

    lmFileList = os.listdir(lmDirectory)
    point=[0,0,0]

    #get order of mirrored sets
    if len(mirrorIndexText) != 0:
      mirrorIndexList=mirrorIndexText.split(",")
      mirrorIndexList=[int(x) for x in mirrorIndexList]
      mirrorIndex=np.asarray(mirrorIndexList)
    else:
      print("Error: no landmark index for mirrored mesh")

    semilandmarkOption = bool(slmDirectory and outputSLMDirectory and (len(mirrorSLMIndexText) != 0 ))
    if semilandmarkOption:
      mirrorSLMIndexList=mirrorSLMIndexText.split(",")
      mirrorSLMIndexList=[int(x) for x in mirrorSLMIndexList]
      mirrorSLMIndex=np.asarray(mirrorSLMIndexList)

    for meshFileName in os.listdir(meshDirectory):
      if(not meshFileName.startswith(".")):
        meshFilePath = os.path.join(meshDirectory, meshFileName)
        currentMeshNode = slicer.util.loadModel(meshFilePath)
        subjectID = os.path.splitext(meshFileName)[0]
        currentLMNode = self.getLandmarkFileByID(lmDirectory, subjectID)
        if currentLMNode:
          lmFilePath = os.path.join(lmDirectory, subjectID)
          targetPoints = vtk.vtkPoints()
          for i in range(currentLMNode.GetNumberOfControlPoints()):
            point = currentLMNode.GetNthControlPointPosition(i)
            targetPoints.InsertNextPoint(point)

          mirrorTransform = vtk.vtkTransform()
          mirrorTransform.SetMatrix(mirrorMatrix)
          mirrorTransformNode=slicer.mrmlScene.AddNewNodeByClass("vtkMRMLTransformNode","Mirror")
          mirrorTransformNode.SetAndObserveTransformToParent(mirrorTransform)

          # apply transform to the current surface mesh and landmarks
          currentMeshNode.SetAndObserveTransformNodeID(mirrorTransformNode.GetID())
          currentLMNode.SetAndObserveTransformNodeID(mirrorTransformNode.GetID())
          slicer.vtkSlicerTransformLogic().hardenTransform(currentMeshNode)
          slicer.vtkSlicerTransformLogic().hardenTransform(currentLMNode)

          # apply rigid transformation
          sourcePoints = vtk.vtkPoints()
          mirrorLMNode =slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode",subjectID)
          for i in range(currentLMNode.GetNumberOfControlPoints()):
            point = currentLMNode.GetNthControlPointPosition(mirrorIndex[i])
            mirrorLMNode.AddControlPoint(point)
            sourcePoints.InsertNextPoint(point)

          rigidTransform = vtk.vtkLandmarkTransform()
          rigidTransform.SetSourceLandmarks( sourcePoints )
          rigidTransform.SetTargetLandmarks( targetPoints )
          rigidTransform.SetModeToRigidBody()
          rigidTransformNode=slicer.mrmlScene.AddNewNodeByClass("vtkMRMLTransformNode","Rigid")
          rigidTransformNode.SetAndObserveTransformToParent(rigidTransform)

          # compute normals
          self.computeNormals(currentMeshNode)

          currentMeshNode.SetAndObserveTransformNodeID(rigidTransformNode.GetID())
          mirrorLMNode.SetAndObserveTransformNodeID(rigidTransformNode.GetID())
          slicer.vtkSlicerTransformLogic().hardenTransform(currentMeshNode)
          slicer.vtkSlicerTransformLogic().hardenTransform(mirrorLMNode)

          # optional semi-landmark alignment
          if semilandmarkOption:
            currentSLMNode = self.getLandmarkFileByID(slmDirectory, subjectID)
            mirrorSLMNode =slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode",subjectID)
            for i in range(currentSLMNode.GetNumberOfControlPoints()):
              point = currentSLMNode.GetNthControlPointPosition(mirrorSLMIndex[i])
              mirrorSLMNode.AddControlPoint(point)
            if currentSLMNode :
              currentSLMNode.SetAndObserveTransformNodeID(mirrorTransformNode.GetID())
              slicer.vtkSlicerTransformLogic().hardenTransform(currentSLMNode)
              currentSLMNode.SetAndObserveTransformNodeID(rigidTransformNode.GetID())
              slicer.vtkSlicerTransformLogic().hardenTransform(currentSLMNode)
              outputSLMName = subjectID + '_mirror.mrk.json'
              outputSLMPath = os.path.join(outputSLMDirectory, outputSLMName)
              slicer.util.saveNode(mirrorSLMNode, outputSLMPath)
              slicer.mrmlScene.RemoveNode(currentSLMNode)
              slicer.mrmlScene.RemoveNode(mirrorSLMNode)

          # save output files
          outputMeshName = subjectID + '_mirror.ply'
          outputMeshPath = os.path.join(mirrorMeshDirectory, outputMeshName)
          slicer.util.saveNode(currentMeshNode, outputMeshPath)
          outputLMName = subjectID + '_mirror.mrk.json'
          outputLMPath = os.path.join(mirrorLMDirectory, outputLMName)
          slicer.util.saveNode(mirrorLMNode, outputLMPath)

          # clean up
          slicer.mrmlScene.RemoveNode(currentLMNode)
          slicer.mrmlScene.RemoveNode(currentMeshNode)
          slicer.mrmlScene.RemoveNode(mirrorTransformNode)
          slicer.mrmlScene.RemoveNode(rigidTransformNode)
          slicer.mrmlScene.RemoveNode(mirrorLMNode)

  def runDCAlign(self, baseMeshPath, baseLMPath, meshDirectory, landmarkDirectory, outputDirectory, optionCPD, optionErrorOutput):
    if optionErrorOutput:
      self.errorCheckPath = os.path.join(outputDirectory, "errorChecking")
      if not os.path.exists(self.errorCheckPath):
        os.mkdir(self.errorCheckPath)
    baseNode = slicer.util.loadModel(baseMeshPath)
    baseMesh = baseNode.GetPolyData()
    baseLandmarks=self.fiducialNodeToPolyData(baseLMPath).GetPoints()
    #base, modelExt = os.path.splitext(baseMeshPath)
    modelExt=['ply','stl','vtp']
    self.modelNames, models = self.importMeshes(meshDirectory, modelExt)
    landmarks = self.importLandmarks(landmarkDirectory)
    self.outputDirectory = outputDirectory
    if not(optionCPD):
      denseCorrespondenceGroup = self.denseCorrespondenceBaseMesh(landmarks, models, baseMesh, baseLandmarks)
    else:
      denseCorrespondenceGroup = self.denseCorrespondenceCPD(landmarks, models, baseMesh, baseLandmarks)

    self.addMagnitudeFeature(denseCorrespondenceGroup, self.modelNames, baseMesh)

    # save results to output directory
    outputModelName = 'decaResultModel.vtp'
    outputModelPath = os.path.join(outputDirectory, outputModelName)
    slicer.util.saveNode(baseNode, outputModelPath)

  def runDCAlignSymmetric(self, baseMeshPath, baseLMPath, meshDir, landmarkDir, mirrorMeshDir, mirrorLandmarkDir, outputDir, optionCPD, optionErrorOutput, optionPointOutput):
    if optionErrorOutput:
      self.errorCheckPath = os.path.join(outputDir, "errorChecking")
      if not os.path.exists(self.errorCheckPath):
        os.mkdir(self.errorCheckPath)
    baseNode = slicer.util.loadModel(baseMeshPath)
    baseMesh = baseNode.GetPolyData()
    baseLandmarks=self.fiducialNodeToPolyData(baseLMPath).GetPoints()
    #base, modelExt = os.path.splitext(baseMeshPath)
    modelExt=['ply','stl','vtp']
    self.modelNames, models = self.importMeshes(meshDir, modelExt)
    landmarks = self.importLandmarks(landmarkDir)
    modelMirrorNames, mirrorModels = self.importMeshes(mirrorMeshDir, modelExt)
    mirrorLandmarks = self.importLandmarks(mirrorLandmarkDir)
    self.outputDirectory = outputDir
    if not(optionCPD):
      denseCorrespondenceGroup = self.denseCorrespondenceBaseMesh(landmarks, models, baseMesh, baseLandmarks)
      denseCorrespondenceGroupMirror = self.denseCorrespondenceBaseMesh(mirrorLandmarks, mirrorModels, baseMesh, baseLandmarks)

    else:
      denseCorrespondenceGroup = self.denseCorrespondenceCPD(landmarks, models, baseMesh, baseLandmarks)
      denseCorrespondenceGroupMirror = self.denseCorrespondenceCPD(mirrorLandmarks, mirrorModels, baseMesh, baseLandmarks)

    self.addMagnitudeFeatureSymmetry(denseCorrespondenceGroup, denseCorrespondenceGroupMirror, self.modelNames, baseMesh)

    # save results to output directory
    outputModelName = 'decaSymmetryResultModel.vtp'
    outputModelPath = os.path.join(outputDir, outputModelName)
    slicer.util.saveNode(baseNode, outputModelPath)

  def runMean(self, landmarkDirectory, meshDirectory, modelExt, outputDirectory):
    modelExt=['ply','stl','vtp']
    self.modelNames, models = self.importMeshes(meshDirectory, modelExt)
    landmarks = self.importLandmarks(landmarkDirectory)
    [denseCorrespondenceGroup, closestToMeanIndex] = self.denseCorrespondence(landmarks, models)
    print("Sample closest to mean: ", closestToMeanIndex)
    # compute mean model
    averagePolyData = self.computeAverageModelFromGroup(denseCorrespondenceGroup, closestToMeanIndex)
    averageModelNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLModelNode', 'meanTemplate')
    averageModelNode.CreateDefaultDisplayNodes()
    averageModelNode.SetAndObservePolyData(averagePolyData)
     # compute mean landmarks
    averageLandmarkNode = self.computeAverageLM(landmarks)
    averageLandmarkNode.GetDisplayNode().SetPointLabelsVisibility(False)

    # save results to output directory
    outputModelName = 'decaMeanModel.ply'
    outputModelPath = os.path.join(outputDirectory, outputModelName)
    slicer.util.saveNode(averageModelNode, outputModelPath)

    outputLMName = 'decaMeanModel.mrk.json'
    outputLMPath = os.path.join(outputDirectory, outputLMName)
    slicer.util.saveNode(averageLandmarkNode, outputLMPath)

  def getLandmarkFileByID(self, directory, subjectID):
    fileList = os.listdir(directory)
    for fileName in fileList:
      fileNameBase = Path(fileName)
      while fileNameBase.suffix in {'.fcsv', '.mrk', '.json'}:
        fileNameBase = fileNameBase.with_suffix('')
      if subjectID == str(fileNameBase):
        #print(fileName + " found matching: " + subjectID)
        # if file with this subject id exists, load into scene
        filePath = os.path.join(directory, fileName)
        currentNode = slicer.util.loadMarkups(filePath)
        return currentNode

  def runAlign(self, baseMeshPath, baseLMPath, meshDirectory, lmDirectory, ouputMeshDirectory, outputLMDirectory, removeScaleOption, slmDirectory, outputSLMDirectory):
    semilandmarkOption = bool(slmDirectory and outputSLMDirectory)
    targetPoints = vtk.vtkPoints()
    point=[0,0,0]
    baseMeshNode = slicer.util.loadModel(baseMeshPath)
    baseLMNode = slicer.util.loadMarkups(baseLMPath)
    # Set up base points for transform
    for i in range(baseLMNode.GetNumberOfControlPoints()):
      point = baseLMNode.GetNthControlPointPosition(i)
      targetPoints.InsertNextPoint(point)

    # Transform each subject to base
    for meshFileName in os.listdir(meshDirectory):
      if(not meshFileName.startswith(".")):
        lmFileList = os.listdir(lmDirectory)
        meshFilePath = os.path.join(meshDirectory, meshFileName)
        subjectID = os.path.splitext(meshFileName)[0]
        #print("working on ", subjectID)
        currentLMNode = self.getLandmarkFileByID(lmDirectory, subjectID)
        if currentLMNode :
          try:
            currentMeshNode = slicer.util.loadModel(meshFilePath)
          except:
            slicer.mrmlScene.RemoveNode(currentLMNode)
            continue
          # set up transform between base lms and current lms
          sourcePoints = vtk.vtkPoints()
          for i in range(currentLMNode.GetNumberOfControlPoints()):
            point = currentLMNode.GetNthControlPointPosition(i)
            sourcePoints.InsertNextPoint(point)
            transform = vtk.vtkLandmarkTransform()
            transform.SetSourceLandmarks(sourcePoints)
            transform.SetTargetLandmarks(targetPoints)
          if not removeScaleOption:
            transform.SetModeToRigidBody()
          else:
            transform.SetModeToSimilarity()

          transformNode=slicer.mrmlScene.AddNewNodeByClass("vtkMRMLTransformNode","Alignment")
          transformNode.SetAndObserveTransformToParent(transform)

          # apply transform to the current surface mesh and landmarks
          currentMeshNode.SetAndObserveTransformNodeID(transformNode.GetID())
          currentLMNode.SetAndObserveTransformNodeID(transformNode.GetID())
          slicer.vtkSlicerTransformLogic().hardenTransform(currentMeshNode)
          slicer.vtkSlicerTransformLogic().hardenTransform(currentLMNode)

          # save output files
          outputMeshName = subjectID + '_align.ply'
          outputMeshPath = os.path.join(ouputMeshDirectory, outputMeshName)
          slicer.util.saveNode(currentMeshNode, outputMeshPath)
          outputLMName = subjectID + '_align.mrk.json'
          outputLMPath = os.path.join(outputLMDirectory, outputLMName)
          slicer.util.saveNode(currentLMNode, outputLMPath)

          # optional semi-landmark alignment
          if semilandmarkOption :
            currentSLMNode = self.getLandmarkFileByID(slmDirectory, subjectID)
            if currentSLMNode :
              currentSLMNode.SetAndObserveTransformNodeID(transformNode.GetID())
              slicer.vtkSlicerTransformLogic().hardenTransform(currentSLMNode)
              outputSLMName = subjectID + '_align.mrk.json'
              outputSLMPath = os.path.join(outputSLMDirectory, outputSLMName)
              slicer.util.saveNode(currentSLMNode, outputSLMPath)
              slicer.mrmlScene.RemoveNode(currentSLMNode)

          # clean up
          try:
            slicer.mrmlScene.RemoveNode(currentLMNode)
            slicer.mrmlScene.RemoveNode(currentMeshNode)
            slicer.mrmlScene.RemoveNode(transformNode)
            slicer.mrmlScene.RemoveNode(baseMeshNode)
            slicer.mrmlScene.RemoveNode(baseLMNode)
          except:
            print(f"could not find nodes to remove for {subjectID}")

  def distanceMatrix(self, a):
    """
    Computes the euclidean distance matrix for n points in a 3D space
    Returns a nXn matrix
     """
    id,jd=a.shape
    fnx = lambda q : q - np.reshape(q, (id, 1))
    dx=fnx(a[:,0])
    dy=fnx(a[:,1])
    dz=fnx(a[:,2])
    return (dx**2.0+dy**2.0+dz**2.0)**0.5

  def numpyToFiducialNode(self, numpyArray, nodeName):
    fiducialNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode',nodeName)
    for point in numpyArray:
      fiducialNode.AddControlPoint(point)
    return fiducialNode

  def computeAverageLM(self, fiducialGroup):
    sampleNumber = fiducialGroup.GetNumberOfBlocks()
    pointNumber = fiducialGroup.GetBlock(0).GetNumberOfPoints()
    groupArray_np = np.empty((pointNumber,3,sampleNumber))
    for i in range(sampleNumber):
      pointData = fiducialGroup.GetBlock(i).GetPoints().GetData()
      pointData_np = vtk_np.vtk_to_numpy(pointData)
      groupArray_np[:,:,i] = pointData_np
    #Calculate mean point positions of aligned group
    averagePoints_np = np.mean(groupArray_np, axis=2)
    averageLMNode = self.numpyToFiducialNode(averagePoints_np, "MeanTemplateLM")
    return averageLMNode

  def fiducialNodeToPolyData(self, path):
    point = [0,0,0]
    polydataPoints = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    [success,fiducialNode] = slicer.util.loadMarkupsFiducialList(path)
    if not success:
      print("Could not load landmarks: ", inputFilePath)
      return
    for i in range(fiducialNode.GetNumberOfControlPoints()):
      point = fiducialNode.GetNthControlPointPosition(i)
      points.InsertNextPoint(point)
    polydataPoints.SetPoints(points)
    slicer.mrmlScene.RemoveNode(fiducialNode)
    return polydataPoints

  def importLandmarks(self, topDir):
    fiducialGroup = vtk.vtkMultiBlockDataGroupFilter()
    for file in sorted(os.listdir(topDir)):
      if file.endswith(".fcsv") or file.endswith(".json"):
        print("reading: ", file)
        inputFilePath = os.path.join(topDir, file)
        # may want to replace with vtk reader
        polydataPoints = self.fiducialNodeToPolyData(inputFilePath)
        fiducialGroup.AddInputData(polydataPoints)
    fiducialGroup.Update()
    return fiducialGroup.GetOutput()

  def importMeshes(self, topDir, extensions):
      modelGroup = vtk.vtkMultiBlockDataGroupFilter()
      fileNameList = []
      for file in sorted(os.listdir(topDir)):
        if file.endswith(tuple(extensions)):
          print("reading: ", file)
          base, ext = os.path.splitext(file)
          fileNameList.append(base)
          inputFilePath = os.path.join(topDir, file)
          # may want to replace with vtk reader
          modelNode = slicer.util.loadModel(inputFilePath)
          modelGroup.AddInputData(modelNode.GetPolyData())
          slicer.mrmlScene.RemoveNode(modelNode)
      modelGroup.Update()
      return fileNameList, modelGroup.GetOutput()

  def procrustesImposition(self, originalLandmarks, sizeOption):
    procrustesFilter = vtk.vtkProcrustesAlignmentFilter()
    if(sizeOption):
      procrustesFilter.GetLandmarkTransform().SetModeToRigidBody()

    procrustesFilter.SetInputData(originalLandmarks)
    procrustesFilter.Update()
    meanShape = procrustesFilter.GetMeanPoints()
    return [meanShape, procrustesFilter.GetOutput()]

  def getClosestToMeanIndex(self, meanShape, alignedPoints):
    import operator
    sampleNumber = alignedPoints.GetNumberOfBlocks()
    procrustesDistances = []

    for i in range(sampleNumber):
      alignedShape = alignedPoints.GetBlock(i)
      meanPoint = [0,0,0]
      alignedPoint = [0,0,0]
      distance = 0
      for j in range(meanShape.GetNumberOfPoints()):
        meanShape.GetPoint(j,meanPoint)
        alignedShape.GetPoint(j,alignedPoint)
        distance += np.sqrt(vtk.vtkMath.Distance2BetweenPoints(meanPoint,alignedPoint))

      procrustesDistances.append(distance)

    try:
      min_index, min_value = min(enumerate(procrustesDistances), key=operator.itemgetter(1))
      return min_index
    except:
      return 0

  def denseCorrespondence(self, originalLandmarks, originalMeshes, writeErrorOption=False):
    meanShape, alignedPoints = self.procrustesImposition(originalLandmarks, False)
    sampleNumber = alignedPoints.GetNumberOfBlocks()
    denseCorrespondenceGroup = vtk.vtkMultiBlockDataGroupFilter()
    # get base mesh as the closest to the mean shape
    baseIndex = self.getClosestToMeanIndex(meanShape, alignedPoints)
    baseMesh = originalMeshes.GetBlock(baseIndex)
    baseLandmarks = originalLandmarks.GetBlock(baseIndex).GetPoints()
    for i in range(sampleNumber):
      correspondingMesh = self.denseSurfaceCorrespondencePair(originalMeshes.GetBlock(i),
      originalLandmarks.GetBlock(i).GetPoints(), alignedPoints.GetBlock(i).GetPoints(),
      baseMesh, baseLandmarks, meanShape, i)
      denseCorrespondenceGroup.AddInputData(correspondingMesh)

    denseCorrespondenceGroup.Update()
    return denseCorrespondenceGroup.GetOutput(), baseIndex

  def denseCorrespondenceCPD(self, originalLandmarks, originalMeshes, baseMesh, baseLandmarks, writeErrorOption=False):
    meanShape, alignedPoints = self.procrustesImposition(originalLandmarks, False)
    sampleNumber = alignedPoints.GetNumberOfBlocks()
    denseCorrespondenceGroup = vtk.vtkMultiBlockDataGroupFilter()

    # assign parameters for CPD
    parameters = {
      "SpacingTolerance": .04,
      "CPDIterations": 100,
      "CPDTolerence": 0.001,
      "alpha": 2,
      "beta": 2,
     }

    for i in range(sampleNumber):
      correspondingPoints = self.runCPDRegistration(originalMeshes.GetBlock(i), baseMesh, parameters)
      # convert to vtkPoints
      correspondingMesh = self.convertPointsToVTK(correspondingPoints)
      correspondingMesh.SetPolys(baseMesh.GetPolys())
      # convert to polydata
      denseCorrespondenceGroup.AddInputData(correspondingMesh)
      # write ouput
      if writeErrorOption:
        plyWriterSubject = vtk.vtkPLYWriter()
        plyWriterSubject.SetFileName("/Users/sararolfe/Dropbox/SlicerWorkspace/SMwSML/Data/UBC/DECAOutCPD/" + str(i) + ".ply")
        plyWriterSubject.SetInputData(correspondingMesh)
        plyWriterSubject.Write()

        plyWriterBase = vtk.vtkPLYWriter()
        plyWriterBase.SetFileName("/Users/sararolfe/Dropbox/SlicerWorkspace/SMwSML/Data/UBC/DECAOutCPD/base.ply")
        plyWriterBase.SetInputData(baseMesh)
        plyWriterBase.Write()

    denseCorrespondenceGroup.Update()
    return denseCorrespondenceGroup.GetOutput()

  def runCPDRegistration(self,sourceData, targetData, parameters):
    from open3d import geometry
    from open3d import utility

    # Downsample meshes
    sourceFilter=vtk.vtkCleanPolyData()
    sourceFilter.SetTolerance(parameters["SpacingTolerance"])
    sourceFilter.SetInputData(sourceData)
    sourceFilter.Update()
    #sourceArray = sourceFilter.GetOutput().GetPoints().GetData()
    sourceArray = sourceData.GetPoints().GetData()

    targetFilter=vtk.vtkCleanPolyData()
    targetFilter.SetTolerance(parameters["SpacingTolerance"])
    targetFilter.SetInputData(targetData)
    targetFilter.Update()
    #targetArray = targetFilter.GetOutput().GetPoints().GetData()
    targetArray = targetData.GetPoints().GetData()

    # Convert to pointcloud for scaling
    sourceArray_np = vtk_np.vtk_to_numpy(sourceArray)
    targetArray_np = vtk_np.vtk_to_numpy(targetArray)
    sourceCloud = geometry.PointCloud()
    sourceCloud.points = utility.Vector3dVector(sourceArray_np)
    targetCloud = geometry.PointCloud()
    targetCloud.points = utility.Vector3dVector(targetArray_np)
    cloudSize = np.max(targetCloud.get_max_bound() - targetCloud.get_min_bound())
    targetCloud.scale(25 / cloudSize, center = False)
    sourceCloud.scale(25 / cloudSize, center = False)
    # Convert back to numpy for cpd
    sourceArray = np.asarray(sourceCloud.points,dtype=np.float32)
    targetArray = np.asarray(targetCloud.points,dtype=np.float32)
    registrationOutput = self.cpd_registration(targetArray, sourceArray, parameters["CPDIterations"], parameters["CPDTolerence"], parameters["alpha"], parameters["beta"])
    deformed_array, _ = registrationOutput.register()
    outputCloud = geometry.PointCloud()
    outputCloud.points = utility.Vector3dVector(deformed_array)
    outputCloud.scale(cloudSize/25, center = False)
    return np.asarray(outputCloud.points)

  def cpd_registration(self, targetArray, sourceArray, CPDIterations, CPDTolerence, alpha_parameter, beta_parameter):
    from pycpd import DeformableRegistration
    output = DeformableRegistration(**{'X': targetArray, 'Y': sourceArray,'max_iterations': CPDIterations, 'tolerance': CPDTolerence}, alpha = alpha_parameter, beta  = beta_parameter)
    return output

  def denseCorrespondenceBaseMesh(self, originalLandmarks, originalMeshes, baseMesh, baseLandmarks):
    meanShape, alignedPoints = self.procrustesImposition(originalLandmarks, False)
    sampleNumber = alignedPoints.GetNumberOfBlocks()
    denseCorrespondenceGroup = vtk.vtkMultiBlockDataGroupFilter()
    for i in range(sampleNumber):
      correspondingMesh = self.denseSurfaceCorrespondencePair(originalMeshes.GetBlock(i),
      originalLandmarks.GetBlock(i).GetPoints(), alignedPoints.GetBlock(i).GetPoints(),
      baseMesh, baseLandmarks, meanShape, i)
      denseCorrespondenceGroup.AddInputData(correspondingMesh)

    denseCorrespondenceGroup.Update()
    return denseCorrespondenceGroup.GetOutput()

  def denseSurfaceCorrespondencePair(self, originalMesh, originalLandmarks, alignedLandmarks, baseMesh, baseLandmarks, meanShape, iteration):
    # TPS warp target and base mesh to meanshape
    meanTransform = vtk.vtkThinPlateSplineTransform()
    meanTransform.SetSourceLandmarks( originalLandmarks )
    meanTransform.SetTargetLandmarks( meanShape )
    meanTransform.SetBasisToR() # for 3D transform

    meanTransformFilter = vtk.vtkTransformPolyDataFilter()
    meanTransformFilter.SetInputData(originalMesh)
    meanTransformFilter.SetTransform(meanTransform)
    meanTransformFilter.Update()
    meanWarpedMesh = meanTransformFilter.GetOutput()

    meanTransformBase = vtk.vtkThinPlateSplineTransform()
    meanTransformBase.SetSourceLandmarks( baseLandmarks )
    meanTransformBase.SetTargetLandmarks( meanShape )
    meanTransformBase.SetBasisToR() # for 3D transform

    meanTransformBaseFilter = vtk.vtkTransformPolyDataFilter()
    meanTransformBaseFilter.SetInputData(baseMesh)
    meanTransformBaseFilter.SetTransform(meanTransformBase)
    meanTransformBaseFilter.Update()
    meanWarpedBase = meanTransformBaseFilter.GetOutput()

    # write ouput
    if hasattr(self,"errorCheckPath"):
      plyWriterSubject = vtk.vtkPLYWriter()
      print(self.modelNames)
      print(iteration)
      plyName = "subject_" + self.modelNames[iteration] + ".ply"
      plyPath = os.path.join(self.errorCheckPath, plyName)
      plyWriterSubject.SetFileName(plyPath)
      plyWriterSubject.SetInputData(meanWarpedMesh)
      plyWriterSubject.Write()

      plyWriterBase = vtk.vtkPLYWriter()
      plyName = "base.ply"
      plyPath = os.path.join(self.errorCheckPath, plyName)
      plyWriterBase.SetFileName(plyPath)
      plyWriterBase.SetInputData(meanWarpedBase)
      plyWriterBase.Write()

    # Dense correspondence
    cellLocator = vtk.vtkCellLocator()
    cellLocator.SetDataSet(meanWarpedMesh)
    cellLocator.BuildLocator()

    point = [0,0,0]
    correspondingPoint = [0,0,0]
    correspondingPoints = vtk.vtkPoints()
    cellId = vtk.reference(0)
    subId = vtk.reference(0)
    distance = vtk.reference(0.0)
    for i in range(meanWarpedBase.GetNumberOfPoints()):
      meanWarpedBase.GetPoint(i,point)
      cellLocator.FindClosestPoint(point,correspondingPoint,cellId, subId, distance)
      correspondingPoints.InsertPoint(i,correspondingPoint)

    #Copy points into mesh with base connectivity
    correspondingMesh = vtk.vtkPolyData()
    correspondingMesh.SetPoints(correspondingPoints)
    correspondingMesh.SetPolys(meanWarpedBase.GetPolys())

    # Apply inverse warping
    inverseTransform = vtk.vtkThinPlateSplineTransform()
    inverseTransform.SetSourceLandmarks( meanShape )
    inverseTransform.SetTargetLandmarks( originalLandmarks )
    inverseTransform.SetBasisToR() # for 3D transform

    inverseTransformFilter = vtk.vtkTransformPolyDataFilter()
    inverseTransformFilter.SetInputData(correspondingMesh)
    inverseTransformFilter.SetTransform(inverseTransform)
    inverseTransformFilter.Update()

    return inverseTransformFilter.GetOutput()

  def convertPointsToVTK(self, points):
    array_vtk = vtk_np.numpy_to_vtk(points, deep=True, array_type=vtk.VTK_FLOAT)
    points_vtk = vtk.vtkPoints()
    points_vtk.SetData(array_vtk)
    polydata_vtk = vtk.vtkPolyData()
    polydata_vtk.SetPoints(points_vtk)
    return polydata_vtk

  def computeAverageModelFromGroup(self, denseCorrespondenceGroup, baseIndex):
    sampleNumber = denseCorrespondenceGroup.GetNumberOfBlocks()
    pointNumber = denseCorrespondenceGroup.GetBlock(0).GetNumberOfPoints()
    groupArray_np = np.empty((pointNumber,3,sampleNumber))

    # get base mesh as closest to the meanshape
    baseMesh = denseCorrespondenceGroup.GetBlock(baseIndex)

     # get points as array
    for i in range(sampleNumber):
      alignedMesh = denseCorrespondenceGroup.GetBlock(i)
      alignedMesh_np = vtk_np.vtk_to_numpy(alignedMesh.GetPoints().GetData())
      groupArray_np[:,:,i] = alignedMesh_np

    #Calculate mean point positions of aligned group
    averagePoints_np = np.mean(groupArray_np, axis=2)
    averagePointsPolydata = self.convertPointsToVTK(averagePoints_np)

    #Copy points into mesh with base connectivity
    averageModel = vtk.vtkPolyData()
    averageModel.SetPoints(averagePointsPolydata.GetPoints())
    averageModel.SetPolys(baseMesh.GetPolys())
    return averageModel

  def addMagnitudeFeature(self, denseCorrespondenceGroup, modelNameArray, model):
    sampleNumber = denseCorrespondenceGroup.GetNumberOfBlocks()
    pointNumber = denseCorrespondenceGroup.GetBlock(0).GetNumberOfPoints()
    statsArray = np.zeros((pointNumber, sampleNumber))
    magnitudeMean = vtk.vtkDoubleArray()
    magnitudeMean.SetNumberOfComponents(1)
    magnitudeMean.SetName("Magnitude Mean")
    magnitudeSD = vtk.vtkDoubleArray()
    magnitudeSD.SetNumberOfComponents(1)
    magnitudeSD.SetName("Magnitude SD")

     # get distance arrays
    for i in range(sampleNumber):
      alignedMesh = denseCorrespondenceGroup.GetBlock(i)
      magnitudes = vtk.vtkDoubleArray()
      magnitudes.SetNumberOfComponents(1)
      magnitudes.SetName(modelNameArray[i])
      for j in range(pointNumber):
        modelPoint = model.GetPoint(j)
        targetPoint = alignedMesh.GetPoint(j)
        distance = np.sqrt(vtk.vtkMath.Distance2BetweenPoints(modelPoint,targetPoint))
        magnitudes.InsertNextValue(distance)
        statsArray[j,i]=distance

      model.GetPointData().AddArray(magnitudes)

    for i in range(pointNumber):
      pointMean = statsArray[i,:].mean()
      magnitudeMean.InsertNextValue(pointMean)
      pointSD = statsArray[i,:].std()
      magnitudeSD.InsertNextValue(pointSD)

    model.GetPointData().AddArray(magnitudeMean)
    model.GetPointData().AddArray(magnitudeSD)

  def addMagnitudeFeatureSymmetry(self, denseCorrespondenceGroup, denseCorrespondenceGroupMirror, modelNameArray, model):
    sampleNumber = denseCorrespondenceGroup.GetNumberOfBlocks()
    pointNumber = denseCorrespondenceGroup.GetBlock(0).GetNumberOfPoints()
    statsArray = np.zeros((pointNumber, sampleNumber))
    magnitudeMean = vtk.vtkDoubleArray()
    magnitudeMean.SetNumberOfComponents(1)
    magnitudeMean.SetName("Magnitude Mean")
    magnitudeSD = vtk.vtkDoubleArray()
    magnitudeSD.SetNumberOfComponents(1)
    magnitudeSD.SetName("Magnitude SD")

     # get distance arrays
    for i in range(sampleNumber):
      alignedMesh = denseCorrespondenceGroup.GetBlock(i)
      mirrorMesh = denseCorrespondenceGroupMirror.GetBlock(i)
      magnitudes = vtk.vtkDoubleArray()
      magnitudes.SetNumberOfComponents(1)
      magnitudes.SetName(modelNameArray[i])
      for j in range(pointNumber):
        modelPoint = model.GetPoint(j)
        targetPoint1 = alignedMesh.GetPoint(j)
        targetPoint2 = mirrorMesh.GetPoint(j)
        distance = np.sqrt(vtk.vtkMath.Distance2BetweenPoints(targetPoint1,targetPoint2))
        magnitudes.InsertNextValue(distance)
        statsArray[j,i]=distance

      model.GetPointData().AddArray(magnitudes)

    for i in range(pointNumber):
      pointMean = statsArray[i,:].mean()
      magnitudeMean.InsertNextValue(pointMean)
      pointSD = statsArray[i,:].std()
      magnitudeSD.InsertNextValue(pointSD)

    model.GetPointData().AddArray(magnitudeMean)
    model.GetPointData().AddArray(magnitudeSD)
