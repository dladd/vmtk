#!/usr/bin/env python

import sys
import vtk
import vtkvmtk
import pypes
import numpy

vmtkpcvwriter = 'vmtkPcvWriter'

class vmtkPcvWriter(pypes.pypeScript):

    def __init__(self):

        pypes.pypeScript.__init__(self)

        self.Image = None
        self.ReferenceImage = None
        self.Output = None
        self.Method = ''
	self.Tolerance = 1E-8
        self.Result = ''
        self.ResultLog = ''

        self.InputFileName = ''
        self.OutputFileName = ''

        self.NumberSlices = 36
        self.ImageScalarType = 4

        self.DataExtent = [-1, -1, -1, -1, -1, -1]
        self.DataSpacing = [1.0, 1.0, 1.0]
        self.DataOrigin = [0.0, 0.0, 0.0]
        self.DataByteOrder = 'littleendian'
        self.DataScalarType = 'float'
        self.DesiredOrientation = 'native'
        self.HeaderSize = 0
        self.FileDimensionality = 3
        self.Flip = [0, 0, 0]
        self.InputDirectoryName = ''
        self.Format = ''

        self.NumberOfTimesteps = 0
        self.TimeIncrement = 0.0
        self.HighThreshold = 0
        self.LowThreshold = -100
        self.AutoOrientDICOMImage = 1
        self.InputFilePrefix = ''
        self.ImageNumberPadding = 5
        self.NumberOfFieldComponents = 3
        self.NumberOfSlices = 36

        self.MatrixCoefficients = []
        self.InvertMatrix = 0
        self.Matrix4x4 = None
        self.Rotation = [0.0,0.0,0.0]
        self.Translation = [0.0,0.0,0.0]
        self.Scaling = [1.0,1.0,1.0]
        self.Flip = [0, 0, 0]

        self.SetScriptName('vmtkpcvwriter')
        self.SetScriptDoc('extracts velocity vector data from a PCV image stack (assumes a particular folder structure- see below)')
        self.SetInputMembers([
            ['Format','f','str',1,'["vtkxml","vtk","dicom","raw","meta","tiff","png"]','file format'],
            ['Image','i','vtkImageData',1,'','the input image','vmtkimagereader'],
            ['InputDirectoryName','d','str',1,'','input directory name - dicom only'],
            ['ReferenceImage','r','vtkImageData',1,'','the reference image to compare against','vmtkimagereader'],
            ['Tolerance','tolerance','float',1,'','tolerance for numerical comparisons'],
            ['NumberOfTimesteps','numbertimesteps','int',1,'','number of timesteps for which there is PCV data'],
            ['NumberSlices','slices','int',1,'','number of slices per stack for the PCV data'],
            ['ImageScalarType','imagescalar','int',1,'','VTK scalar type for image data (choose 5 for magnitude, 4 for flow)'],
            ['TimeIncrement','timeincrement','float',1,'','time increment between PCV timesteps'],
            ['LowThreshold','lowthreshold','float',1,'','lower threshold value for reference image'],
            ['HighThreshold','highthreshold','float',1,'','upper threshold value for reference image'],
            ['Rotation','rotation','float',3,'','rotations around the x-,y- and z-axis for the PCV image stack'],
            ['Translation','translation','float',3,'','translation in the x-,y- and z-directions'],
            ['Scaling','scaling','float',3,'','scaling of the x-,y- and z-directions'],
            ['Flip','flip','bool',3,'','toggle flipping of the corresponding axis'],
            ['OutputFileName','ofile','str',1,'','output file name'],
            ['InputFilePrefix','prefix','str',1,'','input file prefix (e.g. foo_)']
            ])
        self.SetOutputMembers([
            ['Image','o','vtkImageData',1,'','the output image','vmtkimagewriter'],
            ['Result','result','bool',1,'','Output boolean stating if images are equal or not'],
            ['ResultLog','log','str',1,'','Result Log']
            ])


    def ReadITKIO(self):
        if self.InputFileName == '':
            self.PrintError('Error: no InputFileName.')
        reader = vtkvmtk.vtkvmtkITKArchetypeImageSeriesScalarReader()
        reader.SetArchetype(self.InputFileName)
        reader.SetDefaultDataSpacing(self.DataSpacing)
        reader.SetDefaultDataOrigin(self.DataOrigin)
        reader.SetOutputScalarTypeToNative()
        if self.DesiredOrientation == 'native':
            reader.SetDesiredCoordinateOrientationToNative()
        elif self.DesiredOrientation == 'axial':
            reader.SetDesiredCoordinateOrientationToAxial()
        elif self.DesiredOrientation == 'coronal':
            reader.SetDesiredCoordinateOrientationToCoronal()
        elif self.DesiredOrientation == 'sagittal':
            reader.SetDesiredCoordinateOrientationToSagittal()
        reader.SetSingleFile(0)
        reader.Update()
        self.Image = vtk.vtkImageData()
        self.Image.DeepCopy(reader.GetOutput())
        matrix = reader.GetRasToIjkMatrix()
        self.RasToIjkMatrixCoefficients = [
            matrix.GetElement(0,0), matrix.GetElement(0,1), matrix.GetElement(0,2), matrix.GetElement(0,3),
            matrix.GetElement(1,0), matrix.GetElement(1,1), matrix.GetElement(1,2), matrix.GetElement(1,3),
            matrix.GetElement(2,0), matrix.GetElement(2,1), matrix.GetElement(2,2), matrix.GetElement(2,3),
            matrix.GetElement(3,0), matrix.GetElement(3,1), matrix.GetElement(3,2), matrix.GetElement(3,3)]

        matrix.Invert()
        origin = [matrix.GetElement(0,3), matrix.GetElement(1,3), matrix.GetElement(2,3)] 
        translationToOrigin = [-origin[0], -origin[1], -origin[2]] 

        for i in range(3):
            direction = [matrix.GetElement(0,i), matrix.GetElement(1,i), matrix.GetElement(2,i)]
            vtk.vtkMath.Normalize(direction)
            matrix.SetElement(0,i,direction[0])
            matrix.SetElement(1,i,direction[1])
            matrix.SetElement(2,i,direction[2])
        matrix.SetElement(0,3,0.0)
        matrix.SetElement(1,3,0.0)
        matrix.SetElement(2,3,0.0)
 
        transform = vtk.vtkTransform()
        transform.PostMultiply()
        transform.Translate(translationToOrigin)
        transform.Concatenate(matrix)
        transform.Translate(origin)

        matrix = transform.GetMatrix()
        self.XyzToRasMatrixCoefficients = [
            matrix.GetElement(0,0), matrix.GetElement(0,1), matrix.GetElement(0,2), matrix.GetElement(0,3),
            matrix.GetElement(1,0), matrix.GetElement(1,1), matrix.GetElement(1,2), matrix.GetElement(1,3),
            matrix.GetElement(2,0), matrix.GetElement(2,1), matrix.GetElement(2,2), matrix.GetElement(2,3),
            matrix.GetElement(3,0), matrix.GetElement(3,1), matrix.GetElement(3,2), matrix.GetElement(3,3)]


    def filterReferenceImage(self):

        threshold = vtk.vtkImageThreshold()
        threshold.SetInput(self.ReferenceImage)

        self.PrintLog('Filtering reference image...')
        threshold.ThresholdBetween(self.LowThreshold, self.HighThreshold)
        threshold.SetInValue(1)
        threshold.SetOutValue(0)
        threshold.SetOutputScalarType(self.ImageScalarType)
        threshold.Update()

        self.ReferenceImage = threshold.GetOutput()

        array = self.ReferenceImage.GetPointData().GetArray(0)
        arrayName = array.GetName()
        self.PrintLog('Reference array name = ' + str(arrayName))

        if (self.OutputFileName == ''):
            self.PrintError('Error: no OutputFileName.')

        geomFileName=self.OutputFileName + '.geom'
        self.PrintLog('Writing Geometry data to: ' + str(geomFileName))
        f=open(geomFileName, 'w')
        array = self.ReferenceImage.GetPointData().GetArray(0)
        self.numberDataPoints = 0
        for i in range(self.ReferenceImage.GetNumberOfPoints()):
            if array.GetComponent(i,0) == 1 :
                self.numberDataPoints += 1
        firstLine = "%d\n" % self.numberDataPoints
        f.write(firstLine)
        for i in range(self.ReferenceImage.GetNumberOfPoints()):
            if array.GetComponent(i,0) == 1 :
                point = self.ReferenceImage.GetPoint(i)
                line = "%f %f %f\n" % (point[0],point[1],point[2])
                f.write(line)

        self.PrintLog('Number of data points in region = ' + str(self.numberDataPoints))

        f.close()
        self.PrintLog('done!')


    def transformImage(self):

        if (self.Flip[0] == 1) | (self.Flip[1] == 1) | (self.Flip[2] == 1):
            temp0 = self.Image
            if self.Flip[0] == 1:
                flipFilter = vtk.vtkImageFlip()
                flipFilter.SetInput(self.Image)
                flipFilter.SetFilteredAxis(0)
                flipFilter.Update()
                temp0 = flipFilter.GetOutput()
            temp1 = temp0
            if self.Flip[1] == 1:
                flipFilter = vtk.vtkImageFlip()
                flipFilter.SetInput(temp0)
                flipFilter.SetFilteredAxis(1)
                flipFilter.Update()
                temp1 = flipFilter.GetOutput()
            temp2 = temp1
            if self.Flip[2] == 1:
                flipFilter = vtk.vtkImageFlip()
                flipFilter.SetInput(temp1)
                flipFilter.SetFilteredAxis(2)
                flipFilter.Update()
                temp2 = flipFilter.GetOutput()
            self.Image = temp2

        if self.Translation != [0.0,0.0,0.0] or self.Rotation != [0.0,0.0,0.0] or self.Scaling != [1.0,1.0,1.0]:        

            resliceFilter = vtk.vtkImageReslice()
            resliceFilter.SetInput(self.Image)

            self.PrintLog('Setting up transform matrix using specified translation, rotation and/or scaling')
            transform = vtk.vtkTransform()
            transform.RotateX(self.Rotation[0])
            transform.RotateY(self.Rotation[1])
            transform.RotateZ(self.Rotation[2])                       
            transform.Translate(self.Translation[0], self.Translation[1], self.Translation[2])
            transform.Scale(self.Scaling[0], self.Scaling[1], self.Scaling[2])
            self.Matrix4x4 = vtk.vtkMatrix4x4()
            self.Matrix4x4.DeepCopy(transform.GetMatrix())

            transform = vtk.vtkMatrixToLinearTransform()
            transform.SetInput(self.Matrix4x4)
            resliceFilter.SetResliceTransform(transform)
            resliceFilter.Update()
            self.Image = resliceFilter.GetOutput()

    def multiplyImageByReference(self):

        composeFilter = vtk.vtkImageMathematics()
        composeFilter.SetInput1(self.Image)
        composeFilter.SetInput2(self.ReferenceImage)
        composeFilter.SetOperationToMultiply()
        composeFilter.Update()

        self.Image = composeFilter.GetOutput()

    def readImageValues(self):

        if self.Image.GetPointData().GetScalars().GetName() == None:
            self.Image.GetPointData().GetScalars().SetName('__Scalars')
        array = self.Image.GetPointData().GetArray(0)
        arrayName = array.GetName()
        self.PrintLog('array name = ' + str(arrayName))

        for i in range(self.Image.GetNumberOfPoints()):
            self.fieldData[i][self.velComponent]=array.GetComponent(i,0)

    def writeTimestepValues(self):

        if (self.OutputFileName == ''):
            self.PrintError('Error: no OutputFileName.')
        self.PrintLog('Writing PCV data for time: ' + str(self.timestep))

        pcvFileName=self.OutputFileName + str(self.timestep) + '.pcv'
        f=open(pcvFileName, 'w')
        firstLine = "%d\n" % self.numberDataPoints
        f.write(firstLine)

        dataPoint = 0

        array = self.ReferenceImage.GetPointData().GetArray(0)
        for i in range(self.ReferenceImage.GetNumberOfPoints()):
            if array.GetComponent(i,0) == 1 :
                dataPoint += 1
                line = "%f %f %f\n" % (self.fieldData[i][0],self.fieldData[i][1],self.fieldData[i][2])
                f.write(line)
        f.close()
        self.PrintLog('done!')

    def Execute(self):

        self.filterReferenceImage()

        for time in range(self.NumberOfTimesteps):
            self.timestep = time + 1
            self.fieldData = [[0 for k in range(3)] for j in range(self.ReferenceImage.GetNumberOfPoints())]
            for self.velComponent in range(3):
                if self.velComponent ==0: 
                    componentName = "u"
                elif self.velComponent ==1: 
                    componentName = "v"
                elif self.velComponent ==2: 
                    componentName = "w"
                else:
                    exit

                fileNumber = 1 + self.velComponent*self.NumberOfTimesteps*self.NumberSlices + time*self.NumberSlices
                fileFormatted = str(fileNumber).zfill(5)
                self.InputFileName = str(self.timestep) + "/" + componentName + "/image" + fileFormatted + ".dcm"
                self.PrintLog('reading : ' +str(self.InputFileName))

                self.ReadITKIO()

                self.transformImage()
                self.multiplyImageByReference()

                self.readImageValues()

            self.writeTimestepValues()            

        
if __name__=='__main__':
    main = pypes.pypeMain()
    main.Arguments = sys.argv
    main.Execute()


