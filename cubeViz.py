import glob, os
from pylab import*
from os.path import expanduser, join
from mayavi import mlab
import sys
mlab.close(all=True)


" Data Types******************************************************************" 
geometricHeader = dtype([("nCores", float),("origo",float,(3,)),
                     ("xPoints", float),("xLimits",float,(2,)),
                     ("yPoints", float),("yLimits",float,(2,)),
                     ("zPoints", float),("zLimits",float,(2,)),
                     ])
                     
systemHeader = dtype([("id", float),("charge", float), 
                     ("position", float, (3,))])

dataType = dtype([("density", float)])

"Read files*******************************************************************" 
def readFiles():
    rawDataPath = "/home/milad/kurs/qmd/density"
    stateFiles = glob.glob1(rawDataPath,'*.bin')
    nProcs = len(glob.glob1(rawDataPath,'*_cubeFile0000.bin'))
    nStateFiles = len(stateFiles)/nProcs
    print "# state files: ", nStateFiles

    rawData = [0] * nStateFiles
    atomList= [0] * nStateFiles
    geometricData = [0] * nStateFiles
    
    for i in range(0, nStateFiles):
        for j in range(0, nProcs):
            cubeFile = join(rawDataPath,"id" + str(j) + "_cubeFile" 
                        + "%04d" % i  + ".bin")  
            if os.path.exists(cubeFile):
                header, atoms, densityData = loadCubeFile(cubeFile)
                rawData[i] += densityData
                atomList[i] = atoms
                geometricData[i] = header

    return geometricData, atomList,rawData

"File loader *****************************************************************" 
def loadCubeFile(fileName):
    fileName = expanduser(fileName)
    cubeFile = open(fileName, "rb")
    
    #Header with geometric data
    geometricData = fromfile(cubeFile, dtype=geometricHeader, count = 1)
    nCores = int(geometricData[0][0])
    nX = int(geometricData[0][2])
    nY = int(geometricData[0][4])
    nZ = int(geometricData[0][6])
        
    #Header with core data
    atoms = fromfile(cubeFile, dtype=systemHeader, count = nCores)

    #Density data
    densityData = fromfile(cubeFile, dtype=dataType)["density"]
    densityData = array(densityData).reshape(nX,nY,nZ).transpose()
    cubeFile.close()
    
    return geometricData, atoms, densityData
            
"*****************************************************************************"   
class cubeVizualizer:
    """
    Superclass for vizualization of cube files, inlcuding:
        
        1. Core positions
        2. Charge density
    """ 
    def __init__(self, geometricData, atomList, densityData,
                 bgcolor, showCores, showDensity,
                 animate, saveFigure):

        self.geometricData = geometricData
        self.atomList = atomList
        self.densityData  = densityData    
        self.showCores = showCores
        self.showDensity = showDensity
        self.bgcolor = bgcolor
        self.animate = animate
        self.saveFigure = saveFigure
        
        self.vmax = self.densityData[0].max()
        self.vmin = self.densityData[0].min()
        
        self.xLim = array((self.geometricData[0][0][3]))
        self.yLim = array((self.geometricData[0][0][5]))
        self.zLim = array((self.geometricData[0][0][7]))
        self.nX = complex(self.geometricData[0][0][2])
        self.nY = complex(self.geometricData[0][0][4])
        self.nZ = complex(self.geometricData[0][0][6])

        
        self.x, self.y, self.z = mgrid[self.xLim[0]:self.xLim[1]:self.nX, 
                                       self.yLim[0]:self.yLim[1]:self.nY,
                                       self.zLim[0]:self.zLim[1]:self.nZ]

        print "# points: ",int(self.nX.real),int(self.nY.real),int(self.nZ.real)
        print "limits: ", self.xLim, self.yLim, self.zLim
                                       
    
    def vizualize(self):
        # Set background color
        if self.bgcolor == "w":
            bgcolor = (1,1,1)
        elif self.bgcolor == "b":
            bgcolor = (0,0,0)
        else:
            bgcolor = None
        
        mlab.figure("cubeViz", bgcolor= bgcolor, size=(750, 550))
        
        #Show cores?        
        if self.showCores:
            self.coreSource = self.displayCores(self.atomList[0])
        else:
            self.coreSource = None
            
        
        #Show density?    
        if self.showDensity:
            self.densitySource = self.displayDensity(self.densityData[0])
        else:
            self.densitySource = None
            
        
        #Animation?    
        if self.animate:
           self.anim(self)
        
    def displayCores(self,atoms):
        nCores = len(atoms) 
        corePositions = zeros((nCores,3))
                    
            
        for i in range(0,nCores):
            corePositions[i] =  array((atoms[i][2][0], atoms[i][2][1], 
                                        atoms[i][2][2]))
            coreCharge   =  atoms[i][1]
            if i < nCores-4:
                c = (1,0,0)
            else:
                c = (0,1,0)
                
            corePlot = mlab.points3d(corePositions[i,0],corePositions[i,1],
                                 corePositions[i,2],
                                 sqrt(coreCharge),
                                 scale_factor=0.2,
                                 resolution=100,
                                 color = c,
                                 opacity = 1.0)
            
        return corePlot.mlab_source
        
    def displayDensity(self, densityData):
        raise NotImplementedError
        
    def numberOfElectrons(self, densityData):
        dx = linspace(self.xLim[0],self.xLim[1], self.nX.real, retstep=True)[1]
        dy = linspace(self.yLim[0],self.yLim[1], self.nY.real, retstep=True)[1]
        dz = linspace(self.zLim[0],self.zLim[1], self.nZ.real, retstep=True)[1]
        dr = dx * dy * dz
        nElectrons = sum(densityData) * dr   
        print "# electrons: ", nElectrons
        
        
        
    @mlab.animate(delay=10, ui = True)
    def anim(self): 
        f = mlab.gcf()
        i = 0;

        while 1:
         #Update core positions
         if self.coreSource!= None:
             nStep = len(self.atomList)
             atoms = self.atomList[i]
             nCores = len(atoms)
             corePositions = zeros((nCores,3))
             coreCharge = zeros(nCores)
 
             for j in range(0,nCores):
                 corePositions[j] =  array((atoms[j][2][0],
                                     atoms[j][2][1],
                                     atoms[j][2][2]))
                 coreCharge[j]    =  atoms[j][0] 
                 
             self.coreSource.set(x = corePositions[:,0],
                            y = corePositions[:,1],
                            z = corePositions[:,2])
         
         #Update density
         if self.densitySource!= None:
             nStep = len(self.densityData)
             density = self.densityData[i]
             self.numberOfElectrons(density)
             self.densitySource.set(scalars=density)
         
         f.scene.render()

         if self.saveFigure: 
             rawDataPath = "/home/milad/kurs/qmd/density"
             figure = join(rawDataPath, "render/cubeFile" + "%04d" % i+".png")
             mlab.savefig(figure)
         
         i+=1             
         if i >= nStep:
             self.saveFigure = False
             i = 0
         yield

"*****************************************************************************"   
class contourRepresentation(cubeVizualizer):
    """
    Electron density represented by iso-surfaces 
    """ 
    
    def displayDensity(self, densityValues):
        self.numberOfElectrons(densityValues)
        densityPlot = mlab.contour3d(self.x, self.y, self.z, 
                                     densityValues,
                                     vmin = self.vmin,
                                     vmax = self.vmax*0.001,
                                     colormap = "hot",
                                     opacity = 0.1,
                                     line_width = 1.0,
                                     contours=100, transparent=True)
        mlab.outline()
        mlab.axes()
        densityPlot.actor.property.representation = 'surface'
        return densityPlot.mlab_source
"*****************************************************************************" 
class volumeRepresentation(cubeVizualizer):
    """
    Volume rendering to display the electron density.
    """ 
    
    def displayDensity(self, densityValues):
        self.numberOfElectrons(densityValues)
        source = mlab.pipeline.scalar_field(self.x, self.y, self.z, 
                                            densityValues)

        densityPlot = mlab.pipeline.volume(source,
                                    vmin = self.vmin,
                                    vmax = self.vmax*0.005)
        #Set colormap             
        densityPlot = self.setColormap(densityPlot)
    
        return densityPlot.mlab_source
        
    def setColormap(self, densityPlot):
        from tvtk.util.ctf import ColorTransferFunction
        ctf = ColorTransferFunction() 
       
       # Add points to CTF 
        ctf.add_rgb_point(0, 1.0, 1.0, 1.0) 
        ctf.add_rgb_point(0.8, 0.0, 0.0, 1.0) 
        ctf.add_rgb_point(1, 0.0, 0.0, 1.0) 
        
        densityPlot._volume_property.set_color(ctf) 
        densityPlot.update_ctf = True
        
        return densityPlot
"*****************************************************************************" 
class slicer1Representation(cubeVizualizer):
    """
    Electron density represented by slices
    """ 
    
    def displayDensity(self, densityValues):
        self.numberOfElectrons(densityValues)
        source = mlab.pipeline.scalar_field(self.x, self.y, self.z, 
                                            densityValues)    
        densityPlot = mlab.pipeline.image_plane_widget(source,
                        colormap="hot",opacity = 1, transparent=True,                   
                        plane_orientation='x_axes',vmin = self.vmin,
                                vmax = self.vmax*0.01,
                        slice_index=20,)

        densityPlot = mlab.pipeline.image_plane_widget(source,
                    colormap="hot",  opacity = 1,transparent=True,vmin = self.vmin,
                                vmax = self.vmax*0.01,
                    plane_orientation='y_axes',
                    slice_index=20,)         
                    
        densityPlot = mlab.pipeline.image_plane_widget(source,
                colormap="hot",   opacity = 1, transparent=True, vmin = self.vmin,
                                vmax = self.vmax*0.01,                            
                plane_orientation='z_axes',
                slice_index=20,)    
                
        mlab.outline()
#        mlab.axes()
        
        return densityPlot.mlab_source 
            
"*****************************************************************************"
def define_command_line_options(parser=None):
    if parser is None:
        import argparse
        parser = argparse.ArgumentParser()

    parser.add_argument('--showCores', action='store_true', default=1,
                        help='show cores')  
                        
    parser.add_argument('--showDensity', action='store_true', default=True,
                        help='show charge density')  
                        
    parser.add_argument(
        '--bgcolor', type=str, default= "b" , help='background coclor')
                        
    parser.add_argument('--animate', action='store_true', default=0,
                        help='make animation')
    
    parser.add_argument('--saveFigure', action='store_true', default=0,
                        help='save figures for movie')
                        
    parser.add_argument(
        '--vizType', type=int, default=3, 
        help='visual representation of electron density')
                             
    return parser
"*****************************************************************************"      
def main():
    # Read input from the command line
    parser = define_command_line_options()
    args = parser.parse_args()   
    
    #Read files
    geometricData, atomList, densityData = readFiles()
        
    # visual representation:
    if args.vizType ==1:
        vizType = volumeRepresentation
    elif args.vizType ==2:
        vizType = contourRepresentation
    elif args.vizType ==3:
        vizType = slicer1Representation
    else:
        sys.exit("unknwon visual type")
    
     
    cubeViz = vizType(geometricData, atomList, densityData,
                             args.bgcolor, args.showCores,args.showDensity,
                             args.animate, args.saveFigure)   
                                             
    cubeViz.vizualize()

"*****************************************************************************"
if __name__ == '__main__':
    print "--------------------------------------"
    print "             CubeViz                  "
    print "--------------------------------------"
    main()
    







