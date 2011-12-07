import numpy as n
import time
import xml.dom.minidom

class project():
    
    def __init__(self,filename):
        self.dom = xml.dom.minidom.parse(filename)
        self.logs = time.asctime()+'\nLogs for project '+ filename
        self.previews = self.getNode(self.dom,'preview')
    
    # set path to data directory
        self.path = self.getNode(self.dom,'path')
        
    # set image files list
        self.falbedo = self.getNode(self.dom,'albedo')
        self.flst = self.getNode(self.dom,'lst')
        self.falbedo = self.getNode(self.dom,'albedo')
        self.fndvi = self.getNode(self.dom,'ndvi')
        self.femissivity = self.getNode(self.dom,'emissivity')
        self.ffc = self.getNode(self.dom,'fc')
        self.fhv = self.getNode(self.dom,'hv')
        self.flai = self.getNode(self.dom,'lai')
        self.RnDaily = self.getNode(self.dom,'RnDaily')
        self.mask = self.getNode(self.dom,'mask')
        
    # read the header of the firts image to collect image dimensions
        self.ncols = float(self.getAttribute(self.dom,'images','samples'))
        self.nrows = float(self.getAttribute(self.dom,'images','lines'))
        self.xllcorner = 0
        self.yllcorner = 0
        self.cellsize = float(self.getAttribute(self.dom,'images','pixel_size'))
        # NoData policy: a nodata can either be a numeric or a non-numeric (NaN), both being converted to floating point
        self.NODATA_value = float(self.getAttribute(self.dom,'images','no_data'))
        
    # set atmospheric forcing mode
        self.atmmode = self.getAttribute(self.dom,'atmosphere','type')
        
        if self.atmmode == '1D':
            self.logs+= '\nSingle atmospheric forcing'
            #self.pblVars = prjIO.prjReadPBL(self.prjfile)
            self.hg = float(self.getNode(self.dom,'hg'))
            self.hr = float(self.getNode(self.dom,'hr'))
            self.lwdw = float(self.getNode(self.dom,'lwdw'))
            self.pg = float(self.getNode(self.dom,'pg'))
            self.pr = float(self.getNode(self.dom,'pr'))
            self.qg = float(self.getNode(self.dom,'qg'))
            self.qr = float(self.getNode(self.dom,'qr'))
            self.swdw = float(self.getNode(self.dom,'swdw'))
            self.tr = float(self.getNode(self.dom,'tr'))
            self.ur = float(self.getNode(self.dom,'ur'))
            self.avgALS = self.hr *10. # Assumption
        if self.atmmode == '2D':
            self.logs+= '\nMultiple atmospheric forcing'
            self.fhg = self.getNode(self.dom,'hg')
            self.fhr = self.getNode(self.dom,'hr')
            self.flwdw = self.getNode(self.dom,'lwdw')
            self.fpg = self.getNode(self.dom,'pg')
            self.fpr = self.getNode(self.dom,'pr')
            self.fqg = self.getNode(self.dom,'qg')
            self.fqr = self.getNode(self.dom,'qr')
            self.fswdw = self.getNode(self.dom,'swdw')
            self.ftr = self.getNode(self.dom,'tr')
            self.fur = self.getNode(self.dom,'ur')
            self.avgALS = 10000. # Assumption            #self.avgplb = mean(pblheight layer)
            
    # set project prefix (for output filenames)
        self.prefix = self.getNode(self.dom,'prefix')
        self.kbMode = self.getNode(self.dom,'kbMode')
        self.broyden = self.getNode(self.dom,'broyden')
        self.iterate = self.getNode(self.dom,'iterate')
        self.pressureUnit = self.getNode(self.dom,'pressureUnit')
        self.pressureLevel = self.getNode(self.dom,'pressureLevel')
        
    def setGrids(self):
        '''Calculate project scales'''
        self.gridNb = [round(self.nrows * self.cellsize / self.avgALS), round(self.ncols * self.cellsize / self.avgALS)]
        self.pixPerGrid = int(min(self.nrows / self.gridNb[0],self.ncols / self.gridNb[1]))
        self.gridNb = [int(self.nrows / self.pixPerGrid),int(self.ncols / self.pixPerGrid)]
        self.imgDims = [self.gridNb[0] * self.pixPerGrid, self.gridNb[1] * self.pixPerGrid]
        
        self.logs+= '\nALS grids: ' + str(self.gridNb) + '\nGrid size (m): ' + str(self.pixPerGrid * self.cellsize) + '\nImage dimensions (px): ' + str(self.imgDims)
        

        
    #def setAtmLayer(self):
    #    '''Creates atmospheric forcing layers from single column / multiple forcing'''
    #    if self.atmmode == 'single':
            
    def getNode(self,nodes,node):
        """Return the data element of a given node of a dom object"""
        try:
            return nodes.getElementsByTagName(node)[0].firstChild.data
        except IndexError:
            return "undef"
        
    def getAttribute(self,nodes,node,attributeName):
        """Return the attribute of a given node of a dom object"""
        try:
            return nodes.getElementsByTagName(node)[0].attributes[attributeName].value
        except KeyError:
            return "undef"
        
    def read(self,fimage):
        """Return an array of an fimage file reference"""
        if fimage <> "undef":
            img = n.fromfile(self.path+fimage,dtype=n.float32)
            img = img.reshape([self.nrows,self.ncols])
            #img = n.load(self.path+fimage)
            return img

    def writeRaw(self,iContent,iName):
        #mybin=open(self.path+self.prefix+iName+'.bin','w')
        #mybin.write(iContent)
        #mybin.close()
        n.save(self.path+self.prefix+iName,iContent)
