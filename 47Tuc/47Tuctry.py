
# coding: utf-8

# In[1]:


import pyspeckit as ps
from astropy.io import fits
from bokeh.layouts import column

import os

import astropy.units as u
from astropy.io import fits
from spectral_cube import SpectralCube
import matplotlib.pyplot as plt
from bokeh.plotting import output_notebook, figure, show
from bokeh.models import HoverTool, tools
from bokeh.models import Span, Label, Arrow, NormalHead
import numpy as np
#output_notebook()
#get_ipython().magic('matplotlib widget')

import warnings
warnings.filterwarnings('ignore')

#f = fits.open('../Observation1/data/FixedAstrometryDATACUBE_FINAL.fits')

#c = ps.Cube(f[1])
#cube1= SpectralCube.read('../Observation1/data/FixedAstrometryDATACUBE_FINAL.fits',hdu=1)
from astropy.coordinates import SkyCoord, Angle
import astropy.units as u




##Emission Lines
class line(object):
    def __init__(self,name):
        self.name = name

        
        
##Emission Lines
class specto(object):
    def __init__(self,wave,flux,ra,dec,apperture):
        self.wave = wave
        self.flux = flux
        self.ra = ra
        self.dec = dec
        self.apperture = apperture
        
    def plotline(self,title='b'):
        p = self.plot(title,showplot=False)
        linetoplot = self.flux[self.lineindex]
        wavetoplot = self.wave[self.lineindex]
        p.line(wavetoplot,linetoplot,line_color='red')
        show(p)
            
    def plot(self,title='b',showplot=True):
        "Plot with BOkeh"
        x = self.wave
        y = self.flux

        p = figure(plot_width=900, plot_height=500, title=title,active_drag='box_zoom', active_scroll='wheel_zoom')


        p.line(x,y)



        #Tool to get wavelength
        hover2 = HoverTool(
                tooltips=[
                    ("(x,y)", "($x{1}, $y)"),
                ]
            )



        p.add_tools(hover2)
        if showplot:
            show(p)
        return p
        


        
def getspectrafromregion(spectralcube,ra,dec,aperture,xmin=6400.,xmax=6650.):
    """Get a dec,ra in hh:mm:ss and aperture in arcsecond"""
    
    
    sourcename = 'tempspectra'
    region = 'fk5; circle({}, {}, {}")'.format(ra,dec,aperture)

    filename = '{}.fits'.format(sourcename)

    #Spectrum
    subcube = spectralcube.subcube_from_ds9region(region)  
    spectrum = subcube.sum(axis=(1, 2)) 


    if os.path.isfile(filename):
        os.remove(filename)
    spectrum.write(filename)
    
    #print(spectrum)
    pyspec = ps.Spectrum(filename)

    
    os.remove(filename)
    
    
    pyspec.crop(xmin, xmax, unit='angstrom') 
    pyspec.baseline(xmin=xmin, xmax=xmax,exclude=[6520,6600,6660,6700],order=2,subtract=True)

    inds = np.argsort(pyspec.xarr)
    xp2 = pyspec.xarr.value[inds]
    yp2 = pyspec.data[inds]


    x = np.array(xp2)
    y = np.array(yp2)
    
    return specto(x,y,ra,dec,aperture)

def findlines(specto,minwave,maxwave,toleranceinsigma):
    errorindex = np.where( (minwave < specto.wave) &  (maxwave > specto.wave) )[0]
    errorflux = specto.flux[errorindex]
    errorstd = np.std(errorflux)
    lineindex = np.where(specto.flux > toleranceinsigma* errorstd)[0]
    if len(lineindex) > 2:
        #print('Found Line')
        specto.lineindex = lineindex
        specto.line = True
    else:
        #print('Booo')
        specto.line = False
    

def getlatloglist(spectralcube,latsepa=0.5,logsepa=0.5):
    cubelatextrema = spectralcube.latitude_extrema
    cubelongextrema = spectralcube.longitude_extrema
    seplat = Angle(latsepa*u.arcsec)
    seplog = Angle(logsepa*u.arcsec)

    latlist = np.arange(Angle(cubelatextrema[0]).deg,
                        Angle(cubelatextrema[1]).deg,
                        seplat.deg)

    loglist = np.arange(Angle(cubelongextrema[0]).deg,
                        Angle(cubelongextrema[1]).deg,
                        seplog.deg)

    latloglist = []
    for la in latlist:
        templa = Angle(la*u.deg)    
        for log in loglist:
            templog = Angle(log*u.deg)
            latloglist.append([templa.to_string(unit=u.degree, sep=':'),templog.to_string(unit=u.hour, sep=':')])
    return latloglist


def makeregionfromlist(listcoord,nameregfile):
    with open(nameregfile,'w') as fileone:
        for coord in listcoord:#print(line)
                text='fk5; circle({0},{1},{2}") # color=red text={3} \n'.format(coord[0],coord[1],coord[2],coord[3])
                fileone.write(text)

#c = ps.Cube(cube=cube1)


# # mini try first
# 
# Mini datacube:
# #mini = spectralcube[:,135:160,190:220]
# #mini

# In[2]:


#cube = 'mini.fits'
cube = 'ADP.2016-06-21T05:33:40.720.fits'
spectralcube= SpectralCube.read(cube,hdu=1)
#listlatlog = getlatloglist(spectralcube,latsepa=0.5,logsepa=0.5)
listlatlog = getlatloglist(spectralcube,latsepa=0.5,logsepa=0.5)


# In[ ]:


listcoord = []

for number,vaina in enumerate(listlatlog):
    ra = vaina[0]
    dec = vaina[1]
    #print(ra)
    try:
        tryspec = getspectrafromregion(spectralcube,dec, ra, 0.200,xmin=6400.,xmax=6650.)
        findlines(tryspec,minwave=6400,maxwave=6500,toleranceinsigma=5.)
        if tryspec.line:
        #tryspec.plotline(title='{}, {}, {}'.format(tryspec.ra,tryspec.dec,number))
            listcoord.append([tryspec.ra,tryspec.dec,tryspec.apperture,'{'+str(number)+'}'])

    except:
        #print('a')
        pass
    


# In[5]:


len(listlatlog)


# In[ ]:


makeregionfromlist(listcoord,'foundlines.reg')


# In[ ]:


len(listcoord)


# In[ ]:



