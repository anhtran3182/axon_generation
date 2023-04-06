# David Lloyd, 2020-2023
# This script builds a nerve model from a set of contours and a trajectory file
# The contours are used to define the nerve fiber path, and the trajectory file is used to define the nerve fiber diameter

import csv
import json
import os
import sys
import time
from enum import Enum
from pathlib import Path as path
from random import gauss, uniform
import time
import numpy as np
from pathlib import Path as path
import os
from bs4 import BeautifulSoup
import pandas as pd
from XCore import Color
import numpy as np
import pandas as pd
import s4l_v1.analysis as analysis
import s4l_v1.analysis.viewers as viewers
import s4l_v1.document as document
import s4l_v1.model as model
import s4l_v1.simulation.neuron as neuron
import s4l_v1.units as units
import XCoreModeling
from s4l_v1.model import Vec3
from XCoreModeling import (kBodyClashContained_Abuts, kBodyClashInterlock,
                           kBodyClashNone)
from collections import OrderedDict
import seaborn as sns

um2mm = 0.001

_CFILE = os.path.abspath(sys.argv[0] if __name__ == '__main__' else __file__ )
_CDIR = os.path.dirname(_CFILE)


class NeuronFiber(Enum):
    Myelinated = 1
    Unmyelinated = 2
    

class Contour:
    """
    A class representing a contour object.

    Attributes:
    -----------
    label: str
        The label of the contour.
    name: str
        The name of the contour.
    ptslist: list
        A list of points representing the contour.
    xlist: list
        A list of x-coordinates of the contour.
    ylist: list
        A list of y-coordinates of the contour.
    zlist: list
        A list of z-coordinates of the contour.
    color: str
        The color of the contour in hex format.
    pixpermicron: float
        The number of pixels per micron.
    micronsperpix: float
        The number of microns per pixel.
    xoffset: float
        The x-offset of the contour.
    yoffset: float
        The y-offset of the contour.
    centroid: array_like
        The centroid of the contour.
    area: float
        The area of the contour.
    c_dia_area: float
        The circular diameter of the contour.
    """
    def __init__(self, name=None, label=None, ptslist = None, xlist=[], ylist=[], zlist=[], color="#ffffff", pixpermicron=None, micronsperpix=None, xoffset=None, yoffset=None):
        """
        Constructor for Contour class.

        Parameters:
        -----------
        name: str
            The name of the contour.
        label: str
            The label of the contour.
        ptslist: list
            A list of points representing the contour.
        xlist: list
            A list of x-coordinates of the contour.
        ylist: list
            A list of y-coordinates of the contour.
        zlist: list
            A list of z-coordinates of the contour.
        color: str
            The color of the contour in hex format.
        pixpermicron: float
            The number of pixels per micron.
        micronsperpix: float
            The number of microns per pixel.
        xoffset: float
            The x-offset of the contour.
        yoffset: float
            The y-offset of the contour.
        """

        self.label=label
        self.name = name
        self.ptslist = ptslist
        
        if ptslist is not None:
            self.split_ptslist
        else:
            self.xlist = xlist
            self.ylist = ylist
            self.zlist = zlist
        self.color = color
        
        if pixpermicron is not None:
            self.pixpermicron=pixpermicron
            self.micronsperpix=1/pixpermicron
        elif micronsperpix is not None:
            self.micronsperpix=micronsperpix
            self.pixpermicron=1/micronsperpix
            
        self.xoffset = self.micronsperpix*xoffset
        self.yoffset = self.micronsperpix*yoffset
        
    def from_xml_ele(self, xml_cnt):
        self.label = xml_cnt.get('name')
        self.color = xml_cnt.get('color')
    
        pnts = xml_cnt.find_all('point')
        self.ptslist=[]
        for pt in pnts:
            self.ptslist.append([self.micronsperpix*(float(pt.get('x'))-self.xoffset), self.micronsperpix*(float(pt.get('y'))+self.yoffset),self.micronsperpix*float(pt.get('z'))])
        
        try:
            ptsarr = np.asarray(self.ptslist)    
            self.centroid = np.mean(ptsarr, axis=0)
            
            x,y = ptsarr[:,0], ptsarr[:,1]
            self.area = 1/2 * np.sum(x*np.roll(y, 1) - y*np.roll(x, 1))
            self.c_dia_area = 2*np.sqrt(self.area/np.pi)
        except:
            pass
            
    def from_csv(self, csvpth):
        csvpath = path(csvpth)
        self.label = csvpath.stem
        pntsdf = pd.read_csv(csvpath)
        self.ptslist=[]
        
        #I'm iterating here because i just truly don't want to deal with converting two pandas columns into a 3xn list
        for index, row in pntsdf.iterrows():
            self.ptslist.append([row['X'], row['Y'], 0])
            
        #to close the contour
        self.ptslist.append(self.ptslist[0])
        
            
    def split_ptslist(self):
        if self.ptslist is not None:
            self.xlist = [pt[0] for pt in self.ptslist]
            self.ylist = [pt[1] for pt in self.ptslist]
            self.zlist = [pt[2] for pt in self.ptslist]
            
    def hex_to_rgb(self, value):
        """Return (red, green, blue) for the color given as #rrggbb."""
        value = value.lstrip('#')
        lv = len(value)
        return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))

    def rgb_to_hex(self, red, green, blue):
        """Return color as #rrggbb for the given color values."""
        return '#%02x%02x%02x' % (red, green, blue)
    
def xml2dict(xmlpath, xoffset, yoffset):
    with open(xmlpath, 'r') as f:
        data = f.read()
    
        # Passing the stored data inside
        # the beautifulsoup parser, storing
        # the returned object
        Bs_data = BeautifulSoup(data, "xml")
        if Bs_data is None:
            print()
            print(data)
            print()
        imele = Bs_data.find('image')
        if imele is None:
            imele = Bs_data.find('images')
        mpp = float(imele.find('scale').get('x'))
        
        # Finding all instances of tag
        # `contour`
        b_cnt = Bs_data.find_all('contour')

        cnt_dic = dict()

        for idx, cnt in enumerate(b_cnt):
            newcnt = Contour(micronsperpix=mpp, xoffset=xoffset, yoffset=yoffset)
            newcnt.from_xml_ele(cnt)
            
            cnt_dic[idx] = newcnt
            
    return(cnt_dic)

def pts2model(ptsfile, trajectory, group=None, name=None, posifile=None, negifile=None):
    ptsfile = path(ptsfile)
    
    if name is None:
        name = ptsfile.stem
        
    if negifile is not None:
        negifile = path(negifile)
        
        negiseries = pd.Series(pd.read_csv(negifile)['Labels'].values)
        negiseries = negiseries.drop(negiseries.index[0])
        # print(negiseries)
        # negiseries=negiseries.squeeze()
    else:
        negiseries=None
        
    if posifile is not None:
        posifile = path(posifile)
        posiseries = pd.Series(pd.read_csv(posifile)['Labels'].values)
        posiseries = posiseries.drop(posiseries.index[0])
        # print(posiseries)
        # posiseries=posiseries.squeeze()
    else:
        posiseries=None
        
    if group is None:
        shapegroup = model.CreateGroup(name)
        
    pts_df = pd.read_csv(ptsfile)
    tspline = model.CreateSpline(trajectory)
    #conversion into vector list. I wish this were better, but... sigh. Although it is way more efficient than the other method
    # skfirstclm = False
    for column in pts_df.iteritems():
        # if not skfirstclm:
        #     skfirstclm=True
        #     continue
        #grab name and data
        shapename = column[0]
        shapeseries = column[1].dropna()
        shapeseries = shapeseries.drop_duplicates()
        
        #now we'll check our label files and disregard columns in negiseries or not in posiseries, but only if they exist
        
        label = shapename.split('_')[-1]
        
        
        try:
            #needed to remove first stupid line lol, which doesn't convert
            shapepts = [pt.strip('][').split(', ') for pt in shapeseries]
        except Exception as e:
            print(f"{shapename}: {e}")
            continue
        
        # print(label, type(label))
        if posiseries is not None:
            if int(label) not in posiseries.values:
                continue
        elif negiseries is not None:
            if int(label) in negiseries.values:
                continue
            
        
        # shapepts = shapepts[1:]
        for idx, pt in enumerate(shapepts):
            for idy, coord in enumerate(pt):
                try:
                    shapepts[idx][idy] = float(coord) * um2mm
                except:
                    print(f"Can't convert to float: {coord}")
                    del shapepts[idx][idy]
        # print(shapepts)
        
        shape_veclist = [Vec3(pt) for pt in shapepts]
        shape_veclist.append(shape_veclist[0])
        
        #now we make the shape and extrude it
        #TODO: Add label parsing and axon creation
        fiber = model.CreatePolyLine(shape_veclist)
        fiber = XCoreModeling.SweepAlongPath(fiber, tspline, make_solid=True)
        

        if "epi" in shapename.lower():
            fiber.Name = shapename + "_Fascicle"
        elif "blv" in shapename.lower():
            fiber.Name = shapename + "_Blood"
        elif "myl" in shapename.lower():
            fiber.Name = shapename + "_Myelin"
        else:
            fiber.Name = shapename.lower() + "_Nerve"
        shapegroup.Add(fiber)
            
    tspline.Delete()
    return shapegroup
          
def makeNeurons_pts(measfile, traj, name, sim, axon_group=None, axonType=NeuronFiber.Unmyelinated, mylfile=None, neuron_model=None):
    """
    Creates neuron models from point cloud data and returns a list of NeuronFiber objects.

    Parameters:
        measfile (str): Path to the point cloud measurement file.
        traj (list): List of tuples containing the trajectory points of the neuron.
        name (str): Name of the neuron.
        sim (Simulation): Simulation object used to simulate the neuron models.
        axon_group (str): Name of the axon group. Defaults to None.
        axonType (NeuronFiber): Type of axon fiber. Defaults to NeuronFiber.Unmyelinated.
        mylfile (str): Path to the myelination file. Defaults to None.
        neuron_model (str): Name of the neuron model to use. Defaults to None.

    Returns:
        list: List of NeuronFiber objects representing the neuron models.
    """

    if axon_group is None:
        axon_group = model.CreateGroup(name)
        
    measdf = pd.read_csv(measfile)
    cxlist = list(measdf['Centroid_x'])[:-1]
    cxlist = [float(i)*um2mm for i in cxlist]
    cylist = list(measdf['Centroid_y'])[:-1]
    cylist = [float(i)*um2mm for i in cylist]
    czlist = list(measdf['Centroid_z'])[:-1]
    czlist = [float(i)*um2mm for i in czlist]
    dialist = list(measdf['Circular Diameter (Area)'])[:-1]
    dialist = [float(i)*um2mm for i in dialist]
    labelist = list(measdf.iloc[:,0])[:-1]
    labelist = [str(l) for l in labelist]
    
    if mylfile:
        myldf = pd.read_csv(mylfile)
        myldialist = list(myldf['Circular Diameter (Area)'])[:-1]
    
    idx = 0
    for cx, cy, cz, diameter in zip(cxlist, cylist, czlist, dialist):
        diameter = diameter*um2mm
        temp_points = list()
        for point in traj:
            #print("{0},{1}, 0".format(x_offset, y_offset));
            newpoint = Vec3(point)
            newpoint += Vec3(cx, cy, cz)
            temp_points.append(newpoint)
            
        newspline = model.CreateSpline(temp_points);
        newspline.Name = labelist[idx]+"_ax"
        
        if sim is not None:
            if axonType is NeuronFiber.Unmyelinated:
                # new_neuron_settings = neuron.SundtSennSettings();
                # new_neuron_settings = model.SundtNeuronProperties()
                new_neuron_settings = model.TigerholmNeuronProperties()
                # new_neuron_settings = model.SchildNeuronProperties()
                
                if not neuron_model or neuron_model == 'Tigerholm':
                    new_neuron_settings = model.TigerholmNeuronProperties()
                elif neuron_model == 'Sundt':
                    new_neuron_settings = model.SundtNeuronProperties()
                elif neuron_model == 'Schild':
                    new_neuron_settings = model.SchildNeuronProperties()
                else:
                    raise("Unknown Model for Unmyelinated Axon\nValid choices are Sundt, Schild, or Tigerholm (default)")
                diameter = diameter * 1e6
                # print(diameter)
                # expects in raw microns, so we need to multiply by a million
                new_neuron_settings.AxonDiameter = diameter
            else:
                if mylfile:
                    myldia = float(myldialist[idx])
                    if not neuron_model or neuron_model == 'SENN':
                        new_neuron_settings = model.SennNeuronProperties()
                    elif neuron_model == 'Sweeney':
                        new_neuron_settings = model.SweeneyNeuronProperties();
                    elif neuron_model == 'Small_MRG':
                        new_neuron_settings = model.SmallMrgNeuronProperties()
                    elif neuron_model == 'Rat':
                        new_neuron_settings = model.RatNeuronProperties()
                    else:
                        raise("Unknown Model for Myelinated Axon\nValid choices are Sweeney, Small_MRG, Rat, or SENN (default)")
                    
                    # diameter = diameter * 1e6
                    new_neuron_settings.AxonDiameter = myldia
                    # new_neuron_settings.AxonDiameterAtNode = diameter
                    # new_neuron_settings.NodeDiameter = diameter
                else:
                    raise("No myelin measurement file passed for myelinated fiber")

            components = newspline;
                
            new_neuron_settings.Name = newspline.Name + "_neuron"
            newneuron = model.CreateAxonNeuron(components, new_neuron_settings)
            axon_group.Add(newneuron)
            # sim.Add(newneuron)
        
        # axon_group.Add(newspline)
        newspline.Delete()
        idx = idx+1
        
        # Setup Settings
        setup_settings = sim.SetupSettings
        setup_settings.PerformTitration = True
        setup_settings.DepolarizationDetection = setup_settings.DepolarizationDetection.enum.Threshold
        
        # Neuron Settings
        entities = model.AllEntities()
        axons = [ent for ent in entities if '_neuron' in ent.Name]
        # for ax in axon_group.Entities:
        #     neuron_settings = neuron.AutomaticAxonNeuronSettings(ax)
            # sim.Add(ax)
        #     axon_group.Add(ax)
        #     axon_group.Add(neuron_settings)
    return axon_group

def unique_pts_list(ptslist):
    new_list = []
    for pt in ptslist:
        if pt not in new_list:
            new_list.append(pt)
    return new_list

def xml2model(xmlfile, trajectory, name=None, xoffset=0, yoffset=0):
    ptsfile = path(xmlfile)
    cnt_dict = xml2dict(ptsfile, xoffset=xoffset, yoffset=yoffset)
    grp_dict = {}
    # colordict = {}
    # palette = list(sns.color_palette(None, 20))

    if name is None:
        name = ptsfile.stem

    #shared trajectory spline
    tspline = model.CreateSpline(trajectory)

    #for each contour in the xml
    for cntitem in cnt_dict.items():
        try:
            shapename = cntitem[0]
            cnt = cntitem[1]

            ptslist = cnt.ptslist
            # print(len(ptslist))
            grpnm = cnt.label

            #if it doesn't have a group, give it one and store it in the group dict
            if ('axon' in grpnm.lower()) and ('outer' in grpnm.lower()):
                if  grpnm not in grp_dict.keys():
                    grp_dict[grpnm] = model.CreateGroup(grpnm)
                    # colordict[grpnm] = Color(palette[0][0], palette[0][1], palette[0][2], 1)
                    # palette.pop(0)

                #convert ptslist into a closed list of Vec3 points
                # ptslist = np.unique(ptslist, axis=1)
                # plen = len(ptslist)
                # setpt = set(ptslist)
                # setlen = len(setpt)
                # if plen > setlen:
                ptslist = unique_pts_list(ptslist)  
                shape_veclist = [Vec3(pt) for pt in ptslist]
                shape_veclist.append(shape_veclist[0])


                fiber = model.CreatePolyLine(shape_veclist)
                fiber = XCoreModeling.SweepAlongPath(fiber, tspline, make_solid=True)

                fiber.Name = f"{shapename}_{grpnm}"
                cntcolor = cnt.hex_to_rgb(cnt.color)
                fiber.Color = Color(cntcolor[0], cntcolor[1], cntcolor[2], 1)
                
                grp_dict[grpnm].Add(fiber)
                
                if ('axon' in grpnm.lower()) and ('inner' not in grpnm.lower()):
                    
                    temp_points = list()
                    # newtspline = model.CreateSpline(trajectory)
                    for point in trajectory:
                        #print("{0},{1}, 0".format(x_offset, y_offset));
                        newpoint = Vec3(point)
                        newpoint += Vec3(cnt.centroid[0], cnt.centroid[1], cnt.centroid[2])
                        temp_points.append(newpoint)
                        
                    newspline = model.CreateSpline(temp_points)
                    
                    if 'unmyelinated' in grpnm.lower():
                        if  grpnm+'_UMF' not in grp_dict.keys():
                            grp_dict[grpnm+'_UMF'] = model.CreateGroup(grpnm+'_UMF')
                            
                        newspline.Name =  f"{shapename}_{grpnm}_UMF"
                            
                        # new_neuron_settings = neuron.SundtSennSettings();
                        # new_neuron_settings = model.SundtNeuronProperties()
                        new_neuron_settings = model.TigerholmNeuronProperties()
                        # new_neuron_settings = model.SchildNeuronProperties()
                        diameter = float(cnt.c_dia_area) #* 1e6

                        # expects in raw microns, so we need to multiply by a million
                        new_neuron_settings.AxonDiameter = diameter
                        
                        components = newspline
                        
                        new_neuron_settings.Name = newspline.Name + "_Nrn"
                        newneuron = model.CreateAxonNeuron(components, new_neuron_settings)
                        grp_dict[grpnm+'_UMF'].Add(newneuron)
                        newspline.Delete()
                        
                        
                    elif 'myelinated' in grpnm.lower():
                        if  grpnm+'_MYF' not in grp_dict.keys():
                            grp_dict[grpnm+'_MYF'] = model.CreateGroup(grpnm+'_MYF')
                            
                        newspline.Name =  f"{shapename}_{grpnm}_MYF"
                            
                        new_neuron_settings = model.SennNeuronProperties()
                        #new_neuron_settings = model.SweeneyNeuronProperties()
                        #new_neuron_settings = model.SmallMrgNeuronProperties()
                        #new_neuron_settings = model.RatNeuronProperties()
                        new_neuron_settings.AxonDiameter = float(cnt.c_dia_area) * 1e6
                        
                        components = [newspline]
                            
                        new_neuron_settings.Name = newspline.Name + "_Nrn"
                        print(type(cnt.c_dia_area))
                        # print(new_neuron_settings)
                        newneuron = model.CreateAxonNeuron(components, new_neuron_settings)
                        grp_dict[grpnm+'_MYF'].Add(newneuron)
                        newspline.Delete()
                    
        except Exception as e:
            # print(f'Error in {shapename}_{grpnm}, {e}')
            # print(len(ptslist))
            # continue
            raise(e)

    tspline.Delete()
    return grp_dict
    
def makeElectrodes(centerLocation, spacing = 0.5):
    spacingVec = Vec3(spacing)
    length = 3.0
    radius = 0.127/2
    electrode1 = model.CreateSolidCylinder(Vec3(centerLocation + Vec3(0, length/2.0, spacing/2.0)),
                                           centerLocation + Vec3(0, -length/2, spacing/2),
                                           radius)
                              
    electrode2 = model.CreateSolidCylinder(centerLocation + Vec3(0, length/2.0, -spacing/2.0),
                                           centerLocation + Vec3(0, -length/2.0, -spacing/2.0),
                                           radius)
                                           
    electrode1.Name = "Electrode 1"
    electrode2.Name = "Electrode 2"
 

    
