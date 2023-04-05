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

import build_nerve

def main(data_path=None, project_dir=None):
    nameprefixlist = ["BSN"]#, "SP2"]#, "SP3", "SP4"]

    #### GT ####
    # sp = path(r"C:\Users\dalloyd\Documents\Manuscripts\SPN\Sp1\12_gt_meas")
    # epifile = sp / "Epineurium_GT_cnt_pts.csv"

    sp = path(r"C:\Users\dalloyd\Documents\Capstone 2023 Modeling")
    epifile = sp / "22209-20x_Countour_XYZ.csv"
        
    #Make simulation object
    neurosim = neuron.Simulation()
    neurosim.Name = "Neuro"

    traj = list([Vec3(0,0,0.5), Vec3(0,0,0), Vec3(0,0,-0.5)])

    for idx, nameprefix in enumerate(nameprefixlist):
    #make group for this prefix
        prefixgroup = model.CreateGroup(nameprefix)
        grplist = []
        grp = build_nerve.pts2model(epifile, traj, name=f'EPI_{nameprefix}')
        grplist.append(grp)
        
        for elm in grplist:
            prefixgroup.Add(elm)
        
        trajnew = traj
        grplist = []
        # grp = build_nerve.makeNeurons_pts(myfmeas,trajnew, name=f"MYF_ax_{nameprefix}", sim=neurosim, axonType=build_nerve.NeuronFiber.Myelinated, mylfile=mylmeas, neuron_model='SENN')
        # grplist.append(grp)
        
        for elm in grplist:
            prefixgroup.Add(elm)
        
            
        
        
    plate1 = model.CreateSolidBlock( Vec3(-3,-2,-.5), Vec3(1,2,.5) )
    plate1.Name = 'Saline'
    # plate2 = model.CreateSolidBlock( Vec3(-1,-3.81487,-.5), Vec3(1,-0.7,.5) )
    # plate2.ApplyTransform( Rotation(Vec3(1,0,0),math.radians(180)) )
    # plate2.ApplyTransform( Rotation(Vec3(0,0,1),math.radians(-90)) )
    # plate2.ApplyTransform( Translation(Vec3(-3.68474,-1,0)) )

    # plate2.Name = 'Fat'
    
    build_nerve.makeElectrodes(Vec3(.20,0,0))
    # document.AllSimulations.Add(neurosim)
    

if __name__ == '__main__':
    t = time.time()
    main()
    print(f"{time.time()-t} Seconds Elapsed")