# David Lloyd 2020-2023
# This is a simple tool to make blood vessels for a model

import sys, os, csv

import s4l_v1.document as document
import s4l_v1.analysis as analysis
import s4l_v1.analysis.viewers as viewers
import s4l_v1.model as model
from s4l_v1.model import Vec3	
import XCoreModeling
import s4l_v1.simulation.neuron as neuron

import sys, os
from enum import Enum
import s4l_v1.model as model

um2mm = 0.001
length = 1



radius = 100*um2mm
centerLocation = Vec3(-75*um2mm,100*um2mm,0)
p1 = centerLocation + Vec3(0,0,length/2)
p2 = centerLocation + Vec3(0,0,-length/2)
bv = model.CreateSolidCylinder(p1, p2, radius)
bv.Name = "SP1 Vein Wall"

um2mm = 0.001
radius = 75*um2mm
centerLocation = Vec3(-75*um2mm,100*um2mm,0)
p1 = centerLocation + Vec3(0,0,length/2)
p2 = centerLocation + Vec3(0,0,-length/2)
bv = model.CreateSolidCylinder(p1, p2, radius)
bv.Name = "SP1 Vein"

radius = 125*um2mm
centerLocation = Vec3(0,-150*um2mm,0)
p1 = centerLocation + Vec3(0,0,length/2)
p2 = centerLocation + Vec3(0,0,-length/2)
bv = model.CreateSolidCylinder(p1, p2, radius)
bv.Name = "SP1 Artery Wall"

radius = 100*um2mm
centerLocation = Vec3(0,-150*um2mm,0)
p1 = centerLocation + Vec3(0,0,length/2)
p2 = centerLocation + Vec3(0,0,-length/2)
bv = model.CreateSolidCylinder(p1, p2, radius)
bv.Name = "SP1 Artery"

radius = 500*um2mm
centerLocation = Vec3(0,0,0)
p1 = centerLocation + Vec3(0,0,length/2)
p2 = centerLocation + Vec3(0,0,-length/2)
bv = model.CreateSolidCylinder(p1, p2, radius)
bv.Name = "Fat"
