#Wrapper to run AutoBend on multiple models and with repeats from command line
#run from command line using mayapy
#cd "C:/Program Files/Autodesk/Maya2020/bin" #set to maya location
#run mayapy
#mayapy.exe "Z:/lab/NSF synapsid project data share/Maya/Wrapper_vert_bending.py" #set to wrapper location

#####################################################################
#code begins

#initialize so can run directly from mayapy
import maya.standalone
maya.standalone.initialize(name='python')

#Wrapper function to run vertebral bending on a series of maya models
#arrange maya files in a single folder
#source functions
import sys
sys.path.append('Z:\\lab\\NSF synapsid project data share\\Maya')#pathname
import AutoBend as ab#source automated bending functions
import os
import maya.cmds as cmds
from maya.api.OpenMaya import MVector, MMatrix, MPoint
import numpy as np
import maya.api.OpenMaya as om2

#supress error mgs
cmds.scriptEditorInfo(sw=True,se=True)

#set directories
dir='Z:/lab/NSF synapsid project data share/Maya/Validation models'
outdir = 'Z:/lab/NSF synapsid project data share/Maya/Validation results'
models = os.listdir(dir)

for mod in range(0,len(models)):
    #open the file
    modi = dir + '/' + models[mod]
    cmds.file(modi, i=True)
    
    #Count the joints
    joints = cmds.ls('axes*',type='surfaceShape')#NOT working
    joints = list(map(lambda sub:int(''.join( 
          [ele for ele in sub if ele.isnumeric()])), joints))
    njoint = len(joints)
    
    #Measure areas and spacing
    Area_all= []
    spacing_all=[]
    for i in joints:
        ant=str(i)
        post=str(i+1)
        Ant = 'V'+ant+':Mesh'
        Post = 'V'+post+':Mesh'
        COR1 = 'COR_locator'+ant
        
        #vertebral area averaged by joint
        v_area_ant = cmds.polyEvaluate(Ant,a=True)
        v_area_post = cmds.polyEvaluate(Post,a=True)
        v_area = (v_area_ant+ v_area_post)/2
        Area_all.append(v_area)
          
       
        #intervertebral spacing
        ddist = ab.wDist("cent_D"+ant,"Pcent_D"+ant, COR1)
        vdist = ab.wDist("cent_V"+ant,"Pcent_V"+ant, COR1)
        ldist = ab.wDist("cent_L"+ant,"Pcent_L"+ant, COR1)
        rdist = ab.wDist("cent_R"+ant,"Pcent_R"+ant, COR1)
        dis = np.mean([ddist,vdist,ldist,rdist])
        spacing_all.append(dis)
        
    
    v_area = np.mean(Area_all)
    j_space = np.mean(spacing_all)
    
    #table headings
    bendres = [['joint', 'intersection threshold', 'joint spacing',  'strain', 'Translation',
    'starting int',  'TorL','TorL_type', 'TorR','TorR_type', 'LatL','LatL_type', 'LatR','LatR_type', 'DFL','DFL_type',
     'VFL', 'VFLtype']]
    
    #vary joint spacing 
    for js in (0.1,-0.1):
    
        for m in joints:
            ant=str(m)
            post=str(m+1)
            new_axis = 'COR_locator'+ant
            Post = 'V'+post+':Mesh'
            spacing_i=spacing_all[joints.index(m)]  
            space = spacing_i*js
            
            #joint spacing
            ab.adjust_spacing(new_axis, Post, space)
        
        #vary intersection threshold
        for it in (0.0025, 0.005):
          
            #vary centrum and zygapophysis strain
            for strain in (0.45, 0.55):
                
                #set translation
                for Transl in [True, False]:
                     
                        
                    #run experiments
                    for i in joints:
            
                        jointname = str(i) + '_' + str(i+1)
                        test = ab.digital_bending(joint=i, angle=0.5, it=it, strain=strain, v_area=v_area, 
                        ZygTest=False, CentTest=True, Transl=Transl)
                        output = [jointname, it, js, strain, Transl]
                        output = output + test
                        bendres.append(output)
    
    
    #Save results - change directory depending how is mapped
    listresults = bendres
    namei = models[mod]
    namei= namei.replace('.mb','')
    fname = outdir + '/' + namei + '.txt'
    
    with open(fname, 'w') as filehandle:
        for listitem in listresults:
            filehandle.write('%s\n' % listitem)
        
    cmds.file(modified=0)#shut down the file
    cmds.file(new=True)