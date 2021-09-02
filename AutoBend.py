import maya.cmds as cmds
from maya.api.OpenMaya import MVector, MMatrix, MPoint
import numpy as np
import maya.api.OpenMaya as om2

##########################################################################################################################
#Make objects

#Function to set up scene ready for automated bending. Creates objects needed for each joint and sends them to the joint center. Alignment of axes and placement of locators must be done by hand.

#COR - COR object name
#base_axis - axis object
#new_axis - name of new axis to be created
#Post - name of posterior vertebra
#Post1 - name of duplicate
#Ant - name of anterior vertebra
#Ant1 - name of duplicate
#Ant2 - name of second duplicate
#Zyg_ant - name of locator on zyapophyses of anterior vertebra
#Zyg_post - name of locator on Zygapophyses of posterior vertebra
#cent_X - anterior centrum landmark
#Pcent_X - posterior centrum landmark


def make_objects(ant):
    
    #define objects
    ant=str(i)
    COR = 'COR'+ant+':Mesh'#may or may not need mesh here
    base_axis = 'axes'+ant
    Zyg_ant = 'Zyg_ant'+ant
    Zyg_post = 'Zyg_post'+ant
    cent_D = 'cent_D' + ant
    cent_V = 'cent_V' + ant
    cent_L = 'cent_L' + ant
    cent_R = 'cent_R' + ant
    Pcent_D='Pcent_D'+ant
    Pcent_V='Pcent_V'+ant
    Pcent_L='Pcent_L'+ant
    Pcent_R='Pcent_R'+ant
    
    #Import verts/COR/and create giant axes - must do by hand
    cmds.xform(COR, cp= True)#center pivot
    #move axes
    cmds.matchTransform(base_axis,COR)#move to the COR

    #make locators for zygapophyses
    cmds.spaceLocator(n = Zyg_ant)
    cmds.matchTransform(Zyg_ant,COR)
    cmds.spaceLocator(n = Zyg_post)
    cmds.matchTransform(Zyg_post,COR)
    
    #make locators for centrum
    cmds.spaceLocator(n = cent_D)
    cmds.matchTransform(cent_D,COR)
    cmds.spaceLocator(n = cent_V)
    cmds.matchTransform(cent_V,COR)
    cmds.spaceLocator(n = cent_L)
    cmds.matchTransform(cent_L,COR)
    cmds.spaceLocator(n = cent_R)
    cmds.matchTransform(cent_R,COR)
    cmds.spaceLocator(n = Pcent_D)
    cmds.matchTransform(Pcent_D,COR)
    cmds.spaceLocator(n = Pcent_V)
    cmds.matchTransform(Pcent_V,COR)
    cmds.spaceLocator(n = Pcent_L)
    cmds.matchTransform(Pcent_L,COR)
    cmds.spaceLocator(n = Pcent_R)
    cmds.matchTransform(Pcent_R,COR)
    
    #remove COR
    cmds.delete(COR)
    
###########################################################################################
#Function to duplicate vertebrae for the boolean operation

def duplicate_verts(Ant, Post, Ant1, Post1):
    #Duplicate vertebral meshes
    cmds.duplicate(Post, n= Post1)
    cmds.duplicate(Ant, n= Ant1)

    
#############################################################################################    
#Set up one joint

#Function to set up parenting relationships required for automated bending. This function works on either a single joint or the first joint in a series.
#Intersection object is created to measure overlap, then posterior vertebra is parented to anterior for bending. A joint locator at COR is also created.

##COR - COR object name
#base_axis - axis object
#new_axis - name of new axis to be created
#Post - name of posterior vertebra
#Post1 - name of duplicate
#Ant - name of anterior vertebra
#Ant1 - name of duplicate
#Ant2 - name of second duplicate
#Zyg_ant - name of locator on zyapophyses of anterior vertebra
#Zyg_post - name of locator on Zygapophyses of posterior vertebra

def set_up_one_joint(base_axis, new_axis, Post, Post1, Ant, Ant1, Zyg_ant, Zyg_post, inter, cent_D, cent_V, cent_L, cent_R):

    #Set up intersection
    cmds.polyCBoolOp(Ant1, Post1, op=3, n=inter)

    #Set up parenting relationships
    cmds.parent(Ant1,Ant)
    cmds.parent(Zyg_ant,Ant)
    cmds.parent(base_axis,Ant)
    cmds.parent(Post1,Post)
    cmds.parent(Zyg_post,Post)
    
    #centrum objects
    cmds.parent(Cent_D,Ant)
    cmds.parent(Cent_V,Ant)
    cmds.parent(Cent_L,Ant)
    cmds.parent(Cent_R,Ant)

    #Create joint center locator
    cmds.spaceLocator(n = new_axis)
    cmds.matchTransform(new_axis,base_axis)
    cmds.parent(new_axis,base_axis)
    cmds.parent(Post,new_axis)


##########################################################################
#Set up subsequent joints

#Function to set up parenting structure in the case of multiple joints nested hierarchically. Sets up all joints after the first.

##COR - COR object name
#base_axis - axis object
#new_axis - name of new axis to be created
#Post - name of posterior vertebra
#Post1 - name of duplicate
#Ant - name of anterior vertebra
#Ant1 - name of duplicate
#Ant2 - name of second duplicate
#Zyg_ant - name of locator on zyapophyses of anterior vertebra
#Zyg_post - name of locator on Zygapophyses of posterior vertebra

def set_up_sub_joint(base_axis, new_axis, Post, Post1, Ant, Ant1, Zyg_ant, Zyg_post, inter, 
cent_D, cent_V, cent_L, cent_R, Pcent_D, Pcent_V, Pcent_L, Pcent_R):

    #Set up intersection
    cmds.polyCBoolOp(Ant1, Post1, op=3, n=inter)

    #Set up parenting relationships
    cmds.parent(Ant1,Ant)
    cmds.parent(Zyg_ant,Ant)
    
    #centrum objects
    cmds.parent(cent_D,Ant)
    cmds.parent(cent_V,Ant)
    cmds.parent(cent_L,Ant)
    cmds.parent(cent_R,Ant)
    
    #axis
    cmds.parent(base_axis,Ant)
    cmds.parent(Post1,Post)
    cmds.parent(Zyg_post,Post)
    
    #centrum objects
    cmds.parent(Pcent_D,Post)
    cmds.parent(Pcent_V,Post)
    cmds.parent(Pcent_L,Post)
    cmds.parent(Pcent_R,Post)

    #Create joint center locator
    cmds.spaceLocator(n = new_axis)
    cmds.matchTransform(new_axis,base_axis)
    cmds.parent(new_axis,base_axis)
    cmds.parent(Post,new_axis)


##########################################################################
#Automated bending code

#Function to perform the bending analysis in each directions
# Vertebrae are rotated around COR until either the overlap exceeds threshold 'start_' or the zygapophyses disarticulate in ventroflexion.

#distbone - distal vertebra
#boo - intersection object
#new_axis - locator axis object from above
#start_x - overlap threshold for x
#start_y - overlap threshold for y
#start_z - overlap threshold for z
#range_x - increments of movement for x (min is 1 degree)
#range_y - increment of movement for y
#range_z - increment of movement for z
#zyg_ant - locator on zygapophysis of anterior vertebra
#zyg_post - locator on zygapophysis of posterior vertebra
#cent_X - locator on anterior vertebra
#Pcent_X - locator on posterior vertebra
#v_area - total size of vertebra
#disart_fact - how much separation of joint surfaces allowed (1.5 is 50% from starting point)
#joint_tol - tolerance for zyg disarticulation, total=0, 0.5=half way
#ZygTest - True/False - will testing based on zygapophyses be carried out
#CentTest - True/False - will testing based on centrum strain be carried out
#Transl_fact - amount of translation allowed as % vertebral area e.g., 0.005

#wrapper for bending analysis to make it neater

#joint - joint no
#angle - increment for bending
#it - intersection threshold
#strain - cent/zyg strain in both directions
#v_area - vertebra area for scaling

def digital_bending(joint,v_area, angle=0.5, it=0.005, cent_strain=0.5, zyg_strain=0.5,  
ZygTest=True, CentTest=True, Transl=False, Transl_fact=0.005):

    #define objects
    ant=str(joint)
    post = str(joint+1)
    distbone = 'V'+post+':Mesh'
    boo = 'Inter'+ant
    new_axis = 'COR_locator'+ant
    start_x = it
    start_y = it  
    start_z = it
    range_x = angle
    range_y = angle
    range_z = angle
    Zyg_ant = 'Zyg_ant'+ant
    Zyg_post = 'Zyg_post'+ant
    cent_D = 'cent_D' + ant
    cent_V = 'cent_V' + ant
    cent_L = 'cent_L' + ant
    cent_R = 'cent_R' + ant
    Pcent_D='Pcent_D'+ant
    Pcent_V='Pcent_V'+ant
    Pcent_L='Pcent_L'+ant
    Pcent_R='Pcent_R'+ant
    cent_strain = cent_strain
    zyg_strain = zyg_strain

    angles = bending_analysis(distbone, boo, new_axis, start_x, start_y, start_z, range_x, range_y, range_z,
    Zyg_ant, Zyg_post, v_area, cent_D, cent_V, cent_L, cent_R, cent_strain, zyg_strain, 
    Pcent_D, Pcent_V, Pcent_L, Pcent_R, ZygTest, CentTest, Transl, Transl_fact)

    return angles

#####################################################################################################################################

##################################################################################################################################

def bending_analysis(distbone, boo, new_axis, start_x, start_y, start_z, range_x, range_y, range_z,
 Zyg_ant, Zyg_post, v_area, cent_D, cent_V, cent_L, cent_R, cent_strain, zyg_strain, 
 Pcent_D, Pcent_V, Pcent_L, Pcent_R,  ZygTest, CentTest, Transl, Transl_fact):
     
###############################
#Starting objects
     
    #original position of joint
    pos0=cmds.xform(new_axis, t=True, q=True, ws=False)#local space
    pos_ws=cmds.xform(new_axis, t=True, q=True, ws=True)#position in world space

    #Translation
    #Centrum height, aprox half centrum length
    #cent_H = wDist(cent_D, cent_V, new_axis)   
    #Or sqrt(v-area) aprox 3 times CL
    trnsl=np.sqrt(v_area)*Transl_fact #sqaure root intersection threshold, equivalent to ~3% CL  

    #Measure intersection and set threshold
    start_a = cmds.polyEvaluate(boo,a=True)
    start_x = (start_a)+(v_area*start_x)#allow overlap as fraction of vert area
    start_y = (start_a)+(v_area*start_y)
    start_z = (start_a)+(v_area*start_z)
 
    #Distance between centrum markers
    comp_fact = 1-cent_strain
    ten_fact = 1+cent_strain
    cent_dist_D = wDist(cent_D, Pcent_D, new_axis)
    cent_dist_D_thresh_comp= cent_dist_D * comp_fact
    cent_dist_D_thresh_ten= cent_dist_D * ten_fact
    cent_dist_V = wDist(cent_V, Pcent_V, new_axis)
    cent_dist_V_thresh_comp= cent_dist_V * comp_fact
    cent_dist_V_thresh_ten= cent_dist_V * ten_fact
    cent_dist_L = wDist(cent_L, Pcent_L, new_axis)
    cent_dist_L_thresh_comp= cent_dist_L * comp_fact
    cent_dist_L_thresh_ten= cent_dist_L * ten_fact
    cent_dist_R = wDist(cent_R, Pcent_R, new_axis)
    cent_dist_R_thresh_comp = cent_dist_R * comp_fact
    cent_dist_R_thresh_ten = cent_dist_R * ten_fact
    #torsion
    cent_dist_T = np.mean([wDist(cent_D, Pcent_D, new_axis),wDist(cent_V, Pcent_V, new_axis),wDist(cent_L, Pcent_L, new_axis),wDist(cent_R, Pcent_R, new_axis)])
    cent_dist_T_thresh_tor = cent_dist_T * ten_fact
    
    #Disarticulation test threshold
    #Get original zyg marker spacing
    l1 = world_pos(Zyg_ant) * world_matrix(new_axis).inverse()
    l2 = world_pos(Zyg_post) * world_matrix(new_axis).inverse()
    ltest_orig = l1[0]-l2[0] #l1.distanceTo(l2)
    joint_dis = 1-zyg_strain #test for 50% disart
    joint_over = 1+zyg_strain
    ltest_thresh = ltest_orig * joint_dis
    ltest_strain = ltest_orig * joint_over #zygapophyseal strain
    

    
#########################################################################

    #Make angles results objects        
    angles = []      
    angles.append(round(start_a,3))

    #conduct bending experiment
    #######################################################################
    #X=torsion- intersection and joint strain
    #Left
    for i in np.linspace(0,45,num=(45/range_x)+1):
        cmds.rotate(i,0,0,new_axis)
        a = cmds.polyEvaluate(boo,a=True)
        if a > start_x:
            angles.append(i)
            angles.append('intersect')
            break
        
        #joint strain
        if CentTest:    
            cent_dist = np.mean([wDist(cent_D, Pcent_D, new_axis),wDist(cent_V, Pcent_V, new_axis),
            wDist(cent_L, Pcent_L, new_axis),wDist(cent_R, Pcent_R, new_axis)])
            if cent_dist > cent_dist_T_thresh_tor:
                angles.append(i)
                angles.append('cent_ten')
                break
            
        #Maximum    
        if i == 45:
            angles.append(45)
            angles.append('none')
            
    #Right        
    for i in np.linspace(0,-45,num=(45/range_x)+1):
        cmds.rotate(i,0,0,new_axis)
        a = cmds.polyEvaluate(boo,a=True)
        if a > start_x:
            angles.append(i)
            angles.append('intersect')
            break
            
        #joint strain
        if CentTest:    
            cent_dist = np.mean([wDist(cent_D, Pcent_D, new_axis),wDist(cent_V, Pcent_V, new_axis),
            wDist(cent_L, Pcent_L, new_axis),wDist(cent_R, Pcent_R, new_axis)])
            if cent_dist > cent_dist_T_thresh_tor:
                angles.append(i)
                angles.append('cent_ten')
                break  
        
        #Maximum          
        if i == -45:
            angles.append(-45)
            angles.append('none')
            
    #######################################################################################
    #Y = lateroflexion - intersection only
    #Left
    for i in np.linspace(0,45,num=(45/range_y)+1):
        cmds.rotate(0,i,0,new_axis)
        a = cmds.polyEvaluate(boo,a=True)
        if a > start_y:
            if Transl:
                cmds.xform(new_axis,translation=[pos0[0],pos0[1], pos0[2]-trnsl],a=True, objectSpace=True)
                a = cmds.polyEvaluate(boo,a=True)
                if a > start_y:
                    angles.append(i)
                    angles.append('intersect')
                    break    
            else :
                if a > start_y:
                    angles.append(i)
                    angles.append('intersect')
                    break      
            
        #Test for centrum distance
        if CentTest:
        #compression
            cent_dist = wDist(cent_R, Pcent_R, new_axis)
            if cent_dist < cent_dist_R_thresh_comp:
                angles.append(i)
                angles.append('cent_comp')
                break
            #tension
            cent_dist = wDist(cent_L, Pcent_L, new_axis)
            if cent_dist > cent_dist_L_thresh_ten:
                angles.append(i)
                angles.append('cent_ten')
                break  
              
        #Maximum
        if i == 45:
            angles.append(45)
            angles.append('none')
            
    cmds.xform(new_axis,translation=pos_ws,a=True, ws=True)#return to original position        
    
    #right
    for i in np.linspace(0,-45,num=(45/range_y)+1):
        cmds.rotate(0,i,0,new_axis)
        a = cmds.polyEvaluate(boo,a=True)
        if a > start_y:
            if Transl:
                cmds.xform(new_axis,translation=[pos0[0],pos0[1], pos0[2]+trnsl],a=True, objectSpace=True)
                a = cmds.polyEvaluate(boo,a=True)
                if a > start_y:
                    angles.append(i)
                    angles.append('intersect')
                    break    
            else :
                if a > start_y:
                    angles.append(i)
                    angles.append('intersect')
                    break
            
        #Test for centrum distance
        if CentTest:
        #compression
            cent_dist = wDist(cent_L, Pcent_L, new_axis)
            if cent_dist < cent_dist_L_thresh_comp:
                angles.append(i)
                angles.append('cent_comp')
                break  
            #tension
            cent_dist = wDist(cent_R, Pcent_R, new_axis)
            if cent_dist > cent_dist_R_thresh_ten:
                angles.append(i)
                angles.append('cent_ten')
                break    
        
        #Maximum
        if i == -45:
            angles.append(-45)
            angles.append('none')
            
    cmds.xform(new_axis,translation=pos_ws,a=True, ws=True)#return to original position        
        
    ################################################################################
    #Z=Sagittal bending
    #dorsoflexion
    for i in np.linspace(0,45,num=(45/range_z)+1):
        cmds.rotate(0,0,i,new_axis)
        a = cmds.polyEvaluate(boo,a=True)
        if a > start_z:
            if Transl:
                cmds.xform(new_axis,translation=[pos0[0],pos0[1]-trnsl, pos0[2]],a=True, objectSpace=True)
                a = cmds.polyEvaluate(boo,a=True)
                if a > start_z:
                    angles.append(i)
                    angles.append('intersect')
                    break    
            else :
                if a > start_z:
                    angles.append(i)
                    angles.append('intersect')
                    break
            
        #Test for centrum distance
        if CentTest:
        #compression
            cent_dist = wDist(cent_D, Pcent_D, new_axis)
            if cent_dist < cent_dist_D_thresh_comp:
                angles.append(i)
                angles.append('cent_comp')
                break 
            #tension
            cent_dist = wDist(cent_V, Pcent_V, new_axis)
            if cent_dist > cent_dist_V_thresh_ten:
                angles.append(i)
                angles.append('cent_ten')
                break
                
        #Zyg strain
        if ZygTest:
            l1 = world_pos(Zyg_ant) * world_matrix(new_axis).inverse()
            l2 = world_pos(Zyg_post) * world_matrix(new_axis).inverse()
            ltest = l1[0]-l2[0]
            if ltest > ltest_strain:
                angles.append(i)
                angles.append('jstrain')  
                break
            
        #Maximum     
        if i == 45:
            angles.append(45)
            angles.append('none')    
    
    cmds.xform(new_axis,translation=pos_ws,a=True, ws=True)#return to original position        
    
    #ventroflexion
    #intersection, strain
    #test for intersection
    for i in np.linspace(0,-45,num=(45/range_z)+1):
        cmds.rotate(0,0,i,new_axis)
        a = cmds.polyEvaluate(boo,a=True)
        if a > start_z:
            if Transl:
                cmds.xform(new_axis,translation=[pos0[0],pos0[1]+trnsl, pos0[2]],a=True, objectSpace=True)
                a = cmds.polyEvaluate(boo,a=True)
                if a > start_z:
                    angles.append(i)
                    angles.append('intersect')
                    break    
            else :
                if a > start_z:
                    angles.append(i)
                    angles.append('intersect')
                    break
            
        #test for disarticulation
        if ZygTest:
            l1 = world_pos(Zyg_ant) * world_matrix(new_axis).inverse()
            l2 = world_pos(Zyg_post) * world_matrix(new_axis).inverse()
            ltest = l1[0]-l2[0] #l1.distanceTo(l2)
            if ltest < ltest_thresh:
                angles.append(i)
                angles.append('disart')   
                break

        #Test for centrum distance
        if CentTest:
        #compression
            cent_dist = wDist(cent_V, Pcent_V, new_axis)
            if cent_dist < cent_dist_V_thresh_comp:
                angles.append(i)
                angles.append('cent_comp')
                break  
            #tension 
            cent_dist = wDist(cent_D, Pcent_D, new_axis)
            if cent_dist > cent_dist_D_thresh_ten:
                angles.append(i)
                angles.append('cent_ten')
                break 
            
        #Maximum      
        if i == -45:
            angles.append(-45)
            angles.append('none')
    
   
    cmds.xform(new_axis,translation=pos_ws,a=True, ws=True)#return to original position        
             
    print(angles)
    
    #Set back to neutral position
    cmds.rotate(0,0,0,new_axis)
    
    return angles

#############################################################################

#custom functions to get object positions
def world_matrix(obj):
    """'
    convenience method to get the world matrix of <obj> as a matrix object
    """
    return MMatrix( cmds.xform(obj, q=True, matrix=True, ws=True))


def world_pos(obj):
    """'
    convenience method to get the world position of <obj> as an MPoint
    """
    return MPoint( cmds.xform(obj, q=True, t=True, ws=True, s=False)) 

def wDist(obj1,obj2, new_axis):
    #calculate distance between locators
    D1 = world_pos(obj1) * world_matrix(new_axis).inverse()
    dis = D1.distanceTo(world_pos(obj2) * world_matrix(new_axis).inverse())
    return dis

###############################################################
#Set joint spacing

def adjust_spacing(new_axis, Post, space):
    
    mmt=float(space)/2

    #get original position
    pos0 = cmds.getAttr(new_axis+".translate")[0]
    
    #move posterior vert
    #cmds.move(new_axis,moveX=mmt,a=True, objectSpace=True)
    cmds.xform(new_axis,translation=(mmt,0,0),r=True, objectSpace=True)
    
    #get new position and calc rel movement
    pos1 = cmds.getAttr(new_axis+".translate")[0]
    m = pos1[0]-pos0[0]
    
    #adjust COR to stay in the middle
    cmds.xform(Post,translation=(m,0,0), r=True)#, localSpace=True)

##############################################################
#Robs function for finding the center of an object
def getVtxCentroid( shapeNode ) :
    vtxWorldPosition = []    # will contain positions un space of all object vertex
    vtxIndexList = cmds.getAttr( shapeNode+".vrts", multiIndices=True )
    for i in vtxIndexList :
    	curPointPosition = cmds.xform( str(shapeNode)+".pnts["+str(i)+"]", query=True, translation=True, worldSpace=True )
    	vtxWorldPosition.append( curPointPosition )
    vtxX = [x[0] for x in vtxWorldPosition]
    vtxY = [x[1] for x in vtxWorldPosition]
    vtxZ = [x[2] for x in vtxWorldPosition]
    Xcentre = sum(vtxX)/len(vtxX) # You'll want to change these lines to get geometric not arithmetic centre
    Ycentre = sum(vtxY)/len(vtxY)
    Zcentre = sum(vtxZ)/len(vtxZ)
    centroid = [Xcentre,Ycentre,Zcentre]
    return(centroid)

    
######################################################################
#from open maya tutorial from https://pastebin.com/6xtHkt7t
#put a locator on nearest point on mesh or calculate distance

def dagObjFromName(name):
        sel = om2.MSelectionList()
        sel.add(name)
       
        dag =sel.getDagPath(0)
       
        mobj = sel.getDependNode(0)
       
        return mobj, dag
 
 
def locatorToPoint(locator):
        pos = cmds.xform(locator, q=1, rp=1, ws=1)
        return om2.MPoint(pos[0], pos[1], pos[2],1.0)
        

def createLocatorAtClosestPoint(mySphere, myLocator):
        obj,dag = dagObjFromName(mySphere)     
        meshFn = om2.MFnMesh(dag)
 
        point = locatorToPoint(myLocator)
 
        resultPoint = meshFn.getClosestPoint(point,om2.MSpace.kWorld)[0]
       
        result = cmds.spaceLocator()[0]
       
        cmds.setAttr(result+'.t', resultPoint[0], resultPoint[1], resultPoint[2])
        
def distToClosestPoint(mySphere, myLocator):
        obj,dag = dagObjFromName(mySphere)     
        meshFn = om2.MFnMesh(dag)
 
        point = locatorToPoint(myLocator)
        locpos = cmds.xform(myLocator, q=1, rp=1, ws=1)
 
        resultPoint = meshFn.getClosestPoint(point,om2.MSpace.kWorld)[0]
        resultPos = [resultPoint[0], resultPoint[1], resultPoint[2]]
        
        space = cmds.distanceDimension(sp=locpos, ep=resultPos)
        dis = cmds.getAttr("distanceDimension1"+".distance")
        cmds.delete("distanceDimension1")
        if cmds.objExists("locator1"):
            cmds.delete("locator1")
        
        return dis
