import numpy as np
from FEA_MAIN import *

"""This is the file you want to run, 2024 starts at line 506, fea_main does maths"""

"""function to set the force equivalent vector \n
            w = pghd
            coords in element system
           'UDL' = Uniformly Distributed Load. 𝑤̅ = distributed shear load intensity\n
           'LVL' = Linearly Varying Distributed Load. where the line equation is y= W(x/L) and 𝑤̅ = distributed shear load intensity\n
           'PL' = Point load. where a is the distance of the point load from x_e = 0 \n
           'MSPL' = as above but a = 0.5L\n
           'DAL' = Distributed Axial Load.  W = the magnitude of the distributed axial force\n
           'CAL' = Concentrated Axial Load. where a is the distance of the point load from x_e = 0 \n 
"""

def printA_e(structure: Structure):
    """
    Print element assembly matrix for each element\n
    Must be solved already
    """
    print("Element assembly matrix")
    for element in structure.elements:
        print(f"{element.name}: \n {element.AssemblyMatrix} \n")
    print("Done\n")
    
def printf_eq(structure: Structure):
    """Prints each force equivalent terms for the elements in a structure (f_eq)\n

    Args:
        structure (Structure): Structure must be solved
    """
    print("Elements force equivalent vector")
    for element in structure.elements:
        print(f"{element.name}: \n {element.forceEquivalent} \n")
    print("Done\n")
    
def printGlobalForces(structure: Structure):
    """Prints rxn forces of each element in global.

    Args:
        structure (Structure): _description_
    """
    print("Global reaction forces and moments (F)")
    for element in structure.elements:
        print(f"{element.name}: \n {element.GlobalForces} \n")
    print("Done\n")
    
def printforces(structure: Structure):
    """Prints rxn forces of each element in local (f).

    Args:
        structure (Structure): _description_
    """
    print("local element forces and moments (f)")
    for element in structure.elements:
        print(f"{element.name}: \n {element.localforces} \n")
    print("Done\n")
    
def printstress(structure: Structure):
    """print elements stress and strain\n
    Only considers axial stress\n
    A negative value is compressive 

    Args:
        structure (Structure): _description_
    """
    print("Elements stress and strain only axial")
    for element in structure.elements:
        print(f"{element.name}: \n stress: {element.stress} | strain: {element.strain}\n")
    print("Done\n")
    
def printd_e(structure: Structure):
    """print element displacement vector in element co-ordinates (d_e)

    Args:
        structure (Structure): _description_
    """
    print("Element displacement vector in element co-ordinates (d_e) mm")
    for element in structure.elements:
        print(f"{element.name}: \n {element.deflectionlocal * 1e3}\n")
    print("Done\n")
    
def printK_G(structure: Structure):
    """Prints Global Stiffness matrix (K_G) and each elements to the stiffness matrix (K_G_e)

    Args:
        structure (Structure): _description_
    """
    print("Global Stiffness matrix (K_G) and each elements to the stiffness matrix (K_G_e) x10^6")
    print(f"Global Stiffness matrix (K_G): \n {structure.GlobalStiffnessMatrix * 1e-6}\n")
    for element in structure.elements:
        print(f"{element.name} KGe: \n {element.GlobalStiffnessMatrix * 1e-6}\n")
    print("Done\n")

def printInterpolate(frame: FrameElement, point: float, coords= 'local'):
    """prints interpolated values at the given point

    Args:
        frame (FrameElement): _description_
        point (float): _description_
        coords (str, optional): _description_. Defaults to 'local'.
    """
    print(f"{frame.name} displacement at L={point}: {frame.interpolation(0, point, coords)}")

def printReactions(structure: Structure):

    print("Nodal reaction loads (R)")
    for node in structure.nodes:
        print(f"{node.name}:\n {node.Reactions}\n")
    print("Done\n")

def printstandardtestQ(structure: Structure):
    """prints most answers for questions most commonly asked

    Args:
        structure (Structure): _description_
    """
    #print element assembly matrix for each element 
    printA_e(structure)
    
    #print force eq for each element
    printf_eq(structure)
    
    #print Q for structure
    print(f"Global Forcing Vector (Q):\n {structure.externalLoads}\n")
    
    #print deflections of the system (q)
    print(f"Deflections of the system (q) in mm:\n {structure.GlobalDisplacement * 1e3}\n")

    #local element forces and moments
    printReactions(structure)

def distloadangle(alpha):
    """function that returns an x , y multiplier to help with non nodal loads at an angle\n
    Args:
        alpha (float): Angle in degrees. either the angle between Y_e and the load direction or alpha if load is all in Y_G.

    Returns:
        (x, y): A tuple than contants an x and y mulitpler 
    """
    return (np.sin(np.radians(alpha)), np.cos(np.radians(alpha)))

def length(x1, y1, x2, y2):
    
    return np.sqrt((x2-x1)**2+(y2-y1)**2)

def alpha(structure: Structure):

    for element in structure.elements:
        print(f"{element.name} alpha: \n {element.alpha * 180 / np.pi}\n")
    print("Done\n")

def printlegnth(structure: Structure):
    for element in structure.elements:
        print(f"{element.name} alpha: \n {element.L}\n")
    print("Done\n")
    
    
    

def test2020_1():
    E = 10e9
    I = 1e-6
    A = 4e-4
    
    #set up nodes with nodal forces
    node1 = Node("node1", XPos=0, YPos=0, DoFX=False, DoFY=False, ExternalLoadX= 0, ExternalLoadY=0)
    node2 = Node("node2", XPos= 0,YPos= 0.5, DoFX= True, DoFY= True, ExternalLoadX= 0, ExternalLoadY= 0, DoFRotZ=True)
    node3 = Node("node3", XPos=None, YPos= None, DoFX= True, DoFY= True, ExternalLoadX= 0, ExternalLoadY= 0, DoFRotZ=True)
    node4 = Node("node4", XPos=None, YPos=None, DoFX= False, DoFY= False,ExternalLoadX= 0,ExternalLoadY= 0)
    
    #set up of frames and dist loads
    frame1 =FrameElement("frame1", 90, E=E, A=A, I=I, L=0.5, node1=node1, node2=node2)
    frame2 =FrameElement("frame2", 45, E=E, A=A, I=I, L=2, node1=node2, node2=node3,
                         distributedLoadtype=['CAL', 'PL'], distributedLoadpoint=[1.5, 1.5], distributedLoadforce=[694.59, -3939.23])
    frame3 =FrameElement("frame3", -75, E=E, A=A, I=I, L=2, node1=node3, node2=node4)
    
    #set up of structure and added frames and solved
    structure = Structure("structure1")
    structure.add_element(frame1)
    structure.add_element(frame2)
    structure.add_element(frame3)
    structure.solve()
    
    #print element assembly matrix for each element 
    printA_e(structure)
    
    #print force eq for each element
    printf_eq(structure)
    
    #print Q for structure
    print(f"Global Forcing Vector (Q):\n {structure.externalLoads}\n")
    
    #print deflections of the system (q)
    print(f"Deflections of the system (q) in mm:\n {structure.GlobalDisplacement * 1e3}\n")
    
    #element forces and moments
    printGlobalForces(structure)
    
    #local element forces and moments
    printforces(structure)
    
    #print elements stress and strain
    printstress(structure)
    
    #element displacement vector by interpolating
    printd_e(structure)
    print(f"Displacement interpolation\n{frame1.name} @ {frame1.L/2}: {frame1.interpolation(0, point=frame1.L/2, coords='global')}")
    
    #matplotlib structure at 25x deformations and element interpolated for 10 points
    structure.plot(25, 10)
    
def test2020_2():
    E = 30e9
    I = 6e-6
    A = 2e-4
    
    #set up nodes with nodal forces
    node1 = Node("node1", XPos=0, YPos=0, DoFX=False, DoFY=False, ExternalLoadX= 0, ExternalLoadY=0)
    node2 = Node("node2", XPos= None,YPos= None, DoFX= True, DoFY= True, ExternalLoadX= 0, ExternalLoadY= 0, DoFRotZ=True)
    node3 = Node("node3", XPos=None, YPos= None, DoFX= False, DoFY= True, ExternalLoadX= 0, ExternalLoadY= 0)
    
    #set up of frames and dist loads
    frame1 =FrameElement("frame1", 80, E=E, A=A, I=I, L=3, node1=node1, node2=node2)
    frame2 =FrameElement("frame2", 0, E=E, A=A, I=I, L=2, node1=node2, node2=node3,
                         distributedLoadtype=['UDL'], distributedLoadforce=[-10000])
    
    #set up of structure and added frames and solved
    structure = Structure("structure1")
    structure.add_element(frame1)
    structure.add_element(frame2)
    structure.solve()
    
    #print element assembly matrix for each element 
    printA_e(structure)
    
    #print KG and KGe
    printK_G(structure)
    
    #print force eq for each element
    printf_eq(structure)
    
    #print Q for structure
    print(f"Global Forcing Vector (Q):\n {structure.externalLoads}\n")

    #print deflections of the system (q)
    print(f"Deflections of the system (q) in mm:\n {structure.GlobalDisplacement * 1e3}\n")

    #local element forces and moments
    printforces(structure)

    #matplotlib structure at 25x deformations and element interpolated for 10 points
    structure.plot(10, 10)

def test2021_1():
    E = 200e9
    I = 1e-6
    A = 1e-4
    
    #set up nodes with nodal forces
    node1 = Node("node1", XPos=0, YPos=0, DoFX=False, DoFY=False, ExternalLoadX= 0, ExternalLoadY=0)
    node2 = Node("node2", XPos= None,YPos= None, DoFX= True, DoFY= False, ExternalLoadX= 0, ExternalLoadY= 0, DoFRotZ=True)
    node3 = Node("node3", XPos=None, YPos= None, DoFX= True, DoFY= True, ExternalLoadX= 25000, ExternalLoadY= -10000, DoFRotZ=True)
    
    #set up of frames and dist loads
    frame1 =FrameElement("frame1", 0, E=E, A=A, I=I, L=3, node1=node1, node2=node2,
                         distributedLoadtype=['UDL'], distributedLoadforce=[-40000])
    frame2 =FrameElement("frame2", 0, E=E, A=A, I=I, L=2, node1=node2, node2=node3)
    
    #set up of structure and added frames and solved
    structure = Structure("structure1")
    structure.add_element(frame1)
    structure.add_element(frame2)
    structure.solve()
    
    #print element assembly matrix for each element 
    printA_e(structure)
    
    #print KG and KGe
    printK_G(structure)
    
    #print force eq for each element
    printf_eq(structure)
    
    #print Q for structure
    print(f"Global Forcing Vector (Q):\n {structure.externalLoads}\n")

    #print deflections of the system (q)
    print(f"Deflections of the system (q) in mm:\n {structure.GlobalDisplacement * 1e3}\n")

    #local element forces and moments
    printforces(structure)

    #matplotlib structure at 25x deformations and element interpolated for 10 points
    structure.plot(10, 10)
    
def test2021_2():
    E = 200e9
    I = 7.5e-6
    A = 9e-4
    
    #set up nodes with nodal forces
    node2 = Node("node2", XPos=0, YPos= 0, DoFX= False, DoFY= False, ExternalLoadX= 0, ExternalLoadY=0, DoFRotZ= False)
    node1 = Node("node1", XPos= 0,YPos= 2.5, DoFX= True, DoFY= True, ExternalLoadX= 0, ExternalLoadY= 0, DoFRotZ=True)
    node3 = Node("node3", XPos=None, YPos= None, DoFX= True, DoFY= True, ExternalLoadX= 0, ExternalLoadY= 0, DoFRotZ=True)
    node4 = Node("node4", XPos=None, YPos=None, DoFX=False, DoFY=False, ExternalLoadX= 0, ExternalLoadY=0, DoFRotZ= False)

    
    #set up of frames and dist loads
    frame1 =FrameElement("frame1", -90, E=E, A=A, I=I, L=2.5, node1=node1, node2=node2,
                         distributedLoadtype=['LVL'], distributedLoadforce=[73573])
    frame2 =FrameElement("frame2", -36.87, E=E, A=A, I=I, L=np.sqrt(2**2+1.5**2), node1=node1, node2=node3)
    frame3 =FrameElement("frame3", -90, E=E, A=A, I=I, L=1, node1=node3, node2=node4)

    #set up of structure and added frames and solved
    structure = Structure("structure1")
    structure.add_element(frame1)
    structure.add_element(frame2)
    structure.add_element(frame3)
    structure.solve()
    
    #print element assembly matrix for each element 
    printA_e(structure)
    
    #print KG and KGe
    printK_G(structure)
    
    #print force eq for each element
    printf_eq(structure)
    
    #print Q for structure
    print(f"Global Forcing Vector (Q):\n {structure.externalLoads}\n")

    #print deflections of the system (q)
    print(f"Deflections of the system (q) in mm:\n {structure.GlobalDisplacement * 1e3}\n")

    #local element forces and moments
    printGlobalForces(structure)

    #stress and strains
    printstress(structure)

    #interpolate at a point
    printInterpolate(frame2, (frame2.L)/2, 'global')

    #matplotlib structure at 25x deformations and element interpolated for 10 points
    structure.plot(10, 10)

def test2022_1():
    E = 200e9
    I = 5e-6
    A = 2e-4
    L1 = 4
    L2 = np.sqrt(3**2 + 4**2)
    L3 = 3
    
    #set up nodes with nodal forces
    node1 = Node("node1", XPos=0, YPos=-4, DoFX=True, DoFY=False, ExternalLoadX= 0, ExternalLoadY=0, DoFRotZ=True)
    node2 = Node("node2", XPos= None,YPos= None, DoFX= False, DoFY= False, ExternalLoadX= 0, ExternalLoadY= 0, DoFRotZ=False)
    node3 = Node("node3", XPos=None, YPos= None, DoFX= True, DoFY= True, ExternalLoadX= 20000, ExternalLoadY= 40000, DoFRotZ=True)
    
    #set up of frames and dist loads
    frame1 =FrameElement("frame1", 180, E=E, A=A, I=I, L=L1, node1=node1, node2=node2,
                         distributedLoadtype=['UDL'], distributedLoadforce=[25000])
    frame2 =FrameElement("frame2", 36.86, E=E, A=A, I=I, L=L2, node1=node2, node2=node3)
    frame3 =FrameElement("frame3", -90, E=E, A=A, I=I, L=L3, node1=node3, node2=node1)
    
    #set up of structure and added frames and solved
    structure = Structure("structure1")
    structure.add_element(frame1)
    structure.add_element(frame2)
    structure.add_element(frame3)
    structure.solve()
    
    #print element assembly matrix for each element 
    printA_e(structure)
    
    #print KG and KGe
    printK_G(structure)
    
    #print force eq for each element
    printf_eq(structure)
    
    #print Q for structure
    print(f"Global Forcing Vector (Q):\n {structure.externalLoads}\n")

    #print deflections of the system (q)
    print(f"Deflections of the system (q) in mm:\n {structure.GlobalDisplacement * 1e3}\n")

    #local element forces and moments
    printReactions(structure)

    #print element forces in global coords
    printGlobalForces(structure)

    printInterpolate(frame1, frame1.L/2, 'global')

    #matplotlib structure at 25x deformations and element interpolated for 10 points
    structure.plot(25, 10)

def test2022_2():
    E = 200e9
    I = 4.5e-6
    A = 3e-4
    L1 = 1.5
    L2 = 2
    L3 = 2.5
    
    #set up nodes with nodal forces
    node1 = Node("node1", XPos=0, YPos=0, DoFX=False, DoFY=False, ExternalLoadX= 0, ExternalLoadY=0, DoFRotZ=False)
    node2 = Node("node2", XPos= None,YPos= None, DoFX= True, DoFY= True, ExternalLoadX= 0, ExternalLoadY= 0, DoFRotZ=True)
    node3 = Node("node3", XPos=None, YPos= None, DoFX= False, DoFY= False, ExternalLoadX= 0, ExternalLoadY= 0, DoFRotZ=False)
    node4 = Node("node4", XPos= None,YPos= None, DoFX= True, DoFY= True, ExternalLoadX= 0, ExternalLoadY= 0, DoFRotZ=True)

    disloadX, disloadY = distloadangle(15)

    #set up of frames and dist loads
    frame1 =FrameElement("frame1", 0, E=E, A=A, I=I, L=L1, node1=node1, node2=node2,
                         distributedLoadtype=['UDL'], distributedLoadforce=[-5000])
    frame2 =FrameElement("frame2", 15, E=E, A=A, I=I, L=L2, node1=node2, node2=node4,
                         distributedLoadtype=['DAL', 'UDL'], distributedLoadforce=[-5000*disloadX, -5000*disloadY])
    frame3 =FrameElement("frame3", -100, E=E, A=A, I=I, L=L3, node1=node2, node2=node3)
    
    #set up of structure and added frames and solved
    structure = Structure("structure1")
    structure.add_element(frame1)
    structure.add_element(frame2)
    structure.add_element(frame3)
    structure.solve()

    printstandardtestQ(structure)

    printstress(structure)

    printd_e(structure)

    printInterpolate(frame3, frame3.L/2, 'global')

    #matplotlib structure at 10x deformations and element interpolated for 10 points
    structure.plot(30, 10)

def test2023_1():
    E = 200e9
    I = 6e-6
    A = 2e-4
    L1 = 3
    L2 = 5
    L3 = 5
    L4 = 6
    
    #set up nodes with nodal forces
    node1 = Node("node1", XPos=0, YPos=0, DoFX=False, DoFY=False, ExternalLoadX= 0, ExternalLoadY=0, DoFRotZ=True)
    node2 = Node("node2", XPos= 3,YPos= 4, DoFX= True, DoFY= True, ExternalLoadX= 0, ExternalLoadY= 0, DoFRotZ=True)
    node3 = Node("node3", XPos=6, YPos= 0, DoFX= False, DoFY= True, ExternalLoadX= 0, ExternalLoadY= 0, DoFRotZ=False)
    node4 = Node("node4", XPos= 6, YPos= 4, DoFX= False, DoFY= True, ExternalLoadX= 0, ExternalLoadY= 0, DoFRotZ=False)

    disloadX, disloadY = distloadangle(15)

    #set up of frames and dist loads
    frame1 =FrameElement("frame1", 0, E=E, A=A, I=I, L=L1, node1=node2, node2=node4,)
    frame2 =FrameElement("frame2", 53.1, E=E, A=A, I=I, L=L2, node1=node1, node2=node2,)
    frame3 =FrameElement("frame3", -53.1, E=E, A=A, I=I, L=L3, node1=node2, node2=node3)
    frame4 = FrameElement("frame4", 0, E=E, A=A, I=I, L=L4, node1=node1, node2=node3, 
                          distributedLoadtype=['UDL'], distributedLoadforce=[-36000])
    
    #set up of structure and added frames and solved
    structure = Structure("structure1")
    structure.add_element(frame1)
    structure.add_element(frame2)
    structure.add_element(frame3)
    structure.add_element(frame4)
    structure.solve()

    printstandardtestQ(structure)

    printforces(structure)

    printInterpolate(frame4, frame4.L)

    structure.plot(5, 10)\
    
def test2023_2():
    E = 200e9
    I = 1e-4
    A = 2e-4
    L1 = 4
    L2 = 3
    L3 = 4
    
    #set up nodes with nodal forces
    node1 = Node("node1", XPos=0, YPos=0, DoFX=False, DoFY=False, ExternalLoadX= 0, ExternalLoadY=0, DoFRotZ=False)
    node2 = Node("node2", XPos= None,YPos= None, DoFX= True, DoFY= True, ExternalLoadX= 0, ExternalLoadY= 0, DoFRotZ=True)
    node3 = Node("node3", XPos=None, YPos= None, DoFX= False, DoFY= False, ExternalLoadX= 0, ExternalLoadY= 0, DoFRotZ=True)
    node4 = Node("node4", XPos= None,YPos= None, DoFX= False, DoFY= False, ExternalLoadX= 0, ExternalLoadY= 0, DoFRotZ=False)

    #set up of frames and dist loads
    frame1 =FrameElement("frame1", -30, E=E, A=A, I=I, L=L1, node1=node1, node2=node2)
    frame2 =FrameElement("frame2", -90, E=E, A=A, I=I, L=L2, node1=node2, node2=node3,
                         distributedLoadtype=['LVL'], distributedLoadforce=[-117720])
    frame3 =FrameElement("frame3", -180, E=E, A=A, I=I, L=L3, node1=node3, node2=node4,
                         distributedLoadtype=['UDL'], distributedLoadforce=[-117720])
    
    #set up of structure and added frames and solved
    structure = Structure("structure1")
    structure.add_element(frame1)
    structure.add_element(frame2)
    structure.add_element(frame3)
    structure.solve()

    printstandardtestQ(structure)

    printstress(structure)

    printd_e(structure)

    printInterpolate(frame1, 3, 'global')

    #matplotlib structure at 10x deformations and element interpolated for 10 points
    structure.plot(100, 10)

def test2024_1():
    E = 200e9
    I = 5e-6
    A = 3e-4
    L1 = 5
    L2 = 4
    L3 = 3
    
    #set up nodes with nodal forces
    node1 = Node("node1", XPos=0, YPos=2.4, DoFX=False, DoFY=False, ExternalLoadX= 0, ExternalLoadY=0, DoFRotZ=False)
    node2 = Node("node2", XPos= None,YPos= None, DoFX= True, DoFY= True, ExternalLoadX= 15000, ExternalLoadY= -20000, DoFRotZ=True)
    node3 = Node("node3", XPos=None, YPos= None, DoFX= True, DoFY= False, ExternalLoadX= 0, ExternalLoadY= 0, DoFRotZ=True)

    #set up of frames and dist loads
    frame1 =FrameElement("frame1", 0, E=E, A=A, I=I, L=L1, node1=node1, node2=node2,
                          distributedLoadtype=['UDL'], distributedLoadforce=[-30000])
    frame2 =FrameElement("frame2", -36.9, E=E, A=A, I=I, L=L2, node1=node1, node2=node3)
    frame3 =FrameElement("frame3", -(180-53.1), E=E, A=A, I=I, L=L3, node1=node2, node2=node3)
    
    #set up of structure and added frames and solved
    structure = Structure("structure1")
    structure.add_element(frame1)
    structure.add_element(frame2)
    structure.add_element(frame3)
    structure.solve()

    printstandardtestQ(structure)

    printforces(structure)

    printstress(structure)

    printd_e(structure)

    printInterpolate(frame2, 3, 'global')
    printInterpolate(frame2, 3, 'local')

    #matplotlib structure at 20x deformations and element interpolated for 10 points
    structure.plot(20, 10)

def test2024_2():
    E = 200e9
    I = 5e-6
    A = 3e-4
    L1 = 2.5
    L2 = np.sqrt(3**2 + 0.5**2)
    
    #set up nodes with nodal forces
    node1 = Node("node1", XPos=0, YPos=0, DoFX=False, DoFY=False, ExternalLoadX= 0, ExternalLoadY=0, DoFRotZ=False)
    node2 = Node("node2", XPos= None,YPos= None, DoFX= True, DoFY= True, ExternalLoadX= 0, ExternalLoadY= 0, DoFRotZ=True)
    node3 = Node("node3", XPos=None, YPos= None, DoFX= False, DoFY= True, ExternalLoadX= 0, ExternalLoadY= 0, DoFRotZ=False)

    disloadX, disloadY = distloadangle(9.46)
    
    #set up of frames and dist loads
    frame1 =FrameElement("frame1", 90, E=E, A=A, I=I, L=L1, node1=node1, node2=node2)
    frame2 =FrameElement("frame2", -9.46, E=E, A=A, I=I, L=L2, node1=node2, node2=node3,
                         distributedLoadtype=['DAL', 'UDL'], distributedLoadforce=[1000*4*disloadX, -1000*4*disloadY])
    
    print(1000*4*disloadX)
    print(-1000*4*disloadY)

    #set up of structure and added frames and solved
    structure = Structure("structure1")
    structure.add_element(frame1)
    structure.add_element(frame2)
    structure.solve()

    printstandardtestQ(structure)

    printforces(structure)

    #matplotlib structure at 25x deformations and element interpolated for 10 points
    structure.plot(25, 10)

def Assignment1Bar():
    E = 200e9
    A = 8.29e-4
    L1 = 2
    L2 = length(0,0,1.3333,1.2)
    L3 = 1.3333
    L4 = length(1.3333, 1.2, 2, 0)
    L5 = length(1.3333, 1.2, 2.6666, 0.8)
    L6 = length(2, 0, 2.6666, 0.8)
    L7 = 2
    L8 = length(2.6666, 0.8, 4, 0)

    node1 = Node("node1", XPos=0, YPos=0, DoFX=False, DoFY=False, ExternalLoadX= 0, ExternalLoadY=0, DoFRotZ=True)
    node2 = Node("node2", XPos=0, YPos=1.2, DoFX=False, DoFY=False, ExternalLoadX= 0, ExternalLoadY=0, DoFRotZ=True)
    node3 = Node("node3", XPos= None,YPos= None, DoFX= True, DoFY= True, ExternalLoadX= 0, ExternalLoadY= 0, DoFRotZ=True)
    node4 = Node("node4", XPos= None,YPos= None, DoFX= True, DoFY= True, ExternalLoadX= 0, ExternalLoadY= 0, DoFRotZ=True)
    node5 = Node("node5", XPos= None,YPos= None, DoFX= True, DoFY= True, ExternalLoadX= 0, ExternalLoadY= 0, DoFRotZ=True)
    node6 = Node("node6", XPos= None,YPos= None, DoFX= True, DoFY= True, ExternalLoadX= 0, ExternalLoadY= -25000, DoFRotZ=True)

    frame1 =BarElement("frame1", 0,                                           E=E, A=A, L=L1, node1=node1, node2=node4)
    frame2 =BarElement("frame2", np.degrees(np.arctan(1.2/1.3333)),           E=E, A=A, L=L2, node1=node1, node2=node3)
    frame3 =BarElement("frame3", 0,                                           E=E, A=A, L=L3, node1=node2, node2=node3)
    frame4 =BarElement("frame4", -np.degrees(np.arctan(1.2/(2-1.3333))),      E=E, A=A, L=L4, node1=node3, node2=node4)
    frame5 =BarElement("frame5", -np.degrees(np.arctan((1.2-0.8)/1.3333)),    E=E, A=A, L=L5, node1=node3, node2=node5)
    frame6 =BarElement("frame6", np.degrees(np.arctan(0.8/(2.6666-2))),       E=E, A=A, L=L6, node1=node4, node2=node5)
    frame7 =BarElement("frame7", 0,                                           E=E, A=A, L=L7, node1=node4, node2=node6)
    frame8 =BarElement("frame8", -np.degrees(np.arctan(0.8/(4-2.6666))),      E=E, A=A, L=L8, node1=node5, node2=node6)

    structure = Structure("FRAME STRUCTURE")
    structure.add_element(frame1)
    structure.add_element(frame2)
    structure.add_element(frame3)
    structure.add_element(frame4)
    structure.add_element(frame5)
    structure.add_element(frame6)
    structure.add_element(frame7)
    structure.add_element(frame8)
    structure.solve()

    printK_G(structure)

    alpha(structure)

    printstandardtestQ(structure)

    structure.plot(20, 10)
    

def Assignment1Frame():
    E = 200e9
    I = 2.04e-7
    A = 8.29e-4
    L1 = 2
    L2 = length(0,0,1.3333,1.2)
    L3 = 1.3333
    L4 = length(1.3333, 1.2, 2, 0)
    L5 = length(1.3333, 1.2, 2.6666, 0.8)
    L6 = length(2, 0, 2.6666, 0.8)
    L7 = 2
    L8 = length(2.6666, 0.8, 4, 0)

    node1 = Node("node1", XPos=0, YPos=0, DoFX=False, DoFY=False, ExternalLoadX= 0, ExternalLoadY=0, DoFRotZ=False)
    node2 = Node("node2", XPos=0, YPos=1.2, DoFX=False, DoFY=False, ExternalLoadX= 0, ExternalLoadY=0, DoFRotZ=False)
    node3 = Node("node3", XPos= None,YPos= None, DoFX= True, DoFY= True, ExternalLoadX= 0, ExternalLoadY= 0, DoFRotZ=True)
    node4 = Node("node4", XPos= None,YPos= None, DoFX= True, DoFY= True, ExternalLoadX= 0, ExternalLoadY= 0, DoFRotZ=True)
    node5 = Node("node5", XPos= None,YPos= None, DoFX= True, DoFY= True, ExternalLoadX= 0, ExternalLoadY= 0, DoFRotZ=True)
    node6 = Node("node6", XPos= None,YPos= None, DoFX= True, DoFY= True, ExternalLoadX= 0, ExternalLoadY= -25000, DoFRotZ=True)

    frame1 =FrameElement("frame1", 0,                                           E=E, A=A, I=I, L=L1, node1=node1, node2=node4)
    frame2 =FrameElement("frame2", np.degrees(np.arctan(1.2/1.3333)),           E=E, A=A, I=I, L=L2, node1=node1, node2=node3)
    frame3 =FrameElement("frame3", 0,                                           E=E, A=A, I=I, L=L3, node1=node2, node2=node3)
    frame4 =FrameElement("frame4", -np.degrees(np.arctan(1.2/(2-1.3333))),      E=E, A=A, I=I, L=L4, node1=node3, node2=node4)
    frame5 =FrameElement("frame5", -np.degrees(np.arctan((1.2-0.8)/1.3333)),    E=E, A=A, I=I, L=L5, node1=node3, node2=node5)
    frame6 =FrameElement("frame6", np.degrees(np.arctan(0.8/(2.6666-2))),       E=E, A=A, I=I, L=L6, node1=node4, node2=node5)
    frame7 =FrameElement("frame7", 0,                                           E=E, A=A, I=I, L=L7, node1=node4, node2=node6)
    frame8 =FrameElement("frame8", -np.degrees(np.arctan(0.8/(4-2.6666))),      E=E, A=A, I=I, L=L8, node1=node5, node2=node6)

    structure = Structure("FRAME STRUCTURE")
    structure.add_element(frame1)
    structure.add_element(frame2)
    structure.add_element(frame3)
    structure.add_element(frame4)
    structure.add_element(frame5)
    structure.add_element(frame6)
    structure.add_element(frame7)
    structure.add_element(frame8)
    structure.solve()

    printstandardtestQ(structure)
    print(structure.elements[6].GlobalStiffnessMatrix)

    structure.plot(20, 10)
    

def main():
    Assignment1Bar()
    

if __name__ == "__main__":
    main()