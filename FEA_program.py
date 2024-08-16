import numpy as np
from FEA_MAIN import *

"""function to set the force equivalent vector \n
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
    """Prints each force equivalent terms for the elements in a structure\n

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
    print("Element displacement vector in element co-ordinates")
    for element in structure.elements:
        print(f"{element.name}: \n {element.deflectionlocal}\n")
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
    
    #support reaction forces and moments
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
    printforces(structure)

    #matplotlib structure at 25x deformations and element interpolated for 10 points
    structure.plot(10, 10)


def main():
    test2022_1()
    

if __name__ == "__main__":
    main()