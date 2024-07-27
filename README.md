Finite Element Analysis (FEA) Python Implementation
Introduction
This project provides a Python implementation of Finite Element Analysis (FEA) for structural analysis, focusing on bar elements and nodes. The implementation includes classes for nodes, bar elements, and the overall structure. This README serves as a comprehensive guide to understanding and using the provided code.

Classes and Methods
Node Class
The Node class represents the fundamental unit in an FEA model. Nodes are defined by their positions and degrees of freedom (DoF).

Constructor
python
def __init__(self, name, XPos, YPos, DoFX, DoFY, EternalLoadX, EternalLoadY, ZPos=None, DoFZ=None, EternalLoadZ=None):
  name: Identifier for the node.
  XPos, YPos, ZPos: Coordinates of the node in the X, Y, and optionally Z directions.
  DoFX, DoFY, DoFZ: Boolean values indicating if the node has degrees of freedom in the X, Y, and Z directions.
  EternalLoadX, EternalLoadY, EternalLoadZ: External loads applied to the node in the X, Y, and Z directions.

BarElement Class
The BarElement class represents a bar element made up of two nodes.

Constructor
python
def __init__(self, name, alpha, E, A, L, node1: Node, node2: Node):
  name: Identifier for the bar element.
  alpha: Angle of the bar element in degrees.
  E: Young's modulus of the material.
  A: Cross-sectional area of the bar element.
  L: Length of the bar element.
  node1, node2: Nodes that form the endpoints of the bar element.
Methods
setlocalStiffnessMatrix(): Sets the local stiffness matrix for the bar element.
setTransMatrix(): Sets the transformation matrix based on the angle alpha.
setStiffnessMatrix(): Sets the global stiffness matrix for the bar element.
setAssemblyMatrix(): Sets the assembly matrix for the bar element.

Structure Class
The Structure class represents the overall structure made up of elements.

Constructor
python
Copy code
def __init__(self, name):
  name: Identifier for the structure.
  Methods
  add_element(element): Adds an element to the structure.
  addNodes(): Adds all nodes to the structure's node list.
  setGlobalStiffnessMatrix(): Sets the global stiffness matrix for the structure.
  setexeternalLoads(): Sets the external loads for the structure.
  setGlobalDisplacement(): Solves for the global displacements of the structure.
  setElementRxnForces(): Sets the reaction forces for each element.
  setElementDisplacement(): Sets the displacement for each element.

Example Usage
The following example demonstrates how to use the classes to define a structure and perform FEA.

python
Copy code
def main():
    node1 = Node("Node1", 0, 0, False, False, 0, 0)
    node2 = Node("Node2", 1, 0, True, True, 0, -20000)
    node3 = Node("Node3", 0, 0, False, False, 0, 0)

    bar1 = BarElement("bar1", alpha=0, E=200e9, A=400e-6, L=1.1, node1=node1, node2=node2)
    bar2 = BarElement("bar2", alpha=55, E=200e9, A=600e-6, L=0.8, node1=node2, node2=node3)

    structure = Structure("struct")
    structure.add_element(bar1)
    structure.add_element(bar2)

    structure.setGlobalDisplacement()
    structure.setElementRxnForces()
    structure.setElementDisplacement()

    print(bar2.ReactionForces)

if __name__ == "__main__":
    main()
