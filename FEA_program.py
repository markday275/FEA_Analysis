import numpy as np
from FEA_MAIN import *

def main():
    node1 = Node("Node1", 0, 0, False, False, 0, 0)
    node2 = Node("Node2", 0, 0, True, True, 0, -20000)
    node3 = Node("Node3", 0, 0, False, False, 0, 0)

    bar1 = BarElement("bar1",alpha= 0,E= 200e9, A=400e-6, L=1.1, node1= node1, node2= node2)
    bar2 = BarElement("bar2",alpha= 55,E= 200e9, A=600e-6, L=0.8, node1= node2, node2= node3)


    structure = Structure("struct")
    structure.add_element(bar1)
    structure.add_element(bar2)

    structure.setGlobalDisplacement()

    print(structure.GlobalDisplacement)
    

if __name__ == "__main__":
    main()