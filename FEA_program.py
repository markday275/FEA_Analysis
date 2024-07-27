import numpy as np
import FEA_MAIN as fea

def main():
    node1 = fea.Node("Node1", 0, 0, False, False, 0, 0)
    node2 = fea.Node("Node2", 1, 0, True, True, 0, -20000)
    node3 = fea.Node("Node3", 0, 0, False, False, 0, 0)

    bar1 = fea.BarElement("bar1",alpha= 0,E= 200e9, A=400e-6, L=1.1, node1= node1, node2= node2)
    bar2 = fea.BarElement("bar2",alpha= 55,E= 200e9, A=600e-6, L=0.8, node1= node2, node2= node3)


    structure = fea.Structure("struct")
    structure.add_element(bar1)
    structure.add_element(bar2)

    structure.setGlobalDisplacement()

    print(structure.GlobalDisplacement)
    

if __name__ == "__main__":
    main()