import numpy as np
np.set_printoptions(linewidth=200)

class Node:
    """
    Nodes class smallest thing in a FEA model\n
    All other classes the built up by nodes. A node must have a name 
    """
    def __init__(self, name, XPos, YPos, DoFX, DoFY, EternalLoadX, EternalLoadY, ZPos = None, DoFZ = None, EternalLoadZ = None):
        """
        Class for Nodes that make up elements
        """
        self.name = name
        self.XPos = XPos
        self.YPos = YPos
        self.ZPos = ZPos
        self.DofX = DoFX
        self.DofY = DoFY
        self.DofZ = DoFZ
        self.EternalLoadX = EternalLoadX
        self.EternalLoadY = EternalLoadY
        self.EternalLoadZ = EternalLoadZ

    def __repr__(self):
        return (f"Node(name={self.name},"
                f"LocalDeflection={self.LocalDeflection}, GlobalDeflection={self.GlobalDeflection})")

class BarElement:
    """
    Bar Element: Made up of 2 nodes\n
    Alpha is assumed to be in degrees not needs if both element pos are defined
    """
    def __init__(self, name, alpha, E, A, L, node1: Node, node2: Node):
        """
        Class for Bar elements.
        alpha is in degrees
        """
        self.name = name
        self.alpha = alpha * np.pi / 180
        self.E = E
        self.A = A
        self.L = L
        self.node1 = node1
        self.node2 = node2
        self.localStiffnessmatrix = None
        self.StiffnessMatrix = None
        self.AssemblyMatrix = None
        self.GlobalForces = None
        self.TransMatrix = None
        self.deflection = None
        self.forces = None
        self.stress = None
        self.strain = None

    def __repr__(self):
        return (f"BarElement(name={self.name}, "
                f"node1={self.node1}, node2={self.node2}, "
                f"localStiffnessmatrix=\n{self.localStiffnessmatrix}, "
                f"StiffnessMatrix=\n{self.StiffnessMatrix})")
    
    def setlocalStiffnessMatrix(self):
        """
        returns the global stiffness matrix as a numpy nested list ie [[],[]]
        """
        K_e = (self.E * self.A / self.L) * np.array([[1, -1],
                                                     [-1, 1]])
        self.localStiffnessmatrix = K_e

    def setTransMatrix(self):
        if self.alpha is None:
            alpha = np.arctan2(self.node2.YPos - self.node1.YPos, self.node2.XPos - self.node1.XPos)
        else:
            alpha = self.alpha

        self.TransMatrix = np.array([[np.cos(alpha), np.sin(alpha), 0, 0],
                                 [0, 0, np.cos(alpha), np.sin(alpha)]])
        
    def setStiffnessMatrix(self):
        """
        returns the stiffness matrix as a numpy nested list ie [[],[]]\n
        IN GLOBAL CORDANATES 
        """
        if self.localStiffnessmatrix is None:
            self.setlocalStiffnessMatrix()
        if self.TransMatrix is None:
            self.setTransMatrix()
        K_e_global = np.transpose(self.TransMatrix) @ self.localStiffnessmatrix @ self.TransMatrix
        self.StiffnessMatrix = K_e_global

    def setAssemblyMatrix(self, dofnum, dofmap, nodelist): 
        A_e = np.zeros((dofnum,4))
        for node in nodelist:
            if node == self.node1:
                i = nodelist.index(node)
            if node == self.node2:
                i2 = nodelist.index(node)
        if self.node1.DofX :
            A_e[dofmap[i][0] - 1][0] = 1
        if self.node2.DofX:
            A_e[dofmap[i2][0] - 1][2] = 1
        if self.node1.DofY:
            A_e[dofmap[i][1] - 1][1] = 1
        if self.node2.DofY:
            A_e[dofmap[i2][1] - 1][3] = 1
        self.AssemblyMatrix = A_e

class BeamElement:
    def __init__(self, name, alpha, E, I, L, node1: Node, node2: Node):
        self.name = name
        self.alpha = alpha * np.pi / 180
        self.E = E
        self.I = I
        self.L = L
        self.node1 = node1
        self.node2 = node2
        self.localStiffnessmatrix = None
        self.StiffnessMatrix = None
        self.AssemblyMatrix = None
        self.GlobalForces = None
        self.TransMatrix = None
        self.deflection = None
        self.forces = None
        self.stress = None
        self.strain = None

    def __repr__(self):
        return (f"BeamElement(name={self.name}, "
                f"node1={self.node1}, node2={self.node2}, "
                f"localStiffnessmatrix=\n{self.localStiffnessmatrix}, "
                f"StiffnessMatrix=\n{self.StiffnessMatrix})")
    
    def setlocalStiffnessMatrix(self):
        """
        returns the global stiffness matrix as a numpy nested list ie [[],[]]
        """
        K_e = (self.E * self.I / (self.L ** 3)) * np.array([[12,        6*self.L,       -12,        6*self.L],
                                                            [6*self.L,  4*(self.L**2),  -6*self.L,  2*(self.L**2)],
                                                            [-12,       -6*self.L,      12,         -6*self.L],
                                                            [6*self.L,  2*(self.L**2),  -6*self.L,  4*(self.L**2)]])
        self.localStiffnessmatrix = K_e

    def setTransMatrix(self):
        if self.alpha is None:
            alpha = np.arctan2(self.node2.YPos - self.node1.YPos, self.node2.XPos - self.node1.XPos)
        else:
            alpha = self.alpha

        self.TransMatrix = np.array([[np.cos(alpha), np.sin(alpha), 0, 0],
                                     [0, 0, np.cos(alpha), np.sin(alpha)]])
        
    def setStiffnessMatrix(self):
        """
        returns the stiffness matrix as a numpy nested list ie [[],[]]\n
        IN GLOBAL CORDANATES 
        """
        if self.localStiffnessmatrix is None:
            self.setlocalStiffnessMatrix()
        if self.TransMatrix is None:
            self.setTransMatrix()
        K_e_global = np.transpose(self.TransMatrix) @ self.localStiffnessmatrix @ self.TransMatrix
        self.StiffnessMatrix = K_e_global

    def setAssemblyMatrix(self, dofnum, dofmap, nodelist): 
        A_e = np.zeros((dofnum,4))
        for node in nodelist:
            if node == self.node1:
                i = nodelist.index(node)
            if node == self.node2:
                i2 = nodelist.index(node)
        if self.node1.DofX :
            A_e[dofmap[i][0] - 1][0] = 1
        if self.node2.DofX:
            A_e[dofmap[i2][0] - 1][2] = 1
        if self.node1.DofY:
            A_e[dofmap[i][1] - 1][1] = 1
        if self.node2.DofY:
            A_e[dofmap[i2][1] - 1][3] = 1
        self.AssemblyMatrix = A_e

class Structure:
    """
    Is the overall structure that is made up of elements.\n
    Takes elements by add_elements
    """
    def __init__(self, name):
        self.name = name
        self.nodes = []
        self.elements = []
        self.GlobalStiffnessMatrix = None
        self.externalLoads = None
        self.GlobalDisplacement = None
        self.Totdof = None
        self.dofmap = None
        self.ReactionForces = None

    def __repr__(self):
        return (f"Structure(name={self.name},\n nodes=[{', '.join([str(node.name) for node in self.nodes])}], "
                f"\nelements=[{', '.join([str(element.name) for element in self.elements])}])")

    def add_element(self, element):
        self.elements.append(element)
        self.elements.sort(key=lambda s: s.name)
        self.addNodes()

    def addNodes(self):
        """
        adds all the nodes into a list for the assmebly matrix
        should be called by add_element
        """
        for element in self.elements:
            if element.node1 not in self.nodes:
                self.nodes.append(element.node1)
            if element.node2 not in self.nodes:
                self.nodes.append(element.node2)
        self.nodes.sort(key=lambda s: s.name)

    def setTotdof(self):
        self.Totdof = 0
        self.dofmap = []
        for node in self.nodes:
            nodemap = [0,0]
            if node.DofX:
                self.Totdof += 1
                nodemap[0] = self.Totdof
            if node.DofY:
                self.Totdof += 1
                nodemap[1] = self.Totdof
            self.dofmap.append(nodemap)

    def setGlobalStiffnessMatrix(self):
        k_g_e = []
        for element in self.elements:
            if element.AssemblyMatrix is None:
                if self.Totdof is None:
                    self.setTotdof()
                element.setAssemblyMatrix(self.Totdof, self.dofmap, self.nodes)
            if element.StiffnessMatrix is None:
                element.setStiffnessMatrix()
            k_g_e.append(element.AssemblyMatrix @ element.StiffnessMatrix @ element.AssemblyMatrix.T)
        self.GlobalStiffnessMatrix = np.array(sum(k_g_e))

    def setexeternalLoads(self):
        Q = []
        for node in self.nodes:
            if node.DofX :
                Q.append([node.EternalLoadX])
            if node.DofY:
                Q.append([node.EternalLoadY])
        self.externalLoads = np.array(Q)
            
    def setGlobalDisplacement(self):
        if self.GlobalStiffnessMatrix is None:
            self.setGlobalStiffnessMatrix()
        if self.externalLoads is None:
            self.setexeternalLoads()
        self.GlobalDisplacement = np.linalg.solve(self.GlobalStiffnessMatrix, self.externalLoads) 

    def setElementGlobalForces(self):
        if self.GlobalDisplacement is None:
            self.setGlobalDisplacement()
        for element in self.elements:
            if element.StiffnessMatrix is None:
                element.setStiffnessMatrix()
            if element.AssemblyMatrix is None:
                element.setAssemblyMatrix()
            element.GlobalForces = element.StiffnessMatrix @ element.AssemblyMatrix.T @ self.GlobalDisplacement

    def setElementDisplacement(self):
        if self.GlobalDisplacement is None:
            self.setGlobalDisplacement()
        for element in self.elements:
            if element.TransMatrix is None:
                element.setTransMatrix()
            if element.AssemblyMatrix is None:
                element.setAssemblyMatrix()
            element.deflection = element.TransMatrix @ element.AssemblyMatrix.T @ self.GlobalDisplacement

    def setElementForces(self):
        for element in self.elements:
            if element.localStiffnessmatrix is None:
                element.setlocalStiffnessMatrix
            if element.deflection is None:
                self.setElementDisplacement()
            element.forces = element.localStiffnessmatrix @ element.deflection

    def setReactionForces(self):
        rxn = []
        for i in range(len(self.nodes)):
            if (not (self.nodes[i].DofX and self.nodes[i].DofY)):
                elementswithReactions = []
                for element in self.elements:
                    if ((self.nodes[i] == element.node1 ) or (self.nodes[i] == element.node2)):
                        elementswithReactions.append(element)
                rxnx , rxny = (0,0)
                for elerxn in elementswithReactions:
                    if elerxn.GlobalForces is None:
                        self.setElementGlobalForces()
                    rxnx += elerxn.GlobalForces[0][0]
                    rxny += elerxn.GlobalForces[1][0]
                rxn.append([rxnx, rxny])
            else:
                rxn.append([0,0])
        self.ReactionForces = rxn

    def setElementStressStrain(self):
        for element in self.elements:
            if element.deflection is None:
                self.setElementDisplacement()
            element.strain = (element.deflection[1][0] - element.deflection[0][0]) / element.L
            element.stress = element.E * element.strain
        
def main():
    node1 = Node("Node1", 0, 0, False, False, 0, 0)
    node2 = Node("Node2", 0, 0, True, True, 10000, 40000)

    bar1 = BeamElement("bar1",alpha= 90,E= 200e9, I=2e-5, L=6, node1= node1, node2= node2)


    structure = Structure("struct")
    structure.add_element(bar1)

    structure.setElementDisplacement()

    print(bar1.deflection)
    

if __name__ == "__main__":
    main()


        


