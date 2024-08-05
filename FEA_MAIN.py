import numpy as np
import matplotlib.pyplot as plt
np.set_printoptions(linewidth=200, suppress=True)

class Node:
    """
    Nodes class smallest thing in a FEA model\n
    All other classes the built up by nodes. A node must have a name 
    """
    def __init__(self, name, XPos, YPos, DoFX, DoFY, EternalLoadX, EternalLoadY, ZPos = 0, DoFRotZ = False, EternalMomentZ = 0):
        """
        Class for Nodes that make up elements
        """
        self.name = name
        self.XPos = XPos
        self.YPos = YPos
        self.ZPos = ZPos
        self.DofX = DoFX
        self.DofY = DoFY
        self.DofRotZ = DoFRotZ
        self.EternalLoadX = EternalLoadX
        self.EternalLoadY = EternalLoadY
        self.EternalMomentZ = EternalMomentZ

    def __repr__(self):
        return (f"Node(name={self.name}")

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
                                     [0, 0, np.cos(alpha), np.sin(alpha)],
                                     [0,0,0,0],
                                     [0,0,0,1]])
        
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

class FrameElement:
    def __init__(self, name: str, alpha: int, E: float, A: float, I:float, L: float, node1: Node, node2: Node):
        """
        alpha is in degrees
        """
        self.name = name
        self.alpha = alpha * np.pi / 180
        self.E = E
        self.A = A
        self.I = I
        self.L = L
        self.node1 = node1
        self.node2 = node2
        self.localStiffnessmatrix = None
        self.StiffnessMatrix = None
        self.AssemblyMatrix = None
        self.GlobalForces = None
        self.TransMatrix = None
        self.deflectionlocal = None
        self.Globaldeflection = None
        self.forces = None
        self.stress = None
        self.strain = None

    def __repr__(self):
        return (f"FrameElement(name={self.name}, "
                f"node1={self.node1}, node2={self.node2} ")
    
    def setlocalStiffnessMatrix(self):
        """
        returns the global stiffness matrix as a numpy nested list ie [[],[]]
        """
        beta = (self.A * (self.L**2)) / self.I
        Lsquared = self.L**2

        K_e = ((self.E * self.I) / (self.L**3)) * np.array([[beta,  0,          0,          -beta,  0,          0],
                                                            [0,     12,         6*self.L,   0,      -12,        6*self.L],
                                                            [0,     6*self.L,   4*Lsquared, 0,      -6*self.L,  2*Lsquared],
                                                            [-beta, 0,          0,          beta,   0,          0],
                                                            [0,     -12,        -6*self.L,  0,      12,         -6*self.L],
                                                            [0,     6*self.L,   2*Lsquared, 0,      -6*self.L,  4*Lsquared]])
        self.localStiffnessmatrix = K_e

    def setTransMatrix(self):
        if self.alpha is None:
            alpha = np.arctan2(self.node2.YPos - self.node1.YPos, self.node2.XPos - self.node1.XPos)
        else:
            alpha = self.alpha

        submatrix = np.array([[np.cos(alpha),   np.sin(alpha), 0],
                              [-np.sin(alpha),  np.cos(alpha), 0],
                              [0,               0,             1]])
        zeros = np.zeros_like(submatrix)

        temp = np.vstack((submatrix, zeros))
        temp2 = np.vstack((zeros, submatrix))
        temp3 = np.hstack((temp, temp2))
        self.TransMatrix = temp3
        
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
        A_e = np.zeros((dofnum,6))
        for node in nodelist:
            if node == self.node1:
                i = nodelist.index(node)
            if node == self.node2:
                i2 = nodelist.index(node)
        if self.node1.DofX :
            A_e[dofmap[i][0] - 1][0] = 1
        if self.node2.DofX:
            A_e[dofmap[i2][0] - 1][3] = 1

        if self.node1.DofY:
            A_e[dofmap[i][1] - 1][1] = 1
        if self.node2.DofY:
            A_e[dofmap[i2][1] - 1][4] = 1

        if self.node1.DofRotZ:
            A_e[dofmap[i][2] - 1][2] = 1
        if self.node2.DofRotZ:
            A_e[dofmap[i2][2] - 1][5] = 1

        self.AssemblyMatrix = A_e

    def interpolation(self, numpoints, point=None):
        Psi = lambda x: 1 - x/self.L
        Psi2 = lambda x: x/self.L

        N1 = lambda x: (1 -                         ((3*x**2) / (self.L**2)) +  ((2*x**3) / (self.L**3)))
        N2 = lambda x: (((x**3) / (self.L**2)) -    ((2*x**2) / (self.L)) +     x)
        N3 = lambda x: (((3*x**2) / (self.L**2)) -  ((2*x**3) / (self.L**3)))
        N4 = lambda x: (((x**3) / (self.L**2)) -    ((x**2) / (self.L)))

        u = lambda x: Psi(x)*self.deflectionlocal[0][0] + Psi2(x) * self.deflectionlocal[3][0]
        v = lambda x: N1(x)*self.deflectionlocal[1][0] + N2(x)*self.deflectionlocal[2][0] + N3(x)*self.deflectionlocal[4][0] + N4(x)*self.deflectionlocal[5][0]
        if point is None:
            x = np.linspace(0, self.L, numpoints)

            xs_interpolated = u(x) * np.cos(self.alpha) - v(x) * np.sin(self.alpha)
            ys_interpolated = u(x) * np.sin(self.alpha) + v(x) * np.cos(self.alpha)

            return (xs_interpolated, ys_interpolated)
        else:
            return(u(point)+v(point))





class Structure:
    """
    Is the overall structure that is made up of elements.\n
    Takes elements by add_elements
    """
    def __init__(self, name:str):
        self.name = name
        self.nodes = []
        self.elements = []
        self.GlobalStiffnessMatrix = None
        self.externalLoads = None
        self.GlobalDisplacement = None
        self.Totdof = None
        self.dofmap = None
        self.ReactionForces = None
        self.ax = plt.subplots()[1]

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
            nodemap = [0,0,0]
            if node.DofX:
                self.Totdof += 1
                nodemap[0] = self.Totdof
            if node.DofY:
                self.Totdof += 1
                nodemap[1] = self.Totdof
            if node.DofRotZ:
                self.Totdof += 1
                nodemap[2] = self.Totdof
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
            if node.DofRotZ:
                Q.append([node.EternalMomentZ])
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
            element.Globaldeflection = element.AssemblyMatrix.T @ self.GlobalDisplacement
            element.deflectionlocal = element.TransMatrix @ element.Globaldeflection

    def setElementForces(self):
        for element in self.elements:
            if element.localStiffnessmatrix is None:
                element.setlocalStiffnessMatrix
            if element.deflectionlocal is None:
                self.setElementDisplacement()
            element.forces = element.localStiffnessmatrix @ element.deflectionlocal

    def setReactionForces(self):
        rxn = []
        for i in range(len(self.nodes)):
            if (not (self.nodes[i].DofX and self.nodes[i].DofY and self.nodes[i].DofRotZ)):
                elementswithReactions = []
                for element in self.elements:
                    if ((self.nodes[i] == element.node1 ) or (self.nodes[i] == element.node2)):
                        elementswithReactions.append(element)
                rxnx = 0
                rxny = 0
                rxnRotz = 0
                for elerxn in elementswithReactions:
                    if elerxn.GlobalForces is None:
                        self.setElementGlobalForces()
                    rxnx += elerxn.GlobalForces[0][0]
                    rxny += elerxn.GlobalForces[1][0]
                    rxnRotz += elerxn.GlobalForces[2][0]
                rxn.append([rxnx, rxny, rxnRotz])
            else:
                rxn.append([0,0,0])
        self.ReactionForces = np.array(rxn)

    def setElementStressStrain(self):
        for element in self.elements:
            if element.deflectionlocal is None:
                self.setElementDisplacement()
            element.strain = (element.deflectionlocal[1][0] - element.deflectionlocal[0][0]) / element.L
            element.stress = element.E * element.strain

    def plotelements(self, numpoints):
        for element in self.elements:
            xs = np.linspace(element.node1.XPos, element.node2.XPos, numpoints)
            ys = np.linspace(element.node1.YPos, element.node2.YPos, numpoints)
            self.ax.plot(xs, ys, 'b.-')

    def plotdeflectedelements(self, scale, numpoints):
        for element in self.elements:
            if element.Globaldeflection is None:
                self.setElementDisplacement()

            xs = np.linspace(element.node1.XPos, element.node2.XPos, numpoints)
            ys = np.linspace(element.node1.YPos, element.node2.YPos, numpoints)
            xs_inter, ys_inter = element.interpolation(numpoints)
            xs_defle = xs + xs_inter * scale
            ys_defle = ys + ys_inter * scale
            self.ax.plot(xs_defle, ys_defle, 'r.-')

    def plot(self, scale, numpoints = 2):
        self.plotelements(numpoints)
        self.plotdeflectedelements(scale, numpoints)
        self.ax.grid(True)
        plt.show()

    def solve(self):
        self.setGlobalStiffnessMatrix()
        self.setexeternalLoads()
        self.setGlobalDisplacement()
        self.setElementGlobalForces()
        self.setElementDisplacement()
        self.setElementForces()
        self.setReactionForces()
        self.setElementStressStrain()

def main():
    node1 = Node("Node1", 0, 0, False, False, 0, 0)
    node2 = Node("Node2", 0, 3, True, True, 10000, 0, DoFRotZ= True)
    node3 = Node("Node3", 4.5, 3, True, True, 10000, 0, DoFRotZ=True)
    node4 = Node("Node4", 4.5, 0, False, False, 0, 0)

    frame1 = FrameElement("Frame1", alpha= 90, E= 200e9, A= 5e-4, I=1e-5, L=3, node1= node1, node2= node2)
    frame2 = FrameElement("Frame2", alpha= 0, E= 200e9, A= 5e-4, I=1e-5, L=4., node1= node2, node2= node3)
    frame3 = FrameElement("Frame3", alpha=270, E=200e9, A=5e-4, I=1e-5, L=3, node1=node3, node2=node4)


    structure = Structure("struct")
    structure.add_element(frame1)
    structure.add_element(frame2)
    structure.add_element(frame3)

    structure.solve()
    print(structure.GlobalDisplacement* 1e3)
    structure.plot(20, 10)


    #print(frame1.deflection)
    #structure.setElementDisplacement()
    #print(frame1.Globaldeflection)
    

if __name__ == "__main__":
    main()


        


