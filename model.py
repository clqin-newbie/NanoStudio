import math
import numpy as np
from ase import Atoms
from pprint import pprint
from ase.visualize import view

class Nano(object):
    """
    Abstract base class for nano-structure construction
    Base class for graphene, twisted graphene and carbon nanotube
    """
    def __init__(self, n, m, period, pbc, bond_length, atom_type, vacuum):
        """
        Initialize basic parameters for nano-structure
        :param n: Chiral index n
        :param m: Chiral index m
        :param period: Supercell repetition period along x,y,z
        :param pbc: Periodic boundary conditions for x,y,z directions (bool list)
        :param bond_length: C-C bond length (unit: Angstrom)
        :param atom_type: Element symbol list of atoms
        :param vacuum: Vacuum layer thickness for non-periodic direction (unit: Angstrom)
        """
        self.origin_n = n
        self.origin_m = m
        self.bond_length = bond_length
        self.atom_type = atom_type
        self.period = period 
        self.pbc = pbc
        self.vacuum = vacuum

    def rotation_a1a2(self, theta):
        """
        Rotation matrix under oblique coordinate system (a1, a2 basis)
        :param theta: Rotation angle (radian)
        :return: 2D rotation matrix for a1-a2 oblique lattice
        """
        q1 = math.cos(theta)
        q2 = 1 / math.sqrt(3) * math.sin(theta)
        matrics = np.array([[q1 - q2, -2 * q2], [2 * q2, q1 + q2]])  # 变换矩阵
        return matrics
    
    def a1a2_to_ij(self, theta):
        """
        Convert oblique coordinates (a1,a2) to Cartesian coordinates (i,j)
        i axis is parallel to a1 axis
        :param theta: Angle between a1 and a2 (radian)
        :return: Coordinate transformation matrix
        """
        
        q1 = math.cos(theta)
        q2 = math.sin(theta)
        trans_matrics = np.array([[1, q1], [0, q2]])
        return trans_matrics
    
    def calc_angle(theta):
        """
        Calculate equivalent angle for chiral rotation correction
        :param theta: Original angle (degree)
        :return: Corrected angle (degree)
        """
        q = theta // 30
        r = theta % 30
        if q % 2 == 1:
            return 30 - r
        else:
            return r
        
    def rotation_ij(self, theta):
        """
        2D rotation matrix under Cartesian coordinate system (i,j basis)
        :param theta: Rotation angle (radian)
        :return: Standard 2D rotation matrix
        """
        q1 = math.cos(theta)
        q2 = math.sin(theta)
        rotation_matrics = np.array([[q1, -q2], [q2, q1]])  # 旋转矩阵
        return rotation_matrics

    def _get_structures(self, lattice, atom_coord_info, scaled=True):
        """
        Generate ASE Atoms object from lattice and atomic coordinates
        :param lattice: 3x3 cell matrix of simulation box
        :param atom_coord_info: Array contains coordinates and element symbols
        :param scaled: Whether to use scaled fractional coordinates
        :return: ASE Atoms object
        """
        atom_coord = atom_coord_info[:, 0:3].astype('float')
        elements = list(atom_coord_info[:, 3])
        atoms = Atoms(elements, positions=atom_coord, cell=lattice, pbc=self.pbc)
        if scaled:
            scaled_positions = atoms.get_scaled_positions()
            scaled_positions[:, 0] %=  1
            atoms.set_scaled_positions(scaled_positions)
        
        for idx, condi in enumerate(self.pbc):
            if condi == False:
                atoms.center(vacuum=self.vacuum, axis=idx)
        return atoms

    def extended_gcd(self, a, b):
        """
        Extended Euclidean Algorithm
        Solve a*x + b*y = gcd(a, b)
        :param a: Integer a
        :param b: Integer b
        :return: (gcd value, coefficient x, coefficient y)
        """
        if b == 0:
            return (a, 1, 0)
        else:
            g, x1, y1 = self.extended_gcd(b, a % b)
            x = y1
            y = x1 - (a // b) * y1
            return (g, x, y)

    def show_properties(self, layer_num=1):
        """Print structural properties in console"""
        properties = self._get_properties(layer_num)
        properties.pop('lattice')
        pprint(properties)

    def _get_properties(self):
        """
        Calculate basic structural properties
        :return: Dictionary contains chiral angle, pbc, bond length etc.
        """
        chiral_theta = math.acos(((2 * self.origin_n) + self.origin_m) / (
                    2 * math.sqrt(self.origin_m * self.origin_m + self.origin_m * self.origin_n + self.origin_n * self.origin_n))) * 180 / math.pi # 手性角度（度）
        properties = {"pbc": self.pbc, "period": self.period, "chirality (n, m)": (self.origin_n, self.origin_m), "bond_length": self.bond_length, "chiral_theta":  chiral_theta}
        return properties

class Graphene(Nano):
    """Subclass for single/multi-layer graphene structure generation"""
    def __init__(self, n, m, p, q, xy_period=[1, 1], xy_pbc=[True, True], bond_length=1.42, atom_type=['C'], vacuum=7.5):
        """
        Initialize graphene parameters
        :param n: Chiral index n
        :param m: Chiral index m
        :param p: Component of translation vector T
        :param q: Component of translation vector T
        :param xy_period: Supercell repetition along x, y
        :param xy_pbc: PBC for x, y directions
        :param bond_length: C-C bond length (Ang)
        :param atom_type: Atomic element
        :param vacuum: Vacuum layer thickness (Ang)
        """
        xy_pbc.append(False)
        xy_period.append(0)
        super().__init__(n, m, xy_period, xy_pbc, bond_length, atom_type, vacuum)
        self.n = n * xy_period[0]
        self.m = m * xy_period[0]
        self.origin_p = p
        self.origin_q = q
        self.p = p * xy_period[1]
        self.q = q * xy_period[1]

    def get_structures(self, layer_num=1, translations=[(0, 0)], layer_spacing=3.4):
        """
        Build multi-layer graphene structure
        :param layer_num: Number of graphene layers
        :param translations: In-plane translation for each layer (fractional coordinate)
        :param layer_spacing: Interlayer distance (Ang)
        :return: ASE Atoms object of multi-layer graphene
        """
        properties = self._get_properties()
        lattice = properties.pop('lattice')
        atom_info_list = []
        if layer_num != len(translations):
            raise ValueError("num必须和translations长度相同")
        temp_layer_spacing = 0 
        for trans in translations:
            atom_coord_info = self._build_coords(properties, translation=trans)
            atom_coord_info[:, 2] += temp_layer_spacing
            atom_info_list.append(atom_coord_info)
            temp_layer_spacing += layer_spacing
        new_atom_coord_info = np.vstack(atom_info_list)
        struct = self._get_structures(lattice, new_atom_coord_info)
        return struct
   
    def _build_coords(self, properties, translation):
        """
        Generate atomic fractional coordinates for single graphene layer
        :param properties: Pre-calculated structural properties
        :param translation: In-plane translation offset
        :return: Array of atomic coordinates + element symbols
        """
        a = self.bond_length * math.sqrt(3)  # 基矢的长度
        # 手性角

        N = properties['total_atom_num']  # 石墨烯的原子数
    
        C = np.array([[self.n], [self.m]]) # 以a1,a2为基

        A = np.matrix([[1/3+translation[0]], [1/3+translation[1]]])  # 以a1,a2为基
        
        AB = np.matrix([[1 / 3], [1 / 3]])  # AB向量以a1,a2为基
        b_a1_move_length = AB[0, 0]  # A原子沿a1轴移动的距离 -->B
        b_a2_move_length = AB[1, 0]  # A原子沿a2轴移动的距离 -->B
        g2 = math.gcd(self.n, self.m)  # 基础旋转单元的个数
  
        _, z, s = self.extended_gcd(self.m, -self.n)
 
        h = np.array([[z], [s]])  # 以a1,a2为基
        
        # 基础原子组移动的参数
        C_move =  C / g2
        base_a1_move, base_a2_move = C_move[0, 0], C_move[1, 0]
 
        # 基础单元移动的参数
        unit_a1_move = h[0, 0]  # 沿着a1轴移动的距离
        unit_a2_move = h[1, 0]  # 沿着a2轴移动的距离
  
        single_atom_num = int(N / 2)  # 单类原子的个数
        A_atom_coord = np.zeros([single_atom_num, 4], dtype='O')
        B_atom_coord = np.zeros([single_atom_num, 4], dtype='O')

        # 旋转得到所有的原子
        index = 0
        for j in range(g2):  # 基础原子移动的次数，得到基础单元
            for i in range(int(single_atom_num / g2)):  # 基础单元移动的次数
                a1_move = float(j * base_a1_move + i * unit_a1_move)
                a2_move = float(j * base_a2_move + i * unit_a2_move)
                A_temp = A.T + np.array([a1_move, a2_move])
                B_temp = A.T + np.matrix([a1_move + b_a1_move_length, a2_move + b_a2_move_length])
                A_atom_coord[index, 0:2] = A_temp
                B_atom_coord[index, 0:2] = B_temp
                index += 1
        layer = np.vstack((A_atom_coord, B_atom_coord))
        # 判断原子的种类数，并给出相应的元素符号
        num = len(self.atom_type)
        if num == 1:
            layer[:, 3] = self.atom_type[0]
        elif num == 2:
            layer[0:single_atom_num, 3] = self.atom_type[0]
            layer[single_atom_num:, 3] = self.atom_type[1]
     
        data1 = self.a1a2_to_ij(math.pi/3) @  layer[:, 0:2].T  * a  # 以i,j为基
        layer[:, 0:2] = data1.T
    
        return layer

    def _get_properties(self, layer_num=1):
        """
        Calculate lattice, atom number and other properties for graphene
        :param layer_num: Layer number
        :return: Dictionary of structural properties
        """
        properties = super()._get_properties()
        a = self.bond_length * math.sqrt(3)
        a1 = np.array([[1], [0]]) * a  # 以i，j为基
        a2 = np.array([[0.5], [math.sqrt(3) / 2]]) * a  # 以i，j为基
   
        N = int(2*np.abs(self.n*self.q-self.m*self.p)) * layer_num # 石墨烯的原子数
        C = self.n * a1 + self.m * a2  # 手性向量
        T = self.p * a1 + self.q * a2   # 
 
        lattice = list(np.array(C).flatten()) + [0] + list(np.array(T).flatten()) + [0] + [0, 0, 0.1]
        lattice = np.array(lattice).reshape(3, 3)
        properties.update({"(p, q)": (self.p, self.q), "layer_num": layer_num, "bond_length": self.bond_length, 'total_atom_num': N, 'lattice': lattice})
        return properties
 

class TwistGraphene(Nano):
    """Subclass for twisted bilayer / multi-layer graphene (TBG)"""
    def __init__(self, n, m, xy_period=(1, 1), xy_pbc=[True, True], bond_length=1.42, atom_type=['C'], vacuum=7.5):
        """
        Initialize twisted graphene parameters
        :param n: Chiral index n
        :param m: Chiral index m
        :param xy_period: Supercell repetition along x,y
        :param xy_pbc: PBC for x,y directions
        :param bond_length: C-C bond length (Ang)
        :param atom_type: Atomic element
        :param vacuum: Vacuum layer thickness (Ang)
        """
        xy_pbc.append(False)
        xy_period.append(0)
        super().__init__(n, m, xy_period, xy_pbc, bond_length=bond_length, atom_type=atom_type, vacuum=vacuum)
        self.n = n * xy_period[0]
        self.m = m * xy_period[0]
    
    def _get_single_layer(self, n, m, translation=(0, 0)):
        """
        Generate coordinates for single twisted graphene layer
        :param n: Chiral index
        :param m: Chiral index
        :param translation: In-plane translation offset
        :return: Atomic coordinate array of single layer
        """
        T = np.array(self.rotation_a1a2(np.pi/3) @ np.array([[n], [m]]))
        T = T.flatten().tolist()
        p = T[0] * self.period[1]
        q = T[1] * self.period[1]

        graphene = Graphene(n, m, p=p, q=q, xy_period=self.period[:2], xy_pbc=self.pbc[:2],  bond_length=self.bond_length, atom_type=self.atom_type, vacuum=self.vacuum)
        properties = graphene._get_properties()
        layer = graphene._build_coords(properties, translation)
        return layer
    
    def get_structures(self, layer_num=2, translations=[(0, 0), (0, 0)], layer_sequence=[0, 1], layer_spacing=3.4):
        """
        Build twisted multi-layer graphene structure
        :param layer_num: Total layer number
        :param translations: In-plane translation for each layer
        :param layer_sequence: Layer type index (0: original, 1: rotated)
        :param layer_spacing: Interlayer distance (Ang)
        :return: ASE Atoms object of twisted graphene
        """
        properties = self._get_properties(layer_num=1)
        lattice = properties.pop('lattice')
        atom_info_list = []
        if layer_num != len(translations) and layer_num != len(layer_sequence):
            raise ValueError("layer_num, translations以及layer_sequence长度必须相同")
        
        temp_layer_spacing = 0 
        for seq, trans in zip(layer_sequence, translations):
            if seq == 0:             
                n, m = self.n, self.m
            else:
                n, m = self.m, self.n
            atom_coord_info = self._get_single_layer(n, m, translation=trans)
            atom_coord_info[:, 2] += temp_layer_spacing
            if seq != 0:
                temp = self.rotation_ij(math.radians(-properties['twist_gamma'])) @ atom_coord_info[:, :2].T
                atom_coord_info[:, :2] = temp.T
            atom_info_list.append(atom_coord_info)
            temp_layer_spacing += layer_spacing
   
        new_atom_coord_info = np.vstack(atom_info_list)
        struct = self._get_structures(lattice, new_atom_coord_info)
        return struct

    def _get_properties(self, layer_num=2):
        """
        Calculate twist angle, lattice and atom number for twisted graphene
        :param layer_num: Layer number
        :return: Dictionary of structural properties
        """
        properties = super()._get_properties()
        a = self.bond_length * math.sqrt(3)
        a1 = np.matrix([[1], [0]]) * a  # 以i，j为基
        a2 = np.matrix([[0.5], [math.sqrt(3) / 2]]) * a  # 以i，j为基
        twist_gamma = math.acos(0.5 * (self.m * self.m + self.n * self.n + 4 * self.m * self.n) / (
                    self.m * self.m + self.n * self.n + self.m * self.n)) * 180 / math.pi  # 魔角（弧度）
      
        N = int(2 * ((self.n + self.m) * self.n + self.m * self.m)) * layer_num # 单层石墨烯的原子数
        
        C = self.n * a1 + self.m * a2 # 手性向量
        T = self.rotation_ij(math.pi / 3) * C # 

        lattice = list(np.array(C).flatten()) + [0] + list(np.array(T).flatten()) + [0] + [0, 0, 0.1]
        lattice = np.array(lattice).reshape(3, 3)
        properties.update({"layer_num": layer_num, 'lattice': lattice, 'total_atom_num': 2 * N, 'twist_gamma': twist_gamma})
        return properties


class Nanotube(Nano):
    """Subclass for single/multi-wall carbon nanotube (CNT)"""

    def __init__(self, n=[10], m=[0], z_period=[1], z_pbc=True, bond_length=1.42, atom_type=['C'], vacuum=7.5):
        """
        Initialize carbon nanotube parameters
        :param n: List of chiral index n for each tube
        :param m: List of chiral index m for each tube
        :param z_period: Supercell repetition along tube axis (z)
        :param z_pbc: PBC along tube axis
        :param bond_length: C-C bond length (Ang)
        :param atom_type: Atomic element
        :param vacuum: Vacuum layer for radial direction (Ang)
        """
        pbc = [False, False, z_pbc]
        self.n = n
        self.m = m
        # period = [0, 0, z_period]
        self.z_period = z_period
        self.pbc =pbc
        self.bond_length = bond_length
        self.atom_type =atom_type
        self.vacuum = vacuum
        
        # 平移T向量
        self.p, self.q = [], []
        for n, m in zip(n, m):
            g1 = math.gcd(n + 2 * m, m + 2 * n)  # n+2m,m+2n的最大公约数
            self.p.append(-(2 * m + n) / g1)
            self.q.append((m + 2 * n) / g1)
    
    def get_structures(self, translation=[(0, 0)]):
        """
        Build carbon nanotube structure
        :param translation: In-plane offset for each tube
        :return: ASE Atoms object of nanotube(s)
        """
        properties, layer_properties = self._get_properties()
        lattice = properties.pop('lattice')
        atom_coord_info = []
        for i, (n, m, trans) in enumerate(zip(self.n, self.m, translation)):
            N = layer_properties[f"layer{i+1}"]['total_atom_num']  # 石墨烯的原子数
            d = layer_properties[f"layer{i+1}"]['diameter']
            atom_coord = self._build_coords(n, m, N, d, translation=trans)
            atom_coord_info.append(atom_coord)
        atom_coord_info = np.vstack(atom_coord_info)
        struct = self._get_structures(lattice, atom_coord_info, scaled=False)
        return struct
   
    def _build_coords(self, n, m, N, d, translation):
        """
        Generate 3D cylindrical coordinates for single carbon nanotube
        :param n: Chiral index n
        :param m: Chiral index m
        :param N: Total atom number
        :param d: Tube diameter (Ang)
        :param translation: In-plane translation offset
        :return: Atomic coordinate array + element symbols
        """
        a = self.bond_length * math.sqrt(3)  # 基矢的长度
        C = self.a1a2_to_ij(math.pi/3) @ np.array([[n], [m]]) * a # 以i，j为基
        C_norm = np.linalg.norm(C)  # 手性向量的模

        A = self.a1a2_to_ij(math.pi/3) @ np.array([[1/3+translation[0]], [1/3+translation[1]]])  # 以i，j为基
        AB = self.a1a2_to_ij(math.pi/3) @ np.array([[1 / 3], [1 / 3]]) * a # AB向量以C，T为基
        b_rotation_angel = (2 * math.pi * AB.T@C / C_norm ** 2).item()   # A原子沿原点旋转的角度（弧度） -->B
        b_T_move = np.linalg.norm(np.cross(AB.flatten(), C.flatten())) / C_norm  # A原子沿T轴移动的距离 --> B
      
        g2 = math.gcd(n, m)  # 基础旋转单元的个数

        _, z, s = self.extended_gcd(m, -n)

        H = self.a1a2_to_ij(math.pi/3) @ np.matrix([[z], [s]]) * a # 以i，j为基
        # 基础原子旋转的参数
        base_rotation_angel = 2 * math.pi / g2
        
        # 基础旋转单元旋转的参数
        unit_rotation_angel = (2 * math.pi / C_norm ** 2 * H.T @ C).item() # 沿原点旋转的角度
        unit_T_move = np.linalg.norm(np.cross(H.flatten(), C.flatten())) / C_norm
    
        single_atom_num = int(N / 2)  # 单类原子的个数
        A_atom_coord = np.zeros([single_atom_num, 4], dtype='O')
        B_atom_coord = np.zeros([single_atom_num, 4], dtype='O')
  
        origin = A.flatten().tolist()  # 原点
        origin.append(0)
 
        # 旋转得到所有的原子
        index = 0
        for j in range(g2):  # 基础原子移动的次数，得到基础单元
            for i in range(int(single_atom_num / g2)):  # 基础单元移动的次数
                rotation_angel = j * base_rotation_angel + i * unit_rotation_angel
                T_move = float(i * unit_T_move)

                A_atom_coord[index, 0:3] = origin + d / 2 * np.array(
                [math.cos(rotation_angel), math.sin(rotation_angel), 0]) + np.array([0, 0, T_move])
                B_atom_coord[index, 0:3] = origin + d / 2 * np.array(
                [math.cos(rotation_angel + b_rotation_angel), math.sin(rotation_angel + b_rotation_angel),
                    0]) + np.array([0, 0, b_T_move + T_move])
                index += 1
        layer1 = np.vstack((A_atom_coord, B_atom_coord))

        # 判断原子的种类数，并给出相应的元素符号
        num = len(self.atom_type)
        if num == 1:
            layer1[:, 3] = self.atom_type[0]
        elif num == 2:
            layer1[0:single_atom_num, 3] = self.atom_type[0]
            layer1[single_atom_num:, 3] = self.atom_type[1]
        atom_coord_info = layer1
        return atom_coord_info
    
    def _get_properties(self):
        """
        Calculate diameter, length, atom number and lattice for nanotube(s)
        :return: Global properties dict, per-tube properties dict
        """
        layer_properties = {}
        a = self.bond_length * math.sqrt(3)
        layer = 1
        max_d = 0
        max_T = 0
        total_atom_num = 0
        for n, m, p, q, z_num in zip(self.n, self.m, self.p, self.q, self.z_period):
            chiral_theta = math.acos(((2 * n) + m) / (
                    2 * math.sqrt(m * m + m * n + n * n))) * 180 / math.pi # 手性角度（度）
            N = int(2*np.abs(n*q-m*p)) * z_num # 石墨烯的原子数
            C = self.a1a2_to_ij(math.pi/3) @ np.array([[n], [m]]) * a # 以i，j为基
            T = self.a1a2_to_ij(math.pi/3) @ np.array([[p], [q]]) * a *z_num # 以i，j为基
            C_norm = np.linalg.norm(C)  # 手性向量的摸
            d = float(C_norm / math.pi)  # 纳米管的直径
            T_norm = float(np.linalg.norm(T))  # 平移向量即周期的长度
            layer_properties[f"layer{layer}"] = {"chirality (n, m)": (n, m), "period": (0, 0, z_num), "length": T_norm, 'diameter': d, 'total_atom_num': N, "chiral_theta":  chiral_theta}
            layer += 1
            total_atom_num += N
            if d > max_d:
                max_d = d
            if T_norm > max_T:
                max_T = T_norm
        lattice = [max_d, 0, 0, 0, max_d, 0, 0, 0, max_T]
        lattice = np.array(lattice).reshape(3, 3)
        properties = {'total_atom_num': total_atom_num, "pbc": self.pbc, "lattice": lattice, "max_diameter": max_d, "max_length": max_T, "pbc": self.pbc, "bond_length": self.bond_length}
        return properties, layer_properties
    
    def show_properties(self):
        """Print all nanotube structural properties"""
        properties, layer_properties = self._get_properties()
        properties.pop('lattice')
        pprint(properties)
        pprint(layer_properties)

if __name__ == '__main__':
    # 石墨烯
    # model = Graphene(4, 2, -2, 6, xy_period=[1, 1], xy_pbc=[True, True])
    # struct = model.get_structures(layer_num=3, translations=[(0, 0), (1/3, 1/3), (-1/3, 2/3)])
    # # struct = model.get_structures()
    # properties = model.show_properties(layer_num=3)
    # struct.write("./1.vasp", format='vasp')
    # view(struct)

    # # 魔角
    # model = TwistGraphene(4, 2, xy_period=[2, 2], xy_pbc=[True, True])
    # struct = model.get_structures(layer_num=3, translations=[(0, 0), (0, 0), (0, 0)], layer_sequence=[1, 0, 1])
    # properties = model.show_properties(layer_num=3)
    # struct.write("./1.vasp", format='vasp')
    # view(struct)

    # 纳米管
    model = Nanotube([15, 6], [15, 6], z_period=[5, 5], atom_type=['C'], z_pbc=False)
    struct = model.get_structures(translation=[(0, 0), (0, 0)])
    properties = model.show_properties()
    struct.write("./1.vasp", format='vasp')
    view(struct)
    model.show_properties()
