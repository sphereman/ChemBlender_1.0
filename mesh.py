import bpy
import bmesh
import numpy as np
import re
from . import _math, read, node
from .Chem_data import ELEMENTS_DEFAULT

def create_object(coll, name, verts, edges, faces):
    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata(verts, edges, faces)
    object = bpy.data.objects.new(name, mesh)
    coll.objects.link(object)
    return object

# 结构参数
def get_link_vert(vert, edge):
    """
    已知边的一个点，获取另一个点
    """
    if edge.verts[0].index == vert.index:
        return edge.verts[1]
    else:
        return edge.verts[0]


# 属性设置
############################################################################################
attr_types = {
    'FLOAT': 'value',
    'INT': 'value',
    'BOOLEAN': 'value',
    'QUATERNION': 'value',
    'FLOAT_COLOR': 'color',
    'FLOAT_VECTOR': 'vector',
}

def get_attr(obj, attr_name, datatype, domain):
    """
    获取对象特定域上的特定属性
    """
    bpy.ops.object.mode_set(mode='OBJECT')
    bm = bmesh.new()
    bm.from_mesh(obj.data)
    attr = obj.data.attributes.get(attr_name)
    count = len(bm.verts) if domain == 'VERT' else len(bm.edges)
    if datatype == 'FLOAT_COLOR':
        m = 4
    elif datatype == 'FLOAT_VECTOR':
        m = 3
    else:
        m = 1
    tot_count = count*m
    seq = [0]*tot_count
    attr.data.foreach_get(attr_types[datatype],seq)
    attr.data.update()

    new_seq = []
    if datatype == 'FLOAT_COLOR':
        for i in range(count):
            new_seq.append(seq[m*i])
            new_seq.append(seq[m*i+1])
            new_seq.append(seq[m*i+2])
            new_seq.append(seq[m*i+3])
    elif datatype == 'FLOAT_VECTOR':
        for i in range(count):
            new_seq.append(seq[m*i])
            new_seq.append(seq[m*i+1])
            new_seq.append(seq[m*i+2])
    else:
        new_seq = seq
    bm.free()
    return new_seq


def get_layer(bm, datatype, domain, name):
    try:
        layers = {
            ('INT', 'VERT'): bm.verts.layers.int,
            ('INT', 'EDGE'): bm.edges.layers.int,
            ('FLOAT', 'VERT'): bm.verts.layers.float,
            ('FLOAT', 'EDGE'): bm.edges.layers.float,
            ('FLOAT_VECTOR', 'VERT'): bm.verts.layers.float_vector,
            ('FLOAT_VECTOR', 'EDGE'): bm.edges.layers.float_vector,
        }
        return layers[(datatype, domain)].get(name) or layers[(datatype, domain)].new(name)
    except KeyError:
        raise ValueError(f"Invalid datatype '{datatype}' or domain '{domain}'. Must be ('INT', 'FLOAT', 'FLOAT_VECTOR') and ('VERT', 'EDGE').")

def add_attr(mesh_obj, attr_name, datatype, domain, attr_values):
    """
    给MESH对象特定的域添加特定的属性值.
    """
    attr = mesh_obj.data.attributes.get(attr_name)
    if not attr:
        attr = mesh_obj.data.attributes.new(attr_name, datatype, domain)
    attr.data.foreach_set(attr_types[datatype], attr_values)
    return attr

def add_scaffold_attr(scaffold, ATOMS, AtomicNum, BOND_ORDERS, VDW_R, COVA_R, RingNum):
    #STICK_R = np.full(len(AtomicNum), 0.18) # [0.18 for i in AtomicNum]
    #WIRE_R = np.full(len(AtomicNum), 0.03)
    Aromatic = np.full(len(BOND_ORDERS), False)
    atom_colors = []
    for element in ATOMS:
        for value in ELEMENTS_DEFAULT[element.capitalize()][3]:
            atom_colors.append(value)
    attributes = (
        {'name': 'atomic_num',       'type':'INT',      'domain':'POINT',  'values':AtomicNum},
        {'name': 'vdw_radius',       'type':'FLOAT',    'domain':'POINT',  'values':VDW_R},
        {'name': 'covalent_radius',  'type':'FLOAT',    'domain':'POINT',  'values':COVA_R},
        #{'name': 'stick_radius',     'type':'FLOAT',    'domain':'POINT',  'values':STICK_R},
        #{'name': 'wire_radius',      'type':'FLOAT',    'domain':'POINT',  'values':WIRE_R},
        {'name': 'bond_order',       'type':'INT',      'domain':'EDGE',   'values':BOND_ORDERS},
        {'name': 'ring_num',         'type':'INT',      'domain':'EDGE',   'values':RingNum},
        {'name': 'is_aromatic',      'type':'BOOLEAN',   'domain':'EDGE',  'values':Aromatic},
        {'name': 'colour',     'type':'FLOAT_COLOR',    'domain':'POINT',  'values':atom_colors},
    ) 
    for attr in attributes:
        add_attr(scaffold, attr['name'], attr['type'], attr['domain'], attr['values'])

def mesh_to_mol_scaffold(obj):
    atoms_num = len(obj.data.vertices)
    bonds_num = len(obj.data.edges)
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_all(action='SELECT')
    bpy.ops.mesh.delete(type='ONLY_FACE')
    bpy.ops.object.mode_set(mode='OBJECT')
    ATOMS = ['C']*atoms_num
    AtomicNum = [6]*atoms_num
    VDW_R = [1.70]*atoms_num
    COVA_R = [0.76]*atoms_num
    BOND_ORDERS = [1]*bonds_num
    RingNum = [0]*bonds_num
    obj['Type'] = 'scaffold'
    obj['Elements'] = list(set(ATOMS))
    add_scaffold_attr(obj, ATOMS, AtomicNum, BOND_ORDERS, VDW_R, COVA_R, RingNum)


# 设置原子
def set_sel_atoms_attr(obj, atomic_num):
    from rdkit import Chem
    periodic_table = Chem.GetPeriodicTable()
    vdw_radius = periodic_table.GetRvdw(atomic_num)
    covalent_radius = periodic_table.GetRcovalent(atomic_num)
    sel_vert_idxs = get_sel_idxs(obj, 'VERT')
    Atomic_Nums = get_attr(obj, 'atomic_num', 'INT', 'VERT')
    VDW_R = get_attr(obj, 'vdw_radius', 'FLOAT', 'VERT')
    COVA_R = get_attr(obj, 'covalent_radius', 'FLOAT', 'VERT')
    COLORS = get_attr(obj, 'colour', 'FLOAT_COLOR', 'VERT')
    print(COLORS)
    element = list(ELEMENTS_DEFAULT.keys())[atomic_num]
    if atomic_num > 0:
        for idx in sel_vert_idxs:
            Atomic_Nums[idx] = atomic_num
            VDW_R[idx] = vdw_radius
            COVA_R[idx] = covalent_radius
            COLORS[4*idx],COLORS[4*idx+1],COLORS[4*idx+2],COLORS[4*idx+3] = ELEMENTS_DEFAULT[element][3]
    add_attr(obj, 'atomic_num', 'INT', 'VERT', Atomic_Nums)
    add_attr(obj, 'covalent_radius', 'FLOAT', 'VERT', COVA_R)
    add_attr(obj, 'vdw_radius', 'FLOAT', 'VERT', VDW_R)
    add_attr(obj, 'colour', 'FLOAT_COLOR', 'VERT', COLORS)


# 设置键
def set_sel_bonds_attr(obj, bond_order, ring_num):
    sel_edge_idxs = get_sel_idxs(obj, 'EDGE')
    Bond_Orders = get_attr(obj, 'bond_order','INT','EDGE')
    Ring_Nums = get_attr(obj, 'ring_num','INT','EDGE')
    Aromatic = get_attr(obj, 'is_aromatic','BOOLEAN','EDGE')
    if bond_order > 0:
        for idx in sel_edge_idxs:
            Bond_Orders[idx] = bond_order
            if bond_order == 12:
                Aromatic[idx] = True
    if ring_num > 0:
        for idx in sel_edge_idxs: Ring_Nums[idx] =ring_num
    add_attr(obj, 'bond_order', 'INT', 'EDGE', Bond_Orders)
    add_attr(obj, 'ring_num', 'INT', 'EDGE', Ring_Nums)
    add_attr(obj, 'is_aromatic', 'BOOLEAN', 'EDGE', Aromatic)


# 选择
############################################################################################
def select_all(obj, domain):
    bpy.ops.object.mode_set(mode='OBJECT')
    bpy.ops.object.select_all(action='DESELECT')
    obj.select_set(True)
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_mode(type=domain)
    bpy.ops.mesh.select_all(action='SELECT')

def deselect_all(obj):
    bpy.ops.object.mode_set(mode='OBJECT')
    bpy.ops.object.select_all(action='DESELECT')
    obj.select_set(True)
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_all(action='DESELECT')
    bpy.ops.object.mode_set(mode='OBJECT')

def get_sel_idxs(obj, domain):
    """
    获取所选对象特定域的编号
    """
    bpy.ops.object.mode_set(mode='EDIT')
    bm = bmesh.from_edit_mesh(obj.data)
    sel_idxs = []
    if domain == 'VERT':
        sel_idxs = [vert.index for vert in bm.verts if vert.select]
    elif domain == 'EDGE':
        sel_idxs = [edge.index for edge in bm.edges if edge.select]
    bmesh.update_edit_mesh(obj.data)
    bm.free()
    bpy.ops.object.mode_set(mode='OBJECT')
    return sel_idxs

def select_verts(obj, elements):   # 'C', 'H'
    bpy.ops.object.mode_set(mode='EDIT')
    bm = bmesh.from_edit_mesh(obj.data)
    atom_id = get_layer(bm, 'INT', 'VERT', 'atomic_num')
    Elements = list(ELEMENTS_DEFAULT.keys())
    for vert in bm.verts:
        if Elements[vert[atom_id]] in elements:
            vert.select_set(True)

def select_edges(obj, bond_syms):   # ['C','H']
    bpy.ops.object.mode_set(mode='EDIT')
    bm = bmesh.from_edit_mesh(obj.data)
    atom_id = get_layer(bm, 'INT', 'VERT', 'atomic_num')
    Elements = list(ELEMENTS_DEFAULT.keys())
    for edge in bm.edges:
        atom1 = Elements[edge.verts[0][atom_id]]
        atom2 = Elements[edge.verts[1][atom_id]]
        if [atom1,atom2] in bond_syms or [atom2,atom1] in bond_syms:
            edge.select_set(True)

def select_bond_orders(obj, bond_orders):
    bpy.ops.object.mode_set(mode='EDIT')
    bm = bmesh.from_edit_mesh(obj.data)
    bond_order = get_layer(bm, 'INT', 'EDGE', 'bond_order')
    for edge in bm.edges:
        if re.sub('[\W]','',f'{edge[bond_order]}') in bond_orders:
            edge.select_set(True)

# 数值计算
############################################################################################
def calc_length():
    """计算边的长度，需在编辑模式下选中边之后执行"""
    if bpy.context.mode == 'EDIT_MESH':
        obj = bpy.context.edit_object
        bm = bmesh.from_edit_mesh(obj.data)
        locations = []
        for edge in bm.edges:
            if edge.select:
                locations.append(edge.verts[0].co)
                locations.append(edge.verts[1].co)
                break
        if len(locations) == 0: return "N.A."

        pos1 = np.array(locations[0])
        pos2 = np.array(locations[1])
        distance = _math.get_length(pos1-pos2)
        distance = str(float('%.3f' % distance)) + " Å"   # 保留小数点后三位
        return distance

def calc_angle():
    """计算两边的夹角，需在编辑模式下选中两条边后执行"""
    if bpy.context.mode == 'EDIT_MESH':
        obj = bpy.context.edit_object
        bm = bmesh.from_edit_mesh(obj.data)
        bm.edges.ensure_lookup_table()
        locations = []
        for edge in bm.edges:
            if edge.select:
                locations.append(edge.verts[0].co)
                locations.append(edge.verts[1].co)
                if len(locations) == 4: break
        if len(locations) < 4: return "N.A."

        pos1 = np.array(locations[0])
        pos2 = np.array(locations[1])
        pos3 = np.array(locations[2])
        pos4 = np.array(locations[3])
        epsilon = 0.000001
        if _math.get_length(pos2-pos3) < epsilon or _math.get_length(pos1-pos4) < epsilon:
            vec1 = _math.normalize(pos2-pos1)
            vec2 = _math.normalize(pos3-pos4)
        else:
            vec1 = _math.normalize(pos1-pos2)
            vec2 = _math.normalize(pos3-pos4)
        rad = np.arccos(np.dot(vec1,vec2))
        angle = np.rad2deg(rad)
        angle = str(float('%.2f' % angle)) + "°"   # 保留小数点后三位
        return angle
    

# 结构优化 
############################################################################################
def scaffold_to_mol(scaffold):
    """
    从分子骨架对象获取RDKit可读的mol对象
    """
    from rdkit import Chem
    atoms = []
    bonds = []

    # 读取原子和化学键属性
    Atomic_Nums = get_attr(scaffold, 'atomic_num', 'INT', 'VERT')
    Bond_Orders = get_attr(scaffold, 'bond_order', 'INT', 'EDGE')
    Bond_Orders = [1 if x == 0 else x for x in Bond_Orders]
    add_attr(scaffold, 'bond_order', 'INT', 'EDGE', Bond_Orders)
    bond_type_map = {
        1: Chem.BondType.SINGLE,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE,
        12: Chem.BondType.AROMATIC
    }

    # 记录原子和化学键信息
    for vert, atomic_num in zip(scaffold.data.vertices, Atomic_Nums):
        atoms.append((vert.index, vert.co, atomic_num))
    for edge, bond_order in zip(scaffold.data.edges, Bond_Orders):
        bonds.append((edge.vertices[0], edge.vertices[1], bond_order))

    # 创建RDKit分子并添加原子和化学键
    rdmol = Chem.RWMol()
    atom_map = {}

    for idx, pos, atomic_num in atoms:
        atom = Chem.Atom(atomic_num)
        atom_idx = rdmol.AddAtom(atom)
        atom_map[idx] = atom_idx
    
    conformer = Chem.Conformer(rdmol.GetNumAtoms())
    for idx, pos, _ in atoms:
        # 使用 rdkit.Geometry.Point3D 来设置坐标
        conformer.SetAtomPosition(atom_map[idx], (pos[0], pos[1], pos[2]))
    rdmol.AddConformer(conformer, assignId=True)

    for v1, v2, bond_order in bonds:
        bond_type = bond_type_map.get(bond_order, None)
        if bond_type is None:
            raise ValueError(f"未知的 bond_order: {bond_order}")
        rdmol.AddBond(atom_map[v1], atom_map[v2], bond_type)
    
    # 设置分子名称为标题
    molecule_name = scaffold.name if hasattr(scaffold, 'name') else "Unnamed"
    rdmol.SetProp("_Name", molecule_name)
    
    # 生成最终的Mol对象
    mol = rdmol.GetMol()
    Chem.SanitizeMol(mol)
    return mol, atom_map


def mol_optimize(scaffold, addHs):
    """
    对选中的分子骨架进行构象优化，并更新顶点坐标
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    import numpy as np

    mol, atom_map = scaffold_to_mol(scaffold)
    if addHs: mol = Chem.AddHs(mol)

    # **强制更新原子属性，计算隐含价电子**
    mol.UpdatePropertyCache(strict=False)

    # 生成构象并进行MMFF优化
    status = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    #if status == -1:
    #    raise RuntimeError("构象生成失败")
    AllChem.MMFFOptimizeMolecule(mol)

    # 获取优化后的原子坐标
    conf = mol.GetConformer()
    new_positions = {i: np.array(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())}

    # 更新Blender中的顶点坐标
    for vert in scaffold.data.vertices:
        idx = atom_map[vert.index]
        vert.co = new_positions[idx]


# 常见元素的原子价规则
#####################
valence_dict = {
    1: 1, # H
    6: 4, # C
    7: 3, # N
    8: 2, # O
    9: 1, # F
    15: 3, # P
    16: 2, # S
    17: 1 # Cl
    }
#####################

def add_hydrogens(scaffold):
    """
    给分子骨架结构补全氢原子，确保价态饱和，并考虑空间构象
    """
    # 读取原子和化学键属性
    Atomic_Nums = get_attr(scaffold, 'atomic_num', 'INT', 'VERT')
    Bond_Orders = get_attr(scaffold, 'bond_order', 'INT', 'EDGE')

    bpy.ops.object.mode_set(mode='EDIT')
    bm = bmesh.new()
    bm.from_mesh(scaffold.data)
    
    # 计算每个原子的成键情况
    atom_links = {v.index: {
        "bond_count": 0,
        "bond_order_sum": 0,
        "neighbors": []
    } for v in bm.verts}
    
    for edge in bm.edges:
        v1,v2 = edge.verts
        bond_order = Bond_Orders[edge.index]  # 读取键级
        if bond_order == 12: bond_order = 1.5
        atom_links[v1.index]["bond_count"] += 1
        atom_links[v1.index]["bond_order_sum"] += bond_order
        atom_links[v1.index]["neighbors"].append(v2)

        atom_links[v2.index]["bond_count"] += 1
        atom_links[v2.index]["bond_order_sum"] += bond_order
        atom_links[v2.index]["neighbors"].append(v1)
    
    # 计算每个顶点处新添加的氢原子位置
    new_verts = []
    for v in bm.verts:
        atomic_num = Atomic_Nums[v.index]
        if atomic_num not in valence_dict:
            continue   # 跳过不需要补氢的原子
        max_valence = valence_dict[atomic_num]
        current_valence = atom_links[v.index]["bond_order_sum"]
        num_hydrogens = max_valence - current_valence
        if num_hydrogens <= 0:
            continue  # 价键已满，无需补氢

        bond_count = atom_links[v.index]["bond_count"]
        hybridization = 3 + bond_count - atom_links[v.index]["bond_order_sum"]   # sp: 1, sp2: 2, sp3: 3
        new_H_pos = new_hydrogens(v, bond_count, hybridization, Atomic_Nums)
        new_verts.append((v.index, new_H_pos))

    # 添加氢原子顶点和化学键
    bm = bmesh.from_edit_mesh(scaffold.data)
    sum_H = 0
    for new_vert in new_verts:
        index = new_vert[0]
        new_H_pos = new_vert[1]
        for pos in new_H_pos:
            sum_H += 1
            new_vert = bm.verts.new(pos)
            bm.verts.ensure_lookup_table()
            new_edge = bm.edges.new([new_vert, bm.verts[index]])
    bmesh.update_edit_mesh(scaffold.data)
    
    # 设置新的原子和键属性
    Atomic_Nums = get_attr(scaffold, 'atomic_num', 'INT', 'VERT')
    Bond_Orders = get_attr(scaffold, 'bond_order', 'INT', 'EDGE')
    VDW_R = get_attr(scaffold, 'vdw_radius', 'FLOAT', 'VERT')
    COVA_R = get_attr(scaffold, 'covalent_radius', 'FLOAT', 'VERT')
    Atomic_Nums = [1 if x == 0 else x for x in Atomic_Nums]
    Bond_Orders = [1 if x == 0 else x for x in Bond_Orders]
    VDW_R = [1.2 if x == 0 else x for x in VDW_R]
    COVA_R = [0.31 if x == 0 else x for x in COVA_R]
    add_attr(scaffold, 'atomic_num', 'INT', 'VERT', Atomic_Nums)
    add_attr(scaffold, 'vdw_radius', 'FLOAT', 'VERT', VDW_R)
    add_attr(scaffold, 'covalent_radius', 'FLOAT', 'VERT', COVA_R)
    add_attr(scaffold, 'bond_order', 'INT', 'EDGE', Bond_Orders)
    bpy.ops.object.mode_set(mode='OBJECT')


def new_hydrogens(vert, bond_count, hybridization, Atomic_Nums):
    """
    根据已有的键连接和杂化类型，计算剩余的成键方向
    """
    new_directions = []
    new_H_pos = []
    deflection_angle = 0 if bond_count == hybridization else 54.75
    branches = valence_dict[Atomic_Nums[vert.index]]
    if bond_count == 0:  # 单独的点
        if branches == 1:
            new_directions.append((1,0,0))
        elif branches == 2:
            new_directions.append((0.5, 0.866, 0.0))
            new_directions.append((0.5, -0.866, 0.0))
        elif branches == 3:
            new_directions.append((0.5, 0.866, 0.0))
            new_directions.append((0.5, -0.866, 0.0))
            new_directions.append((-1.0, 0.0, 0.0))
        else:
            new_directions.append((0.57735, 0.57735, 0.57735))
            new_directions.append((0.57735, -0.57735, 0.57735))
            new_directions.append((-0.57735, 0.57735, 0.57735))
            new_directions.append((-0.57735, -0.57735, 0.57735))
    else:
        new_directions = _math.branches_dir(vert, hybridization-bond_count+branches-3, deflection_angle)[0]

    for direction in new_directions:
        new_H_pos.append(np.array(vert.co) + 1.09*np.array(direction))

    return new_H_pos
    



