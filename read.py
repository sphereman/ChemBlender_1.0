import bpy
import os
import numpy as np
import requests


def add_bonds_based_on_distance(COORDS, COVA_R, BONDS, BOND_ORDERS):
    for i in range(len(COORDS)):
        for j in range(i+1, len(COORDS)):
            pos_i = np.array(COORDS[i])
            pos_j = np.array(COORDS[j])
            dist = np.linalg.norm(pos_i-pos_j)
            if dist < (COVA_R[i]+COVA_R[j]+0.4):   # 0.4 Å as tolerance distance
                BONDS.append((i,j))
                BOND_ORDERS.append(1)  # 键级默认预设为1

def attr_values_from_mol(mol):
    from rdkit import Chem
    periodic_table = Chem.GetPeriodicTable()
    from rdkit.Chem import AllChem
    """
    给原子和键设置基本属性的列表值.
    """
    ATOMS, AtomicNum, COORDS, BONDS, BOND_ORDERS = [],[],[],[],[]
    VDW_R, COVA_R, RingNum = [],[],[]
    atoms = mol.GetAtoms()
    bonds = mol.GetBonds()
    ring_info = mol.GetRingInfo()
    if not mol.GetNumConformers():
        AllChem.EmbedMolecule(mol)   # Generate a conformer
    conformer = mol.GetConformer()   # Get the first conformer
    cursor = bpy.context.scene.cursor.location
    for atom in atoms:
        ATOMS.append(atom.GetSymbol())
        pos = conformer.GetAtomPosition(atom.GetIdx())
        COORDS.append((pos.x+cursor.x, pos.y+cursor.y, pos.z+cursor.z))
        atomic_num = atom.GetAtomicNum()
        AtomicNum.append(atomic_num)
        VDW_R.append(periodic_table.GetRvdw(atomic_num))
        COVA_R.append(periodic_table.GetRcovalent(atomic_num))
    for i, bond in enumerate(bonds):
        BONDS.append((bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()))
        BOND_ORDERS.append(bond.GetBondType())
        RingNum.append(find_in_list(i, ring_info.BondRings()))
    # print(RingNum)
    if not BONDS:
        add_bonds_based_on_distance(COORDS, COVA_R, BONDS, BOND_ORDERS)
    return ATOMS, AtomicNum, COORDS, BONDS, BOND_ORDERS, VDW_R, COVA_R, RingNum

def find_in_list(num, lst):
    for index, sublist in enumerate(lst, start=1):
        if num in sublist:
            return index
    return 0

def mol_2D_to_3D(mol):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.AddHs(mol)  # 添加氢原子
    AllChem.EmbedMolecule(mol)  # 生成3D构象
    AllChem.UFFOptimizeMolecule(mol)  # 力场优化
    return mol

def is_valid_smiles(smiles_string):
    from rdkit import Chem
    """
    判断给定字符串是否是有效的SMILES字符串.
    """
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if mol:
            print(f"输入内容为SMILES字符串.")
            return True
        else:
            return False
    except Exception as e:
        print(f"解析SMILES时发生错误: {e}")
        return False

def check_type(moltext):
    """
    判断输入字符串文本的类型: XYZ/MOL/SDF/PDB/SMILES/CIF/CML等.
    根据文本格式或后缀名进行判断.
    """
    if os.path.exists(moltext):
        text_type = os.path.basename(moltext).split('.')[-1]
    elif is_valid_smiles(moltext):
        text_type = 'smiles'
    elif moltext.isdigit():
        text_type = 'cid'
    else:
        text_type = None
    return text_type

def download_sdf_from_pubchem(cid):
    """
    从PubChem网站根据CID号载入对应分子的SDF格式文件, 如Aspirin CID: 2244.
    此处下载均为2D分子, 需经过优化得到3D坐标.
    """
    try:
        # 获取SDF文件内容
        sdf_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF"
        response_1 = requests.get(sdf_url)
        if response_1.status_code != 200:
            print(f"Failed to download SDF file for CID {cid}.")
            return None, None, None

        # 获取分子名称和IUPAC名称
        name_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IUPACName,Title/JSON"
        response_2 = requests.get(name_url)
        if response_2.status_code != 200:
            print(f"Failed to get molecular names for CID {cid}.")
            return None, None, None
        
        data = response_2.json()
        properties = data.get("PropertyTable", {}).get("Properties", [{}])[0]
        common_name = properties.get("Title", "N/A")
        iupac_name = properties.get("IUPACName", "N/A")

        # 返回SDF内容和分子名称
        sdf_content = response_1.text
        return sdf_content, common_name, iupac_name
    
    except requests.exceptions.RequestException as e:
        print(f"An error occurred: {e}")
        return None, None, None

def read_MOL(moltext):
    from rdkit import Chem
    text_type = check_type(moltext)
    # get mol object from moltext
    if text_type == 'smiles':
        mol = Chem.MolFromSmiles(moltext)
        mol = mol_2D_to_3D(mol)
    elif text_type == 'xyz':
        mol = Chem.MolFromXYZFile(moltext)
    elif text_type == 'mol' or text_type == 'sdf':
        mol = Chem.MolFromMolFile(moltext, removeHs=False)
    elif text_type == 'pdb':
        mol = Chem.MolFromPDBFile(moltext)
    elif text_type == 'cid':   # from PubChem
        block = download_sdf_from_pubchem(moltext)[0]
        mol = Chem.MolFromMolBlock(block)
        mol = mol_2D_to_3D(mol)   # conformation optimization
    return attr_values_from_mol(mol)