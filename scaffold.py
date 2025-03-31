import bpy
import os
from .read import read_MOL, download_sdf_from_pubchem
from .mesh import create_object, add_scaffold_attr
from .Chem_data import preset_smiles
from . import node
language = 1 if 'zh_HAN' in bpy.context.preferences.view.language else 0

# ------------------------------------------------------------------------------------
class ErrorDialogOperator(bpy.types.Operator):
    bl_idname = "wm.error_dialog"
    bl_label = "错误提示"

    message: bpy.props.StringProperty(default="发生了未知错误") # type: ignore

    def execute(self, context):
        self.report({'INFO'}, "关闭错误提示")
        return {'FINISHED'}

    def invoke(self, context, event):
        return context.window_manager.invoke_props_dialog(self)

    def draw(self, context):
        layout = self.layout
        layout.label(text=self.message)


class MESH_OT_SCAFFOLD_BUILD(bpy.types.Operator):
    bl_idname = "chem.scaffold_build"
    bl_label = "创建骨架" if language else "Build Scaffold"
    bl_desciption = "生成分子骨架的球棍模型" 
    bl_options = {'REGISTER','UNDO'}
    
    def text_input(self, mytool):
        if mytool.choose == 'File':
            moltext = mytool.filetext
        elif mytool.choose == 'SMILES' or mytool.choose == "PubChem":
            moltext = mytool.moltext
        elif mytool.choose == "Saccharides":
            moltext = preset_smiles[mytool.Saccharides][1]
        elif mytool.choose == "Amino_Acids":
            moltext = preset_smiles[mytool.Amino_Acids][1]
        elif mytool.choose == "Polymer_Units":
            moltext = preset_smiles[mytool.Polymer_Units][1]
        return moltext

    def name_input(self, mytool, moltext):
        if os.path.exists(moltext):
            molname = os.path.basename(moltext).split('.')[0]
        elif mytool.choose == "Saccharides":
            molname = preset_smiles[mytool.Saccharides][0]
        elif mytool.choose == "Amino_Acids":
            molname = preset_smiles[mytool.Amino_Acids][0]
        elif mytool.choose == "Polymer_Units":
            molname = preset_smiles[mytool.Polymer_Units][0]
        elif mytool.choose == "PubChem":
            molname = download_sdf_from_pubchem(moltext)[1]
        else:
            molname = moltext
        return molname

    def execute(self, context):
        mytool = context.scene.my_tool
        try:
            moltext = self.text_input(mytool)
            molname = self.name_input(mytool, moltext)
            ATOMS, AtomicNum, COORDS, BONDS, BOND_ORDERS, VDW_R, COVA_R, RingNum = read_MOL(moltext)
            coll = bpy.data.collections.new('Scaffold_'+molname)
            bpy.context.scene.collection.children.link(coll)
            mol_scaffold = create_object(coll, molname, COORDS, BONDS, [])
            mol_scaffold['Type'] = 'scaffold'
            mol_scaffold['Elements'] = list(set(ATOMS))
            add_scaffold_attr(mol_scaffold, ATOMS, AtomicNum, BOND_ORDERS, VDW_R, COVA_R, RingNum)
            bpy.context.view_layer.objects.active = mol_scaffold
            mol_scaffold.select_set(True)
            GN_mol = node.add_geometry_nodetree(mol_scaffold, "GN_"+molname, "NodeTree_"+molname)
            node.Ball_Stick_nodetree(GN_mol)
            return {'FINISHED'}
        except Exception as e:
            print(e)
            text = f"操作失败: 请检查输入内容或第三方Python库配置。{e}" if language else f"Operation failed: Please check the input text or python site-expackages. {e}"
            bpy.ops.wm.error_dialog('INVOKE_DEFAULT', message=text)
            return {'CANCELLED'}
