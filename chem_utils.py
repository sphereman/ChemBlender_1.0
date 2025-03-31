import bpy
from . import mesh, node
from bpy.types import Operator
from bpy.props import IntProperty, BoolProperty, StringProperty
language = 1 if 'zh_HAN' in bpy.context.preferences.view.language else 0

class SelectButton(Operator):
    bl_idname = "mol3d.select"
    bl_label = "选择" if language else "Select"
    bl_description = "根据输入内容选择点或边"
    bl_options = {'REGISTER','UNDO'}

    def simplify_text(self, text):
        for sep in ',;/| ': text = text.replace(sep, ' ')
        elements = text.split()
        elements = [e.capitalize() if e.isalpha() else e for e in elements]
        elements = ['bond' if e == 'Bond' else e for e in elements]
        if 'all' in text.lower() or '*' in text:
            elements += ['atom', 'bond']
        if 'atom' in text.lower() or 'ball' in text.lower():
            elements += ['atom']
        if 'bond' in text.lower() or 'ball' in text.lower():
            elements += ['bond']
        # e.g. ['C','H']
        bond_syms = [[x.split('-')[0].capitalize(),x.split('-')[1].capitalize()] for x in elements if '-' in x]
        return list(set(elements)), bond_syms
    
    def execute(self, context):
        mytool  = context.scene.my_tool
        text = mytool.select_text
        text, bond_syms = self.simplify_text(text)
        ao = bpy.context.object
        if not ao:
            return {'CANCELLED'}
        else:
            mesh.deselect_all(ao)
            mesh.select_verts(ao, text)
            mesh.select_edges(ao, bond_syms)
            mesh.select_bond_orders(ao, text)
            if 'atom' in text: mesh.select_all(ao, 'VERT')
            if 'bond' in text: mesh.select_all(ao, 'EDGE')
            if 'all' in text or '*' in text: mesh.select_all(ao, 'VERT')
            return {'FINISHED'}
        
class DistanceButton(Operator):
    bl_idname = "chem.button_distance"
    bl_label = "Calc Distnace"
    bl_description = "选择一条边计算长度"

    def execute(self, context):
        mytool = context.scene.my_tool
        ao = bpy.context.object   # active object
        if not ao or bpy.context.mode == 'OBJECT':
            return {'CANCELLED'}
        else:
            distance = mesh.calc_length()
            # put the distance value into the string of the output field.
            mytool.distance = distance
            return {'FINISHED'}

class AngleButton(Operator):
    bl_idname = "chem.button_angle"
    bl_label = "Calc Angle"
    bl_description = "选择两条边计算夹角"

    def execute(self, context):
        mytool = context.scene.my_tool
        ao = bpy.context.object
        if not ao or bpy.context.mode == 'OBJECT':
            return {'CANCELLED'}
        else:
            angle = mesh.calc_angle()
            mytool.angle = angle
            return {'FINISHED'}

class SetAtomsButton(Operator):
    bl_idname = 'chem.set_atoms'
    bl_label = "设置原子" if language else "Atom Setting"
    bl_description = "设置选中点的原子属性"
    bl_options = {'REGISTER', 'UNDO'}

    text = "原子序数" if language else "Atomic Number"
    atomic_num: IntProperty(name=text, default=0, min=0, max=118) # type: ignore

    def execute(self, context):
        ao = bpy.context.object
        try:
            mesh.set_sel_atoms_attr(ao, self.atomic_num)
            bpy.ops.object.mode_set(mode='EDIT')
            return {'FINISHED'}
        except Exception as e:
            text = f"请选择正确的对象。" if language else f"Please select right object."
            bpy.ops.wm.error_dialog('INVOKE_DEFAULT', message=text)
            return {'CANCELLED'}
          
class SetBondsButton(Operator):
    bl_idname = 'chem.set_bonds'
    bl_label = "设置键" if language else "Bond Setting"
    bl_description = "设置选中边的键属性"
    bl_options = {'REGISTER', 'UNDO'}

    text = "键级" if language else "Bond Order"
    bond_order: IntProperty(name=text, default=1, min=0, soft_max=3) # type: ignore
    text = "环编号" if language else "Ring Number"
    ring_num: IntProperty(name=text, default=0, min=0, soft_max=10) # type: ignore

    def execute(self, context):
        ao = bpy.context.object
        try:
            mesh.set_sel_bonds_attr(ao, self.bond_order, self.ring_num)
            bpy.ops.object.mode_set(mode='EDIT')
            return {'FINISHED'}
        except Exception as e:
            text = f"请选择正确的对象。" if language else f"Please select right object."
            bpy.ops.wm.error_dialog('INVOKE_DEFAULT', message=text)
            return {'CANCELLED'}

class AddHydrogens(Operator):
    bl_idname = 'chem.add_hydrogens'
    bl_label = 'Add Hydrogens'
    bl_description = "根据饱和度和空间构型补齐氢原子"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        # 获取当前选中的物体
        ao = context.object
        
        # 检查是否选中了有效的物体
        if ao is None or ao.type != 'MESH':
            self.report({'WARNING'}, "请选择一个有效的网格物体！")
            return {'CANCELLED'}
        try:
            # 调用添加氢原子的函数
            mesh.add_hydrogens(ao)
            return {'FINISHED'}
        except Exception as e:
            # 打印异常信息，便于调试
            self.report({'ERROR'}, f"发生错误: {str(e)}")
            return {'CANCELLED'}
            
        
class MMFFUpdateButton(Operator):
    bl_idname = 'chem.mmff_update'
    bl_label = 'Conformation Update'
    bl_description = "使用MMFF力场优化并重新生成分子构象"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        # 获取当前选中的物体
        ao = context.object
        
        # 检查是否选中了有效的网格物体
        if ao is None or ao.type != 'MESH':
            self.report({'WARNING'}, "请选择一个有效的网格物体！")
            return {'CANCELLED'}

        try:
            # 调用mol_optimize函数进行构象更新
            mesh.mol_optimize(ao, addHs=True)
            self.report({'INFO'}, "构象更新成功！")
            return {'FINISHED'}
        except Exception as e:
            # 打印异常信息，便于调试
            self.report({'ERROR'}, f"发生错误: {str(e)}")
            return {'CANCELLED'}
        
class SaveBlockButton(Operator):
    """Convert molecular structure to block file"""
    bl_idname = "chem.export_block"
    bl_label = "保存文件" if language else "Export Block File"
    bl_description = "从分子骨架生成mol文件"
    bl_options = {'REGISTER', 'UNDO'}

    filepath: StringProperty(subtype="FILE_PATH")  # type: ignore 让用户选择保存路径 

    def execute(self, context):
        from rdkit import Chem
        from rdkit.Chem import SDWriter
        ao = context.object
        if ao is None or ao.type != 'MESH':
            self.report({'WARNING'}, "请选择一个有效的网格物体！")
            return {'CANCELLED'}

        try:
            mol = mesh.scaffold_to_mol(ao)[0]
            block = Chem.MolToMolBlock(mol)

            if not self.filepath.lower().endswith('.mol'):
                self.filepath += '.mol'
                
            with open(self.filepath, "w") as f:
                f.write(block)
            self.report({'INFO'}, f"分子已保存为: {self.filepath}")
            return {'FINISHED'}
        except Exception as e:
            self.report({'ERROR'}, f"发生错误: {str(e)}")
            return {'CANCELLED'}
        
    def invoke(self, context, event):
        # 调用文件选择器对话框
        context.window_manager.fileselect_add(self)  # 打开文件保存路径设置窗口
        return {'RUNNING_MODAL'}

class ScaffoldConvertButton(Operator):
    """Convert mesh to molecular scaffold"""
    bl_idname = "chem.mesh2scaffold"
    bl_label = "网格到分子骨架" if language else "Mesh to Mol Scaffold"
    bl_description = "从网格生成分子骨架"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        ao = context.object
        if ao is None or ao.type != 'MESH':
            self.report({'WARNING'}, "请选择一个有效的网格物体！")
            return {'CANCELLED'}

        obj_type = ao.get('Type', None)  # 避免 KeyError
        if obj_type == 'scaffold':
            self.report({'WARNING'}, "该对象已是分子骨架。")
            return {'CANCELLED'}
        
        try:
            mesh.mesh_to_mol_scaffold(ao)
            self.report({'INFO'}, "分子骨架生成成功！")
            GN_mol = node.add_geometry_nodetree(ao, "GN_"+ao.name, "NodeTree_"+ao.name)
            node.Ball_Stick_nodetree(GN_mol)
            return {'FINISHED'}
        except Exception as e:
            self.report({'ERROR'}, f"发生错误: {str(e)}")
            return {'CANCELLED'}