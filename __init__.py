bl_info = {
    "name" : "ChemBlender",
    "author" : "LiHaodong",
    "version" : (1, 0, 0),
    "blender" : (4, 2, 1),
    "location" : "Geometry Nodes Editor > Add > Chem",
    "description" : "Professional Visualization of Molecules and Structures for Scientists and Artiests.",
    "warning" : "",
    "doc_url": "asrc.ustc.edu.cn",
    "category": "Node",
}

import bpy
import os, re, json
from . import auto_load
from .panel import CHEM_texts, CHEM_PT_Build, CHEM_PT_TOOLS
from bpy.types import (
    Operator,
    Menu,
)
from bpy.props import (
    StringProperty,
)


# Get the plugin directory path
dir_path = os.path.dirname(__file__)

geo_node_group = {}
cat_list = []

# 生成分类菜单
def cat_generator():
    global cat_list
    cat_list = []  # 每次调用时重置cat_list
    for item in geo_node_group.items():
        def my_list(self, context):
            layout = self.layout
            for name_group in geo_node_group[self.bl_label]:
                props = layout.operator(
                    NODE_OT_group_add.bl_idname,
                    text=re.sub(r'.*?_', '', name_group), # 去除前缀
                )
                props.name_group = name_group

        menu_type = type("NODE_MT_group_" + item[0], (bpy.types.Menu,), {
            "bl_idname": "NODE_MT_group_" + item[0].replace(" ", "_"),   # Replace spaces with underscores to avoid alpha-numeric suffic warning
            "bl_space_type": 'NODE_EDITOR',
            "bl_label": item[0],
            "draw": my_list,
        })
        if menu_type not in cat_list:
            def generate_menu_draw(name, label): # 强制唯一引用
                def draw_menu(self, context):
                    self.layout.menu(name, text=label)
                return draw_menu
            bpy.utils.register_class(menu_type)
            bpy.types.NODE_MT_chem_GN_menu.append(generate_menu_draw(menu_type.bl_idname, menu_type.bl_label))
            cat_list.append(menu_type)


class NODE_MT_chem_GN_menu(Menu):
    bl_label = "Chem Nodes"
    bl_idname = 'NODE_MT_chem_GN_menu'

    @classmethod
    def poll(cls, context):
        return context.space_data.tree_type == 'GeometryNodeTree'

    def draw(self, context):
        pass


class NODE_OT_group_add(Operator):
    """Add node group"""
    bl_idname = "chem.add_node_group"
    bl_label = "Add node group"
    bl_options = {'REGISTER', 'UNDO'}

    name_group: StringProperty()

    @classmethod
    def poll(cls, context):
        return context.space_data.node_tree is not None  # Ensure there is a node tree in the context

    def execute(self, context):
        old_groups = set(bpy.data.node_groups)

        # 查找 .blend 文件
        for file in os.listdir(dir_path):
            if file.endswith(".blend"):
                filepath = os.path.join(dir_path, file)
                break
        else:
            raise FileNotFoundError(f"未在目录 {dir_path} 中找到 .blend 文件")

        # 加载 .blend 文件中的节点组
        with bpy.data.libraries.load(filepath, link=False) as (data_from, data_to):
            if self.name_group in data_from.node_groups:
                data_to.node_groups.append(self.name_group)
        
        # 确保新导入的节点组在 bpy.data.node_groups 中可用
        new_groups = list(set(bpy.data.node_groups) - old_groups)
        for group in new_groups:
            for node in group.nodes:
                if node.type == "GEOMETRY":
                    name_new = node.node_tree.name.split(".")[0]
                    node.node_tree = bpy.data.node_groups[name_new]
        
        # 删除无效副本节点组
        for group in new_groups:
            if "." in group.name:
                bpy.data.node_groups.remove(group)

        # 添加新的几何节点组到当前节点树
        bpy.ops.node.add_node(type="GeometryNodeGroup")
        node = context.selected_nodes[0]
        node.node_tree = bpy.data.node_groups[self.name_group]
        bpy.ops.transform.translate('INVOKE_DEFAULT')

        return {'FINISHED'}

# Register the plugin's menu
def add_chem_button(self, context):
    if context.area.ui_type == 'GeometryNodeTree':
        self.layout.menu('NODE_MT_chem_GN_menu', text="Chem Nodes", icon='FORWARD')

classes = (
    NODE_OT_group_add,
    # NODE_MT_chem_GN_menu,
)
auto_cls = auto_load.init()
panel_cls = [CHEM_PT_Build,
             CHEM_PT_TOOLS,
             ]
for cls in panel_cls:
    auto_cls.remove(cls)
auto_cls.append(NODE_OT_group_add)

def register():
    global geo_node_group

    # Load GN_menu.json configuration file
    with open(os.path.join(dir_path, "GN_menu.json"), 'r') as f:
        geo_node_group = json.loads(f.read())
    
    if not hasattr(bpy.types, "NODE_MT_chem_GN_menu"):
        bpy.utils.register_class(NODE_MT_chem_GN_menu)
        bpy.types.NODE_MT_add.append(add_chem_button)

    for cls in panel_cls:
        bpy.utils.register_class(cls)
    for cls in auto_cls:
        bpy.utils.register_class(cls)
    bpy.types.Scene.my_tool = bpy.props.PointerProperty(type=CHEM_texts)
    cat_generator()


def unregister():
    if hasattr(bpy.types, "NODE_MT_chem_GN_menu"):
        bpy.types.NODE_MT_add.remove(add_chem_button)
    for cls in panel_cls:
        bpy.utils.unregister_class(cls)
    for cls in auto_cls:
        bpy.utils.unregister_class(cls)
    del bpy.types.Scene.my_tool
    
#if __name__ == "__main__":
#    register()