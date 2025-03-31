import  bpy
import os
# Get the plugin directory path
dir_path = os.path.dirname(__file__)
file = os.path.join(dir_path, "Chem_Nodes.blend")

def add_geometry_nodetree(obj, GN_modifier_name, nodetree_name):
    bpy.context.view_layer.objects.active = obj
    GN_modifier = obj.modifiers.new(GN_modifier_name, 'NODES')
    bpy.ops.node.new_geometry_node_group_assign()
    GN_modifier.node_group.name = nodetree_name
    return GN_modifier

def set_io_nodes(nodetree, in_location, out_location):
    group = nodetree.node_group
    try:
        _input = group.nodes['Group Input']
    except:
        _input = group.nodes['组输入']
    _input.location = in_location
    try:
        _output = group.nodes['Group Output']
    except:
        _output = group.nodes['组输出']
    _output.location = out_location
    _input.label = 'in_node'
    _output.label = 'out_node'
    return (_input, _output)

# add a common geometry node in given nodetree
def add_node(nodetree, name, location, label):
    node = nodetree.node_group.nodes.new(name)
    node.location = location
    node.label = label
    return node

def add_node_group(nodetree, name, location): 
    new_group = nodetree.node_group.nodes.new(type="GeometryNodeGroup")
    new_group.node_tree = bpy.data.node_groups[name]
    new_group.location = location
    return new_group

# create links between two nodes or nodegroups
def nodes_link(nodetree, node_a, socket_out, node_b, socket_in):
    link = nodetree.node_group.links.new
    link(node_a.outputs[socket_out], node_b.inputs[socket_in])

def append(node_group_name, link=False):
    node = bpy.data.node_groups.get(node_group_name)
    if not node or link:
        bpy.ops.wm.append(
            'EXEC_DEFAULT',
            directory = os.path.join(file, 'NodeTree'),
            filename = node_group_name,
            link = link
        )

def Ball_Stick_nodetree(nodetree):
    a,b = (0,0)
    _input, _output = set_io_nodes(nodetree, (a,b), (a+800,b))
    print(dir_path)
    for node_group_name in ["CH_添加分子属性","CH_分子球棍模型","CH_添加分子材质"]:
        append(node_group_name)
    _add_mol_attr = add_node_group(nodetree, "CH_添加分子属性", (a+200,0))
    _ball_stick = add_node_group(nodetree, "CH_分子球棍模型", (a+400,0))
    _ball_stick.inputs[5].default_value = 0.5
    _add_material = add_node_group(nodetree, "CH_添加分子材质", (a+600,0))
    nodes_link(nodetree, _input, 0, _add_mol_attr, 0)
    nodes_link(nodetree, _add_mol_attr, 0, _ball_stick, 0)
    nodes_link(nodetree, _ball_stick, 0, _add_material, 0)
    nodes_link(nodetree, _add_material, 0, _output, 0)