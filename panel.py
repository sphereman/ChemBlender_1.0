import bpy
from bpy.props import StringProperty, EnumProperty
language = 1 if 'zh_HAN' in bpy.context.preferences.view.language else 0

# ------------------------------------------------------------------------------------
# This class lists properties shown in Panel but not in Operator.
class CHEM_texts(bpy.types.PropertyGroup):
    filetext: StringProperty(
        name = "Filepath",
        default = "选择分子文件路径" if language else 'Choose MolFile',
        subtype = 'FILE_PATH'
    ) # type: ignore

    moltext: StringProperty(
        name = "文本输入" if language else "Text Input",
        default = "请输入分子SMILES符或ID号" if language else 'Input Your Text',
        subtype = 'NONE'
    ) # type: ignore

    select_text: StringProperty(
        name = '',
        default = "输入文本" if language else "Input Text Here",
        description = "Text of Element Name for Scaling"
    ) # type:ignore

    choose: EnumProperty(
        name = "选项列表" if language else "Option List",
        default = "File",
        items = [("File","从文件创建",""),
                 ("PubChem","从PubChem创建",""),
                 ("SMILES","从SMILES创建",""),
                 ("Saccharides","糖类",""),
                 ("Amino_Acids","氨基酸",""),
                 ("Polymer_Units", "聚合物单元",""),
        ]
    ) # type: ignore

    Saccharides: EnumProperty(
        name = "糖类" if language else "Saccharides",
        default = "Glc",
        items = [("Glc","葡萄糖",""),
                 ("Fru","果糖",""),
                 ("Gal","半乳糖",""),
                 ("Man","甘露糖",""),
                 ("Xyl","木糖",""),
                 ("Ara","阿拉伯糖",""),
                 ("Rib","核糖",""),
                 ("dRib","脱氧核糖",""),
                 ("Suc","蔗糖",""),
                 ("Lac","乳糖",""),
                 ("Mal","麦芽糖",""),
                 ("I-Mal","异麦芽糖",""),
                 ("Tre","海藻糖",""),
                 ("Cel","纤维二糖",""),
                 ("Mel","蜂蜜三糖",""),
                 ("Raf","棉子糖",""),
                 ("Amy","直链淀粉单元",""),
                 ("Amp","支链淀粉单元",""),
                 ("Cel","纤维素单元",""),
                 ("Cht","几丁质单元",""),
                 ("Gly","糖原单元",""),
        ]
    ) # type: ignore

    Amino_Acids: EnumProperty(
        name = "氨基酸" if language else "Amino acids",
        default = "A",
        items = [("A","丙氨酸","Ala"),
                 ("R","精氨酸","Arg"),
                 ("D","天冬氨酸","Asp"),
                 ("N","天冬酰胺","Asn"),
                 ("C","半胱氨酸","Cys"),
                 ("E","谷氨酸","Glu"),
                 ("Q","谷氨酰胺","Gln"),
                 ("G","甘氨酸","Gly"),
                 ("H","组氨酸","His"),
                 ("I","异亮氨酸","Ile"),
                 ("L","亮氨酸","Leu"),
                 ("K","赖氨酸","Lys"),
                 ("M","甲硫氨酸","Met"),
                 ("F","苯丙氨酸","Phe"),
                 ("P","脯氨酸","Pro"),
                 ("S","丝氨酸","Ser"),
                 ("T","苏氨酸","Thr"),
                 ("W","色氨酸","Trp"),
                 ("Y","酪氨酸","Tyr"),
                 ("V","缬氨酸","Val"),
        ]
    ) # type: ignore

    Polymer_Units: EnumProperty(
        name = "聚合物单元" if language else "Polymer Units",
        default = "PE",
        items = [("PE","聚乙烯单元","Polyethylene"),
                 ("PP","聚丙烯单元","Polypropylene"),
                 ("PS","聚苯乙烯单元","Polystyrene"),
                 ("PVC","聚氯乙烯单元","Polyvinyl chloride"),
                 ("PTFE","聚四氟乙烯单元","Teflon"),
                 ("PVDC","聚偏二氯乙烯单元","Polyvinylidene chloride"),
                 ("PVA","聚乙烯醇单元","Polyvinyl alcohol"),
                 ("PEG","聚乙二醇单元","Polyethylene glycol"),
                 ("PMMA","聚甲基丙烯酸甲酯单元","Polymethyl methacrylate"),
                 ("PAN","聚丙烯腈单元","Polyacrylonitrile"),
                 ("PB","聚丁二烯单元","Polybutadiene"),
                 ("PVAC","聚乙酸乙烯酯单元","Polyvinyl acetate"),
                 ("PLA","聚乳酸单元","Polylactic acid"),
                 ("PET","聚对苯二甲酸乙二醇酯单元","Polyethylene terephthalate"),
                 ("PBT","聚对苯二甲酸丁二醇酯单元","Polybutylene terephthalate"),
                 ("Kevlar","聚对苯二甲酰对苯二胺单元","Poly-p-phenylene terephthamide"),
                 ("PC","聚碳酸酯单元","Polycarbonate"),
                 ("PEEK","聚醚醚酮单元","Polyether ether ketone"),
                 ("PA","聚己内酰胺单元","Polyamide"),
                 ("PI","聚酰亚胺单元","Polyimide"),
                 ("PAA","聚丙烯酸单元","Polyacrylic Acid"),
                 ("PAAm","聚丙烯酰胺单元","Polyacrylamide"),
                 ("PVP","聚乙烯吡咯烷酮单元","Polyvinylpyrrolidone"),
                 ("PDMS","聚二甲基硅氧烷单元","Polydimethylsiloxane"),

        ]
    ) # type: ignore

    distance: StringProperty(
        name = '',
        default = "键长 (Å)" if language else "Distance (Å)",
        description = "Length of an edge in Ångstrom"
    ) # type: ignore

    angle: StringProperty(
        name = '',
        default = "键角 (°)" if language else "Angle (°)",
        description = "Angle of two edges in degree"
    ) # type: ignore

class CHEM_PT_Build(bpy.types.Panel):
    bl_label = "分子骨架创建" if language else "Scaffold Building"
    bl_idname = "CHEM_PT_BUILD"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'ChemBlender'

    def draw(self, context):
        mytool = context.scene.my_tool
        layout = self.layout
        layout.prop(mytool, "choose")
        if mytool.choose == 'Saccharides':
            layout.prop(mytool, "Saccharides")
        elif mytool.choose == 'Amino_Acids':
            layout.prop(mytool, "Amino_Acids")
        elif mytool.choose == 'Polymer_Units':
            layout.prop(mytool, "Polymer_Units")
        elif mytool.choose == 'File':
            layout.prop(mytool, "filetext")
        else:
            layout.prop(mytool, "moltext")

        layout.row()
        text = "创建分子" if language else "Build Molecule"
        layout.operator("chem.scaffold_build", text=text, icon="GREASEPENCIL")

class CHEM_PT_TOOLS(bpy.types.Panel):
    bl_label = "分子工具" if language else "ChenBlender Tools"
    bl_idname = "CHEM_PT_UTILITY"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'ChemBlender'
    bl_options= {'DEFAULT_CLOSED'}

    def draw(self, context):
        layout = self.layout
        mytool = context.scene.my_tool

        row = layout.row(align=True)
        row.prop(mytool, "select_text")
        row.scale_x = 0.75
        text = "选择" if language else "Select"
        row.operator('mol3d.select', text = text)
        row = layout.row(align=True)
        text = "测量长度" if language else "Length Calc."
        row.operator("chem.button_distance", text=text)
        row.prop(mytool,"distance")
        row = layout.row(align=True)
        text = "测量角度" if language else "Angle Calc."
        row.operator("chem.button_angle", text=text)
        row.prop(mytool,"angle")

        row = layout.row()
        text = "设置原子" if language else "Set Atoms"
        row.operator("chem.set_atoms", text=text)
        text = "设置键" if language else "Set Bonds"
        row.operator("chem.set_bonds", text=text)

        row = layout.row()
        text = "补全氢原子" if language else "Add Hydrogens"
        row.operator("chem.add_hydrogens", text=text)
        text = "构象更新" if language else "Struct Optimize"
        row.operator("chem.mmff_update", text=text)
        row = layout.row()
        text = "选中对象导出为分子文件" if language else "Export Scaffold to Block File"
        row.operator("chem.export_block", text=text)
        row = layout.row()
        text = "选中网格转换为分子骨架" if language else "Convert Mesh to Mol Scaffold"
        row.operator("chem.mesh2scaffold", text=text)

