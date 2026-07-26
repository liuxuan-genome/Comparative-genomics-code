#!usr/bin/env python3
# conding: utf-8
# coding:utf-8
'''
############################################################################
#北京组学生物科技有限公司
#author huangls
#date 2023.03.01
#version 1.0
#学习python课程推荐：
#python 入门到精通(生物信息):
#https://zzw.xet.tech/s/2bdd89
################################################################################
'''


import sys, os, argparse, os.path

##################################################

parser = argparse.ArgumentParser(description='This script was used to plot cafe5 result')
parser.add_argument('--nwk', dest='nwk', required=True, help='input tree file with nwk format ')
parser.add_argument('--clade_results', dest='clade_results', required=True, help='clade_results  file from cafe5  ')
parser.add_argument('--gfn', dest='gfn', required=True, type=int,help='total gene family number analysised  ')
parser.add_argument('--scale', dest='scale', required=False,default=10, type=float,help='tree scale,affect branch length default %(default)s ')
parser.add_argument('--pie_with', dest='pie_with', required=False,default=50, type=int,help='pie width  default %(default)s')
parser.add_argument('--fsize', dest='fsize', required=False,default=20, type=float,help='font size   default %(default)s')
#parser.add_argument('--show_branch_length', dest='show_branch_length', required=False,action="store_true",help='whether show_branch_length default %(default)s')



parser.add_argument("-d",'--outdir',dest='outdir',required=False,default=os.getcwd(),help='specify the output file dir,default %(default)s')
parser.add_argument('--prefix', dest='prefix', required=False, default='cafe5_tree_plot',help='file name prefix,default default %(default)s')
#parser.add_argument('--width', dest='width', required=False,default=10, type=int,help='specify the width ,unit is inch default %(default)s ')
parser.add_argument('--height', dest='height',required=False,  default=800,type=int,help='specify the height ,unit is inch default %(default)s')
args = parser.parse_args()
if not os.path.exists(args.outdir): os.mkdir(args.outdir)
dout=os.path.abspath(args.outdir)


from ete3 import Tree,  TextFace,TreeStyle,NodeStyle,faces,AttrFace
import pandas as pd 


tree_data=pd.read_table(args.clade_results,index_col=0)


gf_tatol=args.gfn

pie_data=tree_data.copy()
pie_data.Increase=pie_data.Increase/gf_tatol*100
pie_data.Decrease=pie_data.Decrease/gf_tatol*100

pie_data['noChange']=100-pie_data.Increase-pie_data.Decrease


#nwk="(((Helianthus_annuus<23>:100.197,((Sesamum_indicum<13>:54.5774,Olea_europaea<12>:54.5774)<19>:27.2092,Solanum_lycopersicum<18>:81.7866)<22>:18.4104)<25>:18.4441,(((Populus_trichocarpa<11>:99.6913,(((Acer_yangbiense<1>:5.9205,Acer_truncatum<0>:5.9205)<5>:56.4923,Citrus_sinensis<4>:62.4128)<7>:31.6219,(Arabidopsis_thaliana<3>:75.59,Gossypium_raimondii<2>:75.59)<6>:18.4447)<10>:5.6566)<17>:5.9537,(Juglans_regia<9>:97.9036,Glycine_max<8>:97.9036)<16>:7.7414)<21>:6.3534,(Vitis_vinifera<15>:106.798,Malania_oleifera<14>:106.798)<20>:5.2004)<24>:6.6426)<27>:42.0289,Oryza_sativa<26>:160.67)<28>:20;"

nwk=args.nwk

t = Tree(nwk, format=1)
t.convert_to_ultrametric()

# Basic tree style
ts = TreeStyle()
ts.show_leaf_name = False
ts.show_branch_length = False
ts.show_branch_support = False

ts.scale = args.scale #进化树标尺，值越大进化树越宽

#设置进化树枝显示样式
style = NodeStyle()
style["fgcolor"] = "#0f0f0f"
style["size"] = 0
style["vt_line_color"] = "#000000"
style["hz_line_color"] = "#000000"
style["vt_line_width"] = 2 #进化树枝线条粗细
style["hz_line_width"] = 2
style["vt_line_type"] = 0 # 0 solid, 1 dashed, 2 dotted
style["hz_line_type"] = 0

#为每一个node添加扩张收缩信息
for node in t.traverse():

  node.set_style(style)
  node.img_style["size"] = 0
  node.img_style["shape"] = "square"   # 'circle', 'square' or 'sphere'
  node.img_style["fgcolor"] = "#000000"
  

  #print(node.name)
  if not node.name in tree_data.index:continue
  Increase=tree_data.loc[node.name].Increase
  Decrease=tree_data.loc[node.name].Decrease
  
  exp = TextFace("+"+str(Increase),fsize=args.fsize,fgcolor='red')
  con = TextFace("-"+str(Decrease),fsize=args.fsize,fgcolor='green')
  b = TextFace("/",fsize=args.fsize,fgcolor='blue')
  # Set some attributes
  exp.margin_top = 10
  exp.margin_right = 10
  exp.margin_left = 1
  exp.margin_bottom = 1
  con.margin_top = 10
  con.margin_right = 1
  con.margin_left = 10
  con.margin_bottom = 1
  b.margin_top = 10
  b.margin_right = 0
  b.margin_left = 0
  b.margin_bottom = 1
  #hola.opacity = 0.5 # from 0 to 1
  #hola.inner_border.width = 1 # 1 pixel border
  #hola.inner_border.type = 1  # dashed line
  #hola.border.width = 1
  #hola.background.color = "LightGreen"
  
  # node.add_face(con, column=0, position = "branch-bottom")
  # node.add_face(b, column=1, position = "branch-bottom")
  # node.add_face(exp, column=2, position = "branch-bottom")
  node.add_face(con, column=0, position = "float")
  node.add_face(b, column=1, position = "float")
  node.add_face(exp, column=2, position = "float")
  
  #add pie
  pie=faces.PieChartFace(pie_data.loc[node.name],width=args.pie_with, height=args.pie_with, colors=["red","green","blue"])
  pie.border.width = None
  pie.opacity = 0.5
  pie.margin_left=10
  #faces.add_face_to_node(pie,node, 0, position="branch-bottom")
  
  
  if node.is_leaf(): 
    node.add_face(pie, column=0, position = "branch-top")
   
  else:
    
    node.add_face(pie, column=2, position = "branch-top")
  if node.is_leaf(): 
    
    node.name=node.name.rstrip("<1234567890>")
    N = AttrFace("name", fsize=args.fsize)
    
    node.add_face(N, 1, position="aligned")
    #faces.add_face_to_node(N, node, 0,position="aligned")


  
# 
# for leaf in t.iter_leaves():
#   print(leaf)
#   
# for leaf in t.iter_leaf_names():
#   print(leaf)


root=TextFace(f"MRCA({gf_tatol})",fsize=args.fsize)
root.margin_top = 10
root.margin_right = 10
root.margin_left = 10
root.margin_bottom = 10
t.add_face(root, column=0, position = "branch-top")


#legend add
exp = TextFace("Expansion",fsize=args.fsize,fgcolor='red')
con = TextFace("Contraction",fsize=args.fsize,fgcolor='green')
b = TextFace("/",fsize=args.fsize,fgcolor='blue')
# Set some attributes
exp.margin_top = 10
exp.margin_right = 10
exp.margin_left = 1
exp.margin_bottom = 1
con.margin_top = 10
con.margin_right = 1
con.margin_left = 10
con.margin_bottom = 1
b.margin_top = 10
b.margin_right = 0
b.margin_left = 0
b.margin_bottom = 1


#.legend.add_face(TextFace("Gene families",fsize=10,fgcolor='red'), column=0)

ts.legend.add_face(con, column=0)
ts.legend.add_face(b, column=1)
ts.legend.add_face(exp, column=2)
ts.legend_position=1  #TopLeft corner if 1, TopRight if 2, BottomLeft if 3, BottomRight if 4


t.convert_to_ultrametric()
#t.ladderize(direction=1) #自动排序

os.chdir(dout)
t.render(args.prefix+".pdf", tree_style=ts,h=args.height,dpi=300)
t.render(args.prefix+".png", tree_style=ts,h=args.height, dpi=300)


#t.render(dout+"/"+args.prefix+".pdf", tree_style=ts,h=args.height,dpi=300)
#t.render(dout+"/"+args.prefix+".png", tree_style=ts,h=args.height, dpi=300)



