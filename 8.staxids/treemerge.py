#!/usr/bin/env python3

#imports sets of taxids sorted by loci, generates a tree that notes how it was found, by loci,
# and colors accordingly.

#input: "18S.staxids" and etc for each locus defined in "loci" array below.

#uncomment the "update_taxonomy_database" line for first run, that pulls current lineage info from NCBI.

from ete3 import NCBITaxa
ncbi = NCBITaxa()
# ncbi.update_taxonomy_database()

from ete3 import Tree, TreeStyle, TextFace,add_face_to_node
trees={}
loci=["18S","COI","28S","16S"]
colors={"18S":"pink","COI":"blue","28S":"green","16S":"orange"}
taxids=[]
taxids_by_locus={}
for locus in loci:
    cur_taxids=[]
    cur_filename =  locus+'.staxids'
    print ("Reading",cur_filename+":",)
    with open(cur_filename, 'r') as infile:
        for cur_line in infile:
            cur_taxid = int(cur_line)
            if cur_taxid not in cur_taxids:
                cur_taxids.append(cur_taxid)
    taxids += (cur_taxids)
    taxids_by_locus[locus] = cur_taxids
    print (len(taxids_by_locus[locus]))



tree = ncbi.get_topology(taxids,intermediate_nodes=False)
# for locus in loci:
#     print len(taxids_by_locus[locus])
#print tree.get_ascii(attributes=["sci_name", "rank"])
# #
# for node in tree.iter_descendants("postorder"):
#     # Do some analysis on node
#     node.name = node.sci_name
#     # print node.sci_name,node.taxid

# tree.show(tree_style=ts)

def layout(node):
    if node.is_leaf():
        node.name = node.sci_name
        loci_list=[]
        for locus in loci:

            if node.taxid in taxids_by_locus[locus]:
                loci_list.append(locus)
        if len(loci_list)==1:
            node.img_style["fgcolor"] = colors[loci_list[0]]
            node.img_style["size"] = 10
            node.name = node.sci_name +"[" +loci_list[0] + "]"
        else:
            node.img_style["fgcolor"] = "red"
            node.img_style["size"] = 10* len(loci_list)
            primer_list=""
            for locus in loci_list:
                primer_list+="["+locus+"]"
            node.name = node.sci_name + " " + primer_list
    else:
        text = node.rank +":" +node.sci_name
        F = TextFace(text, tight_text=True)
        add_face_to_node(F, node, column=0, position="branch-right")
        # node.img_style["fgcolor"]="black"
        # node.name="foo"



    

ts = TreeStyle()
ts.layout_fn = layout
ts.force_topology=True
ts.mode="r" # valid modes are c and r
#ts.scale=200


tree.render("mytree.pdf", w=1000, units="mm",tree_style=ts)

# write tree
# tree.write(format=1, outfile="full_tree.nw",quoted_node_names=True)
# write annotations file
