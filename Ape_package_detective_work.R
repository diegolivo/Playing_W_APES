install.packages("ape")
install.packages("phytools")
install.packages("tidytree")
install.packages("ggtree")
install.packages("geiger")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")

library(ggtree) ##ggtree is mainly used to change the aestehtics of a tree
library(ape) ##Builds trees and checks for fundamental characteristics of a tree
library(tidyverse)
library(tidytree) ##Tidytree allows you to manipulate tree data in a tibble format
library(phytools)## Adds more tools like midpoint rooting and defining color schemes
library(geiger) ##Geiger is a package that allows us to check a dataset to our tree to ensure consistency and no data was lost

##To make life easier, remove all spaces/underscores in tip labels on a tree

##For testing I am using my cluster from my wigeon fecal paper cluster 17 which is unclassified cressdnaviruses

midroot_cl17 <-  read.tree(file = "mpr_Cl17.nwk" ) ## read.tree tells you general information about tree file including numer of nodes and tips with all tip labels and converts to a new object defined as phylo

##Important note, most commands need a defined object as phylo to run in ape package

ur_cl17 <- read.tree(file = "ur_Cl_17.nwk" ) ##Comparing against unrooted tree (ur)

is.rooted(midroot_cl17) ##Will inform you if your tree is rooted with true or false 
          
is.rooted(ur_cl17) ##compare against unrooted tree

comparePhylo(ur_cl17, midroot_cl17) ##Will inform you and compare number of nodes, tips/tip labels, the rooting differences, and splits

checkValidPhylo(ur_cl17) ##CheckValidPHylo will inform you the elements of object phy and inform you if errors when converted to the phy object

checkValidPhylo(midroot_cl17)

plot(ur_cl17) ##When you have defined phy objects you can plot the trees in RStudio and visualize while you manipulate

plot(midroot_cl17)

## Side note: when plotting, if you recieve "Error in plot.new() : figure margins too large" just adjust the plot pane and rerun command


labeled_nodes_mpr_17 <- midroot_cl17 |> makeNodeLabel(method = "number") ##This command makes lables for you nodes with parameters defined in method = "", in this case I used numbers to define my nodes

##WARNING: MakeNodeLabels will NOT actually correspond to the actual nodes, this WILL THROW OFF YOUR PARSING OF A TREE

nodelabels(labeled_nodes_mpr_17) ##Allows you to see how the nodes are defined in seperate list, this is still confusinf when the command is run



rotated_node5_mpr_17 <- rotate(labeled_nodes_mpr_17, node = c("'MH973752 C17'", "'KY487927 C17'")) ##You can rotate the tips in a clade with the rotate function. define the two tips you want to rotate with node = c(). IN this case I rotated the two tips in Node5


##Compare the roated tree with the unrotated tree. 

plot(rotated_node5_mpr_17,  show.node.label = TRUE)

plot(labeled_nodes_mpr_17,  show.node.label = TRUE)


##Removing tip labels

plot(dup_cl9_mpr,show.node.label = TRUE)  


dropped_tip_cl9 <- drop.tip(dup_cl9_mpr, tip= "'MK032746_1'", trim.internal = TRUE, subtree = FALSE,
         root.edge = 0, rooted = is.rooted(dup_cl9_mpr), collapse.singles = FALSE,
         interactive = FALSE)

##For dropping one tip, use the drop.tip() command. Define tip to drop with tip="TipName" 

##trim.internal will delete the internal branches from the tip 

##root.edge = INTEGER will root use the number of defined internal branches to build a new root edge


#Comparing dropped tip from duplicate MK032746

plot(dropped_tip_cl9)  
plot(dup_cl9_mpr) 


## We can turn using the Package tidytree, we can turn our phylogenetic data into a tibble so we can refer and understand how we want to manipulate our trees 

dup_df_cl9 <- as_data_frame(dup_cl9_mpr) ##turing duplicate cl9 to tibble
df_cluster_9 <- as_data_frame(dropped_tip_cl9) ##turning dropped tip cl9 into a tibble 

view(df_cluster_9)
view(dup_df_cl9)

##This is useful because you can see node position and the corresponding tip label. When trying to define and manipulate tree you can see what variable a command is using 

checkLabel(dropped_tip_cl9)
zoom(dropped_tip_cl9, focus = 3:4, subtree = FALSE, col = rainbow,show.node.label = TRUE) ##Focus zooms into node number defined in df_cluster_9
zoom(dup_cl9_mpr, focus = 3:4, subtree = FALSE, col = rainbow,show.node.label = TRUE) ##Focus zooms into node number in dup_df_cl9


##Midpoint rooting via phytools package

r_mpr_cl17 <- midpoint.root(ur_cl17)

##Compare the two 
plot(r_mpr_cl17) 
plot(midroot_cl17)

##Based on visual comparison the two trees and midpoint rooted the same but the clades are in different orientations

df_cluster_17 <- as_data_frame(r_mpr_cl17)

##Aesthetics of tree
plot.phylo(r_mpr_cl17, edge.color = "red", tip.color = "blue", align.tip.label = TRUE) ##Set the standard for colors on trees and align our tip labels
           
  
##Edge/branch coloring

sample <- which.edge(r_mpr_cl17, "'Bird K19 MH 385 lt 2802 cir 2661 C17'") ##To color specific branches from our tree, we first need to define the edges by using which.edge() command to creat a vector of the tips of interest
sample  ##Our vector with edges of our samples 
co <- rep("black", Nedge(r_mpr_cl17)) ##WE set the default edge (branch) color to black for all edges 
co[sample] <- "red" ##We specify the values of our vector to be red
par(xpd = TRUE) ##Open a graphical device (or reset it if already open) and draw an empty plot with the limit found at the previous step.
plot(r_mpr_cl17, edge.col = co) ##Now we plot our tree with the colors we defined with edge.col = co
edgelabels() ##Shows the numbering of edges for easier visualization 


mycol2 <- def(r_mpr_cl17$tip.label, "'Bird K19 MH 385 lt 2802 cir 2661 C17'" = "red", default = "blue") ##Here we can define which tips we want to color and set a default by using def(phylo$tip.label, x= "color1", default = "color2")

plot(r_mpr_cl17, tip.color = mycol2) ##Add our defined color scheme to tip.color in plot()

par(xpd = TRUE) ##Open a graphical device (or reset it if already open) and draw an empty plot with the limit found at the previous step.

plot(r_mpr_cl17, tip.color = mycol2, edge.color = co, show.node.label = TRUE, align.tip.label = TRUE) ## All together



gtr80_cluster17 <- read.tree("Cluster_17t_gtr80.nwk") |>
  midpoint.root() ##Import tree with defined aLRT branch support and midpoint root it

df_gtr80 <- as_data_frame(gtr80_cluster17) ##Turning tree to data frame helps alot


plot(gtr80_cluster17)



rename.tips()


##Do the edge color and node color step by step, R doesnt like it when you try to do it all once

sample <- which.edge(gtr80_cluster17, "'Bird K19 MH 385 lt 2802 cir 2661 C17'")
co <- rep("blue", Nedge(gtr80_cluster17))
co[sample] <- "red"
mycol2 <- def(gtr80_cluster17$tip.label, "'Bird K19 MH 385 lt 2802 cir 2661 C17'" = "red", default = "blue")
plot(gtr80_cluster17, edge.color = co, tip.color = mycol2)

##For ALRT support we need to ensure that our tree already has it, node labels will show you the aLRT branch support

par(xpd = TRUE)

plot(gtr80_cluster17, tip.color = mycol2, edge.color = co, align.tip.label = TRUE, show.node.label = TRUE) |> ##Showing node labels as they show aLRT Branch support
  edgelabels(text = c(100, 95.5, 100), edge = c(1,2,8), adj = c(0.5, 0.5), frame = "circle", ##So far this is the best way to show aLRT by defining the edge as a vector of labels with the corresponding edge number to get the label at the halfway mark of the branch 
             pch = NULL, thermo = NULL, pie = NULL, piecol = NULL, ##Need these to input a shape for the edge label, will test these parameters
                col = "black", bg = "white")  ##Color the labels


##The method above I found to be easier then trying to manipulate nodes from below


##Testing other methods for aLRT branch support using node labels and to be honest, it sucks

plot(gtr80_cluster17, tip.color = mycol2, edge.color = co, align.tip.label = TRUE, show.node.label = TRUE) 
  nodelabels( c(100,100,95), node = c(8,9,11), adj = c(5, 0.5), frame = "rect",
           pch = NULL, thermo = NULL, pie = NULL, piecol = NULL,
           col = "black", bg = "lightblue", horiz = FALSE,
           width = NULL, height = NULL)
  

##Import data of names and amino acid sequences for cluster 17

  
fxd_gtr80_cluster17 <- read.tree("fxd_Cluster_17t_gtr80.nwk") ##Adjust names that have names of tips without underscores and notes

df_gtr80 <- as_data_frame(fxd_gtr80_cluster17) ##Create df to compare and troubleshoot during testing phase
  
names_cl17 <- read.csv("Cluster_17_names_AA.csv",  header = TRUE, row.names = 1) ##Import  CSV with row.names =1 to indicate that geiger is looking at row one for comparative analysis for tip names, header = TRUE keeps the names of coluns thesame

##Testing name.check function from geiger

##The name,check

name.check(fxd_gtr80_cluster17, names_cl17, data.names = NULL) ##For name.check, the command wants an input of a tree, dataset

##data.names = NULL will look for the names of the tips from both tree without looking at order

name.check(fxd_gtr80_cluster17, names_cl17, data.names = TRUE) ##Changing  data.names will tell you that it did not find the order of tip names based on the order of the ree


##Playing around with different inputs and outputs with geiger name.check()

##in this example I used cluster 9 that has two sequences with the same accesion number(MK03247) with one note as genono and one as unclass

fxd_gtr80_cluster9 <- read.tree("fxd_Cluster_9_alignment_gtr80.nwk") ##Adjust names that have names of tips without underscores and notes

plot(fxd_gtr80_cluster9)

mpr_cl9 <- midpoint.root(fxd_gtr80_cluster9) ##Midpoint root tree

plot(mpr_cl9)


names_cl9 <- read.csv("Cluster_9_names.csv",  header = TRUE, row.names = 1) ##Import  CSV with row.names =1 to indicate that geiger is looking at row one for comparative analysis for tip names, header = TRUE keeps the names of coluns thesame

name.check(fxd_gtr80_cluster9, names_cl9, data.names = NULL) ##Without editing the tree or dataset the output will say OK meaning that all the data from the tree is in the dataset we chose


dropped_tip_cl9 <- drop.tip(mpr_cl9, tip= "MK032746unclass", trim.internal = TRUE, subtree = FALSE,
                            root.edge = 0, rooted = is.rooted(mpr_cl9), collapse.singles = FALSE,
                            interactive = FALSE) ##dropping the MK03247unclass sequence tip

plot(dropped_tip_cl9)

##Now when we check it, data_not_tree will inform us of sequences in the dataset not present in the tree
##IF we didn't removed the name from our dataset instead then we will recieve tree_not_data meaning sequence in tree is not found in data

name.check(dropped_tip_cl9,names_cl9, data.names = NULL)

## IMPORTANT NOTE: IF FORMATTING OF NAMES IS INCONSISTENT THEN IT WILL BE DIFFICULT TO KNOW WHERE THE ERRORS ARE!!!!!




