# Olfactory-Bulb---AR---2024
R code and data for the publication "But how does it smell? An investigation of olfactory bulb size among living and fossil primates and other euarchontoglirans"
File:
EG-Fossil-Sept-Treemaking-OB --> R code for making a time calibrated tree using the paleotree function "bin_cal3TimePaleoPhy". Requires the following files: int_times_OB.csv, taxon_times_OB.csv, Node-Legend.csv, mytree.nex. Outputs file: timecaltree-OB.nex.
int_times_OB.csv --> bin ages needed to "runbin_cal3TimePaleoPhy". 
taxon_times_OB.csv --> designation of what bin each taxon blelongs in needed to "runbin_cal3TimePaleoPhy". 
Node-Legend.csv --> the ages of specific nodes in the tree using in the "runbin_cal3TimePaleoPhy" function. Node ages obtained from original extant only tree downloaded from vertlife.org. Ages are also seen in mytree.nex file.
mytree.nex --> nexus tree file with extant taxa (original tee from vertlife.org) and fossils added in mesquite. Also used to run analyses in script: OB-script_May.
timecaltree-OB.nex --> the resulting time calibrated tree used in analyses in OB-script_May R code.
OB-Classifier-log.csv --> the datafile with OB measurements and clade groupings used to run analyses in OB-script_May R code.
OB-script-May --> R code to perfom analyses on OB data. 
