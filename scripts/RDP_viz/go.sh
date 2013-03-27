# Running this script will use the various Javascript programs to generate visualizable data 
# from RDP training set data and RDP classifier output (augmented with abundance data from ref_select)
# After running this script, the results can be visualized by loading the .html files in this directory 
# a modern (e.g. HTML 5 compliant) browser.
# Running this script requires that Coffeescript and Node.js be installed on this system:
#
# Tested with Coffeescript v1.2.0: http://coffeescript.org/
# Tested with Node.js v0.6.13: http://nodejs.org/
#
# Under OS X, these can be installed through MacPorts:
# e.g. "sudo port install coffee-script"

cat ./test_data/trainset6_db_taxid.txt | coffee ../RDP_train_to_tree.coffee > RDP_expand.json

cat ./test_data/GG2.d40_filtered_class.out | coffee ../RDP_tree_dev.coffee ./RDP_expand.json | tee GG2_RDP_new.json | gawk 'BEGIN { print "var tree_struct = "; } (1); END { print ";"}' > GG2_data_new.json
cat ./test_data/GG3.d50_filtered_class.out | coffee ../RDP_tree_dev.coffee ./RDP_expand.json | tee GG3_RDP_new.json | gawk 'BEGIN { print "var tree_struct = "; } (1); END { print ";"}' > GG3_data_new.json

#cat ./test_data/GG2.d40_filtered_class.out | js -e 'RDP_heir_fn = "./RDP_viz/RDP_expand.json";' -f ../RDP_tree.js | tee GG2_RDP.json | gawk 'BEGIN { print "var tree_struct = "; } (1); END { print ";"}' > GG2_data.json
#cat ./test_data/GG3.d50_filtered_class.out | js -e 'RDP_heir_fn = "./RDP_viz/RDP_expand.json";' -f ../RDP_tree.js | tee GG3_RDP.json | gawk 'BEGIN { print "var tree_struct = "; } (1); END { print ";"}' > GG3_data.json

cat GG2_RDP_new.json | coffee ../RDP_tree_nexus_dev.coffee > GG2_RDP_new.nex
cat GG3_RDP_new.json | coffee ../RDP_tree_nexus_dev.coffee > GG3_RDP_new.nex

#cat GG2_RDP.json | js -f ../RDP_tree_nexus_dev.coffee > GG2_RDP.nex
#cat GG3_RDP.json | js -f ../RDP_tree_nexus_dev.coffee > GG3_RDP.nex

coffee ../RDP_tree_merge_dev.coffee ./GG2_RDP.json ./GG3_RDP.json | gawk 'BEGIN { print "var tree_struct = "; } (1); END { print ";"}' > GG2_GG3_merged_data_new.json

# js -e 'input_files = { GG2 : "./RDP_viz/GG2_RDP.json", GG3 : "./RDP_viz/GG3_RDP.json" };' -f ../RDP_tree_merge.js | gawk 'BEGIN { print "var tree_struct = "; } (1); END { print ";"}' > GG2_GG3_merged_data.json


