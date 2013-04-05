# Make updated pipeline diagram
dot -Tpng SEAStAR_Pipeline_Stages.gv >SEAStAR_Pipeline_Stages.png

marked SEAStAR_User_Guide.md > SEAStAR_User_Guide.html
marked ../README.md > ../README.html
marked RDP_vignette.md > RDP_vignette.html

# Make updated pdf user guide document
# This is run twice to get the table of contents and back reference correct.
cp SEAStAR_User_Guide.html ..
cp SEAStAR_User_Guide.md ..
cp style.css ..
cp SEAStAR_Pipeline_Stages.png ..
cp RDP_vignette.html ../vignettes/RDP
cp 16S_abundance.png ../vignettes/RDP

