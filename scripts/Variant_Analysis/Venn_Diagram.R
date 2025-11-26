#This script was made to make a venn diagram showing the overlapp between the the different datasets

 setwd("/Users/stella/Desktop/Master_Projekt/ABraOM_SNV_data")

 #install.packages("VennDiagram")
 library(VennDiagram)
 

grid.draw(venn.plot);
grid.newpage

#these numbers were retrieved using the python script "vcf_intersections.py"
 venn.plot <- draw.triple.venn(
   area1 = 1209,
   area2 = 12141,
   area3 = 29638,
   n12 = (9+950),
   n23 = (8642+950),
   n13 = (172+950),
   n123 = 950,
   category = c("ABraOM", "ClinVar", "gnomAD"),
   fill = c("#d2009e", "#f15249", "#ef8f00"),
   cex = 2.5,
   col = "black",
   cat.cex = 1.65,
   cat.pos = c(285, 105, 190),
   cat.dist = 0.09,
   cat.just = list(c(0.05, 0), c(1, 1), c(0,-4)),
   ext.pos = 30,
   ext.dist = -0.05,
   ext.length = 0.85,
   ext.line.lwd = 2,
   ext.line.lty = "dashed",
   fontfamily = "Arial",
   cat.fontfamily = "Arial",
   lwd = 1,
   ity = dash
 );
 grid.newpage();
 #change backround colour to match other figure
 grid.draw(rectGrob(gp = gpar(fill = "#f6f6f6")))
 
 grid.draw(venn.plot);
 
