###### The Big and Tragic Rearrangement of Mouse Biogeography ---
### Date: 2.19.2023
### Figure Number: 4
### Determine intersecting features between datasets

setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/")

## Numbers of Genera that overlap --- 

ucla_o_lum_genera <- readRDS(file = "regionalluminalgenera.RDS")
length(ucla_o_lum_genera) # 44
ucla_o_muc_genera <- readRDS(file = "regionalmucosalgenera.RDS")
length(ucla_o_muc_genera) # 44

ucla_v_muc_genera <- readRDS(file = "immdefgenera.RDS")
length(ucla_v_muc_genera) # 39

cs_lum_genera <- readRDS(file = "CS_Facility_Luminal_Genera.RDS")
length(cs_lum_genera) # 43 
cs_muc_genera <- readRDS(file="CS_Facility_Mucosal_Genera.RDS")
length(cs_muc_genera) # 39

length(union(cs_lum_genera, ucla_o_lum_genera)) #58 
muc1 <- union(ucla_o_muc_genera, ucla_v_muc_genera)
length(muc1) #58 
muc2 <- union(cs_muc_genera, ucla_v_muc_genera)
length(muc2) #58
muc3 <- union(ucla_o_muc_genera, cs_muc_genera)
length(muc3) #54
length(union(muc1,cs_muc_genera)) #65
17/65

ucla_o_v_muc_together <- intersect(ucla_o_muc_genera, ucla_v_muc_genera)
length(ucla_o_v_muc_together) #25
25/58
ucla_o_cs_muc <- intersect(ucla_o_muc_genera, cs_muc_genera)
length(ucla_o_cs_muc) #29 
29/54
ucla_v_cs_muc <- intersect(ucla_v_muc_genera, cs_muc_genera)
length(ucla_v_cs_muc) #20
20/58
ucla_cs_muc_all <- intersect(ucla_o_v_muc_together, cs_muc_genera)
length(ucla_cs_muc_all) #17 
