library(flowCore)
library(ggcyto)
library(ggplot2)
library(flowWorkspace)
library(flowStats)

cnt = read.FCS('~/Desktop/Translation_Screen/FlowJo_Data/May3_HeLa_sam_Main/Specimen_001_Ctr_hela_004.fcs', transformation = "linearize-with-PnG-scaling", alter.names = T )
sam = read.FCS('~/Desktop/Translation_Screen/FlowJo_Data/May3_HeLa_sam_Main/Specimen_001_HeLaSAM_005.fcs', transformation = "linearize-with-PnG-scaling", alter.names = T )

# SCALE TRANSFORM
# cnt = read.FCS('~/Desktop/Translation_Screen/FlowJo_Data/May3_HeLa_sam_Main/Specimen_001_Ctr_hela_004.fcs', transformation = "scale", alter.names = T )
# sam = read.FCS('~/Desktop/Translation_Screen/FlowJo_Data/May3_HeLa_sam_Main/Specimen_001_HeLaSAM_005.fcs', transformation = "scale", alter.names = T )

# INITIAL TEST EXPERIMENT
# cnt = read.FCS('~/Desktop/Translation_Screen/FlowJo_Data/0429_HeLa_TestRun/Specimen_001_Unstained_003.fcs', transformation = "linearize", alter.names = T )
# sam = read.FCS('~/Desktop/Translation_Screen/FlowJo_Data/0429_HeLa_TestRun/Specimen_001_HeLa_RNA_004.fcs', transformation = "linearize", alter.names = T )

autoplot(cnt, "FSC.A", "SSC.A")
autoplot(sam, "FSC.A", "SSC.A")
autoplot(cnt, "FSC.A", "FSC.W")
autoplot(sam, "FSC.A", "FSC.W")
autoplot(cnt, "SSC.H", "SSC.W")
autoplot(sam, "SSC.H", "SSC.W")

experiment = as(list(cnt, sam), "flowSet")

FSC_A_FSC_W = rectangleGate(filterId = "AW", 
                            "FSC.A" = c(100000, 300000), "FSC.W" = c(100000, 150000))
SSC_A_FSC_A = rectangleGate(filterId = "scatter", 
                            "FSC.A" = c(100000, 250000), "SSC.A" = c(100, 100000))
SSC_H_SSC_W = rectangleGate(filterId = "HW", 
                            "SSC.H" = c(1000, 30000), "SSC.W" = c(100000, 150000))

# # Scale Gates
# FSC_A_FSC_W = rectangleGate(filterId = "AW", 
#                             "FSC.A" = c(0.25, 0.9), "FSC.W" = c(0.25, 0.8))
# SSC_A_FSC_A = rectangleGate(filterId = "scatter", 
#                             "FSC.A" = c(0.25, 0.9), "SSC.A" = c(0, 0.5))
# SSC_H_SSC_W = rectangleGate(filterId = "HW", 
#                             "SSC.H" = c(0, .95), "SSC.W" = c(.3, .7))

morphology_gates = FSC_A_FSC_W & SSC_A_FSC_A &SSC_H_SSC_W
summary(filter(experiment,morphology_gates)) 

cells = Subset(experiment, morphology_gates)
autoplot(cells, "FSC.A", "SSC.A")
autoplot(cells, "FSC.A", "FSC.W")
autoplot(cells, "SSC.H", "SSC.W")
autoplot(cells, "APC.A", "FITC.A")


# Modify Cells 

p <- ggcyto(transform(cells[[2]], 
                      `APC.A` = log10(`APC.A`), 'FITC.A' = log10(`FITC.A`) ) , 
            aes(x = `APC.A`, y =  `FITC.A`))
p <- p + geom_hex(bins = 128) + ggtitle("CRISPR-based Overexpression Screen") + ylim(c(2.5, 5.5) ) + xlim(c(2, 5))
p

p_cnt <- ggcyto(transform(cells[[1]],
                `APC.A` = log10(`APC.A`), `FITC.A` = log10(`FITC.A`  / 2) ) , 
            aes(x = `APC.A`, y =  `FITC.A`))
p_cnt <- p_cnt + geom_hex(bins = 128) + ggtitle("Control Cells") + 
  ylim(c(2.5, 5.5) )  + xlim(c(2, 5))
p_cnt


# ggplot( transform ( cells, 'Ratio' = log10 ( `FITC.A` /  `APC.A`) ) , aes( x= 'Ratio' ) ) + geom_histogram(alpha = .2, fill = "red", bins = 128)
summary ( transform ( cells[[1]], 'Ratio' = `FITC.A` /  ( 2.154 * `APC.A`)) )
summary ( transform ( cells[[2]], 'Ratio' = `FITC.A` / (1.116 * `APC.A`)) )

cnt_ratio = head(transform ( cells[[1]], 'Ratio' = `FITC.A` /  ( 2.1549 * `APC.A`)) , 5847 ) 
sam_ratio = head(transform ( cells[[2]], 'Ratio' = `FITC.A` / (1.1155 * `APC.A`)) , 26611 ) 

sam_ratio = log10(sam_ratio[,11]) 
cnt_ratio = log10(cnt_ratio[,11])

round (quantile(cnt_ratio , seq ( 0,1 , .01)) , 3 ) 
round (quantile(sam_ratio , seq ( 0,1 , .01)) , 3 ) 

sam_ratio[sam_ratio < log10 ( 0.5 ) ] = log10 ( 0.5 ) 
cnt_ratio[cnt_ratio < log10 ( 0.5 ) ] = log10 ( 0.5 ) 
sam_ratio[sam_ratio > log10 ( 3) ] = log10 ( 3 ) 
cnt_ratio[cnt_ratio > log10 ( 3) ] = log10 ( 3 ) 


dat <- data.frame( ratio = c(cnt_ratio,sam_ratio), 
          type = c( rep("cnt", length(cnt_ratio)), rep("sam", length(sam_ratio))  ) )
ggplot(dat, aes(x=ratio, fill=type)) + geom_density(alpha=0.2, position="identity")


# Test GatingSet
gs = GatingSet(experiment)
#gs = workFlow(experiment)

add (gs, SSC_A_FSC_A, parent = "root", name = "debris")
add (gs, SSC_H_SSC_W, parent = "debris")
add (gs, FSC_A_FSC_W , parent = "HW")
getNodes(gs)
# Gate the data
recompute(gs)
autoplot(gs, "AW")
getProp(gs[[1]], "AW" )

# It is possible to do biexp-transformation with ?flowJo_biexp_trans

autoplot(gs[[1]])
autoplot(gs[[2]])

getData(gs[[2]], "HW")

pars = c("APC.A", "FITC.A")
norm <- normalization(normFun=function(x, parameters, ...)
  warpSet(x, parameters, ...),
  parameters=pars,
  arguments=list(grouping="GroupID", monwrd=TRUE),
  normalizationId="Warping")
add (gs, norm, parent = "AW")

# t1 =  estimateLogicle(gs[[1]], c("APC.A", "FITC.A"))
# gs = transform (gs, t1)
# autoplot(flowData(gs)[[2]], "AW", "APC.A" )



## Mar 28 HPG OPP K562 Test
# HPG, OPP stained with FITC; RNA with APC
# FITC Signal was saturated so reduced PMT to 313 and 247. 
# There was a big loss of OPP cells. 
cnt = read.FCS('~/Desktop/Translation_Screen/FlowJo_Data/Mar28_HPG_OPP_K562_Test/file=2-5.fcs',  transformation = "linearize-with-PnG-scaling", alter.names = T)
hpg_247 = read.FCS('~/Desktop/Translation_Screen/FlowJo_Data/Mar28_HPG_OPP_K562_Test/file=2-9.fcs', transformation = "linearize-with-PnG-scaling", alter.names = T)
hpg_313 = read.FCS('~/Desktop/Translation_Screen/FlowJo_Data/Mar28_HPG_OPP_K562_Test/file=2-10.fcs', transformation = "linearize-with-PnG-scaling", alter.names = T)
opp_247 = read.FCS('~/Desktop/Translation_Screen/FlowJo_Data/Mar28_HPG_OPP_K562_Test/file=2-8.fcs', transformation = "linearize-with-PnG-scaling", alter.names = T)
opp_313 = read.FCS('~/Desktop/Translation_Screen/FlowJo_Data/Mar28_HPG_OPP_K562_Test/file=2-7.fcs', transformation = "linearize-with-PnG-scaling", alter.names = T)

autoplot(cnt, "FSC.A", "SSC.A")
autoplot(cnt, "FSC.A", "FSC.W")
autoplot(cnt, "SSC.H", "SSC.W")
autoplot(opp_247, "FSC.A", "SSC.A")
autoplot(opp_247, "FSC.A", "FSC.W")
autoplot(opp_247, "SSC.H", "SSC.W")

FSC_A_FSC_W = rectangleGate(filterId = "AW", 
                            "FSC.A" = c(100000, 250000), "FSC.W" = c(80000, 100000))
SSC_A_FSC_A = rectangleGate(filterId = "scatter", 
                            "FSC.A" = c(100000, 250000), "SSC.A" = c(100, 30000))
SSC_H_SSC_W = rectangleGate(filterId = "HW", 
                            "SSC.H" = c(1000, 25000), "SSC.W" = c(10000, 90000))



morphology_gates = FSC_A_FSC_W & SSC_A_FSC_A &SSC_H_SSC_W

experiment = as(list(cnt, hpg_247, hpg_313, opp_247, opp_313), "flowSet")
summary(filter(experiment,morphology_gates)) 
cells = Subset(experiment, morphology_gates)
autoplot(cells, "FSC.A", "SSC.A")
autoplot(cells, "FSC.A", "FSC.W")
autoplot(cells, "SSC.H", "SSC.W")
autoplot(cells, "APC.A", "FITC.A")


