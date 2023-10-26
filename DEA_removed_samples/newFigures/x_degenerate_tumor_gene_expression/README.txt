READ ME for X Degenerated Tumor Expression Plots
-----------------------------------------------------------

Last updated: 10/12/2023


PURPOSE OF FOLDER: Gene expression boxplots from a series of X-degenerate genes on ICGC Tumor samples. 

DATA USED:
This folder contains gene expression plots from ICGC read count data of all tumor tumor samples used in subsequent analyses. The expression data is raw counts. The only quality parameters that have been applied are keeping genes based on a mean FPKM of at least 0.5 in one group. Also, a TMM of at least 6 reads in at least 10 samples. 

SCRIPT: The script used to generate plot is Sex_check_automation_removed_samples.

HOW TO READ FIGURE: The X-axis is the sample number (three digit identifier) from all tumor-adjacent samples. The Y-axis is the read count expression data (no logarithm was applied). The samples are color coded and ordered by sex and viral etiology.

0 equals male HBV sample
1 equals male HCV sample
2 equals female HBV sample
3 equals female HCV sample

FOLDER CONTENTS: 
---------------------------------------------------

Boxplots were generated for the following genes:

	-- AR
	-- DDX3Y
	-- EIF1AY
	-- KDM5D
	-- TSPY1 
	-- TXLNGY
	-- XIST

Two files were generated for each gene. One file has the y axis capped at 500 to more easily see delineate between low to no gene expression.
 
FILES:
----------------------------------

AR_tumor_expression_500Axis_data.pdf -- Boxplot of AR gene expression data from ICGC tumor samples with y axis capped at 500 

AR_tumor_expression_data.pdf -- Boxplot of AR gene expression data from ICGC tumor samples 

DDX3Y_tumor_expression_500Axis_data.pdf --  Boxplot of DDX3Y gene expression data from ICGC tumor samples with y axis capped at 500 

DDX3Y_tumor_expression_data.pdf -- Boxplot of DDX3Y gene expression data from ICGC tumor samples 

EIF1AY_tumor_expression_500Axis_data.pdf --  Boxplot of EIF1AY gene expression data from ICGC tumor samples with y axis capped at 500 

EIF1AY_tumor_expression_data.pdf -- Boxplot of EIF1AY gene expression data from ICGC tumor samples 

KDM5D_tumor_expression_500Axis_data.pdf --  Boxplot of KDM5D gene expression data from ICGC tumor samples with y axis capped at 500 

KDM5D_tumor_expression_data.pdf -- Boxplot of KDM5D gene expression data from ICGC tumor samples 

TSPY1_tumor_expression_500Axis_data.pdf --  Boxplot of TSPY1 gene expression data from ICGC tumor samples with y axis capped at 500 

TSPY1_tumor_expression_data.pdf -- Boxplot of TSPY1 gene expression data from ICGC tumor samples 

TXLNGY_tumor_expression_500Axis_data.pdf --  Boxplot of TXLNGY gene expression data from ICGC tumor samples with y axis capped at 500 

TXLNGY_tumor_expression_data.pdf -- Boxplot of TXLNGY gene expression data from ICGC tumor samples 

XIST_tumor_expression_500Axis_data.pdf --  Boxplot of XIST gene expression data from ICGC tumor samples with y axis capped at 500 

XIST_tumor_expression_data.pdf --  Boxplot of XIST gene expression data from ICGC tumor samples 
