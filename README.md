# [Machine learning-driven network biomarker discovery and risk scoring system construction for breast cancer prognosis](https://github.com/zpliulab/CNet-Cox)

<!-- CNet-Cox: Network prognostic biomarkers identification by Cox proportional hazards model with connectivity network-regularized constraints -->

![Screenshot](Figure/Framwork_Github.png)


## CNet-Cox
<!--START_SECTION:news-->
* A **C**nnected **Net**work-regularized **Cox** roportional hazards model (**CNet-Cox**) is proposed to perform **feature selection**. 
* In real-world **breast cancer (BRCA) dataset**, we validated the CNet-Cox model is efficient to identify the **connected-network-structured features** that can serve as **prognostic biomarkers**.
* In the comparison study, we also proved the proposed **CNet-Cox** results in better classification performance and feature interpretability than other seven method named **ENet-Cox**, **Lasso-Cox**, **L0-Cox**, **L1/2-Cox**, **SCAD-Cox**, **MCP-Cox** and **Ridge-Cox**.
* If you have any questions about **CNet-Cox**, please directly contact the corresponding author [Prof. Zhi-Ping Liu](https://scholar.google.com/citations?user=zkBXb_kAAAAJ&hl=zh-CN&oi=ao) with the E-mail: zpliu@sdu.edu.cn
<!--END_SECTION:news-->


## Citation
Lingyu Li, Wai-Ki Ching and Zhi-Ping Liu*. "**Machine learning-driven network biomarker discovery and risk scoring system construction for breast cancer prognosis**." Submited to [npj Digital Medicine](https://www.nature.com/npjdigitalmed/).
<!-- [Expert Systems with Applications](https://www.journals.elsevier.com/expert-systems-with-applications/).   -->

## R packages
* [glmSparseNet](https://bioconductor.org/packages/release/bioc/html/glmSparseNet.html) (v1.8.1). 
* [curatedTCGAData](https://www.bioconductor.org/packages/release/data/experiment/html/curatedTCGAData.html) (v1.12.1). 
* [TCGAutils](https://bioconductor.org/packages/release/bioc/html/TCGAutils.html) (v1.10.1). 
* [dplyr](https://cran.r-project.org/web/packages/dtplyr/index.html) (v1.0.8). 
* [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) (v1.30.1). 

## Data
<!--START_SECTION:news-->
* **Code** provides all **R/Matlab/Python** scripts required to reproduce and utilize **CNet-Cox**. 
* **Data** gives necessary input/output files by the **R/Matlab/Python** codes. *Note: Some input files only give the first few lines due to upload size limitation, but this does not affect the results of **CNet-Cox**.*
* **Supp** presents additional supporting files that contained in **CNet-Cox**. 
* **Validation** provides the original supporting files (high-resolution images) used for validating the PRS genes induced by **CNet-Cox** for breast cancer. 
<!--END_SECTION:news-->


## R codes for Data 
The **serial number (1), (2), ..., (16)** represents the order in which the program runs in our work. 
<!--START_SECTION:news-->
* (1) ``TCGA_pro_clin_DE.R``  --  Get data of all samples, select 112 Tumor + 112 Normal samples to and get DEGs.
* (2) ``thetaSelectGEDFN.R``  --  Use GCWs method get top 1% genes, repeat 10 times, make union.
* (3) ``malacards_GEDFN_mamaprint_KEGG.R``  --  Integrate data from MalaCards, KEGG, Mamaprint, GCWs, DEGs to union gene and corresponding expression data.
* (4) ``network_match_union.R``  --  Get the network of union gene in RegNetwork, extract the expression data of TCGA corresponding to union gene, and scale them.
* (5) ``data_splitnew.R``  --  According to the random seeds of other methods, the scaled data of the union gene of TCGA is divided into training data and testing data.
* (6) ``adj_union.R`` ---- Adjacency matrix and its eigenvalues.
* (7) ``cut_union.R`` ---- Diameters and cut-nodes of component of DEGs in RegNetwork.
<!--END_SECTION:news-->


## R codes for Result 
The **serial number (1), (2), ..., (4)** represents the order in which the program runs in our work. 
<!--START_SECTION:news-->
* (1) ``feature_select_all_new.R`` -- Extract the common genes of TCGA and GEO, using the identified 32 genes. 
* (2) ``class_net_Cox.R`` -- Train on TCGA data, predict on GEO data, apply linear Cox classifier for classification, observe results.
* (3) ``network_match_all_new.R`` -- Extract the net information of the biomarkers identified by each method.
* (4) ``ROCplot.R`` -- Plot ROC curves on independent datasets.
<!--END_SECTION:news-->


## Matlab codes
<!--START_SECTION:news-->
**Note**: If you encounter the error **Unrecognized function or variable 'cvx_begin'.**, it means the ``cvx`` package has not been installed.

- It needs download ``cvx`` package from [CVX Home Page](http://cvxr.com/cvx/download/), where *Standard bundles* for Linux, Mac and Windoes can be chosen. 

- Then, unzip the downloaded file to your chosen directory (e.g., D:\Applications\MATLAB\R2020a\toolbox). To install in Matlab, simply run the script at **Command Window** of Matlab:
```bash
cd 'D:\Applications\MATLAB\R2020a\toolbox\cvx'  # my mathlab path
cvx_setup
```

To evaluate CNet-Cox on *Data_train* and *Data_test*, simply run the provided script:
```bash
cd CNet-Cox/Code/Matlab
matlab  # activate matlab R2020a or R2020b
matlab -r Coxmain.m
```
* (1) ``Coxmain.m`` -- main function.
* (2) ``costFunctionCox.m`` -- Objective function.
* (3) ``cvCox.m`` -- Cross validation to select optimal parameters.
* (4) ``cvCoxexample.m`` -- Cross validation to select optimal parameters - an example.
* (5) ``Laplcian_Matrix.m`` -- Laplacian matrix according to the adjacency matrix.
* (6) ``LogitisLapCox.m`` -- LogitisLap function for CV.
* (7) ``LogitisLapCoxexample.m`` -- LogitisLap function for CV - an example.
* (8) ``SGNLR.m`` -- SGNLR function for LogitisLap.
* (9) ``getLambMaxCox.m`` -- getLambMax function for cv.
* (10) ``PredictSVM.m`` -- Predict function on test dataset.
* (11) ``plotROC.m`` -- Roc curve function on test dataset.
* (12) ``confusion.m`` -- Confusion matrix.
* (13) ``printConMat.m`` -- Output confusion matrix.
* (14) ``RiskMatrix.m`` -- Risk mmatrix constuction.
* (15) ``RiskMatrixexample.m`` -- Risk mmatrix constuction - an example.
<!--END_SECTION:news-->


## Python codes
<!--START_SECTION:news-->
* (1) ``Fig1A_Radar_figure.ipynb`` -- Plot the radar figure.
* (2) ``Spatial_validation_externaldata.ipynb`` -- Validate PRS gene on spatial dataset and compare PRS score with MammaPrint and OncotypeDX.
<!--END_SECTION:news-->