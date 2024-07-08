# STPathway -- A Robust Statistical Approach for Finding Informative Spatially Associated Pathways
Spatial transcriptomics offers deep insights into cellular functional localization and communication by mapping gene expression to spatial locations. Traditional methods focus on selecting spatially variable genes, which often miss the complexity of biological pathways and network dynamics. Here, we introduce a novel framework, STPathway, that shifts the focus towards directly identifying functional pathways associated with spatial variability. This method adapts the Brownian distance covariance test in an innovative manner to explore the heterogeneity of biological functions over space. Unlike most other methods, this statistical testing approach is free of gene selection and parameter selection, allowing for the detection of nonlinear and complex dependencies.

![workflow](https://github.com/tianlq-prog/STpathway/blob/main/figure/frame_loc.pdf)

# Tutorial 

For the step by step tutoral, please refer to the notebook:

https://github.com/tianlq-prog/STpathway/blob/main/PDAC_example.ipynb

In the example we utilize the PDAC data as example which contains

```pythonscript
data_preprocessing(pos=pos, neg=neg, 
                       pos_adductlist=["M+H","M+NH4","M+Na","M+ACN+H","M+ACN+Na","M+2ACN+H","2M+H","2M+Na","2M+ACN+H"], 
                       neg_adductlist = ["M-H","M-2H","M-2H+Na","M-2H+K","M-2H+NH4","M-H2O-H","M-H+Cl","M+Cl","M+2Cl"], 
                       idx_feature = 4, match_tol_ppm=5, zero_threshold=0.75, log_transform=True, scale=1000)
```
