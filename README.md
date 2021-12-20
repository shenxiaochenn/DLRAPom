# DLRAPom
## Abstract
In terms of conventional bioinformatics analysis for screening disease-related differentially expressed non-coding RNAs, large number of results are usually obtained, and most of them are difficult to be verified by molecular biological experiments. Finding a reliable target for follow-up disease mechanism research is just like finding a needle in a haystack. Although many advanced bioinformatics analysis methods have been reported in recent years, it is very difficult for most basic medical researchers to understand and use these methods proficiently. Most of these reported methods focused on one or two types of non-coding RNA, and there were few predictions combined with experimental validation to identify the lncRNA-miRNA-mRNA network as a whole functional module. The lack of a reliable and easy-to-operate screening pipeline for disease-related non-coding RNA regulatory axis is a problem that needs to be solved urgently. To address this, we designed a hybrid pipeline, disease-related lncRNA-miRNA-mRNA regulatory axis prediction from multiomics (DLRAPom), to identify risk biomarkers and disease-related lncRNA-miRNA-mRNA regulatory axes by adding a novel machine learning model on the basis of conventional analysis and combining experimental validation. The pipeline consists of four parts, including selecting hub biomarkers by conventional bioinformatic analysis, discovering the most essential protein-coding biomarkers by a novel machine learning model, extracting the key lncRNA-miRNA-mRNA axis, and validating experimentally. In the pipeline, a novel machine learning algorithm, Optimized-XGBoost, was developed and experimental validation was utilized as the key component to ensure the accuracy of extracting regulatory features from the ncRNA-gene-disease association network. Compared with the methods reported previously, our pipeline is distinguished from those methods, not only by adopting the above-mentioned machine learning model but also by verifying the predicted results based on biological experiments. In short, as a flexible pipeline, DLRAPom can contribute to molecular pathogenesis research of diseases, effectively predicting potential disease-related noncoding RNA regulatory networks and providing promising candidates for functional research on disease pathogenesis.
![image](https://github.com/shenxiaochenn/DLRAPom/blob/main/summary_figure.jpg)
### As indicated in the above summary figure
The pipeline consists of four steps, namely, selecting hub biomarkers by conventional bioinformatic analysis, discovering the most essential protein-coding biomarkers by a novel machine learning model, extracting the key lncRNA-miRNA-mRNA axis, and validating experimentally. In the Step1, hub biomarkers are selected by conventional bioinformatics analyses. These results are used as the inputs of the Step2 to discover essential protein coding biomarkers and obtain the importance of each protein coding biomarker. In the Step3, the competing endogenous network is constructed based on the obtained information of lncRNA-miRNA and miRNA-target gene by databases. Among all the constructed regulatory axes, the regulatory axes containing the predicted risk protein coding biomarkers in the novel Optimized-XGBoost model are selected as the main outcomes of our pipeline and would be used for subsequent experimental verification. If there are multiple regulatory axes, the criticality of the regulatory axes is ranked in descending order according to the importance of the predicted protein coding biomarker included in each axis. After the significant expression change of each RNA molecule in the predicted regulatory axes is confirmed, further supportive evidence for the pairwise targeting relationships within the predicted regulatory axes were required. If these targeting relationships have not been reported before, we need to determine whether these targeting relationships exist through the dual-luciferase reporter assay. For a predicted regulatory axis, only when the biological targeting relationships of lncRNA-miRNA and miRNA-mRNA have been both experimentally verified, can the predicted regulatory axis be considered to be targetable and reliable.
