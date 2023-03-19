# Prevalence Mapping for Cancer in WA State
Name: **Walker Azam**

*STAT/CSSS 554*

*Instructor: Dr. John Wakefield (UW Dept. of Biostatistics)*



# Project Background

## Introduction

In this project I will be mapping and predicting Lung Cancer mortality in WA state in male populations. Lung cancer remains one of the most common forms of cancer to develop in both men and women. According to the Washington State Cancer Registry, Lung and Bronchus cancer had the highest mortality count, making up 21.7% of all cancer related mortalities in 2018, despite a 10.5% incidence count. It is also shown that men have higher age-specific rates of lung cancer, often due to difference in smoking habits. Given this I believe performing disease mapping for Lung and Bronchus cancer in Washington State can find regions of high risk, and can help plan public health interventions, such as encouraging early screening for Lung Cancer in men. 

I will specifically be looking at cancer deaths in males, per county, aggregate counts from 2015-2019. The specific site group for cancer as defined by SEER (The National Cancer Institute’s Surveillance Epidemiology and End Results) is Lung & Bronchus (Recode 22030). The data was collected from the Washington Tracking Network (more information in the **Data Background** section). The primary goal of this project is to find suitable models for disease mapping, and to provide discussion on regions of high risk, the validity of the models selected, and covariate analysis with health accessibility metrics. As aforementioned, health outcome's of smoking have been well established in regards to lung cancer, but I was particularly interested in whether health accessibility also had associations with lung cancer morbidity. In particular, this could inform whether Washington state resources could be allocated better to regions of high risk.

## Methodology

To calculate the expected cases of death from Lung & Bronchus Cancer in male populations, the Age Adjusted Rate per 100,000 was utilized, which was also provided by the Washington State Department of Health's WTN. The age-adjusted rate was derived from The National Cancer Institute’s Surveillance Epidemiology and End Results (SEER). 4 model variations will be compared for mapping prevalence of male lung cancer mortality in Washington state, and evaluated for appropriateness.These models will be used to indicate counties that appear to have the highest risk for male lung cancer, for which resources and public health planning can allocate resources or possible interventions. The models include Standard Morbidity Ratio (SMR), IID models, spatial (BYM2) models, and spatial models with covariates. 

For the full report, results, and data background please refer to the `docs` folder, where a pdf copy of the paper and supplementary materials can be found.
