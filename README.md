## Prognostic Significance of Cyclin A1 Expression in pT1 Urothelial Carcinoma of Bladder: A Tissue Microarray Study of 149 Patients Treated by Transurethral Resection

### Description
This repository contains the data analysis report, scripts, and figures of a study evaluating the prognostic significance of cyclin A1 expression in pT1 urothelial carcinomas of the urinaty bladder.

**Main Files**  
* [CA1Results.md](https://github.com/alcideschaux/CyclinA1-pT1Bladder/blob/master/Results/CA1Results.md): this file contains the results of the data analysis.  
* Figures (in .png format) are located in the [figure](https://github.com/alcideschaux/CyclinA1-pT1Bladder/tree/master/Results/figure) folder.  
* [CodeBook.md](https://github.com/alcideschaux/CyclinA1-pT1Bladder/blob/master/CodeBook.md) contains the codebook naming the variables with their labels and descriptions.  
* [ca1.csv](https://github.com/alcideschaux/CyclinA1-pT1Bladder/blob/master/ca1.csv) corresponds to the dataset that was used for data analysis.  

### Material & Method
**Selection of Patients and Building of Tissue Microarrays**  
The study was approved by the appropriate institutional review boards. Specimens from 149 consecutive patients with pT1 bladder urothelial carcinoma treated by TURB followed by BCG and or intravesical mitomycin therapy between 2002 to 2009. Clinical charts and pathology reports were retrieved and glass slides were reviewed in all cases. Patients had no previous history of urothelial lesions elsewhere. Representative areas from formalin-fixed, paraffin-embedded (FFPE) tissue blocks, were cored and spotted in triplicate to build 5 tissue microarrays (TMA).

**Predictor and Outcome Variables and Patients' Follow-Up**  
Predictor variables included the following: 1) pT stage at biopsy defined as tumor confined to muscolaris mucosae (pT1a) or tumor invading muscolaris mucosae or beyond (pT1b); 2) carcinoma in situ (CIS) at biopsy (present, absent); 3) initial treatment following biopsy (TURB alone, TURB + BCG, TURB + mitomycin C); 4) Immunohistochemical expression of cyclin A1 (see below).

Outcome variables included tumor recurrence and tumor progression as the primary endpoints. Tumor recurrence was defined as subsequent development of histologically proven bladder urothelial neoplasm(s) of ≥ pTa stage. Tumor progression was defined as the development of any ≥ pT2 bladder urothelial lesion and or systemic metastases. In cases with more than one recurrence, the highest pT stage was used as the endpoint.

For outcome analysis, “time to event” refers to the time length in months from the initial treatment till the event (i.e., tumor recurrence/tumor progression) or the end of the study. Patients in whom the event was observed were followed-up afterward to monitor the appearance of new lesions or progression of the previous ones. Patients in whom the event was not observed by the end of the study were considered as “censored”. Patients were followed-up for a mean of 23.6 months, a median of 18 months, and a range of 3 to 108 months. No patients were lost at follow-up.

**Immunohistochemistry for Cyclin A1 and Scoring System** 
Immunohistochemistry was performed on 5-micron formalin-fixed paraffin-embedded tissue microarray sections. Antigen retrieval was performed using sodium citrate buffer (pH 6.0) for 30 min, followed by endogeneous peroxidase block with 0.3% H2O2 in phosphate-buffer saline and primary antibody incubation overnight at 4˚C with anti-cyclin A antibody (Abcam, clone 6E6, dilution 1:100, overnight at 4°C). The primary antibody was later detected using the PowerVision+ Poly-HRP IHC detection system (Leica Biosystems), as per manufacturer’s instructions.  After incubation with detection reagents, Diaminobenzidine (DAB) served as a chromogen, and Harris hematoxylin was used as a counterstain. FFPE tonsil tissue was used as positive control.

Cyclin A1 expression was evaluated as the percentage of tumor cells with positive nuclear staining estimated at each TMA spot. Then, cyclin A1 expression was estimated for each case using the pooled arithmetic mean of cyclin A1 expression of the 3 TMA spots. Cases with >10% of positive tumor cells were classified as “Positive cyclin A1 expression” and with <= 10% as “Negative cyclin A1 expression”. Seven cases were excluded from the study due to technical failures in Cyclin A1 staining.

**Statistical Analysis**  
Associations between variables were evaluated using the Mann-Whitney U test or the Fisher's exact test, as appropriate. Time-to-event analyses were done using the Kaplan-Meier method, and survival curves were compared using the log-rank (Mantel-Cox) test. Proportional hazard risks for selected endpoints were estimated using Cox regression analyses. For comparison of predictive models the log-likelihood (LL) for each model was estimated and the log-likelihood ratio test was used for the comparison. For statistical inference, alpha was set at 0.05. Significance levels of P values were adjusted using Hommel's correction in order to control the family-wise error rate. Data were analyzed using R version 3.1.0 "Spring Dance" (R Foundation for Statistical Computing, Vienna, Austria).

### Authors
Enrico Munari (1); Alcides Chaux (1,4); Leonel Maldonado (6); Eva Compérat (5); Trinity J. Bivalacqua (2,3); Mohammad O. Hoque (2,3,6); and George J. Netto (1,2,3)

(1,2,3) Departments of Pathology, Urology & Oncology, Johns Hopkins Medical Institutions, Baltimore, MD  
(4) Office of Scientific Research, Norte University, Asunción, Paraguay  
(5) Pathology Department, Hôpital Pitié-Salpêtrière, Paris, France  
(6) Department of Otolaringology and Head & Neck Surgery, Johns Hopkins Medical Institutions, Baltimore, MD.

Enrico Munari and Alcides Chaux contributed equally to this work.