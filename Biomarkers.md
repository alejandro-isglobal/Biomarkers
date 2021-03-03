Statistical Analysis of Biomarkers
========================================================
author: Alejandro C�ceres, PhD
date: 
autosize: true

Biomarkers: Statistical session
 

<p style="font-size: 70px"> Allmirall 15th May 2021 </p>

Presentation
========================================================
  
- <b>Senior statistician</b> in the Bioinformatics group at Barcelona Institute of Global Health <https://www.isglobal.org/>
  
- Adjunct <b>lecturer</b> in statustics at the Universitat Politectica de Catalunya <https://eebe.upc.edu/es>
  
- Over 13 years of experience analysing biomedical data
- Develop novel analysis methods for <b>biomarker discovery</b>
- High dimensional data including imaging and omic data: genomic, transcriptomic, exposomic, etc.
- I write scientific articles and implement methods in software packages (R/Matlab).


Presentation
========================================================
You can find me at: 
  
- [linkedin](https://es.linkedin.com/in/alejandro-caceres-dominguez-7449aa176)

- [google scholar](https://scholar.google.es/citations?user=s1D-6WAAAAAJ&hl=es)

- [gitHub](https://github.com/alejandro-isglobal)

- [my blog](https://alejandro-isglobal.github.io/)



Presentation
========================================================

Analytical validation of a biomarker: Functional magnetic resonance imaging

<img src="./Biomarkers-figure/fmri.png" style="width:75%"  align="center">
  
  
  
Presentation
========================================================

Discovery of severity Biomarker: Chromosome Y function and risk of disease in men 

<img src="./Biomarkers-figure/jnci.png" style="width:75%"  align="center">
  
  
Aim of the talk
========================================================
  
- What is the utility of surveying thousand/million biomarkers in drug development?
  
- How is biomarker utility assessed in terms of statistical evidence?
  
  
Summary
========================================================
  
1. Introduction:
  - Biomarker definition
- Types of biomarkers (context of use)
- Guidance for qualification and reporting (FDA)
2. Statistical evidence of biomarker utility in assessing treatment efficacy:
  - Biomarker's sensitivity and specificity
  - Regression analyses of biomarker levels on efficacy (stratified and with interactions by treatment)
   - Multiple/Composite Biomarkers
3. Examples:
   - Assessment of HIV antiviral resistance with a composite genomic biomarker (ROC curve, pubmed: https://pubmed.ncbi.nlm.nih.gov/12060770/)
   - Prediction of Brodalumab treatment on psoriasis area severity index using gene expression data. (https://pubmed.ncbi.nlm.nih.gov/31883845/)


Introduction: Biomarker definition
========================================================

Definition:

- A biomarker is a <b>biological measurement</b> that has the potential to inform **decision-making** in relation to a clinical treatment or intervention. 

Properties:

- Source material or matrix 

- Method of measurement

- Porpuse

Introduction: Biomarker definition
========================================================

Source material or matrix:

- specific analyte (e.g., cholesterol)
- anatomic feature (e.g., joint angle)
- physiological characteristic (e.g., blood pressure) 

when they are composite

- how the compoenents are interrelated (e.g. genomic data)
- how is the process of obtaining the biomarker (e.g. algorithm, score)


Introduction: Biomarker definition
========================================================

Method of measurement:

- molecular
- histologic
- radiographic
- physiologic characteristic
- composite (algorithmic)

New biomarkers are continously created as novel methods of measurements and their analysis are constantly developed for adressing outstanding disease-related or treatment-related <b>needs</b>.


Introduction: Types of biomarkers 
========================================================

How biomarkers inform decision making during drug development?

Context of use (COU) is a description of the biomarker's specific use in drug development:
  
  1. Category:
  
  - Disease related: (Who should be studied?)
- Diagnostic (selection)
- Prognostic (stratification)
- Susceptibility/risk (enrich)

- Treatment related: 
  - Predictive (Who should be treated?)
- Safety (Should we stop treaement?)
- Monitoring (Should we continue with treatment?)
- Pharmacodynamic/respose biomarker (What is the outcome of treatment? surrogate endpoint)

Introduction: Types of biomarkers 
========================================================
  | Role        | Description   | 
  | ----------  | ------------- | 
  | Diagnosis of a disease | To make a diagnosis more reliably, more rapidly, or more inexpensively than available methods |
  | Severity assessment | To identify subgroup of patients with a severe form of a disease associated with an increased probability of death or severe outcome
| Risk assessment | To identify subgroup of patients who may experience better (or worse) outcome when expose to an intervention
| Prediction of drug effects | To identify the pharmacological response of a patient exposed to a drug (efficacy, toxicity, and pharmacokinetics) 
| Monitoring | To assess the response to a therapeutic intervention


Introduction: Types of biomarkers
========================================================
  
Context of use (COU) is a description of the biomarker's specific use in drug development:

2. Use:
  - Porpuse of use in drug developmnet (safety biomarker to evaluate drug-induced injury)
  - Stage of use (phase 1 clinical trial)
  - Population (healthy adults, psoriasis patients)
  - Therapeutic mechanisms of action for which the biomarkers offer information.


Introduction: Guidance for qualification (FDA)
========================================================

Qualification of a biomarker is a determination the biomarker can be relied on to have a specific interpretation and application in drug development and regulatory review. 

Biomarker Qualification: Evidentiary Framework [doc](https://alejandro-isglobal.github.io/teaching/docs/fdabiomarker.pdf)

For using a biomarker:

- what has a biomarker been qualified for?  

For developing abiomarker:

- What are the evidentiary requirements for demonstrating the ulitily of the biomarker?



Introduction: Guidance for qualification (FDA)
========================================================

FDA Guidance for biomarker qualification (2018) for section 507 of the Federal Food, and Cosmetic Act (FD&C Act):

A. Define the **need**

 - i.e disease needs and added value of the biomarker to drug development)

B. Define the context of use **COU** 

 - i.e Prediction: identify patients that respond to treatment

C. Assessment of **benefits** and **risks** 

 - i.e. Benefit: Don't treat a patient who won�t benefit (specificity), treat a patient who will benefit (sensitivity). 
- i.e. Risk: not treat a patient who could benefit, treat a patient who won�t benefit. 



D. Determining evidence that is **statistically** sufficient to support COU


Statistical evidence 
========================================================
  
D. Determining evidence that is **statistically** sufficient to support COU

1. Analytical considerations: Is the test reliable?
  
  
- validation of the Biomarkers test�s technical performance
- cost-effectives, feasability
- assessment of measurement error 

Gene expression biomarkers:
  
- RNAseq experiments have are experimentally validated, ready to include in a clinical trial 

- Composite biomarkers derived from them (algorithm) need to be validated

Statistical evidence 
========================================================
  
2. Medical considerations: Is the biomarker clinically useful biomarker?
  
Establish the relationship between a biomarker and an outcome of interest:
  
- Randomized controlled trial
- Single-arm/historical control trial
- Cohort study
- Case-control study (including nested)
- Cross-sectional study
- Case series or case reports
- Registry information
- Meta-analysis

Strongest evidence comes from <b>prospective studies</b> that are specifically designed but data from studies conducted for <b>other purposes</b> can be used to support biomarker utility.


Statistical evidence 
========================================================
  
The aim is to provide statistical evidence for the correlation between the biomarker and the outcome according to the COU:
  
- Predictive biomarker: correlation between the biomarker and efficacy of treatment.  


Sensitivity and specificity 
========================================================
  
Suppose we have a collection of **treated** individuals and for each subject 

- we have measured treatment efficacy as a response variable
- we have a run a test to detect the presence of a biomarker

The Response measurement is dichotomous and has the events:
  
- yes (the patient responded to treatment)
- no (the patient did not respond to treatment)

The Biomarker measurement is dichotomous (dichotomized by a cutoff) and has the events:
  
- positive (the biomarker was detected)
- negative (the biomarker was not detection detected)

Sensitivity and specificity 
========================================================
  
Each individual has two measurements: (Response, Biomarker)

|  Subject  |  Response  |  Biomarker  |
| ------------- | ------------- | ---------- |
| $s_1$         |   yes        | positive |
| $s_2$         |   no         | negative |
| $s_3$         |   yes        | positive |
|...            |   ...        | ...      |
| $s_i$         |   no         | positive* |
|...            |   ...        | ...      |
|...            |   ...        | ...      |
| $s_3$         |   yes        | negative* |
|...            |   ...        | ...      |
  
  
Sensitivity and specificity 
========================================================
  
Let's think first in terms of the response

Within those who responded to treatment (yes), how many detected with the biomarler (positive)?

-<b>Sensitivity</b> (true positive rate)

$$fr(positive|yes)=\frac{n_{positive|yes}}{n_{negative|yes}+n_{negative|yes}}$$


Sensitivity and specificity 
========================================================

Let's think first in terms of the response

Within those who did not respond to treatment (no), how many were not detected with the biomarker (negative)?
  
- <b>specificity</b> (True negative rate)

$$fr(negative|no)=\frac{n_{negative|no}}{n_{positive|no}+n_{negative|no}}$$
  
  
Sensitivity and specificity 
========================================================
  
|  | Response: Yes |  Response: No |
| --------- | --------- | -------- |
| <b>Biomarker: positive</b> | $fr(positive|yes)$ | $fr(positive|no)$ | 
| <b>Biomarker: negative</b> | $fr(negative|yes)$ | $fr(negative|no)$ | 
| <b>sum</b>      | 1                | 1               |
  
  
|  | Response: Yes | Response: No |
| --------- | --------- | -------- |
| <b>Biomarker: positive</b> |True positive rate (<b>sensitivity</b>) | False positive rate| 
| <b>Biomarker: negative</b> | False negative rate| True negative rate (<b>specificity</b>)| 
| <b>sum</b>      | 1                | 1               |
  
Trade-off between sensitivity and specificity needs to be evaluated in the context of use and usefullness of the biomarker's test 

Sensitivity and specificity 
========================================================

Let's think first in terms of the test

<br />
Within those whose biomarker was detected (positive), how many responded to treatment (yes)?
  
- <b>Positive predictive value</b> 
  
$$fr(yes|positive)=\frac{n_{yes|positive}}{n_{yes|positive}+n_{no|positive}}$$
  
Sensitivity and specificity 
========================================================
  
Let's think first in terms of the biomarkers test

<br />
Within those whose biomarker was not detected (negative), how many did not respond to treatment (no)?

- <b>Negative predictive value</b> 

$$fr(yes|negative)=\frac{n_{yes|negative}}{n_{yes|negative}+n_{no|negative}}$$


Sensitivity and specificity 
========================================================

|  | Response: Yes  |  Response: No |  sum  |
| --------- | --------- | -------- | ------ |
| <b>Biomarker: positive</b> | $PPV: fr(yes|positive)$ | $fr(no|positive)$ | 1 |
| <b>Biomarker: negative</b> | $fr(yes|negative)$ | $NPV: fr(no|negative)$ | 1 |


PPV: positive predicted value
NPV: negative predicted value

These are really the values that we want to know. They depend on the probability of response to treatment. 


Sensitivity and specificity 
========================================================

There is a way to convert from sentitivity to positive predicted value (Baye's rule)

$$fr(yes|possitive)=\frac{fr(positive|yes)}{fr(possitive)}fr(yes)$$
  
which can be rewritten 

$$\frac{fr(yes|possitive)}{fr(no|possitive)}=\frac{fr(positive|yes)}{1-fr(negative|no)} \frac{fr(yes)}{1-fr(yes)}$$
  
Sensitivity and specificity 
========================================================
  
$$ODD_{posttest}=LHR*ODD_{pretest}$$
    
$LHR+=\frac{Sensitivity}{1-specificity}$
      
|  | LHR+ |
| --- | --- |
| Excellent Efficacy | >10 |
| Good Efficacy | 5-10 |
| Poor Efficacy  | 1-5 |
| No Efficacy  | 1 |

      
Regression analyses
========================================================
      
      
When the levels of the biomaker are continous then a regression analysis can be used to determine the association with the outcome (response to treatment).
    
The type of correlation depends on the COU.
    
Let�s consider the following [study](https://www.nature.com/articles/s41398-019-0521-7):
      
**Biomarkers for response in major depression: comparing paroxetine and venlafaxine from two randomised placebo-controlled clinical studies.** Carboni et al. *Translational psychiatry*. 2019
    
    
Regression analyses
========================================================
      
Features of the study 
    
Two placebo-controlled studies evaluating the efficacy and tolerability of a novel drug candidates.
    
- Two drug treatments: paroxetine or venlafaxine as active comparators
    
- panel of peripheral biomarkers (including IL-6, IL-10, TNF-a, TNFRII, BDNF, CRP, MMP9 and PAI1) in depressed patients receiving paroxetine, venlafaxine, or placebo
    
Aim: assess the correlation between **biomarker levels** and response outcome: 17 item scale of depression symptoms; responders >50% in reduction from baseline (reduction from 2 to 1 = reduction from 10 to 5)
    
Regression analyses
========================================================
      
Demographics:
      
<img src="./Biomarkers-figure/tab1.JPG" style="width:75%"  align="center">
      
Regression analyses
========================================================
      
Analysis 1.
    
Associations between biomarker levels and depression severity at **base-line**: which are sate biomarkers (CUO: diagnosis)? 
      
$$D_{base} \rightarrow B_{base}$$
      
- paroxetine sutudy: IL-6 (r=0.23, p=0.018), IL-10 (r=0.19, p=0.045), stratifying by sex no significant associations were found for females. 
    
    
- veroxine study: No signfinificnat correlations were found. 
The biomarkers do not show diagnostic capacity.
    
    
    
Regression analyses
========================================================
      
Analysis 2.
    
Associations between biomarkers' changes and changes in depression symtoms: Which are biomakers of treatment efficacy  (CUO: surrogate endpoints)? 

$$\Delta D=D_{w10}-  D_{base} \rightarrow \Delta B= B_{w10} - B_{base}$$

Adjusting for sex and $B_{base}$ and $D_{base}$ in **full** population



Regression analyses
========================================================

<img src="./Biomarkers-figure/tab2.JPG" style="width:75%"  align="center">

- TNF-a, IL-6, IL-10 and CRP significanlty reduced with $\Delta D$ in the paroxetine  study, none in the venlafaxine.

- IL-10 reduced with $\Delta D$ in males in both **studies** 

Regression analyses
========================================================

Analisis 3. 

<img src="./Biomarkers-figure/tab3.JPG" style="width:75%"  align="center">


Regression analyses
========================================================

Analysis 3.

Associations between changes in symptoms and biomarkers' levels at  baseline: Which biomarkers predict improvement in sypmtoms when treated (CUO: prognostic under treatment biomarker)?
      
a. Similar to the sensibility and specificity analysis we can condition (stratify) on treated only
    
$$B_{base|treated} \rightarrow \Delta D$$
      
Adjusting by $D_{base|treated}$ and sex
    
    
- For those treated with paroxetine: IL-10 and TNF-a are at baseline were significantly associated changes in depression symptoms at week 10. 
    
IL-10 and TNF-a showed predictive capacity under paroxetine treatment
    
    
Regression analyses
========================================================
      
Analysis 3.
    
b. Associations between changes in symptoms and biomarkers' levels at  baseline: Which biomarkers predict response when not treated (CUO: Prognosis biomarkers)?

$$B_{base|placebo} \rightarrow \Delta D}$$

Adjusting by $D_{base|treated}$ and sex

- CPR associated with improvement of symptoms in the veroxine study.

This indicates that improvement of symptoms when treated with paroxetine may not be due to placebo effects. 


Regression analyses
========================================================

Analysis 3.

c. Associations of symptoms changes and the interaction between baseline biomarkers' levels and treatment (CUO: Predictive biomarkers) 
    
$$B_{base|placebo}\times T \rightarrow \Delta D$$
      
Adjusting by $D_{base|treated}$ and sex
    
- For those treated with paroxetine: treatment intectiond with IL-10 and TNF-a showed trend to significance (P=0.054, P=0.085).
    
While testing for interactions requires more power, this suggests that individuals with low values of IL-10 will respond better to treatment than those with high values.
    
Where to set the threshold?
      
      
Multiple/Composite Biomarkers
========================================================
      
A main problem when analysing multiple biomarkers, all intependenlty is **multiplicity**
      
- Take a biomarker with no correlation with efficacy and test the correlation in 100 clinical trials: 5% studies with find significant results. 
    
- Take a 100 biomarkers with no correlation with efficacy in one clinical trial: 5% of biomakers will be declared significant. 
    
In any suh clinical trial is almost sure that will have at least one significant biomarker. 
    
- You want that only 5% of null trials report a significant biomarker. 
    
    
    

Multiple/Composite Biomarkers
========================================================
      
      
Correct threshold of significance (correction for multiple comparisons):
      
- Bonferroni: divide the P value by the number of biomarkers. In the Depression study then $P < 0.05/8 = 0.0062$: None of the results are significant!
      
- False discovery rate: Order the 8 Pvalues from lower to higher: $P_i$ ($i=1...8$) and select $i$ such that $P_i \leq i/8*0.05$. All P values between 0 and i are declared significant.      
    
Both methods are widely implemented in statistical software. Bonferrroni is more conservative than FDR, and FDR is most commonly used in omic studies.  
    
    
Multiple/Composite Biomarkers
========================================================
      
      
Another alternative is to **construct** a composite Biomarker. 
    
- Computational processed and/or algorithms using machine leaning and AI to discover a subset of  individuals where treatment effect is maximum. 
    
    
- Use the biomarkers to measure a new biological quantity that could be in the disease pathway. 
    
    
Multiple/Composite Biomarkers
========================================================
Genetic mosaicisms 
    
<img src="./Biomarkers-figure/mosaicism_med.jpeg" style="width:50%"  align="center">
      
One most common somatic mutation in man is the loss of chromosome Y.
    
In a sample we can measure the lost of RNA transciption from genes in chromosome Y from cells that do not produce RNA from chromosme Y becase they lost it.
    
- We convert 1000 of biomarkers into 1 with biological sense. 
    
- We have shown that the ost of trancription of chromosome Y is associates with cancer, BMI, lower immune cell count in blood.   
    
    
    
Examples
========================================================
      
If we have detected one biomarker with continuous levels that is associated with efficacy, how do we select the threshold for clinical applications?
      
[Diversity and complexity of HIV-1 drug resistance: a bioinformatics approach to predicting phenotype from genotype](https://pubmed.ncbi.nlm.nih.gov/12060770)
    
- Response:  Antiretroviral drug resistance
    
- Biomarker: Score that stratifies patients with genomic data using a machine learning method
    
ROC curve
========================================================
      

```r
library(RCurl)
hiv <- read.delim("https://alejandro-isglobal.github.io/data/hiv.txt")
head(hiv)
```

```
  response      test
1        1 -0.438185
2        1 -0.766791
3        1  0.695282
4        1 -0.689079
5        1  0.325977
6        1  0.704040
```
    
    
ROC curve
========================================================
      

```r
table(hiv$response)
```

```

 -1   1 
267  78 
```
    
- response: <b>no</b>, resistance to drug treatment: <b>1</b>
- response: <b>yes</b>, no resistance to drug treatment: <b>-1</b>
      
ROC curve
========================================================
      

```r
br <- seq(-2,2,0.25)
    
hist(hiv$test[hiv$response==-1], 
br=br, freq=F,xlab="RF", main="")
    
hist(hiv$test[hiv$response==1], 
br=br, freq=F, add=T, col="blue")
    
legend("toprigh",
       legend=c("yes=no resis.", "no=resis."),          
       col=c(1,2), lty=1)
```
    
    
ROC curve
========================================================
      
![plot of chunk unnamed-chunk-4](Biomarkers-figure/unnamed-chunk-4-1.png)
cutoff at $-0.5$
      
- Biomarker: <b>positive</b>: $test < -0.5$ 
- Biomarker: <b>negative</b>: $test > -0.5$ 
      
      
      
 ROC curve
 ========================================================
      
 $Biomarker: positive$ when there was no resistance (yes)
    
 Sensitivity: $fr_{[cutoff=-0.5]}(positive|yes)$
      

```r
  mean(hiv$test[hiv$response==-1] < -0.5 )
```

```
[1] 0.9438202
```
    
![plot of chunk unnamed-chunk-6](Biomarkers-figure/unnamed-chunk-6-1.png)
    
ROC curve
========================================================
      
$Biomarker: positive$ when there was no resistance (no)
    
False positive rate, $1 - specificity$: $fr_{[cutoff=-0.5]}(positive|no)$
      

```r
mean(hiv$test[hiv$response==1] < -0.5 )
```

```
[1] 0.2692308
```
    
    
    
![plot of chunk unnamed-chunk-8](Biomarkers-figure/unnamed-chunk-8-1.png)
    
    
ROC curve
========================================================
(1 - specificity, Sensitivity)$_{cutoff}$
</br>=(false positive rate,true positive rate)$_{cutoff}$=$(0.269, 0.943)_{-0.05}$
      

```r
library(cvAUC)
out <- cvAUC(-hiv$test, -hiv$response) #calcular ROC
plot(out$perf, col="blue", main="ROC") #plot
lines(c(0,1),c(0,1)); points(0.269, 0.943, pch=16) #cutoff=-0.5
```

![plot of chunk unnamed-chunk-9](Biomarkers-figure/unnamed-chunk-9-1.png)
    
ROC curve
========================================================
Area under the curve
    
$AUC=Pr(X2 < X1)$
      
Where $X1$ is the outcome of a positive test  and $X2$ the outcome of a negative test 
    

```r
ci.cvAUC(-hiv$test, -hiv$response)
```

```
$cvAUC
[1] 0.9047825

$se
[1] 0.02276935

$ci
[1] 0.8601554 0.9494096

$confidence
[1] 0.95
```
    
ROC curve
========================================================
      

```r
library(pROC)
rocobj <- roc(-hiv$response, -hiv$test)
coords(rocobj, "best", best.method="youden")
```

```
  threshold specificity sensitivity
1  0.700003   0.7948718   0.9138577
```
    
It optimizes $sensitivity-(1-specificity)$
      
ROC curve
========================================================
      
      

```r
library(cvAUC)
out <- cvAUC(-hiv$test, -hiv$response) #calcular ROC
plot(out$perf, col="blue", main="ROC") #plot
lines(c(0,1),c(0,1)); points(0.269, 0.943, pch=16)#cutoff=-0.5
points(1-0.794, 0.9138, pch=16, col="red")  #optimal
```

![plot of chunk unnamed-chunk-12](Biomarkers-figure/unnamed-chunk-12-1.png)
    
    
Prediction of treatment effects
========================================================
      
[Short-term transcriptional response to IL-17 receptor-A antagonism in the treatment of psoriasis](https://pubmed.ncbi.nlm.nih.gov/31883845/]). JACI. 2020
    
<img src="./Biomarkers-figure/abs.JPG" style="width:75%"  align="center">
      
      
      
Prediction of treatment effects
========================================================
      
- They used transcription data from a panel of genes associated with psoriasis. 
    
- They show that the improvement in  of psoriasis transcriptome with Brodalumbad treatment by responders. 
    
COU: supprogate end point. Improvement of psoriasis transcriptome showing causal action of the drug in a biological pathway. 
    
    
<img src="./Biomarkers-figure/paper.JPG" style="width:75%"  align="center">
      
      
Prediction of treatment effects
========================================================
      
The question of whether the biomarkers can be used to predict Brodalumbad remains. 
    
I downloaded the data from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117468) to find predictors of the efficacy. 
    



























```
Error in load(data, envir = envir, ...) : 
  n�mero m�gico de archivo de restauraci�n inv�lido (el archivo puede estar da�ado) -- ning�n dato recargado
```
