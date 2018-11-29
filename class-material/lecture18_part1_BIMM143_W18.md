---
title: "BIMM-143, Lecture 18"
subtitle: "Investigating cancer genomics datasets"
output:
    html_document:
        keep_md: true
        theme: cosmo
        highlight: pygments
        toc: true
        toc_float: true
        toc_depth: 2
        number_sections: false
        df_print: paged
---


<style> 
h1, .h1 {
    margin-top: 50px;
}
h2, .h2, h3, .h3 {
    margin-top: 30px;
}
body{
  font-size: 13pt;
}
.q_box { 
  display: block;
  border: 1px solid Gray; 
  border-radius: 8px;
    width: 90%;
    left: 6px;
    padding: 5px 5px 5px 15px;
    color: Red;
    margin-bottom: 25px;
}

</style> 





**BIMM-143 Lecture 18:**  
Barry Grant &lt; <http://thegrantlab.org/bimm143/> &gt;  
2018-11-28   (13:43:34 PST on Wed, Nov 28)  
  

# Overview

Cancer is fundamentally a disease of the genome, caused by changes in the DNA, RNA, and proteins of a cell that push cell growth into overdrive. Identifying the genomic alterations that arise in a given cancer can help researchers decode how a particular cancer develops and improve upon the diagnosis and treatment of cancers based on their distinct molecular abnormalities.  

<img align="right" src="https://bioboot.github.io/bggn213_S18/class-material/head_up_cancer.png">

With the ability to sequence whole genomes and exomes, attention has turned to trying to understand the full spectrum of genetic mutations that underlie cancer.  

The genomes or exomes of tens of thousands of cancers have now been sequenced. Analyzing this data can yield important new insights into cancer biology.  This is important because it is estimated that cancer will strike 40% of people at some point in their lifetime with frequently devastating effects.  


# 1. The NCI Genomic Data Commons
The National Cancer Institute (NCI) in the US has established the [**Genomic Data Commons**](https://gdc.cancer.gov/about-gdc) (or **GDC** for short) for sharing cancer genomics data-sets.  

This includes data from the large scale projects such as **Cancer Genome Atlas** (TCGA) and other projects.  The TGCA project aims to generate comprehensive, multi-dimensional maps of the key genomic changes in major types and sub-types of cancer. As of writing, TCGA has analyzed matched tumor and normal tissues from over 11,000 patients covering 33 cancer types and sub-types.

You can get a feel for the types of cancer data contained in the NCI-GDC by visiting their web portal: [https://portal.gdc.cancer.gov](https://portal.gdc.cancer.gov).

## Exploring the GDC online
Visit the NCI-GDC web portal and enter p53 into the search box.

<div class="q_box">
> **Q1**. How many *Cases* (i.e. patient samples) have been found to have p53 mutations?

> **Q2**. What are the top 6 misssense mutations found in this gene? <br> **HINT:** Scroll down to the 'TP53 - Protein' section and mouse over the displayed plot. For example **R175H** is found in 156 cases. 

> **Q3**. Which domain of the protein (as annotated by PFAM) do these mutations reside in?

> **Q4**. What are the top 6 *primary sites* (i.e. cancer locations such as Lung, Brain, etc.) with p53 mutations and how many *primary sites* have p53 mutations been found in? <br> **HINT:** Clicking on the number links in the *Cancer Distribution* section will take you to a summary of available data accross *cases*, *genes*, and *mutations* for p53. Looking at the *cases* data will give you a ranked listing of *primary sites*.
</div>



Return to the NCI-GDC homepage and using a similar search and explore strategy answer the following questions:

<div class="q_box">
> **Q5**. What is the most frequentely mutated position associated with cancer in the **KRas** protein (i.e. the amino acid with the most mutations)?  

> **Q6**. Are KRas mutations common in Pancreatic Adenocarcinoma (i.e. is the Pancreas a common '*primary site*' for KRas mutations?).

> **Q6**. What is the '*TGCA project*' with the most KRas mutations? 

> **Q7**. What precent of cases for this '*TGCA project*' have KRas mutations and what precent of cases have p53 mutations? <br> **HINT:** Placing your mouse over the project bar in the **Cancer Distribution** panel will bring up a tooltip with useful summary data.  

> **Q8**. How many TGCA Pancreatic Adenocarcinoma *cases* (i.e. patients from the TCGA-PAAD project) have RNA-Seq data available?
</div>

<br>

<img align="right" src="ras_mutants_pov.png">


> **Side-Note**: If Barry forgets please remind him to show the top miss-sense mutation sites on the protein structure (and discuss mechanism) as well as demo how to find RNA-Seq FASTQ files and other data such as biopsy images etc.   



By now it should be clear that the NCI-GDC is a rich source of both genomic and clinical data for a wide range of cancers. For example, at the time of writing there are 5,306 files associated with Pancreatic Adenocarcinoma and 14,278 for Colon Adenocarcinoma. These include RNA-Seq, WXS (whole exome sequencing), Methylation and Genotyping arrays as well as well as rich metadata associated with each case, file and biospecimen.  


# 2. The GenomicDataCommons R package

The [GenomicDataCommons](https://bioconductor.org/packages/release/bioc/html/GenomicDataCommons.html) Bioconductor package provides functions for querying, accessing, and mining the NCI-GDC in R. Using this package allows us to couple large cancer genomics data sets (for example the actual RNA-Seq, WXS or SNP data) directly to the plethora of state-of-the-art bioinformatics methods available in R. This is important because it greatly facilitates both targeted and exploratory analysis of molecular cancer data well beyond that accessible via a web portal.

This section highlights how one can couple the [GenomicDataCommons](https://bioconductor.org/packages/release/bioc/html/GenomicDataCommons.html) and [maftools](https://bioconductor.org/packages/release/bioc/html/maftools.html) bioconductor packages to quickly gain insight into public cancer genomics data-sets.

We will first use functions from the `GenomicDataCommons` package to identify and then fetch somatic variant results from the NCI-GDC and then provide a high-level assessment of those variants using the `maftools` package. The later package works with [Mutation Annotation Format](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification) or **MAF** format files used by GDC and others to store somatic variants. 


The workflow will be:

- Install packages if not already installed
- Load libraries
- Identify and download somatic variants for a representative TCGA dataset, in this case pancreatic adenocarcinoma.
- Use maftools to provide rich summaries of the data.


```r
source("https://bioconductor.org/biocLite.R")
biocLite(c("GenomicDataCommons", "maftools"))
```

Once installed, load the packages, as usual.


```r
library(GenomicDataCommons)
library(maftools)
```

Now lets check on GDC status:


```r
GenomicDataCommons::status()
```

```
## $commit
## [1] "acaf65369f6cea8337c6b59f0888761c9ed34654"
## 
## $data_release
## [1] "Data Release 13.0 - September 27, 2018"
## 
## $status
## [1] "OK"
## 
## $tag
## [1] "1.17.1"
## 
## $version
## [1] 1
```

If this statement results in an error such as `SSL connect error`, then please see the [troubleshooting section here](https://bioconductor.org/packages/release/bioc/vignettes/GenomicDataCommons/inst/doc/overview.html#ssl-connection-errors).  


# 3. Querying the GDC from R

We will typically start our interaction with the GDC by searching the resource to find data that we are interested in investigating further. In GDC speak this is called *"Querying GDC metadata"*. Metadata here refers to the extra descriptive information associated with the actual patient data (i.e. 'cases') in the GDC. 

> **For example**: Our query might be '**find how many patients were studied for each major project**' or '**find and download all gene expression quantification data files for all pancreatic cancer patients**'.  We will answer both of these questions below.


The are four main sets of metadata that we can query with this package, namely `cases()`, `projects()`, `files()`, and `annotations()`. We will start with `cases()` and use an example from the package associated [publication](https://www.biorxiv.org/content/biorxiv/early/2017/04/04/117200.full.pdf) to answer our first question above (i.e. find the number of cases/patients across different projects within the GDC): 


```r
cases_by_project <- cases() %>%
  facet("project.project_id") %>%
  aggregations()
head(cases_by_project)
```

```
## $project.project_id
##               key doc_count
## 1           FM-AD     18004
## 2      TARGET-NBL      1127
## 3       TCGA-BRCA      1098
## 4      TARGET-AML       988
## 5       TARGET-WT       652
## 6        TCGA-GBM       617
## 7         TCGA-OV       608
## 8       TCGA-LUAD       585
## 9       TCGA-UCEC       560
## 10      TCGA-KIRC       537
## 11      TCGA-HNSC       528
## 12       TCGA-LGG       516
## 13      TCGA-THCA       507
## 14      TCGA-LUSC       504
## 15      TCGA-PRAD       500
## 16   NCICCR-DLBCL       489
## 17      TCGA-SKCM       470
## 18      TCGA-COAD       461
## 19      TCGA-STAD       443
## 20      TCGA-BLCA       412
## 21      TARGET-OS       381
## 22      TCGA-LIHC       377
## 23      TCGA-CESC       307
## 24      TCGA-KIRP       291
## 25      TCGA-SARC       261
## 26      TCGA-LAML       200
## 27      TCGA-ESCA       185
## 28      TCGA-PAAD       185
## 29      TCGA-PCPG       179
## 30      TCGA-READ       172
## 31      TCGA-TGCT       150
## 32      TCGA-THYM       124
## 33      TCGA-KICH       113
## 34       TCGA-ACC        92
## 35      TCGA-MESO        87
## 36       TCGA-UVM        80
## 37      TARGET-RT        75
## 38      TCGA-DLBC        58
## 39       TCGA-UCS        57
## 40      TCGA-CHOL        51
## 41    CTSP-DLBCL1        45
## 42    TARGET-CCSK        13
## 43 VAREPOP-APOLLO         7
```

Note that the **facet()** and **aggregations()** functions here are from the `GenomicDataCommons` package and act to group all cases by the project id and then count them up. 

If you use the **View()** function on our new `cases_by_project` object you will find that the data we are after is accessible via `cases_by_project$project.project_id`. 

<div class="q_box">
> **Q9**. Write the R code to make a barplot of the cases per project. Lets plot this data with a log scale for the y axis (`log="y"`), rotated axis labels (`las=2`) and color the bar coresponding to the TCGA-PAAD project.


```r
x <- cases_by_project$project.project_id

# Make a custom color vector for our plot
colvec <- rep("lightblue", nrow(x))
colvec[___] <- "red"

# Plot with 'log' for y axis and rotate labels with 'las'
#par(___)  
barplot(___, names.arg=___, log="y", col=colvec, las=2)
```


![](lecture18_part1_BIMM143_W18_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

</div>

Lets take another snippet of code from their package vignette and adapt it to answer our second question from above - namely '*find all gene expression data files for all pancreatic cancer patients*':



```r
## This code snipet is taken from the package vignette
file_records <- files() %>%
  filter(~ cases.project.project_id == "TCGA-PAAD" &
    data_type == "Gene Expression Quantification" &
    analysis.workflow_type == "HTSeq - Counts") %>%
  response_all()
```

The above R code was copied from the package documentation and all we have changed is the addition of our 'TCGA-PAAD' project to focus on the TGCA pancreatic cancer project related data. We will learn more about were the various names and options such as `analysis.workflow_type` etc. come from below (basically they are from the CDG database itself!).

In RStudio we can now use the **View()** function to get a feel for the data organization and values in the returned `file_records` object.  


```r
View(file_records)
```

We should see that `file_records$results` contains a row for every RNA-Seq data file from the 'TCGA-PAAD' project. At the time of writing this was 182 RNA-Seq data files.


```r
nrow(file_records$results)
```

```
## [1] 182
```


```r
head(file_records$results)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["data_type"],"name":[1],"type":["chr"],"align":["left"]},{"label":["updated_datetime"],"name":[2],"type":["chr"],"align":["left"]},{"label":["file_name"],"name":[3],"type":["chr"],"align":["left"]},{"label":["submitter_id"],"name":[4],"type":["chr"],"align":["left"]},{"label":["file_id"],"name":[5],"type":["chr"],"align":["left"]},{"label":["file_size"],"name":[6],"type":["int"],"align":["right"]},{"label":["id"],"name":[7],"type":["chr"],"align":["left"]},{"label":["created_datetime"],"name":[8],"type":["chr"],"align":["left"]},{"label":["md5sum"],"name":[9],"type":["chr"],"align":["left"]},{"label":["data_format"],"name":[10],"type":["chr"],"align":["left"]},{"label":["acl"],"name":[11],"type":["list"],"align":["right"]},{"label":["access"],"name":[12],"type":["chr"],"align":["left"]},{"label":["state"],"name":[13],"type":["chr"],"align":["left"]},{"label":["data_category"],"name":[14],"type":["chr"],"align":["left"]},{"label":["type"],"name":[15],"type":["chr"],"align":["left"]},{"label":["experimental_strategy"],"name":[16],"type":["chr"],"align":["left"]}],"data":[{"1":"Gene Expression Quantification","2":"2018-09-11T22:53:33.639292+00:00","3":"49895f4a-72ac-4d5e-ba56-8c8bb5de4758.htseq.counts.gz","4":"49895f4a-72ac-4d5e-ba56-8c8bb5de4758_count","5":"d257277b-072f-4b6c-bead-07332de2a533","6":"253547","7":"d257277b-072f-4b6c-bead-07332de2a533","8":"2016-05-29T10:41:56.268330-05:00","9":"661f22b698d45a10f6c00e420c6a2fbd","10":"TXT","11":"<chr [1]>","12":"open","13":"released","14":"Transcriptome Profiling","15":"gene_expression","16":"RNA-Seq","_rn_":"1"},{"1":"Gene Expression Quantification","2":"2018-09-11T22:53:33.639292+00:00","3":"8a799dfa-c1b5-4b13-9c91-6cbfe2abbc9f.htseq.counts.gz","4":"8a799dfa-c1b5-4b13-9c91-6cbfe2abbc9f_count","5":"167aef29-9e90-4bd1-ab3c-49bdb9866939","6":"260246","7":"167aef29-9e90-4bd1-ab3c-49bdb9866939","8":"2016-05-26T21:10:39.562381-05:00","9":"fd8ed974721299c7fce17d0722d6e8e2","10":"TXT","11":"<chr [1]>","12":"open","13":"released","14":"Transcriptome Profiling","15":"gene_expression","16":"RNA-Seq","_rn_":"2"},{"1":"Gene Expression Quantification","2":"2018-09-11T22:53:33.639292+00:00","3":"b78b6f49-3cb2-452a-a41c-6dfa90e631db.htseq.counts.gz","4":"b78b6f49-3cb2-452a-a41c-6dfa90e631db_count","5":"0c931ae0-0169-4084-be7f-e50a330baa99","6":"253906","7":"0c931ae0-0169-4084-be7f-e50a330baa99","8":"2016-05-26T21:26:50.741787-05:00","9":"c10791f045d9d3e02747a12c1716baae","10":"TXT","11":"<chr [1]>","12":"open","13":"released","14":"Transcriptome Profiling","15":"gene_expression","16":"RNA-Seq","_rn_":"3"},{"1":"Gene Expression Quantification","2":"2018-09-11T22:53:33.639292+00:00","3":"aec2e0c7-4792-41af-873c-3f3a53ec6d38.htseq.counts.gz","4":"aec2e0c7-4792-41af-873c-3f3a53ec6d38_count","5":"fdf73b53-a45b-4f06-8418-19896fc3d076","6":"255095","7":"fdf73b53-a45b-4f06-8418-19896fc3d076","8":"2016-05-29T10:30:41.561524-05:00","9":"8332437c278e6a16f8af95d13bb24ab4","10":"TXT","11":"<chr [1]>","12":"open","13":"released","14":"Transcriptome Profiling","15":"gene_expression","16":"RNA-Seq","_rn_":"4"},{"1":"Gene Expression Quantification","2":"2018-09-11T22:53:33.639292+00:00","3":"657e19a6-e481-4d06-8613-1a93677f3425.htseq.counts.gz","4":"657e19a6-e481-4d06-8613-1a93677f3425_count","5":"52d2e6bb-80f3-42a7-b2d4-cc72e5fd83f1","6":"254576","7":"52d2e6bb-80f3-42a7-b2d4-cc72e5fd83f1","8":"2016-05-29T11:08:43.811369-05:00","9":"b165dcb355976a361d16fd3a98e39783","10":"TXT","11":"<chr [1]>","12":"open","13":"released","14":"Transcriptome Profiling","15":"gene_expression","16":"RNA-Seq","_rn_":"5"},{"1":"Gene Expression Quantification","2":"2018-09-11T22:53:33.639292+00:00","3":"4172e3f8-3578-4f33-9168-6f8c2b8d0783.htseq.counts.gz","4":"4172e3f8-3578-4f33-9168-6f8c2b8d0783_count","5":"7e374f79-6f9b-4034-b4cb-d71b7404682a","6":"255804","7":"7e374f79-6f9b-4034-b4cb-d71b7404682a","8":"2016-05-30T18:32:45.805450-05:00","9":"ea2325dbb6d75ebd2fd8013d986ada4c","10":"TXT","11":"<chr [1]>","12":"open","13":"released","14":"Transcriptome Profiling","15":"gene_expression","16":"RNA-Seq","_rn_":"6"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

We could download these with standard R tools, or for larger data-sets such as this one, use the packages **transfer()** function, which uses the GDC transfer client (a separate command-line tool) to perform more robust data downloads.


# 4. Variant analysis with R

Note we could go to the NCI-GDC web portal and enter the [Advanced Search page](https://portal.gdc.cancer.gov/query) and then construct a search query to find MAF format somatic mutation files for our 'TCGA-PAAD' project.  

After some exploration of the website I came up with the following query: "`cases.project.project_id in ["TCGA-PAAD"] and files.data_type in ["Masked Somatic Mutation"] and files.data_format in ["MAF"]`".

<div class="q_box">
> **Q9**. How many MAF files for the TCGA-PAAD project were found from this advanced web search? 
</div>

Lets do the same search in R with the **files()** function from the `GenomicDataCommons` package. We will be searching based on `data_type`, `data_format`, and `analysis.workflow_type`. The last term here will focus on only one of the MAF files for this project in GDC, namely the MuTect2 workflow variant calls.


```r
maf.files = files() %>%
    filter(~ cases.project.project_id == 'TCGA-PAAD' &
        data_type == 'Masked Somatic Mutation' &
        data_format == "MAF" &
        analysis.workflow_type == "MuTect2 Variant Aggregation and Masking"
    ) %>%
    response_all()
```

```r
#View(maf.files)
attributes(maf.files)
```

```
## $names
## [1] "results"      "query"        "pages"        "aggregations"
## 
## $class
## [1] "GDCfilesResponse" "GDCResponse"      "list"
```

```r
head(maf.files$results)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["data_type"],"name":[1],"type":["chr"],"align":["left"]},{"label":["updated_datetime"],"name":[2],"type":["chr"],"align":["left"]},{"label":["file_name"],"name":[3],"type":["chr"],"align":["left"]},{"label":["submitter_id"],"name":[4],"type":["chr"],"align":["left"]},{"label":["file_id"],"name":[5],"type":["chr"],"align":["left"]},{"label":["file_size"],"name":[6],"type":["int"],"align":["right"]},{"label":["id"],"name":[7],"type":["chr"],"align":["left"]},{"label":["created_datetime"],"name":[8],"type":["chr"],"align":["left"]},{"label":["md5sum"],"name":[9],"type":["chr"],"align":["left"]},{"label":["data_format"],"name":[10],"type":["chr"],"align":["left"]},{"label":["acl"],"name":[11],"type":["list"],"align":["right"]},{"label":["access"],"name":[12],"type":["chr"],"align":["left"]},{"label":["state"],"name":[13],"type":["chr"],"align":["left"]},{"label":["data_category"],"name":[14],"type":["chr"],"align":["left"]},{"label":["type"],"name":[15],"type":["chr"],"align":["left"]},{"label":["experimental_strategy"],"name":[16],"type":["chr"],"align":["left"]}],"data":[{"1":"Masked Somatic Mutation","2":"2018-09-11T22:53:33.639292+00:00","3":"TCGA.PAAD.mutect.fea333b5-78e0-43c8-bf76-4c78dd3fac92.DR-10.0.somatic.maf.gz","4":"TCGA-PAAD-mutect-public_10.0_new_maf","5":"fea333b5-78e0-43c8-bf76-4c78dd3fac92","6":"6991687","7":"fea333b5-78e0-43c8-bf76-4c78dd3fac92","8":"2017-12-01T17:52:47.832941-06:00","9":"cdddbf7bc36774e85a5033ad1be223ba","10":"MAF","11":"<chr [1]>","12":"open","13":"released","14":"Simple Nucleotide Variation","15":"masked_somatic_mutation","16":"WXS","_rn_":"1"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

<div class="q_box">
> **Q10**. What line in the above code would you modify to return all MAF files for the TCGA-PAAD project? 
</div>

We will use the **ids()** function to pull out the unique identifier for our MAF file


```r
uid <- ids(maf.files)
uid
```

```
## [1] "fea333b5-78e0-43c8-bf76-4c78dd3fac92"
```

Once we have the unique identifier(s) (in this case, only fea333b5-78e0-43c8-bf76-4c78dd3fac92), the **gdcdata()** function downloads the associated files to a cache directory on your computer and returns a file-name for each identifier.


```r
#maffile = gdcdata(uid, destination_dir =".")
maffile = gdcdata(uid)

maffile
```

```
##                                                                                                                               fea333b5-78e0-43c8-bf76-4c78dd3fac92 
## "/Users/barry/Library/Caches/GenomicDataCommons/fea333b5-78e0-43c8-bf76-4c78dd3fac92/TCGA.PAAD.mutect.fea333b5-78e0-43c8-bf76-4c78dd3fac92.DR-10.0.somatic.maf.gz"
```


### MAF analysis
The MAF file is now stored locally and the maftools package workflow, which starts with a MAF file, can proceed, starting with reading the pancreatic cancer MAF file.


```r
vars = read.maf(maf = maffile, verbose = FALSE)
```

With the data now available as a **maftools** MAF object, a lot of functionality is available with little code. While the maftools package offers quite a few functions, here are a few highlights. Cancer genomics and bioinformatics researchers will recognize these plots:


##  Plotting MAF summary.

We can use `plotmafSummary()` function to plot a summary of the maf object, which displays number of variants in each sample as a stacked barplot and variant types as a boxplot summarized by Variant_Classification. We can add either mean or median line to the stacked barplot to display average/median number of variants across the cohort.


```r
plotmafSummary(maf =vars, rmOutlier = TRUE,
               addStat = 'median', dashboard = TRUE,
               titvRaw = FALSE)
```

![](lecture18_part1_BIMM143_W18_files/figure-html/unnamed-chunk-17-1.png)<!-- -->


## Drawing oncoplots

A very useful summary representation of this data can be obtained via so-called *oncoplots*, also known as *waterfall plots*. 


```r
oncoplot(maf = vars, top = 10)
```

![](lecture18_part1_BIMM143_W18_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

You might to run the **oncoplot()** command in the R Console and then zoom the display to see the full plot (as it is rather large and may not appear initially in your Rmarkdown document before Knitting. Another option is to send your plot to a PNG or PDF plot device directly, for example:


```r
# Oncoplot for our top 10 most frequently mutated genes
pdf("oncoplot_panc.pdf")
oncoplot(maf = vars, top = 10, fontSize = 12)
dev.off()
```

```
## quartz_off_screen 
##                 2
```

> **NOTE**: The **oncoplot()** function is a wrapper around ComplexHeatmap's `OncoPrint()` function and there are lots and lots of possible customization options as usual with R graphics.

> **NOTE**: Variants annotated as Multi_Hit are those genes which are mutated more than once in the same sample.

## Oncostrip

We can visualize any set of genes using the **oncostrip()** function, which draws mutations in each sample similar to the graphic on the NCI-GDC web portal. Note that **oncostrip()** can be used to draw any number of genes using the input `top` or `genes` arguments


```r
oncostrip(maf=vars, genes=c("KRAS", "TP53"))
```

![](lecture18_part1_BIMM143_W18_files/figure-html/unnamed-chunk-20-1.png)<!-- -->


Another plot focusing on KRAS in our particular dataset.


```r
lollipopPlot(maf = vars, gene = 'KRAS', 
                         showMutationRate = TRUE, domainLabelSize = 3)
```

```
##    HGNC refseq.ID protein.ID aa.length
## 1: KRAS NM_004985  NP_004976       188
## 2: KRAS NM_033360  NP_203524       189
```

![](lecture18_part1_BIMM143_W18_files/figure-html/unnamed-chunk-21-1.png)<!-- -->

Lets do one for p53 also


```r
lollipopPlot(maf = vars, gene = 'TP53')
```

```
##    HGNC    refseq.ID   protein.ID aa.length
## 1: TP53    NM_000546    NP_000537       393
## 2: TP53 NM_001126112 NP_001119584       393
## 3: TP53 NM_001126118 NP_001119590       354
## 4: TP53 NM_001126115 NP_001119587       261
## 5: TP53 NM_001126113 NP_001119585       346
## 6: TP53 NM_001126117 NP_001119589       214
## 7: TP53 NM_001126114 NP_001119586       341
## 8: TP53 NM_001126116 NP_001119588       209
```

![](lecture18_part1_BIMM143_W18_files/figure-html/unnamed-chunk-22-1.png)<!-- -->


# Summary

Additional functionality is available for both the `GenomicDataCommons` and `maftools` packages not to mention the 100's of other bioinformatics R packages that can now work with this type of data in both exploratory and targeted analysis modes. 

The purpose of this hands-on session was to highlight how one can leverage two such packages to quickly gain insight into rapidly expanding public cancer genomics data-sets. Hopefully this will inspire your further exploration of these and other bioinformatics R packages.


