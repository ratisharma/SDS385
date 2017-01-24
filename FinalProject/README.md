# Replication of Methods in Houser et al. (2015) using R
##Spencer Woody, Final Project for SDS 385 Fall 2016

**December 5, 2016**

This is my final project for SDS 385 fall semester 2016. Here I present a walkthrough of the 2015 paper ["Controlled Measurement and Comparative Analysis of Cellular Components in *E. coli* Reveals Broad Regulatory Changes in Response to Glucose Starvation"][1] by Houser et al. written in R and presented in [R Markdown](http://rmarkdown.rstudio.com)

### How to view this project

UPDATE 01-23-2016

[Got it hosted on github now! Just click here.](https://spencerwoody.github.io/SDS385/FinalProject/) 

~~*For best results*, download the `project.html` file (click [here][2] and hit `Cmd + S`) and open it in your browser. On my laptop, the TeX equations and plots only render in Safari. I have not figured out how to get github to host my html page yet, but hopefully I will have that up soon. Oh, and you'll need to have an Internet connection to render the equations. Alternatively, you can also look at `project.pdf` , (click [here]) but I have formatted my report with the intent of it being viewed as a webpage. well it loks like RStudio is not compiling my Rmd into TeX, so looks like HTML will have to do for now. Sorry about that!~~

### Code

* `project.R` is the R-script used for the analysis.
* `project.Rmd` is the R Markdown-file used to produce the report.

### Folders
 
 * `Data` contains the raw and processed data. `rna.csv` is the raw RNA counts, `protein.csv` is the raw protein counts. `dictionary.csv` is used to map proteins to their transcripts, the "b"-name synonym for each gene, and the operon of the each gene. `entrez_dictionary.csv` contains is used to map a gene to an Entrez ID for the purposes of running DAVID API queries. `sig_rnas.csv` is a list of RNAs determined to be significant by DESeq. All other files contain the processed RNA and protein data. See `project.R` or the report for details. 
 * `DESeq` contains the R-script for running DESeq, and the results of DESeq on the RNA data.
 * `Results` contains the results from fitting the time course function for RNAs (`rna_parameters.csv`) and proteins (`pro_parameters.csv`), and the results of functional annotation clusering for up-regulated RNAs, down-regulated RNAs, up-regulated (`FuncAnnotChart_rna_up`, `FuncAnnotChart_rna_do`, `FuncAnnotChart_pro_up`, and `FuncAnnotChart_pro_do`, respectively).
 * `Prospectus` simply contains my proposal for the project which I turned in to Professor James Scott. 
 
### Shout-outs 
 
A big thank you to Dr. John R. Houser of the Marcotte Lab at UT-Austin, who authored the paper and ran the original analysis, and was kind enough to answer my questions about his techniques and let me use his data, and of course to [Professor James Scott][4] who helped me with the optimization of the time course model and taught a great class this semester. 

[1]: http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004400

[2]: https://raw.githubusercontent.com/spencerwoody/SDS385/master/FinalProject/project.html

[3]: https://github.com/spencerwoody/SDS385/raw/master/FinalProject/example.pdf

[4]: https://github.com/jgscott

