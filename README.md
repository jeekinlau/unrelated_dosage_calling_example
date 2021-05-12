# unrelated_dosage_calling_example
Dosage calling for unrelated genotypes like a cultivar evaluation GWAS panel

The R_script uses the output of Axiom Analysis Suite and assigns dosage calls to each separate probe as separate markers using R package ['fitpoly 3.0.0'](https://cran.r-project.org/web/packages/fitPoly/index.html) Voorrips RE, Gort G, Vosman B, 2011 

Following dosage calling, a custom R script combines the forward and backward strand probes of each SNP and produces two files. The first file contains the dosage calls and the second contains more information as to whether that SNP was called using two agreeing probes, just one probe, or whether the marker was discarded due to conflicting calls.

## UPDATE: tool available
The custom R script has been packaged into a R package hosted on github

To install the R package use the following code:
```
install.packages('devtools')
devtools::install_github("jeekinlau/RoseArrayTools")
```


## Links

* [R_script](https://raw.githubusercontent.com/jeekinlau/unrelated_dosage_calling_example/main/docs/fitPoly_and_custom_script.R)
* [Vignette](https://jeekinlau.github.io/unrelated_dosage_calling_example/fitpoly_tutorial_example.html)