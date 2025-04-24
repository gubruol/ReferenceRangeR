# ReferenceRangeR

A digital tool designed for the estimation and verification of reference intervals using real-world laboratory data.
Built in R with the Shiny framework, the tool enables users to analyze up to 200,000 test results via simple copy-and-paste input - eliminating the need for file uploads or coding. The tool integrates multiple statistical methods, including refineR, TMC, TML, kosmic, and reflimR, and supports sex- and age-based stratification tools and dynamic drift detection algorithm. 
It provides automated data cleaning, visualization, and analysis. Sex-based differences are visualized using violin and boxplots, while statistical significance is assessed via ANOVA. Age stratification is guided by an iterative algorithm that merges age groups based on minimal median deviation within a user-defined number of bins, compared against dataset-derived uncertainty thresholds. The application is accessible via a public web interface and offers a Docker-based installation for secure, local deployment.

# Homepage
https://kc.uol.de/referenceranger/

## Installation (Docker)
```
docker run -d -p 80:3838 --name referenceranger gubruol/referenceranger:latest
```
