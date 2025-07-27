# ReferenceRangeR

A digital tool designed for the estimation and verification of reference intervals using real-world laboratory data.
Built in R with the Shiny framework, the tool enables users to analyze up to 200,000 test results via simple copy-and-paste input - eliminating the need for file uploads or coding. The tool integrates multiple statistical methods, including refineR, TMC, TML, kosmic, and reflimR, and supports sex- and age-based stratification tools and dynamic drift detection algorithm. 
It provides automated data cleaning, visualization, and analysis. Sex-based differences are visualized using violin and boxplots, while statistical significance is assessed by non-parametric tests. Age stratification is guided by an iterative algorithm that merges age groups based on minimal median deviation within a user-defined number of bins, compared against dataset-derived uncertainty thresholds. The application is accessible via a public web interface and offers a Docker-based installation for secure, local deployment.

### Screenshot
![grafik](https://github.com/user-attachments/assets/a199916e-dea0-4bd4-8651-87ec6b57a31f)

## Homepage (Demo)
https://kc.uol.de/referenceranger/

## Docker Installation
You can run ReferenceRangeR locally using Docker. This setup is recommended for secure, isolated deployment. Note: Only supported on x64 (Intel/AMD) platforms. Not compatible with ARM-based systems like Raspberry Pi.
### Prerequisites
[Install Docker](https://docs.docker.com/get-started/) for your operating system.

### Start the App
Open a terminal or command prompt and run:
```
docker run -d -p 80:3838 --name referenceranger gubruol/referenceranger:latest
```
This will download the latest Docker image und map the container to local port 80.

### Access the App
Open you Browser an access the tool at http://localhost
