# ReferenceRangeR

ReferenceRangeR is a browser-based application for the indirect estimation of
reference intervals from routine laboratory data. It implements five published
estimation procedures behind a common interface and provides the pre-analytical
checks required to decide whether an interval should be stratified by sex or age.

The application is intended for research purposes only. It does not constitute an
in-vitro diagnostic medical device, and the intervals it produces must not be
adopted for patient care without independent verification.

---

## 1. Principle

Direct reference interval estimation requires a prospectively recruited cohort of
healthy individuals. Indirect estimation instead uses the results already held in
a laboratory information system. Such a collection is a mixture of results from a
non-pathological and a pathological subpopulation. All methods implemented here
model that mixture and recover the parameters of the non-pathological component,
from which the 2.5th and 97.5th percentiles are derived.

The quality of an indirect estimate depends on the composition of the input data.
The proportion of pathological results, the sample size, the skewness of the
distribution and the presence of values below the limit of quantification all
influence which method is appropriate. No single procedure is superior across all
data sets; Section 5.2 summarises the properties relevant to that choice.

---

## 2. Workflow

The application is organised in three tabs, to be worked through in order:

1. **Input Data** — load the data, assign the columns, inspect the import.
2. **Check Data** — test for sex differences and age dependence.
3. **Calculate Reference Interval** — restrict the population, choose a method,
   estimate the interval and export the report.

Step 2 is optional but recommended: an interval estimated over a population that
should have been stratified is biased regardless of the method used.

A **Proceed to data check** and a **Proceed to calculation** button at the foot of
the Input Data tab advance the workflow without using the navigation bar.

---

## 3. Data input

### 3.1 Accepted sources

Three input paths are available on the Input Data tab.

- **Upload** — files in CSV, TSV, TXT, XLSX or XLS format. The field separator
  and the decimal sign are detected automatically. The upload limit is 200 MB.
- **Paste Data from Clipboard** — a header line followed by at least one data row,
  as copied from a spreadsheet or a query result. Tab, semicolon, comma and
  whitespace are recognised as separators; the separator occurring most often in
  the header line is used.
- **Demo** — 50,000 simulated results with an age-dependent median and a
  pathological fraction of approximately 10 %, provided so that the application
  can be evaluated without transferring real data.

### 3.2 Variables

| Variable  | Status   | Purpose                                              |
|-----------|----------|------------------------------------------------------|
| Result    | required | the measured quantity                                 |
| Age       | optional | age stratification and age drift analysis             |
| Sex       | optional | sex stratification and sex difference testing         |
| Trimester | optional | restriction to a trimester of pregnancy (advanced)    |

Columns are assigned in the **Columns** panel. On import the application proposes
an assignment based on the column names, including the German equivalents
(*Alter*, *Geschlecht*, *Trimenon*). Any proposal can be overridden.

### 3.3 Treatment of individual values

Result values are interpreted as follows.

- A comma is accepted as the decimal separator.
- A value prefixed with `<` is treated as being below the limit of quantification.
  Such values are retained as text and are used only by the TMC method, which
  models left-censored data explicitly. They do not enter the other procedures.
- Empty fields, non-numeric text and values of zero or below are discarded. The
  number discarded per category are reported separately.

Age values are read numerically; entries that cannot be parsed are counted as
unrecognised. Trimester values outside 1 to 3 are likewise counted as
unrecognised and are otherwise ignored.

### 3.4 Harmonisation of sex labels

The labels `m`, `male`, `M`, `männlich` and `Mann` are mapped to male; `f`, `w`,
`female`, `F`, `W`, `weiblich` and `Frau` to female; `d`, `D`, `diverse` to
diverse. Any other non-empty label is unrecognised, and the corresponding results
are excluded from all sex-stratified analyses.

When unrecognised labels are present, an assignment panel appears below the column
selection. Each unrecognised label can be assigned to the female or the male
group; **Reset labels** discards all manual assignments and returns every label to
the unrecognised pool. Reassignment takes effect immediately and is reflected in
the counters and the preview.

### 3.5 Verification of the import

Five counters summarise the import: total records, usable results, and the number
of records carrying an age, a sex and — in advanced mode, where a trimester column
has been assigned — a trimester. Each counter reports the number of unrecognised
entries beneath it whenever that number exceeds zero.

The table below the counters shows the first 100 records. **View original data**
switches it to the columns exactly as they were read, including the records that
were later discarded and the sex labels in their original spelling; **View cleaned
data** returns to the processed values. Comparing the two views is the recommended
way to confirm that the import was interpreted as intended.

---

## 4. Pre-analytical checks

Both checks are performed on the Check Data tab and operate on the population
selected in the sidebar of that tab.

### 4.1 Sex differences

Groups containing fewer than 100 results are excluded, and the number of results
removed on that basis is reported. Calculation is only performed when a minimum of
500 samples are present. Where two groups remain, a Wilcoxon rank-sum
test is applied; where more than two remain, a Kruskal-Wallis test is used.

Statistical significance alone is a weak criterion in this setting, because the
sample sizes typical of indirect estimation will render clinically irrelevant
differences significant. The largest relative difference between group medians is
therefore reported alongside the p-value and compared with the permissible
uncertainty (pU) estimated from the data.

### 4.2 Age drift

The age stratification algorithm aims to detect age-related
trends. It works by iteratively grouping age intervals based
on similarity of their medians. The procedure is initiated
with the allocation of 1,000 equidistant age bins. Subsequently,
the median result for each age bin is calculated, and
the differences between these medians and the median of
the adjacent bin are determined. Bins with the smallest
difference in the median result are then merged iteratively
until either the desired number of groups is attained or
the differences surpass a predetermined threshold. The
threshold is derived from the pU, which represents the
maximum allowable total measurement error and is estimated
from the dataset, as previously outlined. Groups with
fewer than 50 observations are excluded from the analysis.
The resulting strata are subsequently summarized visually
and in tabular form. A median trend, inclusive of the 95 %
confidence intervals, is incorporated through the implementation
of additive quantile regression.

Each row of the age group table carries a **use** button that transfers the from
and to values of that group to the age range fields of both the Check Data and the
Calculate Reference Interval tab, so that a group can be carried directly into an
estimation.

---

## 5. Reference interval estimation

### 5.1 Analysis population

The sidebar of the Calculate Reference Interval tab restricts the data before
estimation. Available restrictions are sex (all, male, female), an age range and,
in advanced mode where a trimester column has been assigned, a trimester. The
restrictions are applied jointly. The selection is recorded in the report.

At least 500 usable results are required. Estimates from smaller collections are
refused rather than reported with a caveat.

### 5.2 Methods

| Method  | Principle                                              | Notable requirement                       |
|---------|--------------------------------------------------------|-------------------------------------------|
| refineR | Box-Cox transformed model of the non-pathological fraction | —                                    |
| TMC     | Chi-square minimization, left-censored values modelled explicitly | uses `<` values        |
| TML     | Truncated maximum likelihood         | —                                         |
| kosmic  | estimation from the central, non-pathological part of the distribution | requires non-integer results |
| reflimR | quantile-based estimation with detection of lognormality | —                                        |

Guidance for method selection, in the absence of a general recommendation:

- pathological fraction above 20 %: refineR, kosmic or TMC;
- high proportion of values below the limit of quantification: TMC;
- strongly skewed distributions: TMC.

kosmic requires that at least one selected result carries a decimal place. Where
the selection consists exclusively of whole numbers the method declines to run and
reports the reason. TMC operates on the character representation of the results
and therefore counts left-censored values in its sample size; the other methods
count only results that could be converted to a number.

### 5.3 Method-specific options

The following options are available in advanced mode.

- **refineR** — a modified Box-Cox transformation, and 0 to 50 bootstrap
  iterations. Bootstrapping yields confidence intervals for both limits at a
  considerable increase in computation time; 0 iterations suppresses them.
- **TML** — a fast mode that rounds results to three significant figures. The
  full-precision mode is substantially slower.

### 5.4 Comparison with established limits

Where a reference interval is available for comparison, for instance one supplied by the assay
manufacturer, its limits can be entered after activating **Compare with known
limits**. The estimated limits are then assessed against equivalence bands derived
from the permissible uncertainty of the entered limits. The bands are shaded in
the plot, and each estimated limit is reported as lying within or outside the
permissible uncertainty of its counterpart. The comparison is available for the
refineR, TMC, TML and reflimR methods.

---

## 6. Reporting

**Export as PDF** produces a single-page report containing the estimated reference interval,
the confidence intervals where these were computed, the transformation parameter,
the data source, the selection applied, the method with its settings and its
literature reference, the comparison with established limits where one was
requested, and the diagnostic plot. The file is named after the method and the
date of export.

---

## 7. Limitations

Indirect estimation assumes that the non-pathological subpopulation dominates the
collection and that its distribution can be described by the model underlying the
chosen method. Where a collection is drawn from a strongly pre-selected population
— a specialist clinic, or a test ordered only on clinical suspicion — that
assumption may not hold, and the estimate will not describe a healthy reference
population.

Estimates should be compared across at least two methods before adoption.
Divergence between methods indicates that the data violate the assumptions of at
least one of them.

---

## 8. Citation

When reporting results obtained with this application, please cite:

Brandhorst G, Voß M, Wosniok W, Arzideh F, Haeckel R, Rosenkranz D, Streichert T,
Petersmann A. ReferenceRangeR: a novel tool designed to facilitate reference
interval estimation and verification. *Clin Chem Lab Med* 2025;64(3):680–687.
doi:10.1515/cclm-2025-0637

---

## 9. References

Each method is documented in the publication linked from its entry in the method
list of the application. The same references are reproduced here.

- refineR: *Scientific Reports*.
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8346497/
- TMC: *Clinical Chemistry and Laboratory Medicine*.
  https://doi.org/10.1515/cclm-2018-1341
- TML: *Clinical Chemistry and Laboratory Medicine*.
  https://doi.org/10.1515/CCLM.2007.249
- kosmic: *Scientific Reports*.
  https://www.nature.com/articles/s41598-020-58749-2
- reflimR: *Journal of Laboratory Medicine*.
  https://doi.org/10.1515/labmed-2023-0042

Source code: https://github.com/gubruol/ReferenceRangeR
Disclaimer: https://kc.uol.de/disclaimer/
