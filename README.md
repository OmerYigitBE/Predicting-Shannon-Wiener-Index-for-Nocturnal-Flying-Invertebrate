# Predicting Shannon Wiener Index (SWI) for Nocturnal Flying Invertebrate
In this project, I analyzed a biodiversity data to build a parametric regression model for SWI - a measure of biodiversity. Then, I investigated a possible association between duration and temperature, using non-parametric techniques.

Codes are provided, and the analysis results with graphs and tables are included in the report. Below is a summary of the analysis.

## Data
The data holds 400 records of flying nocturnal invertebrate (e.g., insects such as moths, beetles, mosquitos, etc.) biodiversity, based on 400 sampling events in 2017 in the Belgian province of Limburg. Each sampling event was carried out on a patch, being a small nature area that has been subject to nature management for at least 3 months. Each patch was visited only once (so the 400 data points reflect 400 different locations). We assume spatially uncorrelated data, as is indicated by an exploratory test for spatial autocorrelation (not shown). Sampling was carried out by catching the invertebrates in a net in the evening. Sampling events were carried out from the beginning of Spring until the middle of Summer, when sampling stopped because of extreme heat. Note that multiple sampling events could be carried out on the same day by different persons. Also note that the dataset is ordered chronologically. However, information about the exact dates is not provided.

The outcome of interest is the Shannon-Wiener index (SWI), which is a measure of biodiversity. SWI is a non-negative metric that is usually in practice smaller than 4.5. Low values denote low diversity, while higher values denote higher diversity.

- SWI: the Shannon-Wiener index for (flying nocturnal) invertebrate diversity
on the patch (non-negative, larger values denote higher diversity)
- SWF: an (adjusted) Shannon-Wiener index for floristic diversity on the patch.
The interpretation of this metric is the same as for SWI (non-negative,
larger values denote higher diversity)
- temperature: temperature at the sampling event (in degrees Celsius)
- size: the size of the sampling patch (in m2)
- management: the number of years that the patch has been subject to nature management
- duration: the duration of a sampling event (in minutes)

## Parametric Analysis
- In the very first step, the data is divided into training and validation datasets.
- Then, an exploratory analysis is done on the training dataset - for checking missingness, outliers, correlations, distributions etc.
- After having a better understanding of the dataset, a simple linear regression model is fitted.
- Then, Gauss-Markov conditions, multicollinearity, studentized residuals and heteroscedasticity are checked.
- In order to solve the heteroscedasticity problem, stepwise regression models and variable transformations are tried. Also, weighted least square regression approach is tried for comparison.
- Two candidate models are compared based on their performances on the validation dataset and the ultimate model is found.
- Then, this model is fitted to the entirety of the dataset to gain insight about the predictors.
  - As temperature increases, biodiversity increases, because the environmental conditions will
be suitable to more species.
  - Similarly, being subjected to nature management longer increases biodiversity.
  - Temperature is found to be more influential than management.

## Non-Parametric Analysis
- The relationship between the variables of interest (duration, temperature) is not strictly linear.
- Quadratic regression models with different orders and non-parametric lowess regression models with different degrees and spans are fitted and compared.
- Ultimately, non-parametric model is found to perform better in different criteria.
