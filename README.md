# Discovering Structural Biomarkers for Retinal Diseases

## Repository structure
```text
discovering-structural-biomarkers-for-retinal-diseases/
│
├── data/
│   └── PatientsData.xlsx              # raw dataset; not uploaded
│
├── R/
│   ├── packages.R                     
│   ├── fda_utils.R                    
│   ├── eval_utils.R                  
│   └── model_defs.R                  
│
├── scripts/
│   ├── 01_preprocessing.R             
│   ├── 02_eda_multivariate.R          
│   ├── 03_covariance_analysis.R      
│   ├── 04_baseline_models.R           
│   ├── 05_eda_functional.R            
│   ├── 06_fda_plots.R                 
│   └── 07_logreg_scores.R             
│
├── datasets/                          # processed .rds datasets; generated locally, not committed
│
├── results/                           # model outputs and .rds results; generated locally, not committed
│
├── figures/                           
│
├── .gitignore                         
├── .Rproj
└── README.md
```


## Dataset presentation

## Analysis outline 

## Multivariate exploratory data analysis

## Classification models

## Functional data analysis

## Functional Data Analysis

### Smoothing

Each patient's six retinal layer profiles are reconstructed as smooth functions
from the discrete OCT thickness measurements. We use a cubic B-spline basis
(`norder = 4`) with $K = 60$ basis functions on the normalized domain
$[0, 1]$, and a second-derivative roughness penalty:

$$
\hat{c} = \arg\min_{c} \sum_{j} \bigl(y_j - x(t_j)\bigr)^2 + \lambda \int \bigl(x''(t)\bigr)^2 \, dt
$$

The smoothing parameter $\lambda$ is selected by Generalized Cross-Validation
(GCV) over a logarithmic grid $\lambda \in 10^{\text{seq}(-4, 4, 0.5)}$. Because
the downstream analyses (MFPCA, functional distances) require all patients to
live in a common functional space, we do **not** select a per-patient
$\lambda$. Instead, for each layer we compute the per-patient GCV across the
grid and choose the single $\lambda$ that minimizes the **average GCV across
all patients** in that layer. This yields one common smoothing parameter per
layer, applied uniformly to every patient, ensuring cross-patient
comparability.

The selected smoothing parameters per layer are stored in
`df_smooth_by_lyr[[ell]]$lambda_opt` and can be extracted with:

```r
sapply(df_smooth_by_lyr, function(layer) layer$lambda_opt)
```
The following are the reconstructions of the OCT scans previously reported.

<p align="center">
  <img src="figures/functional_data_analysis/smooth_healthy5.png" alt="healthy5" width="750">
</p> 

<p align="center">
  <img src="figures/functional_data_analysis/smooth_AMD12.png" alt="AMD5" width="750">
</p> 


### Mean functions per class

Per-class mean thickness functions. For class $i$ and layer $\ell$, the mean
function is defined as

$$
\mu_i^{(\ell)}(t) = \frac{1}{n_i} \sum_{p \in \text{class } i} X_p^{(\ell)}(t),
$$

the average thickness profile over all patients of class $i$ in layer $\ell$.

<p align="center">
  <img src="figures/functional_data_analysis/mean_funs.png" alt="mean functions per class" width="750">
</p> 

### Variance functions per class

Per-class variance functions. For class $i$ and layer $\ell$,

$$
\sigma_i^{(\ell)\,2}(t) = \frac{1}{n_i - 1} \sum_{p \in \text{class } i} \bigl(X_p^{(\ell)}(t) - \mu_i^{(\ell)}(t)\bigr)^2,
$$

the within-class variability of thickness at each retinal position. Inflated
variance functions indicate diseases with heterogeneous structural
presentation.

<p align="center">
  <img src="figures/functional_data_analysis/var_funs.png" alt="variance functions per class" width="750">
</p> 



### Multivariate Functional PCA (MFPCA)

Each patient is described not by one function but by a **vector-valued
function** stacking the six layer profiles:

$$
\mathbf{X}_i(t) = \bigl(X_i^{(1)}(t), \dots, X_i^{(6)}(t)\bigr), \qquad t \in \mathcal{T}=[0,1].
$$

The natural space for these objects is the product Hilbert space
$\mathcal{H} = [L^2(\mathcal{T})]^6$, equipped with the inner product

$$
\langle \mathbf{f}, \mathbf{g} \rangle_{\mathcal{H}} = \sum_{\ell=1}^{6} \int_{\mathcal{T}} f^{(\ell)}(t) \, g^{(\ell)}(t) \, dt.
$$

MFPCA finds the eigenfunctions of the covariance operator $\mathcal{C}$ on
$\mathcal{H}$. These eigenfunctions are themselves vector-valued,
$\boldsymbol{\psi}_k(t) = (\psi_k^{(1)}(t), \dots, \psi_k^{(6)}(t))$,
satisfying

$$
\mathcal{C} \boldsymbol{\psi}_k = \lambda_k \boldsymbol{\psi}_k.
$$

Each patient is then represented by a set of scalar scores

$$
\xi_{ik} = \langle \mathbf{X}_i - \boldsymbol{\mu}, \boldsymbol{\psi}_k \rangle_{\mathcal{H}} = \sum_{\ell=1}^{6} \int_{\mathcal{T}} \bigl(X_i^{(\ell)}(t) - \mu^{(\ell)}(t)\bigr) \psi_k^{(\ell)}(t) \, dt,
$$

with the reconstruction

$$
\mathbf{X}_i(t) \approx \boldsymbol{\mu}(t) + \sum_{k=1}^{K} \xi_{ik} \, \boldsymbol{\psi}_k(t).
$$

The **same score** $\xi_{ik}$ multiplies all six components of
$\boldsymbol{\psi}_k$: a single scalar controls how all six layers jointly
deviate from their means, so the layers are not treated independently — each
eigenfunction encodes a joint spatial pattern across all layers.

We use the `MFPCA` R package (Happ & Greven, 2018).

The 8 first multivariate eigenfunctions $\psi_k^{(\ell)}(t)$ account for 88% of the total variance, 
Each eigenfunction is the joint spatial pattern across the six layers that
captures the $k$-th largest mode of patient-to-patient variation.

<p align="center">
  <img src="figures/functional_data_analysis/mfpca_fun_plot.png" alt="mfpca funs" width="750">
</p> 

The following are perturbation plots: the mean function $\mu^{(\ell)}(t)$ is perturbed by
$+ \bar{\xi}_{ij} \cdot \psi_k^{(\ell)}(t)$, visualizing how a positive or negative score
along component $k$ reshapes the six-layer profile relative to the mean.

<p align="center">
  <img src="figures/functional_data_analysis/mfpca_lyr_perturbations.png" alt="mfpca perturbations" width="750">
</p> 

We can clearly see that only for the layer5 component of the first principal component there is a significant difference among classes. 

### Functional separation analysis

To quantify how much of the total functional variability is attributable to
disease class, we decompose the dispersion in the product space
$\mathcal{H} = [L^2(\mathcal{T})]^6$ into within-class and between-class
components:

$$
W_k = \sum_{i \in k} \sum_{\ell=1}^{6} \int_{\mathcal{T}} \bigl(X_i^{(\ell)}(t) - \hat{\mu}_k^{(\ell)}(t)\bigr)^2 \, dt,
$$

$$
B = \sum_{k=1}^{5} n_k \sum_{\ell=1}^{6} \int_{\mathcal{T}} \bigl(\hat{\mu}_k^{(\ell)}(t) - \hat{\mu}^{(\ell)}(t)\bigr)^2 \, dt.
$$

The functional analog of the ANOVA effect size is the ratio

$$
\eta^2 = \frac{B}{B + \sum_k W_k},
$$

the fraction of total functional variability explained by class membership.
This is a trace-based descriptive measure (it aggregates across layers and
positions without modeling inter-layer covariance) and requires no
distributional assumptions. We compute it via distance-based PERMANOVA
(`vegan::adonis2`) on the pairwise functional $L^2$ distance matrix; the
reported $R^2$ by the function is exactly $\eta^2$.  
The value observed is $\eta^2 = 0.17$,
which suggests a strong heterogeneity even inter-class, most of the variability can't be explained 
just by the membership to a group. 

To localize where along the retina the classes separate, we compute a
pointwise multivariate statistic. At each position $t$, the six layer values
form a multivariate observation in $\mathbb{R}^6$, and a one-way MANOVA across
classes yields Wilks' Lambda $\Lambda(t) \in (0, 1]$, accounting for
inter-layer covariance at that position:

$$
\Lambda(t) = \frac{\det \mathbf{W}(t)}{\det\bigl(\mathbf{W}(t) + \mathbf{B}(t)\bigr)},
$$

where $\mathbf{B}(t)$ and $\mathbf{W}(t)$ are the $6 \times 6$ between- and
within-class scatter matrices at position $t$.

We've  plot pointwise $1 - \Lambda(t)$ along the retina (rather
than $\Lambda(t)$ so that large values correspond to strong separation, 
since small Wilks' Lambda indicates strong class separation). The pronounced
peak in the central foveal region identifies it as the locus of strongest
multivariate class separation, consistent with the foveal concentration of
discriminative signal.

![Pointwise multivariate separation](figures/pointwise_separation.png)


## Results

### References
- Ramsay, J. O. & Silverman, B. W. (2005). *Functional Data Analysis* (2nd ed.).
  Springer.
- Happ, C. & Greven, S. (2018). *Multivariate Functional Principal Component
  Analysis for Data Observed on Different (Dimensional) Domains.* Journal of
  the American Statistical Association, 113(522), 649–659.

### Authors
Adelaide Carnevale, Riccardo Cicuttin, Giulio Dalla Costa, Francesca Elefante      
Supervisor: Dr. Lara Cavinato
