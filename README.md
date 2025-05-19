# The Break in the Cycle: Inositol Pyrophosphate Fluxomics Disentangled via Mathematical Modelling

**Repository for analysis scripts and model files associated with the publication:**

> **"The Break in the Cycle: Inositol Pyrophosphate Fluxomics Disentangled via Mathematical Modelling"**
> Jacques Hermes et al., 2025
> [bioRxiv preprint](https://doi.org/10.1101/2025.03.27.645712)

---

## 🧪 Abstract

This study investigates the metabolic pathways of inositol pyrophosphates (IPPs) in the yeast cell line ΔSPX and the human tumor cell line HCT116. Utilizing pulse-labelling experiments with ¹⁸O water and ordinary differential equation (ODE) models, we explore the synthesis and turnover of the highly phosphorylated IPP, 1,5-InsP₈. Our findings challenge the notion that 1,5-InsP₈ can be synthesized through distinct routes, revealing a linear reaction sequence in both systems.

Using model reduction via the profile likelihood method, we achieved statistically concise identifiability analysis that led to significant biological insights. In yeast, 1,5-InsP₈ production primarily occurs through phosphorylation of 5-InsP⁷, while in HCT116 cells, it is driven mainly by 1-InsP⁷, with variations observed under different phosphate conditions.

---

## 📁 Repository Structure

```
├── Yeast_Low_to_High_Pi/         # Yeast ΔSPX strain, Pi increase
│   ├── Nonreduced_model/
│   ├── First_reduction_step/
│   └── ...
├── HCT116_Normal_to_Normal_Pi/   # HCT116 control condition
│   ├── ...
├── HCT116_Low_to_High_Pi/        # HCT116 Pi increase
│   ├── ...
└── README.md
```

Each top-level folder represents one of the three experimental conditions described in the manuscript. Subfolders contain scripts for modeling, reduction, profiling, and visualization.

---

## 📦 Required R Packages

The following R packages are used throughout the analysis:

```r
install.packages(c(
  "deSolve", "trust", "parallel", "ggplot2", "ggthemes",
  "cowplot", "magrittr", "tidyr", "fBasics",
  "data.table", "purrr"
))

# These packages may be installed from CRAN or GitHub:
install.packages("cOde")
install.packages("dMod")
```

---

## ▶️ How to Use

1. Clone this repository:

   ```bash
   git clone https://github.com/yourusername/The_break_in_the_cycle.git
   ```

2. Open RStudio and load one of the main scripts for a specific condition (e.g. `model_full.R` inside a condition folder).

3. Each script sets its working directory automatically and includes comments guiding execution.

4. Output includes:

   * Parameter estimates
   * Identifiability profiles
   * Plots of model fits and reductions

---

## 📌 Citation

If you use this repository, please cite:

> Hermes, J. et al. (2025). *The Break in the Cycle: Inositol Pyrophosphate Fluxomics Disentangled via Mathematical Modelling*. bioRxiv. [https://doi.org/10.1101/2025.03.27.645712](https://doi.org/10.1101/2025.03.27.645712)


---

## 👤 Author

**Jacques Hermes**
Department of Data Management
University of Freiburg
[jacques.hermes@fdm.uni-freiburg.de](mailto:jacques.hermes@fdm.uni-freiburg.de)

---

## ⚖️ License

This project is licensed under the [MIT License](./LICENSE).
