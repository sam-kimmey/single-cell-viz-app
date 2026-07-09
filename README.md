# CellViz

**An R Shiny app for interactive visualization of spatial single-cell data.**

CellViz makes it easy to explore spatial single-cell data across one or multiple ROIs (regions of interest).

---

## Features

- Visualize spatial single-cell data for a single ROI or across multiple ROIs
- Overlay cell phenotype and neighborhood assignments directly on spatial plots
- Explore biaxial or spatial expression plots of biomarkers
- Interactive, zoomable plots powered by `plotly`
- **Brush to subset cells** — draw a selection directly on the graph to isolate cells in a region and investigate areas of interest in greater detail

---

## Input Data

CellViz accepts a single CSV file containing single-cell data (upload limit: **100 MB**).

| Column          | Required? | Description                     |
|-----------------|-----------|----------------------------------|
| `roi_id`        | Required  | Unique identifier for each ROI  |
| `phenotype`     | Optional  | Cell phenotype assignment - up to 36 phenotypes    |
| `neigh_kmeans`  | Optional  | Cell neighborhood assignment    |
| `sample_group1` | Optional  | Scientific sample groupings     |
| `sample_group2` | Optional  | Scientific sample groupings     |

---

## Running Locally

### Requirements

- R version 4.5.3
- R packages:
  - `shiny`
  - `data.table`
  - `ggplot2`
  - `plotly`
  - `shinymanager`
  - `tidyverse`
  - `pals`
  - `paletteer`
  - `Polychrome`

Install the required packages in R:

```r
install.packages(c(
  "shiny", "data.table", "ggplot2", "plotly", "shinymanager",
  "tidyverse", "pals", "paletteer", "Polychrome"
))
```

### 1. Clone the repository

```bash
git clone https://github.com/sam-kimmey/single-cell-viz-app.git
```

### 2. Navigate to the repository

```bash
cd single-cell-viz-app
```

### 3. Launch the app

```r
shiny::runApp("app.R")
```

---

## Authors

Developed by **Josh Kramer**, **Vini Karumuru**, and **Sam Kimmey**

---
### License and Trademark

The source code of this project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

However, the license does strictly apply to the underlying software code only. It does NOT grant permission to
use the Oregon Physics LLC or MIBIscope name, logos, trademarks, or official branding assets.
All rights to corporate branding, design assets, and trademarks are explicitly reserved by Oregon Physics LLC.
