# Data and code for Bauer et al. (2023) bioRxiv

*Markus Bauer* <a href="https://orcid.org/0000-0001-5372-4174"><img src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16"/></a>, Jakob K. Huber, and Johannes Kollmann <a href="https://orcid.org/0000-0002-4990-3636"><img src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16"/></a>

Data and code for:

Bauer M, Huber JK, Kollmann J (2023) __Fit by design: Developing substrate-specific seed mixtures for functional dike grasslands.__ &ndash; *bioRxiv*

[![DOI:10.1101/2023.03.09.530576](http://img.shields.io/badge/DOI-10.1101/2023.03.09.530576-informational.svg)](https://doi.org/10.1101/2023.03.09.530576)

**Study region**: [Experiment at River Danube](https://www.openstreetmap.org/#map=17/48.83977/12.88445) <br> <br> \## Content of the repository

1.  **Data**: the folder `data` contains
    -   `Raw` and `processed` data of the sites variables (.csv)
    -   `Raw` and `processed` data of the species' abundances (.csv)
    -   `Raw` and `processed` data of the species' traits (.csv)
    -   Raw and processed `raws/spatial` data (.shp)
    -   `photos` of the plots (.jpg)
2.  **Outputs**: the folder `outputs` contains
    -   The figures generated (.tiff)
    -   The tables generated (.html/.png)
    -   The models calculated (.Rdata)
    -   The summary statistics of the models (.csv)
3.  **R**: the folder `R` contains
    -   Scripts to calculate all models (.R)
    -   Scripts to generate all figures and tables (.R)
    -   Metadata script for creating EML file
    -   Folder for calculating habitat types (ESY)
4.  **Markdown**: the folder `markdown` contains
    -   Markdown documents of the analyses with model evaluations and comparisons (.md)

#### Package versioning

The used versions of R and the packages are saved in `2023_danube_dike_experiment/renv.lock`.

You can restore this state by executing `renv::restore()` in the console.

## Citation

[![CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](http://creativecommons.org/licenses/by/4.0/)

This work is licensed under a [Creative Commons Attribution 4.0 International License](http://creativecommons.org/licenses/by/4.0/).

When using the **data available** in this repository, please cite the original publication and the dataset.

**Publication**

> Bauer M, Huber JK, & Kollmann J (2023) Fit by design: Developing substrate-specific seed mixtures for functional dike grasslands. &ndash; *bioRxiv*. <https://doi.org/10.1101/2023.03.09.530576>

**Dataset**

> Bauer M, Huber JK & Kollmann J (2022) Data and code for Bauer et al. (2023) bioRxiv (v1.0.0) [Data set]. &ndash; *Zenodo*. <https://doi.org/10.5281/zenodo.7713396>

Contact [markus1.bauer\@tum.de](mailto:markus1.bauer@tum.de){.email} for any further information.
