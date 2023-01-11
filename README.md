# Danube dike experiment

_Markus Bauer <a href="https://orcid.org/0000-0001-5372-4174"><img src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height = "16"></a>, Jakob Huber, and Johannes Kollmann <a href="https://orcid.org/0000-0002-4990-3636"><img src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height = "16"></a>_  

***

**Journal**: [XXX](https://www.???.??)

[![DOI:10.XXX](http://img.shields.io/badge/DOI-10.14471/2018.38.006-informational.svg)](https://doi.org/10.XXX)

**Study region**: [Experiment at River Danube](https://www.openstreetmap.org/#map=17/48.83977/12.88445)

## Content of the repository

1. __Data__: the folder `data` contains  
    * `Raw` and `processed` data of the sites variables (.csv) 
    * `Raw` and `processed` data of the species' abundances (.csv) 
    * `Raw` and `processed` data of the species' traits (.csv)
    * Raw and processed `raws/spatial` data (.shp)
    * `photos` of the plots (.jpg)
 
2. __Outputs__: the folder `outputs` contains  
    * The figures generated (.tiff)
    * The tables generated (.html/.png)
    * The models calculated (.Rdata)
    * The summary statistics of the models (.csv)
    
3. __R__: the folder `R` contains  
     * Scripts to calculate all models (.R)
    * Scripts to generate all figures and tables (.R)
    * Metadata script for creating EML file
    * Folder for calculating habitat types (ESY)
    
4. __Markdown__: the folder `markdown` contains 
    * Markdown documents of the analyses with model evaluations and comparisons (.md)

***

__Package versioning__

The used versions of R and the packages are saved in `2022_waste_bricks_trees/renv.lock`.

You can restore this state by executing `renv::restore()` in the console.
    
***

[![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg


When using the __data available__ in this repository, please cite the original publication and the dataset.  

__Publication__

> Bauer M, Huber J, & Kollmann J (Under preparation) XXX

__Dataset__

> Bauer M, Huber J & Kollmann J (2022) Data and code for Bauer et al. (submitted): Dike grassland experiment (v1.0.0) [Data set]. – *Zenodo*. https://doi.org/10.xxx

Contact markus1.bauer@tum.de for any further information.  

