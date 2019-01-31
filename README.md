## Japanese national spatial conservation prioritization using Zonation



## Introduction

This repository contains code and configuration files that help to understand how the spatial conservation prioritization analyses for and the result of the article [Spatial conservation prioritization for the East Asian islands: A balanced representation of multitaxon biogeography in a protected area network](https://doi.org/10.1111/ddi.12869). In other words, the repository contains two different types of content:

1. **Zonation setting files**

   This configuration files demonstrate what types of Zonation spatial conservation prioritization analyses were used. 
   

2. **R source code files**

   Code for creating the Zonation configuration files as well as analyzing the results.
   

The content of the repository links to [a sample data set hosted in Zenodo](https://doi.org/10.5281/zenodo.1462723). 

It is important to note that even when coupled with the data available in Zenodo, this repository **cannot be used to replicate** the original analyses of [the published study](https://doi.org/10.1111/ddi.12869). The reason is that the data release excludes a significant number of threatened species and completely omits all the species names because of species conservation concerns (see “Data availability” below for more information). Even so, the content of this repository will help you to understand how the Zonation analyses were constructed and how the results of the prioritization analyses were analyzed. In case you are interested in the full original analyses setups, please contact the following persons:

- Buntarou Kusumoto (<kusumoto.buntarou@gmail.com>)
- Joona Lehtomäki (<joona.lehtomaki@iki.fi>)



## Structure of the repository



#### R

R source files needed to 1) create the setting files, 2) analyzing the results and 3) creating figures.

#### templates

Template files needed for the creation of the Zontion settings files.

#### zsetup_release

The directory contains the Zonation settings files used in the original study. The only difference is, that they are based on [the sample data](https://doi.org/10.5281/zenodo.1462723) set, not the full data set. All the setting files are automatically created by the `R/01_zproject/01_create_zproject_for_release.R` script and the data files in Zenodo.



## Data availability

the data release, we have excluded threatened species (critically endangered, endangered, and  vulnerable species; 1,510 species) and equivalently rare species (distribution range < 2,000 km^2; 1,124 species)  from the viewpoint of rare species conservation. In the release, we also used species ID instead of species scientific name to avoid any risk of systematic poaching of common species which became a serious concern in Japan. 

For scientific uses of the original data including species names,  please contact Prof. Yasuhiro Kubota (<kubota.yasuhiro@gmail.com>).



## License

All code is licensed under [the MIT license](https://opensource.org/licenses/MIT). All other content under [the CC BY 4.0 license](https://creativecommons.org/licenses/by/4.0/).



## References

Lehtomäki J, Kusumoto B, Shiono T, Tanaka T, Kubota Y, Moilanen A. 2018. Spatial conservation prioritization for the East Asian islands: A balanced representation of multitaxon biogeography in a protected area network. Diversity & distributions **2**:18. [DOI: 10.1111/ddi.12869](https://doi.org/10.1111/ddi.12869).