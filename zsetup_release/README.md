## Japanese national spatial conservation prioritization using Zonation

We modeled the distributions of 6235 species (amphibians, birds, freshwater fish, mammals, plants and reptiles) using 4 389 489 occurrence data points. We then applied the Zonation software for spatial conservation prioritization (SCP). First, we analyzed each taxon individually to understand baseline priority patterns. Second, we combined all taxa into an inclusive analysis to identify most important PAN expansions. Human influence was used as a proxy for both disturbance and potential opportunity cost of PAN expansion. 

## Zonation setup structure

Each taxon (amphibians, birds, freshwater fish, mammals, plants, reptiles) has 
the following variants:

* `YY_XXX_caz`
* `YY_XXX_abf`
* `YY_XXX_caz_wgt`
* `YY_XXX_abf_wgt`
* `YY_XXX_caz_wgt_con`
* `YY_XXX_abf_wgt_con`
* `YY_XXX_abf_wgt_con_hm2`
* `YY_XXX_caz_wgt_con_hm2`
* `YY_XXX_abf_wgt_con_hm3`
* `YY_XXX_caz_wgt_con_hm3`


`YY` is a running ID number and `XXX` stands for the taxon code, which are the 
following:

| Code | Description               |
|------|---------------------------|
| amp  | Amphibians                |
| bir  | Birds                     |
| frf  | Freshwater fish           |
| mam  | Mammals                   |
| pla  | Plants                    |
| rep  | Reptiles                  |

The other codes stand for the following:

| Code | Description               |
| ---- | ------------------------- |
| caz  | Core-area Zonation        |
| abf  | Additive benefit function |
| wgt  | Weights                   |
| con  | Condition layer           |
| hm3  | PA mask with 3 levels     |


All variants have groups enabled, grouping based on the IUCN Red-List category:

| Code | Group                     |
|------|---------------------------|
| LC   | 1                         |
| NT   | 2                         |
| VU   | 3                         |
| EN   | 4                         |
| CR   | 5                         |
| EW   | 6                         |
| EX   | 7                         |
| DD   | 8                         |

