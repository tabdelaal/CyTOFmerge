# CyTOFmerge
## Integrating mass cytometry data across multiple panels

The implementation is done in R and Matlab, provided in 2 functions:

1) MarkersSelection

MarkerSelection function can be used to provide a reduced set of markers that can describe the cellular composition of a Cytof dataset, this reduced set of markers can be then used as shared markers while designing your 2nd Cytof panel.

For full description, please check the MarkersSelection function description

2) CombineFCS

This function can be used to combine 2 FCS files having a set of shared markers and return one FCS file (matrix) with the total number of cells is equal to the summation of cells in both FCS files, with each cell has an extended number of measured markers.

For full description, please check the CombineFCS function description

Also, we provide a documentation (vignette) explaning the steps to reproduce the results obtained for the Vortex dataset, as an example, using R.

For citation and further information please refer to:
"CyTOFmerge: Integrating mass cytometry data across multiple panels"
