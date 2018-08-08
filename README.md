# CyTOFmerge
## Merging of two CyTOF datasets with shared markers

The implementation is done in Matlab, provided in 2 functions:

1) function MarkersList = MarkersSelection( FCSFolder,subsets, threshold )

MarkerSelection function can be used to provide a reduced set of markers that can describe the cellular composition of a Cytof dataset, this reduced set of markers can be then used as shared markers while designing your 2nd Cytof panel.

For full description, check MarkersSelection.m

2) function CombineFCS(FCS1,FCS2, outputfilename)
This function can be used to combine 2 FCS files having a set of shared markers and return one FCS file with the total number of cells is equal to the summation of cells in both FCS files, with each cell has an extended number of measured markers.

For full description, please check CombineFCS.m

For citation and further information please refer to this publication:
"CyTOFmerge: Integrating mass cytometry data across multiple panels"
