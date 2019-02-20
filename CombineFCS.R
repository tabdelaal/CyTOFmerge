CombineFCS <- function (FCSfile1,RelevantMarkers1,FCSfile2,RelevantMarkers2,arcsinhTrans=FALSE)
{
  # This function can be used to combine 2 FCS files having a set of shared
  # markers and return one dataframe with the total number of cells is equal
  # to the summation of cells in both FCS files, with each cell has an
  # extended number of measured markers.
  #
  # Input description:
  # FCSfile1: Name (full name with directory) of the first FCS file
  #
  # RelevantMarkers1: integers array indicating the markers (columns) indices to
  #                   be used from FCSfile1.
  #
  # FCSfile2: Name (full name with directory) of the second FCS file
  #
  # RelevantMarkers2: integers array indicating the markers (columns) indices to
  #                   be used from FCSfile2.
  #
  # arcsinhTrans: True = apply arcsinh transformation with cofactor of 5
  #               prior to files combination.
  #               False (default) = no transformation applied.
  #
  # Important Notes: The shared markers short names (PnN or name2 or parameters@data$desc) 
  #                  must be the same in both FCS files, in order to the function to
  #                  identify them and use them for combination.
  #
  # For citation and further information please refer to this publication:
  # "CyTOFmerge: Integrating mass cytometry data across multiple panels"
  
  # Load data and get the Markers names of both files
  library(flowCore)
  FCS1 = read.FCS(FCSfile1,transformation = FALSE,truncate_max_range = FALSE)
  colnames(FCS1@exprs) <- FCS1@parameters@data$desc
  FCS1.data = as.data.frame(FCS1@exprs)[,RelevantMarkers1]
  VarNames1 = colnames(FCS1.data)
  
  FCS2 = read.FCS(FCSfile2,transformation = FALSE,truncate_max_range = FALSE)
  colnames(FCS2@exprs) <- FCS2@parameters@data$desc
  FCS2.data = as.data.frame(FCS2@exprs)[,RelevantMarkers2]
  VarNames2 = colnames(FCS2.data)
  
  # apply arcsinh transformation
  if(arcsinhTrans){
    FCS1.data = asinh(FCS1.data/5)
    FCS2.data = asinh(FCS2.data/5)
  }
  
  # Find shared Markers
  Matches = as.vector(matrix(0,ncol = length(VarNames1)))
  for (i in c(1:length(VarNames1))){
    if(any(VarNames1[i]==VarNames2))
      Matches[i] = which(VarNames1[i]==VarNames2)
  }
  
  # Reorder data
  Shared_Index = which(Matches>0)
  Data1.nonshared = FCS1.data[,-Shared_Index]
  VarNames1.nonshared = VarNames1[-Shared_Index]
  Data2.nonshared = FCS2.data[,-Matches[Shared_Index]]
  VarNames2.nonshared = VarNames2[-Matches[Shared_Index]]
  
  # This is the order data {Shared  Non_Shared}
  FCS1.data = cbind(FCS1.data[,Shared_Index], Data1.nonshared)
  VarNames1 = c(VarNames1[Shared_Index], VarNames1.nonshared)
  FCS2.data = cbind(FCS2.data[,Matches[Shared_Index]], Data2.nonshared)
  VarNames2 = c(VarNames2[Matches[Shared_Index]], VarNames2.nonshared)
  
  # Combine data
  m = length(Shared_Index)           # Number of shared markers
  
  IDX1 = FNN::get.knnx(FCS2.data[,1:m],FCS1.data[,1:m], k = 50, algorithm = "kd_tree")
  IDX1 = IDX1$nn.index
  IDX2 = FNN::get.knnx(FCS1.data[,1:m],FCS2.data[,1:m], k = 50, algorithm = "kd_tree")
  IDX2 = IDX2$nn.index
  
  Data.combine.1 = matrix(0,nrow = dim(FCS1.data)[1],ncol = dim(FCS1.data)[2]+dim(FCS2.data)[2]-m)
  Data.combine.2 = matrix(0,nrow = dim(FCS2.data)[1],ncol = dim(FCS1.data)[2]+dim(FCS2.data)[2]-m)
  
  for (i in c(1:dim(FCS1.data)[1])){
    Data.combine.1[i,1:m] = as.matrix(FCS1.data[i,1:m])
    Data.combine.1[i,(m+1):dim(FCS1.data)[2]] = as.matrix(FCS1.data[i,(m+1):dim(FCS1.data)[2]])
    Data.combine.1[i,(dim(FCS1.data)[2]+1):dim(Data.combine.1)[2]] = as.matrix(apply(FCS2.data[IDX1[i,],(m+1):dim(FCS2.data)[2]],2,median))
  }
  
  for (i in c(1:dim(FCS2.data)[1])){
    Data.combine.2[i,1:m] = as.matrix(FCS2.data[i,1:m])
    Data.combine.2[i,(m+1):dim(FCS1.data)[2]] = as.matrix(apply(FCS1.data[IDX2[i,],(m+1):dim(FCS1.data)[2]],2,median))
    Data.combine.2[i,(dim(FCS1.data)[2]+1):dim(Data.combine.2)[2]] = as.matrix(FCS2.data[i,(m+1):dim(FCS2.data)[2]])
  }
  
  Data.combine = rbind(Data.combine.1,Data.combine.2)
  VarNames.combine = c(VarNames1,VarNames2[(m+1):length(VarNames2)])
  
  Data.combine = as.data.frame(Data.combine)
  colnames(Data.combine) = VarNames.combine
  return (Data.combine)
}