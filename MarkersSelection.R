MarkersSelection <- function (FCSfolder,RelevantMarkers,subsets=FALSE,arcsinhTrans=FALSE,SelectionThreshold=0.8,CorrThreshold=0.7,DownSampleSize=100000)
{
  # MarkerSelection function can be used to provide a reduced set of markers
  # that can describe the cellular composition of a Cytof dataset, this
  # reduced set of markers can be then used as shared markers while designing
  # your 2nd Cytof panel.
  #
  # Input description
  #
  # FCSFolder: extension of the folder having the FCS files of the Cytof data
  #            from which you want to make the selection
  #
  # RelevantMarkers: integers array indicating the markers (columns) indices to
  #                  be used in the analysis.
  #
  # subsets:  False (default) = FCS files are the analysed samples output 
  #           from the Cytof (no clustering involved), which means that the 
  #           selection will only be done using the Distance and the Nearest 
  #           Neighbor scores (excluding the Cluster score).
  #
  #           True = FCS files are the subsets (clusters) defining the 
  #           different cell populations (if your data is already analysed 
  #           and clustered), now all three scores will be included in the
  #           selection.
  # 
  # arcsinhTrans: True = apply arcsinh transformation with cofactor of 5
  #               prior to any analysis.
  #               False (default) = no transformation applied.
  #
  # Selectionthreshold: It can be either a value from [0,1[ which sets the 
  #                     threshold for the used scores to decide the number of
  #                     markers, the default is (0.8).
  #                     Or it can be an integer value from 2 to "the total 
  #                     number of markers" - 1. In this case the function 
  #                     will return the desired number of reduced markers 
  #                     regardless of the performance scores.
  #
  # CorrThreshold: A value from [0,1[ which sets the threshold used to 
  #                exclude highly correlated markers, the default is (0.7).
  #
  # DownSampleSize: The number of cells randomly selected from the data, to
  #                 speed up the calculation of the pairwise average 
  #                 Euclidean distance, the default is (100,000).
  #                 False = no downsampling applied, the full data is used to
  #                 calculate the the pairwise average Euclidean distance.
  #
  # Markers = MarkersSelection(FCSfolder,c(9:47))
  #           This will return the minimal markers list that can acheive a
  #           minimal Selectionthreshold of 0.8 for the Distance and the 
  #           Nearest Neighbor scores (exluding the Cluster score). 
  #
  # Markers = MarkersSelection(FCSfolder,c(9:47),subsets=TRUE)
  #           This will return the minimal markers list that can acheive a
  #           minimal Selectionthreshold of 0.8 for all three scores 
  #           (including Cluster Score).
  #
  # Markers = MarkersSelection(FCSfolder,c(9:47),SelectionThreshold = 0.7)
  #           This will return the minimal markers list that can acheive a
  #           minimal Selectionthreshold of 0.7 for the Distance and the 
  #           Nearest Neighbor scores (exluding the Cluster score).
  #
  # Markers = MarkersSelection(FCSfolder,c(9:47),subsets=TRUE,SelectionThreshold = 0.7)
  #           This will return the minimal markers list that can acheive a
  #           minimal Selectionthreshold of 0.7 for all three scores 
  #           (including Cluster Score).
  #
  # Markers = MarkersSelection(FCSfolder,c(9:47),SelectionThreshold = 10)
  #           This will return the top 10 markers that can describe the
  #           dataset structure regadless of their performance scores.
  #
  # For citation and further information please refer to this publication:
  # "CyTOFmerge: Integrating mass cytometry data across multiple panels"
  
  # check SelectionThreshold validity
  if((SelectionThreshold >= 1 && SelectionThreshold < 2)|| (SelectionThreshold >=2 && SelectionThreshold != floor(SelectionThreshold))){
    stop('SelectionThreshold must be a decimal value between [0,1[, or an integer >= 2')}
  
  # load data
  CyTOF.data = data.frame()
  Subset.labels = data.frame()
  count = 0
  files = list.files(path = FCSfolder, pattern = '.fcs',full.names = TRUE)
    
  library(flowCore)
  for (i in files){
    Temp <- read.FCS(i,transformation = FALSE, truncate_max_range = FALSE)
    colnames(Temp@exprs) <- Temp@parameters@data$desc
    CyTOF.data = rbind(CyTOF.data,as.data.frame(Temp@exprs)[,RelevantMarkers])
    if(subsets){
      count = count + 1
      Subset.labels = rbind(Subset.labels,array(1*count,dim = c(dim(Temp@exprs)[1],1)))
    }
  }
  Subset.labels = as.vector(as.matrix(Subset.labels))
  VarNames = colnames(CyTOF.data)
  
  # apply arcsinh transformation
  if(arcsinhTrans)
    CyTOF.data = asinh(CyTOF.data/5)
  
  # Filter highly correlated markers
  R = cor(CyTOF.data,method = "pearson")
  del = vector()
  for (i in c(1:(dim(R)[1]-1))){
    for (j in c((i+1):dim(R)[1])){
      if(abs(R[i,j]) > CorrThreshold){
        Vi = var(CyTOF.data[,i])
        Vj = var(CyTOF.data[,j])
        if(Vi > Vj)
          del = c(del,j)
        else
          del = c(del, i)
      }
    }
  }
  del = unique(del)
  Deleted.data = CyTOF.data[,del]
  if(!pracma::isempty(del)){
    Deleted.markers = VarNames[del]
    CyTOF.data = CyTOF.data[,-del]
    VarNames = VarNames[-del]
    print('Remove by preprocessing')
    Deleted.markers
  }
  
  # Check if the Selectionthreshold is a desired number of markers
  if(SelectionThreshold >= 2 && SelectionThreshold == floor(SelectionThreshold)){
    PCA = princomp(CyTOF.data)
    Markers.importance = ((PCA$loadings^2)[,1:SelectionThreshold]) %*% ((PCA$sdev^2)[1:SelectionThreshold])
    Sorted.importance = sort(Markers.importance,decreasing = TRUE, index.return = TRUE)
    Markers.list = VarNames[Sorted.importance$ix[1:SelectionThreshold]]
    return(Markers.list)
  }
  
  # Continue here if the SelectionThreshold is between [0,1[
  rank = matrix(0,nrow = dim(CyTOF.data)[2], ncol = dim(CyTOF.data)[2])
  PCA = princomp(CyTOF.data)
  for (i in c(2:dim(CyTOF.data)[2])){
    Markers.importance = ((PCA$loadings^2)[,1:i]) %*% ((PCA$sdev^2)[1:i])
    Sorted.importance = sort(Markers.importance,decreasing = TRUE, index.return = TRUE)
    rank[,i] = Sorted.importance$ix
  }
  
  Cumulative.Variance = cumsum(PCA$sdev^2 / sum(PCA$sdev^2))
  
  # Evaluation to find m
  Cluster_Score = 0
  Distance_Score = 0
  NN_Score = 0
  
  X = sample(dim(CyTOF.data)[1],dim(CyTOF.data)[1])
  if((dim(CyTOF.data)[1]%%2)==0){
    X1 = X[seq(1,length(X)-1,2)]
    X2 = X[seq(2,length(X),2)]
  } else {
    X1 = X[seq(1,length(X),2)]
    X2 = X[seq(2,length(X)-1,2)]
  }
  
  if(DownSampleSize == FALSE || DownSampleSize > dim(CyTOF.data)[1])
    DownSampleSize = dim(CyTOF.data)[1]
  
  Data = cbind(CyTOF.data,Deleted.data)
  Data = Data[sample(dim(Data)[1],DownSampleSize),]
  
  # Calculate the average pairwise Euclidean distance (takes long time)
  Euc_Dist = as.vector(matrix(0,nrow = dim(Data)[1]))
  for (i in c(1:dim(Data)[1])){
    for (j in c(i:dim(Data)[1])){
      Euc_Dist[i]=Euc_Dist[i]+sqrt(sum((Data[i,] - Data[j,])^2))
      Euc_Dist[j]=Euc_Dist[j]+sqrt(sum((Data[i,] - Data[j,])^2))
    }
  }
  Euc_Dist_Avg = Euc_Dist/dim(Data)[1]
  
  m = min(which(Cumulative.Variance > SelectionThreshold)) - 1
  
  while(Distance_Score < SelectionThreshold || NN_Score < SelectionThreshold || (subsets && Cluster_Score < SelectionThreshold)){
    m = m + 1
    # Simulate the two overlapping datasets
    Cut_Index = floor((dim(CyTOF.data)[2]+dim(Deleted.data)[2]-m)/2)+m
    Data1 = as.matrix(CyTOF.data[X1,c(rank[1:m,m],rank[(m+1):Cut_Index,m])])
    Data2 = as.matrix(cbind(CyTOF.data[X2,c(rank[1:m,m],rank[(Cut_Index+1):dim(CyTOF.data)[2],m])],Deleted.data[X2,]))
    Data.Sorted = as.matrix(cbind(CyTOF.data[c(X1,X2),rank[,m]],Deleted.data[c(X1,X2),]))
    if(subsets){
      Subset.labels.1 = Subset.labels[X1]
      Subset.labels.2 = Subset.labels[X2]
      Subset.labels.Sorted = Subset.labels[c(X1,X2)]
    }
    
    # Find the 1st NN in the original dataset
    NN = FNN::get.knn(Data.Sorted, k = 1, algorithm = "kd_tree")
    NN_dist = NN$nn.dist
    
    # Find the 50 neighbors from one dataset to the other
    IDX1 = FNN::get.knnx(Data2[,1:m],Data1[,1:m], k = 50, algorithm = "kd_tree")
    IDX1 = IDX1$nn.index
    IDX2 = FNN::get.knnx(Data1[,1:m],Data2[,1:m], k = 50, algorithm = "kd_tree")
    IDX2 = IDX2$nn.index
    
    Data.combine.1 = matrix(0,nrow = dim(Data1)[1],ncol = dim(Data.Sorted)[2])
    Data.combine.2 = matrix(0,nrow = dim(Data2)[1],ncol = dim(Data.Sorted)[2])
    if(subsets){
      Subset.labels.combine.1 = as.vector(matrix(0,nrow = dim(Data1)[1]))
      Subset.labels.combine.2 = as.vector(matrix(0,nrow = dim(Data2)[1]))
    }
    
    getmode <- function(x) {
      uniqx <- unique(x)
      uniqx[which.max(tabulate(match(x, uniqx)))]
    }
    
    # Combine datasets
    for (i in c(1:dim(Data1)[1])){
      Data.combine.1[i,1:m] = Data1[i,1:m]
      Data.combine.1[i,(m+1):Cut_Index] = Data1[i,(m+1):dim(Data1)[2]]
      Data.combine.1[i,(Cut_Index+1):dim(Data.combine.1)[2]] = apply(Data2[IDX1[i,],(m+1):dim(Data2)[2]],2,median)
      if(subsets){
        Subset.labels.combine.1[i] = getmode(Subset.labels.2[IDX1[i,]])
      }
    }
    
    for (i in c(1:dim(Data2)[1])){
      Data.combine.2[i,1:m] = Data2[i,1:m]
      Data.combine.2[i,(m+1):Cut_Index] = apply(Data1[IDX2[i,],(m+1):dim(Data1)[2]],2,median)
      Data.combine.2[i,(Cut_Index+1):dim(Data.combine.2)[2]] = Data2[i,(m+1):dim(Data2)[2]]
      if(subsets){
        Subset.labels.combine.2[i] = getmode(Subset.labels.1[IDX2[i,]])
      }
    }
    
    Data.combine = rbind(Data.combine.1,Data.combine.2)
    Subset.labels.combine = c(Subset.labels.combine.1,Subset.labels.combine.2)
    
    # Calculate evaluation scores
    if(subsets){
      Cluster_Score = sum(Subset.labels.combine == Subset.labels.Sorted)/length(Subset.labels.combine)
    }
    
    euc_dist = sqrt(rowSums((Data.Sorted-Data.combine)^2))
    
    Distance_Score = (median(Euc_Dist_Avg)-median(euc_dist))/median(Euc_Dist_Avg)
    
    NN_Score = sum(NN_dist > euc_dist)/length(euc_dist)
  }
  
  Markers.list = VarNames[rank[1:m,m]]
  return(Markers.list)
}