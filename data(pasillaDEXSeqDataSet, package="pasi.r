data(pasillaDEXSeqDataSet, package="pasilla")
dxd <- estimateSizeFactors( dxd )
dxd <- estimateDispersions( dxd )

dxd2 <- dxd
assay(dxd2, "counts") <- assay(dxd2, "counts") + 1L


dxd <- testForDEU( dxd )
dxd <- estimateExonFoldChanges( dxd )


       dxd2 <- testForDEU( dxd2 )
       dxd2 <- estimateExonFoldChanges( dxd2 )
