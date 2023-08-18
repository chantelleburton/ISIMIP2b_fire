cut_results <- function(pdata,limits) {
	nlimits=length(limits)
	if (class(pdata)=="RasterBrick" && nlayers(pdata)==1) pdata=pdata[[1]]
	odata=pdata
	odata[]=nlimits+1
	for (i in nlimits:1) odata[pdata<=limits[i]]=i
	odata[is.na(pdata)]=NaN
	return(odata)
}
