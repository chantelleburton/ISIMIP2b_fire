plot_raster_from_raster <-function(z,x_range=NULL,y_range=NULL,
				   limits=seq(min.raster(z),max.raster(z),length.out = 5),
                                   cols=rainbow(length(limits)),
				   coastline=NULL,coast.lwd=par("lwd"), coast.col = "black", 
				   add_legend=TRUE,legend_type='add_raster_legend2',
                                       legend.pos='bottomleft',
				   smooth_image=TRUE,smooth_factor=5,
				   projection=NULL,prj_parameters = NULL,orientation=NULL,
				   e=NULL,e_polygon=TRUE,regions = ".",fill.ocean=FALSE,
				   quick=FALSE,add=FALSE,...) {


    if (is.null(x_range)) x_range=c(xmin(z),xmax(z))
    if (is.null(y_range)) y_range=c(ymin(z),ymax(z))

    if (!is.null(limits) && length(cols)!=(length(limits)+1))
	cols=make_col_vector(cols,ncols=length(limits)+1)

    if (quick) {
	smooth_image=FALSE
	fill.ocean=FALSE
	e_polygon=FALSE
	map_db="world"
    } else map_db='worldHires'

    mar=par("mar")

    if (!is.null(coast.lwd))
	add = addMap(map_db,projection, prj_parameters = prj_parameters,x_range,y_range,regions,
                     orientation,add, lwd = coast.lwd, col = coast.col)
    if (smooth_image) z=disaggregate(z,smooth_factor,method="bilinear")

    plot_raster_map(z,x_range,y_range,limits,cols,coastline,coast.lwd, coast.col,
                    add_legend,legend_type,legend.pos,add=add,projection=projection,
		    orientation=orientation, prj_parameters = prj_parameters,
		    libs_path=libs_path,e=e,regions=regions,fill.ocean=fill.ocean,
		    e_polygon=e_polygon,
		    map_db=map_db,...)
    
    par(mar=mar)
}

addMap <- function(map_db, projection, prj_parameters = NULL, x_range,y_range,regions,
		   orientation, add, lwd, col) {

    if (is.null(projection)) {
    	test = try(map(map_db,interior=FALSE, xlim=range(x_range),ylim=range(y_range),
                       regions=regions,
		       mar=par("mar"), add =add, lwd = lwd, col = col),silent = TRUE)
    } else {
        if (projection == "robinson") {
            if (!exists("coast_lines_proj_robinson")) {
                coast_shapefile = "/data/users/dkelley/mapping_Data/ne_10m_coastline/ne_10m_coastline.shp"
                coast_lines = readOGR(coast_shapefile, layer = ogrListLayers(coast_shapefile))
            #unproj_proj4string = proj4string(coast_lines)
                robin_crs = CRS("+proj=robin +lon_0=0w")
                coast_lines_proj_robinson <<- spTransform(coast_lines, robin_crs)
            }
            xylim = project_roboinson(cbind(x = rep(c(0,x_range), 3), y = rep(c(0,y_range), each = 3)))
            xlim = range(xylim[,1]); ylim = range(xylim[,2])
            coast = crop(coast_lines_proj_robinson, c(xlim, ylim))
            
            test =(plot(coast, add =add, lwd = lwd, col = col))

            
            if (!add) {
                nseq = 100
                xseq = seq(x_range[1], x_range[2], length.out = nseq)
                yseq = seq(y_range[1], y_range[2], length.out = nseq)
                
                xseq = xseq[c(1:nseq, rep(nseq, nseq), nseq:1, rep(1, nseq))]
                yseq = yseq[c(rep(1, nseq), 1:nseq, rep(nseq, nseq), nseq:1)]
                xyseq = project_roboinson(cbind(x = xseq, y = yseq))
                polygon(xyseq[,1], xyseq[,2], border = NA, col = '#EEEEEE')
                plot(coast, add =TRUE, lwd = lwd, col = col)
                addLine <- function(x1, x2 = x1, y1, y2 = y1) {
                    xy = cbind(x = seq(x1, x2, length.out = 100),
                               y = seq(y1, y2, length.out = 100))
                    xy_p = project_roboinson(xy)
                    lines(xy_p[,1], xy_p[,2], lty = 2, col = 'grey', lwd = lwd)
                }
                xs = seq(x_range[1], x_range[2], by = 30)
                lapply(xs, addLine, y1 = y_range[1], y2 = y_range[2])

                ys = seq(y_range[1], y_range[2], by = 30)
                lapply(ys, addLine, x1 = x_range[1], x2 = x_range[2])
                
            }
        } else {
            
	    test = try(map(map_db,projection=projection ,interior=FALSE, 
                           orientation=orientation,
                           parameters = prj_parameters,
		           xlim=range(x_range),ylim=range(y_range),regions=regions, 
                           add = add, lwd = lwd, col = col), silent = TRUE)
        }
    }
    if (add) return(TRUE) else return(!class(test) == "try-error")
}
