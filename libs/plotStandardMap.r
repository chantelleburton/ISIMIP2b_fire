source("libs/return_multiple_from_functions.r")
library(mapproj)

plotStandardMap <- function(x, txt = '', limits, cols, e = NULL, recrop_e = TRUE, 
                            y_range = c(-60, 90), limits_error = c(0.5),
                            ePatternRes = 30,  ePatternThick = 0.6, add = FALSE,
                            projection = NULL,...) {
    if (nlayers(x) == 1) x = x[[1]]
    mask = raster('data/seamask.nc')
    mask = raster::resample(mask, x)
    x[mask != 2] = NaN
    if(is.null(e) && nlayers(x) > 1) {
        if (nlayers(x) == 2) {
            e = 1-x[[2]]
            x = x[[1]]
        } else e = 1-sd.raster(x)
    }
    if(!is.null(e)) e[mask != 2] = NaN
    
    FUN <- function(...) {
        plot_raster_from_raster(x, y_range = y_range, limits = limits, cols = cols,
                                transpose = FALSE, srt = 0, add_legend = FALSE,
                                quick = TRUE, e = e, interior = FALSE,
                                ePatternRes = ePatternRes, ePatternThick = ePatternThick,
                                limits_error = limits_error, projection = projection,...)
    }
    
    FUN(add = add,...)
    addCoastlineAndIce2map()    
    if (!is.null(projection)) {
        FUN(add = TRUE, ...)
        plot_raster_from_raster(mask, y_range = c(-60, 90), limits = c(0.5),
                                cols = c("white", make.transparent("white", 1)),
                                quick = TRUE, interior = FALSE, add_legend = FALSE,
                                coast.lwd = NULL, add = TRUE, projection = projection,...)
    #contour(mask, add = TRUE, drawlabels = FALSE, lwd = 0.5)
    }
    
    mtext(txt, side = 2, line = -2, adj = 0.1)
}

add_icemask <- function() {
	icemask = raster('data/icemask.nc')
	plot_raster_from_raster(icemask, add = TRUE, cols = c('#FFFFFFFF', 'grey'), y_range = c(-60, 90),
						    limits = c(-0.5, 0.5), add_legend = FALSE, interior = FALSE, coast.lwd = 0.67)#, coast.lwd = NULL)
}

addCoastlineAndIce2map <- function() {
    add_icemask()
    
    mask = raster('data/seamask.nc')
    plot_raster_from_raster(mask, y_range = c(-60, 90),
                            add = TRUE, col = c("white", "transparent"), 
                            limits = c(0.5), quick = TRUE, interior = FALSE,
                            coast.lwd = 2, add_legend = FALSE)
    #
    #contour(mask, add = TRUE, drawlabels = FALSE, lwd = 0.5)  

    ployBox <- function(x, y)
        polygon(c(x[1], x[2], x[2], x[1]), c(y[1], y[1], y[2], y[2]), col = "white", border = "white")
    
    ployBox(c(-180, -90), c(-60, 0))
    ployBox(c(-180, -120), c(-60, 25))
    ployBox(c(-50, -19), c(10, 25))
    ployBox(c(-50, -13.5), c(27.5, 34))
    ployBox(c(115, 125), c(-8, -7))
    ployBox(c(104, 111), c(2.5, 8))
    ployBox(c(122, 128), c(2.5, 5))  
    ployBox(c(52.5, 74), c(-25, 12))  
    ployBox(c(89, 95), c(6.5, 15.5))  
    ployBox(c(90, 107.5), c(-17.5, -10))  
    ployBox(c(-32.5, 7.5), c(-40, 2.5)) 
    
    ployBox(c(-61, -23), c(22.5, 44))  
    ployBox(c(-67.5, -57.5), c(31, 33))  
    ployBox(c(130, 180), c(1, 30))  
    ployBox(c(155, 180), c(-5, 5))  
    ployBox(c(-62.5, 142.5), c(-42.5, -60))  
}

StandardLegend <- function(cols, limits, dat, extend_max = TRUE, oneSideLabels = TRUE, transpose = FALSE,
                           plot_loc = NULL, srt = 0, add = FALSE, ...) { 
        
        if (is.null(plot_loc)) {
            if (add) plot_loc = c(0.32, 0.87, 0.025, 0.07)
            else plot_loc = c(0.05, 0.95, 0.6, 0.7)
        }
        add_raster_legend2(cols, limits, dat = dat, add = add,
                           transpose = transpose, srt = srt, oneSideLabels= oneSideLabels,
                           plot_loc = plot_loc,
                           ylabposScling = 1, extend_max = extend_max, ...)
}
