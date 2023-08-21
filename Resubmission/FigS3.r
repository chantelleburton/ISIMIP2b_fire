library(ncdf4)
library(raster)
source("libs/sourceAllLibs.r")
sourceAllLibs("libs/")
graphics.off()


dir = "zenodo_dat/supplement/"
 
models = c("HADGEM2-ES", "GFDL-ESM2M","MIROC5", "IPSL-CM5A-LR")

hist_dir = "historic_on"
futr_dir = "RCP6.0_on"

hist_layers = list(1657:1728, 1:168)
hist_layers = list(1:240, NULL)
futr_layers = list(list(NULL, (2018-2006)*12 + (1:240)),
                   list(NULL, (2047-2006)*12 + (1:240)),
                   list(NULL, (2043-2006)*12 + (1:240)),
                   list((1728-35):1728, 1:204)) 
 #(2047-2006)*12 + (1:120), 

#(2018-2038)	(2047-2067)	(2043-2063)	(2003-2023)


index = c(1, 6, 5, 7, 8, 9, 4, 3, 2, 10)

variables = list("burnt_area.nc", "temp.nc", "precip.nc", "smc.nc", "cveg.nc", "trees.nc", "totalVeg.nc", "crop.nc", 
                 "pas.nc", "specific_humid.nc") # "humid.nc"
names = c("Burnt\nArea (%)", "Temperature\n ~DEG~C", "Precipitation\n(mm~yr-1~)", 
          "Top layer soil\nmoisture (%)", "Vegetation\nCarbon(kgC~m2~)",
          "Tree\nCover(%)",
          "Total Vegetation\nCover (%)", "Cropland\n(%)", "Pasture\n(%)", 
          "Specific\nHumidity (%)")


cols = list(rev(c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695')),
            c("#74a9cf", "#d0d1e6", '#fff7ec','#fdd49e','#fc8d59','#ef6548','#d7301f','#b30000','#7f0000'),
            c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061'),
            c('#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3','#f5f5f5','#c7eae5','#80cdc1','#35978f','#01665e','#003c30'),
            c('#c51b7d','#de77ae','#f1b6da','#fde0ef','#f7f7f7','#e6f5d0','#b8e186','#7fbc41','#4d9221'),
            c('#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3','#f5f5f5','#c7eae5','#80cdc1','#35978f','#01665e','#003c30'),
            c('#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3','#f5f5f5','#c7eae5','#80cdc1','#35978f','#01665e','#003c30'),
            c('#7f3b08','#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7','#d8daeb','#b2abd2','#8073ac','#542788','#2d004b'),
            c('#7f3b08','#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7','#d8daeb','#b2abd2','#8073ac','#542788','#2d004b'),
            c('#8c510a','#bf812d','#dfc27d','#f6e8c3','#f5f5f5','#c7eae5','#80cdc1','#35978f','#01665e'))

levels = list(c(-20, -10, -5, -2, -1,  1, 2, 5, 10, 20),
              c(-0.25, 0, 0.25, 0.5, 1.0, 1.25, 1.5, 1.75)*2, 
              c(-400, -200, -100, -50, -20, -10, 0, 10, 20, 50, 100, 200, 400), 
              c(-0.5, -0.2, -0.1, -0.05, -0.02,  0.02, 0.05, 0.1, 0.2, 0.5),
              c(-8, -6, -4, -2, -1, 0, 1, 2, 4, 7, 8),#c(-4, -2, -1, -0.5, -0.02, 0, 0.02, 0.5, 1, 2, 4),
              c(-40, -20, -10, -5, -2, -1, 0, 1, 2, 5, 10, 20, 40),
              c(-40, -20, -10, -5, -2, -1, 0, 1, 2, 5, 10, 20, 40),
              c(-0.5, -0.2, -0.1, -0.05, -0.02, -0.01,0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)*100,
              c(-0.5, -0.2, -0.1, -0.05, -0.02, -0.01,0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)*100,
              c(-200, -150, -100, -50, 0, 50, 100, 150, 200)/25)
        


forVar <- function(variable, name, col, level) {
    print(variable)
    forModel <- function(model, futr_layer) {
        grabLayer <- function(rs, layers) {
            if (is.null(layers)) return(NULL)
            return(rs[[layers]])
        }

        addLayer.null <- function(layers) {
            r1 = grabLayer(his, layers[[1]])
            r2 = grabLayer(rcp, layers[[2]])
            if (is.null(r1)) return(r2)
            if (is.null(r2)) return(r1)
            return(addLayer(r1, r2))
        }
        his = brick(paste(dir, hist_dir, model, variable, sep = '/'))        
        rcp = brick(paste(dir, futr_dir, model, variable, sep = '/'))
        hist = addLayer.null(hist_layers)
        futr = addLayer.null(futr_layer )
#addLayer(grabLayer(his, hist_layers[[1]]), grabLayer(his, hist_layers[[2]]))
 #       futr = addLayer(grabLayer(his, futr_layer [[1]]), grabLayer(his, futr_layer [[2]]))
         
        
        dat = mean(futr - hist)
        if (variable == "precip.nc") dat = dat*60*60*24*365
        if (variable == "humid.nc") dat = -1000*dat
        if (variable == "cveg.nc") dat = dat
        if (variable == "burnt_area.nc") dat = dat*60*60*24*365*100
        if (variable == "smc.nc") dat = dat/1000
        if (any(variable == c("totalVeg.nc", "trees.nc", "crop.nc",  "pas.nc"))) dat = dat * 100
        
        plotStandardMap(dat, '', level, col)
        if (model == models[1]) mtext.units(side = 2, name, line = -0.5)
        
        if (variable == variables[1]) mtext(model, side = 3)
    }
    
    mapply(forModel, models, futr_layers)
    legendColBar(c(0.3, 0.75), c(0.1, 0.9), col, level, F, T, T, oneSideLabels=F)

}
graphics.off()


lmat = cbind(0, t(matrix(1:(5*length(variables)), nrow = 5)), 0)
lmat = rbind(0, lmat, 0)
widths = c(0.2, 1, 1, 1, 1, 0.3, 0.01)
heights = c(0.2, rep(1, length(variables)), 0.1)

png("figs/bioclimate_maps.png", height = 2*sum(heights)*0.6, width = 3*sum(widths),
    res = 300, units = 'in')
layout(lmat, widths = widths, heights = heights)
par(mar = rep(0, 4))

mapply(forVar, variables[index], names[index], cols[index], levels[index])
dev.off()
