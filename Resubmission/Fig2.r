## set to data downloaded from zenodo
dir = 'zenodo_dat/'

library(ncdf4)
library(raster)
source("libs/sourceAllLibs.r")
sourceAllLibs("libs/")
graphics.off()


gwts = 'data/Tas_vs_year-ISIMIP.csv'

experiments = c("_nofire", "_fire")
models = c(HADGEM2 = "HADGEM2-ES", GFDL = "GFDL-ESM2M", MIROC = "MIROC5", 
            IPSL = "IPSL-CM5A-LR")
variable = "frac"

tfile0 = 'temp/'
futuresID = 'RCP6.0'

gfed = bound = brick("data/GFEDregions.nc")
bound[] = 0
bound[,2:ncol(gfed)] = bound[,2:ncol(gfed)] + (gfed[,2:ncol(gfed)] != gfed[,1:(ncol(gfed)-1)])
bound[2:nrow(gfed),] = bound[2:nrow(gfed),] + (gfed[2:nrow(gfed),] != gfed[1:(nrow(gfed)-1),])
bound[,3:ncol(gfed)] = bound[,3:ncol(gfed)] + (gfed[,3:ncol(gfed)] != gfed[,1:(ncol(gfed)-2)])
bound[3:nrow(gfed),] = bound[3:nrow(gfed),] + (gfed[3:nrow(gfed),] != gfed[1:(nrow(gfed)-2),])
bound = bound != 0


files = list.files(dir, full.names = TRUE)
openDat <- function(experiment, model, variable)   {
    file = files[grepl(experiment, files) & grepl(variable, files) & grepl(substr(model, 1, 4), files, ignore.case=TRUE)]
    if (length(file) != 1) browser()
    brick(file)
}

aa <- function(i, dat, ..., ncount = 11) {
    print(i)
    tfile = paste0(c(tfile0, ...,  ncount, i, '.nc'), collapse = '-')
    if (file.exists(tfile)) return(raster(tfile))
    
    print(tfile)
    out = mean(dat[[(i-ncount):i]])
    writeRaster(out, file = tfile, overwrite = TRUE)
}

openAllP <- function(...) {
    tfile = paste0(c(tfile0, 'running21all',  ..., '.nc'), collapse = '-')
    if (file.exists(tfile)) return(brick(tfile))
    
    datYr = openDat(...)
    dat = layer.apply(21:nlayers(datYr), aa, datYr, 'running21', ..., ncount = 20)
    dat = dat/mean(dat[[c(1, 20)]])
    dat = writeRaster(dat, file = tfile, overwrite = TRUE)
}

openMod <- function(...) 
    dat = lapply(experiments, openAllP, ...)


dats = lapply(models, openMod, variable = variable)


analyseCell <- function(i, dat, gwt, Temp = 1.5, ID) {
    ts1 = as.vector(dat[[1]][i]); ts2 = as.vector(dat[[2]][i])
    if (all(ts1 == 0) && all(ts2 == 0) ) return(rep(NaN, 4))
    
    tgrad  <- function(id, ts) diff(ts[(id-1):id]) > 0
    grad = tgrad(ID, ts1)
    impact = ts1[ID]
    tt2 = (ts2 - impact) > 0
    f_ID = 1+which(tt2[-1] != head(tt2, -1))
    #f_ID = which.min(abs(ts2 - impact))
    if (length(f_ID) == 0) f_ID = length(ts2)

    f_grad = sapply(f_ID, tgrad, ts2)
    
    if (length(f_ID) > 1) {
        gtest = which(f_grad == grad)
        if (length(gtest) > 0) {
            f_ID = f_ID[gtest]
            f_grad = f_grad[gtest]
        } 
        if (length(f_ID)>0) {
            gtest = which.min(abs(f_ID-ID))
            f_grad = f_grad[gtest]
            f_ID = f_ID[gtest]
        }
    }
    if (f_ID == 1) browser()
    return(c(impact, grad, f_ID, f_grad))
}

gwts = read.csv(gwts)[-1,]

analyseModel <- function(dat, model, Temp = 1.5, switch = FALSE) {
    if (switch) tfile = paste('outputs/degreeEquiv2-switch', model, Temp, '.nc', sep = '-')
    else tfile = paste('outputs/degreeEquiv2', model, Temp, '.nc', sep = '-')
    if (file.exists(tfile)) return(brick(tfile))
    mask = !is.na(dat[[1]][[1]])
    gwt = gwts[, grepl(substr(futuresID, 4, 6), names(gwts)) & 
                 grepl(substr(model, 1, 4), names(gwts), ignore.case = TRUE)]
    if (switch) dat = dat[2:1]
    index = which(mask[])#[1:42000]
    ncells = length(index)
    ntest = 1000
    analyseGroup <- function(i) {
        print(i/ncells)
        tfile = paste(tfile0, 'analyseGroup', model, Temp, i, ntest, '.Rd', sep = '-')
        
        if (file.exists(tfile) && Temp < 1.1) {load(tfile); return(out)}
        ids = index[i:(min(i+ntest-1, length(index)))]
        out = sapply(ids,analyseCell, dat, 
                  gwt, Temp = Temp, ID = which.min(abs(gwt - Temp)))
        
        save(out, file = tfile)
        out
    }
    vout = lapply(seq(ntest*0+1, length(index), by = ntest), analyseGroup)
    
    voutAll = do.call(cbind, vout)
    out = dat[[1]][[1:4]]
    out[] = NaN
    fgwt = state = out[[1]]
    out[index] = t(voutAll)
    fgwt[] = gwt[out[[3]][]]
    fgwt[out[[3]][] == length(gwt)] = -1 
    
    stateFinder <- function(id) {
        
        above = out[[1]] > 1; below = !above; up = out[[id]]; down = !up
        state[above & up] = 3
        state[above & down] = 4
        state[below & up] = 2
        state[below & down] = 1
         
        state
    }
    
    state = addLayer(stateFinder(2), stateFinder(4))
    outs  = addLayer(out[[1]], state, fgwt)
    
    outs = writeRaster(outs, file = tfile, overwrite = TRUE)
    return(outs)
}

forTemp <- function(Temp, ...)
    out = mapply(analyseModel, dats, models, Temp, ...)

Temps = c(1, 1.5, 2)
outs = lapply(Temps, forTemp)


##########
## plot ##
##########
colsI = c('#8e0152','#c51b7d','#de77ae','#f1b6da','#fde0ef','#f7f7f7','#e6f5d0','#b8e186','#7fbc41','#4d9221','#276419')
limsI = c(0.1,  0.2, 0.4, 0.6, 0.8, 1, 1.25, 1.6, 2.5, 5, 10)
colsS = c("#8c510a", "#c7eae5", "#01665e", "#f6e8c3")
nlimT = 8
colsT = list(c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7'),
             c('#d1e5f0','#92c5de','#4393c3','#2166ac','#053061'))

out = outs[[3]]
Temp = 2

limsT = c(-1.2, -1, -0.8, -0.6, -0.4, -0.2, -0.1, 0, 0.1, 0.2, 0.4, 0.6, 0.8, 1, 1.2)
maxT = min(apply(gwts[,grepl('6.0', names(gwts))], 2, max))

limsT = limsT[limsT <= (maxT-Temp)]
colsT = c(make_col_vector(colsT[[1]], limits = limsT[limsT <=0]),
          make_col_vector(colsT[[2]], limits = limsT[limsT >=0]))
plotForTemp <- function(Temp, out, outF) {
png(paste0("figs/degreeEquil-", Temp, ".png"), height = 6.5, width = 7.2, units = 'in', res = 300)
layout(rbind(1:4, 5:8, 9, 10:13, 14:17, 18, 18+matrix(1:20, ncol = 4), 39, 40:43, 44), 
       heights = c(1, 1, 0.3, 1, 1, 0.3, rep(1, 5), 0.3, 1, 0.3))

layout(cbind(1:5,6:10), heights = c(1, 1, 1, 1, 0.3))
par(mar = rep(0, 4), oma = c(2, 2, 2, 0))

plotImpact <- function(res, name, mttext = "Impact", addModName = TRUE){
    plotStandardMap(res[[1]], colsI, limsI)
    if (addModName) mtext(name, side = 3, line = 0)
    if (name == models[1]) mtext(mttext, side = 2)
}

plotState <- function(res, name, mttext, id = 2:3, addModName = TRUE, letter = '') {
    r = res[[id[1]]]
    unknown = (res[[id[1]]] ==1 & (res[[id[2]]] == 2 | res[[id[2]]] ==3)) |
              (res[[id[1]]] ==2 & (res[[id[2]]] == 1 | res[[id[2]]] ==4)) |
              (res[[id[1]]] ==3 & (res[[id[2]]] == 1 | res[[id[2]]] ==4)) |
              (res[[id[1]]] ==4 & (res[[id[2]]] == 2 | res[[id[2]]] ==3))
    r[unknown] = NaN
    plotStandardMap(r, colsS, seq(1.5, 3.5))
    plotStandardMap(unknown, cols = c("transparent", "#888888"), limits = c(0.5), add = TRUE)
    mtext(side = 3, adj = 0.1, letter)       
    
    if (name == names(models)[1]) mtext(mttext, side = 3)
    if (addModName) mtext(name, side = 2, line = 0)
    return(unknown)
}

unknown = mapply(plotState, out, names(models), "State", letter = letters[1:length(out)])

plot.new()
legend('center', ncol = 2, pch = 15, col = colsS, pt.cex = 2, cex = 1.15,
      c('reducing', 'recovering', 'increasing', 'diminishing'), bty = 'n')

plotTemp <- function(res, unknown, name, letter = '') {  
    
    gwt = res[[4]] 
    mask = gwt == -1
    gwt = gwt - Temp

    state = res[[2]]

    plotTempState <- function(i) {
        if (i > 0) {
            if (sum(state[] == i, na.rm = TRUE) == 0) {plot.new(); return()}
            gwt[state != i] = NaN
        }
        
        plotStandardMap(gwt, colsT, limsT)
        plotStandardMap(mask, cols = c("transparent", "#000099"), limits = c(0.5), add = TRUE)
        plotStandardMap(unknown, cols = c("transparent", "#888888"), limits = c(0.5),
                        add = TRUE)
        if (name == models[1]) mtext("Equivalent Temperature", side = 3)
        mtext(side = 3, adj = 0.1, letter)  
        #if (name == models[1]) mtext(c('All', 'reducing', 'recovering', 
        #                               'increasing', 'diminishing')[i+1], side = 2)
    }
    lapply(0, plotTempState)
}

mapply( plotTemp, out, unknown, models, 
                 letter = letters[(1+length(out)):(2*length(out))])

StandardLegend(colsT, limsT, out[[1]][[1]], extend_min = TRUE, units = '~DEG~C', oneSideLabels = FALSE, rightx = 0.8)
points(x = 0.9, y = 0.42, pch = 15, col = '#000099', cex = 2)
text(x = 0.9, y = 0.42, adj = c(0.5, 1.2), 'Beyond\nrun', xpd = NA, cex = 1.15)
dev.off()
}

mapply(plotForTemp, Temps, outs, outs)



