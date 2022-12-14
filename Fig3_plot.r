###########
## setup ##
###########
library(raster)
source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs("../gitProjectExtras/gitBasedProjects/R/")
sourceAllLibs("../rasterextrafuns/rasterExtras/R/")
#library(rasterExtras)
#library(gitBasedProjects)
library(ncdf4)
sourceAllLibs("libs/")
graphics.off()

dat = read.csv("data/Tas_vs_year.xlsx - NBP.csv", header = FALSE, stringsAsFactors = FALSE)#[,-1]

models = c("HadGEM2-ES" = "H", "GFDL" = "G", "MIROC" = "M", "IPSL" = "I")

regions = c("Global", "BONA", "TENA", "CEAM", "NHSA", "SHSA", "EURO", "MIDE", 
                      "NHAF", "SHAF", "BOAS", "CEAS", "SEAS", "EQAS", "AUST")

accumulative = FALSE

start_year = 2000
nyears = 20
modTemp <- function(model) {
    id = which(dat[,2] == model) + 1
    if (model == 'H') id = id + 1
    out =  as.numeric(dat[id,-1])
    out[1] = 0
    out
}
temps = lapply(models, modTemp)

cols1 = c('#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3','#f5f5f5')
cols2 = c('#c7eae5','#80cdc1','#35978f','#01665e','#003c30')

years = as.numeric(dat[2,-1])
years[1] = years[2] - 1
forRegions <- function(region, Rid, yearsTest = TRUE, stripes = FALSE, letter = '') {
    
    if (region == "Global") idR = which(dat[,1] == "GlobTemp")
    else idR = which(dat[,2] == region)
    
    if (yearsTest) xs = lapply(temps, function(i) years)
    else xs = temps#lapply(temps, function(i)  i[tyrID,-1])
    xrange = range(unlist(xs))
    if (letter == '') mt = region else mt = paste0(letter, ') ', region)

    modelNBP <- function(model, accTest) {
        idM = which(as.character(dat[,2]) == model)
        id = idR[idR > idM][1]
        
        out = apply(dat[id+(2:1), -1], 2, as.numeric)
        if (accTest) {
            out = t(apply(out, 1, cumsum))
            for (i in 1:2) out[i,] = out[i,] - out[i, years == 2008]
        }
        return(out)
    }

    nbps = lapply(models, modelNBP, accumulative)
    tnbps = lapply(models, modelNBP, FALSE)

    dnps = lapply(nbps, function(i) i[2,] - i[1,])
    tdnps =  lapply(tnbps, function(i) i[2,] - i[1,])
    testKicker <- function(x, temp, syrs = start_year, offset = 10, xs = tdnps, id) {
        
        syr = which(years == syrs)
        x0 = x[syr - 1 +  1:nyears]
        doTest <- function(i) 
            wilcox.test(x0, xs[[1]][i + 1:nyears], paired = FALSE)[[3]]

        pvs = sapply(syr:(length(x) -nyears), doTest)
        #if (id == 3) browser()
        sig = which(pvs < 0.1)[1] + offset
        c(sig + syrs + nyears/2, temp[sig + syr])
    }
    
    fireSig = mapply(testKicker, tdnps, temps, id = 1:3, SIMPLIFY = FALSE)
    
    ddnps = lapply(tdnps, function(x) sapply(20:length(x), function(i) mean(x[(i-19):i])))
    fireSig_an = mapply(testKicker, ddnps, temps, MoreArgs = list(xs = ddnps), SIMPLIFY = FALSE)
    

    if (!stripes) {
        
        test_temps = seq(floor(min(unlist(temps))*10)/10,
                         min(sapply(temps, max))-0.05, by = 0.1)

        binExp <- function(id) {
            binT <- function(t1, dt = 0.1) {
                t1 = t1 - dt/2
                t2 = t1 + dt
                binMod <- function(nbp, temp) {                 
                    test = temp > t1 & temp <= t2
                    c(mean(nbp[id, test]), sum(test))
                }
                out = mapply(binMod, nbps, temps)
                
                c(mean(t1, t2), sum(out[2,]),
                  range(out[1,], na.rm = TRUE), mean(out[1,], na.rm = TRUE))
            }
            out = mapply(binT, test_temps)

        }
        bins = lapply( 1:2, binExp)
        
        addPoly <- function(dat, col) {
            polygon(c(dat[1,], rev(dat[1,])), c(dat[4,], rev(dat[3,])), 
                   col = make.transparent(col, 0.95), border =NA)# make.transparent(col, 0.95))
            
            lineLWD <- function(x, y, lwd, ...) 
                lapply(1:length(x),
                       function(i) lines(x[(i-1):i], y[(i-1):i], lwd = lwd[(i-1):i], ...))
            
            lineLWD(dat[1,], dat[5,], col = make.transparent(col, 0.6), 
                    lwd = 3*sqrt(dat[2,]/max(dat[2,])))
        }
        yrange = unlist(lapply(bins, function(i) i[-(1:2),]))
        yrange = yrange[!is.na(yrange) & !is.infinite(yrange)]
        yrange = range(yrange)
        plot(range(test_temps), yrange, 
              axes = FALSE, xlab = '', ylab = '', type = 'n')   
         
        axis(2)
        axis(1)
        grid()
        xpol = range(sapply(fireSig, function(i) i[2]))
        polygon(xpol[c(1, 2, 2, 1, 1)], yrange[c(1,1,2,2,1)], col = '#00000011', border = NA)
        lapply(xpol, function(x) lines(c(x, x), yrange, lty = 2))
        for (i in 1:3) {
            mapply(addPoly, bins, c('blue', 'red'))
            mapply(addPoly, rev(bins), rev(c('blue', 'red')))
        }
        lines(c(-9E9, 9E9), c(0, 0), lty = 2, lwd = 2)
       
        mtext(side = 3, adj = 0.1, mt, line = -1.7)
        
        return()
    }
    
    plot(xrange, c(0, 14), axes = FALSE, 
        ylim = c(13, 0), xlab = '', ylab = '', type = 'n')

    #axis(1)
    
    
    findColsLims <- function(dats, symm = FALSE) {
        
       
        if (symm) {
            limits = sort(c(0, quantile(dats, seq(0.1, 0.9, 0.2))))
            #test = limits < 0
            limits = unique(limits)
            #browser()
            limits = sort(c(-limits, limits))
            #if (abs(limits[1]) > head(limits, 1)) 
            #    limits = c(limits[test], -rev(limits[test]))
            #else
            #    limits = c(-rev(limits[!test]), limits[!test])
        } else limits = sort(c(0, quantile(dats, seq(0.1, 0.9, 0.1))))
        
        limits = unique(signif(limits, 1))
        limits[limits == -0.09] = -0.1
        limits = limits[limits != -0.0003]
        if (symm) limits = limits[limits != 0]
        colsA = colsB = c()
        if (any(limits < 0))
            colsA = make_col_vector(cols1, ncols = 1+sum(limits <= 0))
        if (any(limits > 0))
            colsB = make_col_vector(cols2, ncols = 1+sum(limits >= 0))
        cols = c(colsA, colsB)
        if ((length(cols) -2) == length(limits)) cols = c(colsA, colsB[-1])
        if ((length(cols) -3) == length(limits)) cols = c(head(colsA, -1), colsB[-1])
        if ((length(cols) -1) != length(limits)) browser()
        
        return(list(limits, cols))
    }
    c(limits, cols) := findColsLims(unlist(nbps))
    c(dlimits, dcols) := findColsLims(unlist(dnps), TRUE)
    
    forModel <- function(model, mi, x, mark1, mark2) {
        idM = which(dat[,2] == model)
        id = idR[idR > idM][1]
        nbp = nbps[[which(models == model)]]#dat[id+(2:1), ]
        #if (accumulative) nbp = cumsum(as.numeric(nbp))
        forExp <- function(i, cols, limits) {
            y = mi + (i-1) * (length(models)+0.5) 
            y = rep(y, length(x))
            if (i == 3) col = as.numeric(nbp[2,-1]) - as.numeric(nbp[1,-1])
            else col = as.numeric(nbp[i,-1])
            if (length(col) == (length(x)-1)) col =c(0, col)
            col = cut_results(col, limits) 
            
            plotLines <- function(x, y, col, lwd = 3, ...) 
                lines(c(x, x), y + c(-0.5, 0.5), col = col, lend = 1, lwd = lwd, ...)
            cols.pt = cols[col] 
            
            mapply(plotLines, x, y, cols.pt)
            if (region == regions[1] && i == 1) text(x = x[7], y = y[7], model, cex = 1.2, font = 2)
            
            if (i == 3) {
                addMark <- function(mark, lty = 2, col = "black", adj = -0.2, srt = 270) {
                    plotLines(mark[1], y[1], col, lwd = 1.5, lty = lty)
                    text(x = mark[1], y = y[1], col = col,
                        adj = c(0.5, adj), round(mark[2], 2), srt = srt)
                }
                
                addMark(mark1, 2)#; addMark(mark2, 3, "white", adj = -0.2, srt = 90)
            }
        }
        
        forExp(1, cols, limits)        
        forExp(2, cols, limits)      
        forExp(3, dcols, dlimits)
        
    }
    
    mapply(forModel, models, 1:length(models), xs, fireSig, fireSig_an)
    
    xx = xrange[2] + diff(xrange) * c(0.017,0.05)
    legendColBar(xx, c(0.5, 8.5), 10, cols, limits, TRUE)
    legendColBar(xx, c(9, 14), 10, dcols, dlimits, TRUE)
    mtext(side = 3, adj = 0.1, mt, line = -1.7)
    if (Rid %% 4 == 1) {
        textSt <- function(y, txt) 
            text(x = xrange[1] - diff(xrange)*0.04, y = y , txt, srt = 90, cex = 1.3, xpd = NA)
        
        textSt(02.5 , 'Without Fire')
        textSt(07.0 , 'With Fire')
        textSt(11.5, 'Difference')
    }
    axis(1)

    out = matrix(unlist(fireSig), nrow = 2)
    colnames(out) = models
    rownames(out) = paste0(region, '-', c('Year', 'Temp'))
    outi = apply(out, 1, range)
    if (any(is.na(outi))) outi = apply(out, 1, function(i) c(min(i, na.rm = TRUE), NA))
    outi = as.vector(outi)
    names(outi) = c('year-min', 'year-max', 'tas-min', 'tas-max')
    outi[3:4] = round(outi[3:4], 2)
    
    if (region == "Global") {
        sigTemp = rbind(out, 
                        mapply(function(yr, nbp) diff(nbp[,which(years == yr)]), out[1,], nbps))

        emssionsAtTemp <- function(tas) {
            forModel <- function(temp, nbp) {
                i = which.min(abs(temp - tas))
                diff(nbp[,i])
            }
            out = mapply(forModel, temps, nbps)
            return(round(out))
        }
        tass = seq(1.1, 2.0, by = 0.1)
        ipccTemp = sapply(tass, emssionsAtTemp)
        colnames(ipccTemp) = tass
        write.csv(ipccTemp, file = 'outputs/allowableEmissins_change.csv')
    }
    return(outi)

    #if (Rid > (length(regions) -4)) axis(1)
}

plotTS <- function(xadj = 0.5) {
    mapply(forRegions, regions, 1:length(regions), TRUE, letter = letters[1:length(regions)])
    mtext.units(outer = TRUE, side = 1, '~DELTA~ T (~DEG~C)', line = 1.25, adj = xadj)
    mtext.units(outer = TRUE, side = 2, 'NBP (PgC ~yr-1~)', line = 0.75)
    plot.new()
    col = c("blue", "red")
    par(mar = c(1, 4, 7, 1))
    for (i in 1:6) {
        legend("topleft", c("Without fire", "With fire"), lwd = 1.5, bty = 'n',
              col = make.transparent(col, 0.6))
        legend("topleft", c("Without fire", "With fire"), lwd = 10, bty = 'n',
               col = make.transparent(col, 0.95))
    }
}
figname <- function(nme) paste0("figs/", nme, "-accu_", accumulative, ".png")



png(figname("NBPstripes-combined"), res = 300, units = 'in', height = 2*6.69291*0.67, width = 2*7.08661)
    par(mar = c(1.5, 1.5, 0.5, 0.25), oma = c(2.2, 2.7, 1, 1.5))
    layout(cbind(t(matrix(1:16, ncol = 4)), rbind(17:18, 17:18, 19:20, 19:20)),
           widths = c(1,1,1,1, 1.5, 1.5))
    plotTS(xadj = 2.5/7)
    subRegions = c(1,3, 5, 6)
    par(mar = c(1.5, 1.5, 0.5, 0.25))
    accumulative = TRUE
    mapply(forRegions, regions[subRegions], subRegions, TRUE, TRUE, letter = letters[16:19])
dev.off()

accumulative = FALSE
png(figname("NBP_vs_temp"), res = 300, units = 'in', height = 12, width = 12)
    par(mfrow = c(4, 4), mar = c(1.5, 1.5, 0, 1.5), oma = c(2.75, 3.0, 0.25, 0.0))
    plotTS()
dev.off()

accumulative = TRUE
png(figname("NBPstripes"), res = 300, units = 'in', height = 14, width = 12)
    par(mfrow = c(4, 4), mar = c(1.5, 1.5, 0, 1.5), oma = c(1.5, 1, 1, 1.5))

    outs = mapply(forRegions, regions, 1:length(regions), TRUE, TRUE, 
                  letter = letters[1:length(regions)])
    write.csv(t(outs), file = 'temp_yrs.csv')

dev.off()

