legendColBar <- function(xx = c(0.3, 0.7), yy = c(0.1, 0.9), cols, limits, switch = FALSE,
                             extend_max = TRUE, extend_min = TRUE, minLab = '', maxLab = '',
                             units = '', transpose = FALSE, oneSideLabels = TRUE, adj = 0.5,
                             add = FALSE) {
        if (!add) plot(c(0, 1), c(0, 1), axes = FALSE, type = 'n')
        cols = make_col_vector(cols, ncols = length(limits)+1)
        ys = seq(yy[1], yy[2], length.out = length(cols) +1)
        addBox <- function(y1, y2, yi, col) {
            print(yi)
            polyFun <- function(x, y) {
                if (transpose) {xi = x; x = y; y = xi}
                    polygon(x, y, col = col, lwd = 2, xpd = NA)
            }
            F1 <- function() polyFun(c(xx, mean(xx), xx[1]), c(y2, y2, y1, y2))
            F2 <- function() polyFun(c(xx, mean(xx), xx[1]), c(y1, y1, y2, y1))
            if (yi == 1 && extend_min) 
                if (switch) F1() else F2()
            else if (yi == length(cols) && extend_max) 
                if (switch) F2() else F1()
            else polyFun(c(xx, rev(xx), xx[1]), c(y2, y2, y1, y1, y2))
                 
        }
        id = 1:length(cols)
        if (switch) {
            id = rev(id)
            cols = rev(cols)
            limits = rev(limits)
        }
        mapply(addBox, ys[-1], head(ys, -1), id, cols)
        
        y = ys
        
        if (!extend_min && minLab == "" && limits[1] == 0) y[2] = mean(y[1:2])
        limits = c(minLab, limits, maxLab)
        limits[limits != ''] = paste0(limits[limits != ''], units)
        x = rep(xx[2] + diff(xx) *0.3, length(y))
        if (is.na(oneSideLabels)) xx = xx[1] - diff(xx) *0.3
        else if (!oneSideLabels) x[seq(1, length(x), by = 2)] = xx[1] - diff(xx) *0.3
        
        if (transpose) { xi = x; x = y; y = xi}
        text(x, y, limits, xpd = NA, adj = adj)
    }
