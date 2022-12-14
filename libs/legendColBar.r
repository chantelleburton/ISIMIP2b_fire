legendColBar <- function(xx, yy, dxl, cols, limits, switch = FALSE,
                             extend_max = TRUE, extend_min = TRUE) {

        ys = seq(yy[1], yy[2], length.out = length(cols) +1)
        addBox <- function(y1, y2, yi, col) {
            print(yi)
            polyFun <- function(x, y) polygon(x, y, col = col, lwd = 2, xpd = TRUE)
            if (yi == 1 && extend_min) 
                polyFun(c(xx, mean(xx), xx[1]), c(y2, y2, y1, y2))
            else if (yi == length(cols) && extend_max) 
                polyFun(c(xx, mean(xx), xx[1]), c(y1, y1, y2, y1))
            else polyFun(c(xx, rev(xx), xx[1]), c(y2, y2, y1, y1, y2))
                 
        }
        id = 1:length(cols)
        if (switch) {
            id = rev(id)
            cols = rev(cols)
            limits = rev(limits)
        }
        mapply(addBox, ys[-1], head(ys, -1), id, cols)
        text(xx[2] + diff(xx) *0.3, head(ys[-1], -1), limits, xpd = NA, adj = 0)
    }
