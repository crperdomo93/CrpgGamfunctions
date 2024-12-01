#' Calculate Differences Relative to a Reference Year and Generate a Plot
#'
#' This function calculates differences in predictions for each year relative 
#' to a reference year (e.g., 1990) using a fitted Generalized Additive Model (GAM). 
#' It computes confidence intervals and generates a plot of these differences by group.
#'
#' @param modelo A fitted GAM object from the `mgcv` package.
#' @param x0 A numeric value specifying the reference year (default: 1990).
#' @param name A string specifying the title of the plot (default: 
#'        "Difference from 1990").
#' @param yetiqueta A string specifying the label for the y-axis of the plot (default: "dif").
#' @param xlab A string specifying the label for the x-axis of the plot (default: "Year").
#'
#' @return A list containing:
#' \item{grafico}{A `ggplot2` object representing the plot of differences and their confidence intervals.}
#' \item{datos}{A `data.frame` containing:
#'   \itemize{
#'     \item `yhatdiffint`: The predicted differences relative to the reference year.
#'     \item `se`: The standard error of the differences.
#'     \item `upr`: The upper bound of the confidence interval.
#'     \item `lwr`: The lower bound of the confidence interval.
#'     \item `Impacted`: The group (e.g., Impacted or Non-Impacted).
#'     \item `Year`: The year.
#'   }}
#'
#' @details
#' The function computes differences in predicted values for each year relative to 
#' a specified reference year (`x0`). It constructs confidence intervals for these 
#' differences and generates a grouped plot using `ggplot2`.
#'
#' @examples
#' \dontrun{
#' # Assuming `m` is a fitted GAM object:
#' result <- fintra2(modelo = m5, x0 = 1990)
#' print(result$grafico)
#' head(result$datos)
#' }
#'
#' @import ggplot2
#' @import mgcv
#' @export




plotIntra = function(modelo = m, x0 = 1990, name = "Difference from 1990", 
                   yetiqueta = "diff", xlab = "Year") {
  
  # Extracting data from the model
  datos = modelo$model
  impsub = levels(datos$Impacted) # Identifying the levels of the variable "Impacted"
  
  # Creating lists to store intermediate results
  r0 = list()       # To store rows for reference year
  lxp0 = list()     # To store design matrices for the reference year
  r = list()        # To store rows for other years
  lxp = list()      # To store design matrices for other years
  xd = list()       # To store differences in design matrices
  yhatdiffint = list() # To store predictions and intervals
  
  # Creating data for the reference year (x0)
  pdat0 <- expand.grid(Year = rep(x0, length(seq(min(datos$Year), max(datos$Year)))), 
                       Impacted = levels(datos$Impacted),
                       Stnumber = c("50010500"))
  
  # Generating the design matrix for the reference year
  xp0 = predict(modelo, newdata = pdat0, type = "lpmatrix", 
                exclude = list('s(Stnumber)', 's(Year,Stnumber)'))
  xp0[, grepl('Stnumber', colnames(xp0))] = 0 # Zeroing out random effects for Stnumber
  
  # Extracting rows for each level of "Impacted" in the reference year
  for (i in 1:length(impsub)) {
    r0[[paste(impsub[i])]] = with(pdat0, Impacted == impsub[i])
    lxp0[[paste("xp0", impsub[i], sep = "-")]] = xp0[r0[[paste(impsub[i])]],]
  }
  
  # Creating data for all years
  pdat <- expand.grid(Year = seq(min(datos$Year), max(datos$Year)), 
                      Impacted = levels(datos$Impacted),
                      Stnumber = c("50010500"))
  
  # Generating the design matrix for all years
  xp = predict(modelo, newdata = pdat, type = "lpmatrix", 
               exclude = list('s(Stnumber)', 's(Year,Stnumber)'))
  xp[, grepl('Stnumber', colnames(xp))] = 0
  
  # Extracting rows for each level of "Impacted" for all years
  for (i in 1:length(impsub)) {
    r[[paste(impsub[i])]] = with(pdat, Impacted == impsub[i])
    lxp[[paste("xp", impsub[i], sep = "-")]] = xp[r[[paste(impsub[i])]],]
  }
  
  # Calculating differences between the current year and the reference year (x0)
  for (i in 1:length(impsub)) {
    xd[[paste("xd", impsub[i], sep = "-")]] = lxp[[i]] - lxp0[[i]]
  }
  
  # Calculating predictions, confidence intervals, and preparing data for plotting
  for (i in 1:length(impsub)) {
    yhatdiffint[[i]] = xd[[i]] %*% coef(modelo)
    yhatdiffint[[i]] = cbind(yhatdiffint[[i]], 
                             sqrt(rowSums((xd[[i]] %*% vcov(modelo)) * xd[[i]])))
    crit = qt(.975, df.residual(modelo)) # Critical value for 95% confidence intervals
    yhatdiffint[[i]] = cbind(yhatdiffint[[i]], 
                             upr = yhatdiffint[[i]][, 1] + (crit * yhatdiffint[[i]][, 2]),
                             lwr = yhatdiffint[[i]][, 1] - (crit * yhatdiffint[[i]][, 2]))
    yhatdiffint[[i]] = as.data.frame(yhatdiffint[[i]])
    yhatdiffint[[i]] = cbind(yhatdiffint[[i]], impsub[i]) # Adding group info
    yhatdiffint[[i]] = cbind(yhatdiffint[[i]], seq(min(datos$Year), max(datos$Year))) # Adding years
    colnames(yhatdiffint[[i]]) = c("yhatdiffint", "se", "upr", "lwr", "Impacted", "Year")
  }
  
  # Combining data for all groups into a single data frame
  difintra = do.call(rbind.data.frame, yhatdiffint)
  
  # Creating the plot
  grafico = ggplot(difintra, aes(x = Year, y = yhatdiffint, group = Impacted)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
    geom_line() +
    facet_wrap(~Impacted) +
    ggtitle(name) +
    ylab(yetiqueta) +
    xlab(xlab) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          plot.title = element_text(hjust = 0.5)) +
    geom_hline(yintercept = 0, linetype = 3, color = "red", lwd = 1)
  
  # Returning the plot and data
  flist = list(grafico = grafico, datos = difintra)
}
