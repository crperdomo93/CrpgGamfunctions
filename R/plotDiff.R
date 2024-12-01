#' Compute Differences in Predicted Values Between Groups and Generate a Plot
#'
#' This function calculates the differences in predicted values between two groups 
#' (e.g., Impacted vs. Non-Impacted) using a fitted Generalized Additive Model (GAM). 
#' It generates confidence intervals and creates a plot of these differences over time. 
#' This function is specifically created to work with the variables of the data frame "data".
#' When using a database with different names, it is necessary to make these changes in the code.
#'
#' @param modelo A fitted GAM object from the `mgcv` package.
#' @param variable A string specifying the variable of interest (default: "Oxygen").
#' @param name A string specifying the title of the plot (default: 
#'        "Difference in dissolved oxygen between impacted and non-impacted").
#' @param exclude A string or vector of strings specifying the terms to exclude 
#'        when calculating predictions (default: "s(Year,Stnumber)").
#' @param ylab A string specifying the label for the y-axis of the plot (default: "DO").
#' @param xlab A string specifying the label for the x-axis of the plot (default: "Year").
#'
#' @return A list containing:
#' \item{data}{A `data.frame` with the following columns:
#'   \itemize{
#'     \item `Year`: The year.
#'     \item `yhat`: The predicted differences between groups.
#'     \item `lower`: The lower bound of the confidence interval.
#'     \item `upper`: The upper bound of the confidence interval.
#'     \item `Diferencia`: A label indicating the comparison (e.g., "Impactado-NO Impactado").
#'   }}
#' \item{grafico}{A `ggplot2` object representing the plot of differences and their confidence intervals.}
#'
#' @details
#' The function generates a new dataset by expanding the variable `Year` over its observed range,
#' along with group labels and a fixed level for random effects (e.g., `Stnumber`). It computes 
#' the differences in predictions between two groups (e.g., Impacted and Non-Impacted) and their 
#' confidence intervals. The output includes a plot showing these differences over time.
#'
#' @examples
#' \dontrun{
#' # Assuming `m` is a fitted GAM object:
#' result <- fdiff2(modelo = m, variable = "Oxygen")
#' print(result$grafico)
#' head(result$data)
#' }
#'
#' @import ggplot2
#' @import mgcv
#' @export
plotDiff = function(modelo = m,
                  variable = "Oxygen", 
                  name = "Difference in dissolved oxygen between impacted and non-impacted",
                  exclude = "s(Year,Stnumber)", ylab = "DO", xlab="Año") {
  
  # Check if more than one exclusion term is provided
  if (length(exclude) > 1) {
    # Extract the dataset used in the GAM model
    basedatos = modelo$model
    
    # Create a new dataset for generating the predictor matrix
    pdat <- expand.grid(
      Year = seq(min(basedatos$Year), max(basedatos$Year), by = 1),
      Impacted = levels(basedatos$Impacted),
      #This is the default, it does not make any changes if it were another station. 
      #It is only so that the column is generated
      Stnumber = c("50010500") 
    )
    
    # Generate the linear predictor matrix excluding specific terms
    xp = predict(modelo, newdata = pdat, type = "lpmatrix",
                 exclude = list(paste(exclude)))
    
    # Identify columns for the two groups (Impacted and Not Impacted)
    c1 = grepl('I', colnames(xp))
    c2 = grepl('N', colnames(xp))
    
    # Identify rows for the two groups in the new dataset
    r1 <- with(pdat, Impacted == 'I')
    r2 <- with(pdat, Impacted == 'N')
    
    # Compute the difference matrix between the two groups
    xd = (xp[r1, ] - xp[r2, ])
    # Calculate the predicted differences
    yhdif = xd %*% coef(modelo)
    
    # Compute confidence intervals
    se <- sqrt(rowSums((xd %*% vcov(modelo)) * xd))
    crit <- qt(.975, df.residual(modelo))
    upr <- yhdif + (crit * se)
    lwr <- yhdif - (crit * se)
    
    # Prepare the data for output and plotting
    difdata = data.frame(
      Year = rep(seq(min(basedatos$Year), max(basedatos$Year))),
      yhat = yhdif,
      lower = lwr,
      upper = upr,
      Diferencia = paste("Impactado", "NO Impactado", sep = "-")
    )
    
    # Create the plot
    grafico = ggplot(difdata, aes(x = Year, y = yhat)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
      geom_line() +
      ggtitle(name) +
      ylab(paste(ylab)) +
      xlab(paste(xlab)) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)
      ) +
      geom_hline(yintercept = 0, linetype = 3, color = "red", lwd = 1)
    
    # Return the data and plot as a list
    l = list(data = difdata, grafico = grafico)
    
  } else {
    # Case where only one exclusion term is provided
    
    # Extract the dataset used in the GAM model
    basedatos = modelo$model
    
    # Create a new dataset for generating the predictor matrix
    pdat <- expand.grid(
      Year = seq(min(basedatos$Year), max(basedatos$Year), by = 1),
      Impacted = levels(basedatos$Impacted),
      Stnumber = c("50010500")
    )
    
    # Generate the linear predictor matrix excluding a single term
    xp = predict(modelo, newdata = pdat, type = "lpmatrix",
                 exclude = list(exclude[1], exclude[2]))
    
    # Identify columns for the two groups (Impacted and Not Impacted)
    c1 = grepl('I', colnames(xp))
    c2 = grepl('N', colnames(xp))
    
    # Identify rows for the two groups in the new dataset
    r1 <- with(pdat, Impacted == 'I')
    r2 <- with(pdat, Impacted == 'N')
    
    # Compute the difference matrix between the two groups
    xd = (xp[r1, ] - xp[r2, ])
    # Calculate the predicted differences
    yhdif = xd %*% coef(modelo)
    
    # Compute confidence intervals
    se <- sqrt(rowSums((xd %*% vcov(modelo)) * xd))
    crit <- qt(.975, df.residual(modelo))
    upr <- yhdif + (crit * se)
    lwr <- yhdif - (crit * se)
    
    # Prepare the data for output and plotting
    difdata = data.frame(
      Year = rep(seq(min(basedatos$Year), max(basedatos$Year))),
      yhat = yhdif,
      lower = lwr,
      upper = upr,
      Diferencia = paste("Impactado", "NO Impactado", sep = "-")
    )
    
    # Create the plot
    grafico = ggplot(difdata, aes(x = Year, y = yhat)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
      geom_line() +
      ggtitle(name) +
      ylab(paste(ylab)) +
      xlab(paste(xlab)) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)
      ) +
      geom_hline(yintercept = 0, linetype = 3, color = "red", lwd = 1)
    
    # Return the data and plot as a list
    l = list(data = difdata, grafico = grafico)
  }
  
  # Return the final list containing data and plot
  l
}
