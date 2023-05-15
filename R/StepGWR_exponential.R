#' StepGWR: a hybrid spatial model that combines the variable selection capabilities of stepwise regression with the predictive power of Geographically Weighted Regression (GWR) model
#'
#' @param data_sp A dataframe containing a response variable and the predictor variables, as well as the coordinates of the locations. In the dataframe, first column is the response variable (y), last two columns are coordinates i.e., Latitude and Longitudes and in between them is the set of predictor variables(X's).
#' @param bw A numeric value specifying the bandwidth parameter for the GWR model. It can be noted that, optimum bandwidth parameter value can vary and depends on the spatial pattern of the dataset.
#' @param split_value  Splitting value for dividing the dataset into training and testing set, e.g. 0.8 or 0.7
#' @param exponential_kernel Spatial weight function of the GWR model, e.g. exponential_kernel
#' @return A list with the following components:
#'   - `Important_vars `: Selected important variables based on stepwise regression
#'   - `GWR_y_pred_test`: The predicted values based on GWR at test locations
#'   - `rrmse`: root means square error
#'   - `R_squared`: R square value
#'   - `mse`: mean squared error
#'   - `mae`: mean absolute error
#' @examples
#' n<- 100
#' p<- 7
#' m<-sqrt(n)
#' id<-seq(1:n)
#' x<-matrix(runif(n*p), ncol=p)
#' e<-rnorm(n, mean=0, sd=1)
#' xy_grid<-expand.grid(c(1:m),c(1:m))
#' Latitude<-xy_grid[,1]
#' Longitude<-xy_grid[,2]
#' B0<-(Latitude+Longitude)/6
#' B1<-(Latitude/3)
#' B2<-(Longitude/3)
#' B3<-(2*Longitude)
#' B4<-2*(Latitude+Longitude)/6
#' B5<-(4*Longitude/3)
#' B6<-2*(Latitude+Longitude)/18
#' B7<-(4*Longitude/18)
#' y<-B0+(B1*x[,1])+(B2*x[,2])+(B3*x[,3])+(B4*x[,4])+(B5*x[,5])+(B6*x[,6])+(B7*x[,7])+e
#' data_sp<-data.frame(y,x,Latitude,Longitude)
#' StepGWR_exp<-StepGWR_exponential(data_sp,0.5,0.8,exponential_kernel)
#' @references
#' 1. Leung, Y., Mei, C. L. and Zhang, W. X. (2000). Statistical tests for spatial non-stationarity based on the geographically weighted regression model. Environment and Planning A, 32(1),9-32.<DOI:10.1068/a3162>.
#' 2. Brunsdon, C., Fotheringham, A.S. and Charlton, M,E. (1996).Geographically weighted regression: a method for exploring spatial non-stationarity. Geogr Anal.28(4),281-298.<DOI:10.1111/j.1538-4632.1996.tb00936.x>.
#' @export
#' @import qpdf
#' @import numbers
#' @import MASS
#' @importFrom stats cor dist coef

StepGWR_exponential<- function(data_sp,bw,split_value,exponential_kernel) {

  # Step1: Generation of simulated spatial population with spatial coordinate points
  data_sp<-as.data.frame(data_sp)
  # Step2: Splitting of data into training and testing sets
  train_idx <- sample(nrow(data_sp), split_value * nrow(data_sp))
  train_data <- data_sp[train_idx,]
  nrow(train_data)
  test_data <- data_sp[-train_idx, ]
  nrow(test_data)
  x_train <-subset(train_data[,-c(1,ncol(train_data)-1,ncol(train_data))])
  y_train <-train_data[,1]
  x_test <- subset(test_data,select = -c(1,ncol(test_data)-1,ncol(test_data)))
  y_test <- test_data[,1]

  # Step3: Important variable selection using Stepwise Variable Selection Approach

  train_data_Step<-subset(train_data,select = -c(ncol(train_data)-1,ncol(train_data)))

  # Perform stepwise variable selection
  step_model <- stepAIC(lm(train_data_Step[,1] ~ ., data = train_data_Step[,-1]), direction = "both")
  # Extract the selected variables
  selected_variables_step <- names(coef(step_model)[-1])
  # Extracting important variables
  imp_var_Step <- x_train[, selected_variables_step]

  # Step5 : Fitting of the GWR model on training data

  coords_train<-cbind(train_data[,ncol(train_data)-1],train_data[,ncol(train_data)])
  dists<- as.matrix(dist(coords_train))

  # Define the Exponential kernel function
  exponential_kernel <- function(dists, bw) {
    exp(-(dists/bw))
  }
  weights <- exponential_kernel(dists, bw)
  dim(weights)
  Xt <- as.matrix(cbind(rep(1, nrow(x_train)),imp_var_Step))
  dim(Xt)
  y_tr<- matrix(y_train)
  dim(y_tr)
  W.train<-list()
  for (i in 1:ncol(weights)){
    t <- diag(weights[,i],nrow=nrow(x_train), ncol=nrow(x_train))
    W.train[[i]]<-t
  }
  W.train

  Beta.train<-list()
  for (i in 1:length(W.train)){
    lm<- solve((t(Xt)%*%W.train[[i]]%*%Xt))%*%t(Xt)%*%W.train[[i]]%*%y_tr
    Beta.train[[i]]<-lm
  }
  Beta.train

  X.tr_row<-list()
  for(i in 1:nrow(Xt)){

    X.tr_row[[i]]<-(Xt[i,])
  }


  y_hat.train<-mapply("%*%", X.tr_row,Beta.train)
  y_hat.train


  # Step6 : Make predictions at the test locations

  Step_imp_X_test<-x_test[, selected_variables_step]

  coords_test<-cbind(test_data[,ncol(test_data)-1],test_data[,ncol(train_data)])
  dists_test <- as.matrix(dist(coords_test))
  weights_test <- exponential_kernel(dists_test, bw)
  dim(weights_test)

  Xtest <- as.matrix(cbind(rep(1, nrow(x_test)),Step_imp_X_test))
  dim(Xtest)
  ytest<- matrix(y_test)
  dim(ytest)

  W.test<-list()
  for (i in 1:ncol(weights_test)){
    test <- diag(weights_test[,i],nrow=nrow(x_test), ncol=nrow(x_test))
    W.test[[i]]<-test
  }
  W.test

  Beta.test<-list()
  for (i in 1:length(W.test)){
    lm_test<- solve((t(Xtest)%*%W.test[[i]]%*%Xtest))%*%t(Xtest)%*%W.test[[i]]%*%ytest
    Beta.test[[i]]<-lm_test
  }
  Beta.test

  X.test_row<-list()
  for(i in 1:nrow(Xtest)){

    X.test_row[[i]]<-(Xtest[i,])
  }


  y_hat.test<-mapply("%*%", X.test_row,Beta.test)
  y_hat.test

  # Step7 : Create summary output
  # Compute accuracy measures

  rrmse_test <- sqrt(mean((y_hat.test - y_test)^2)) / mean(y_test)
  mae_test <- mean(abs(y_hat.test - y_test))
  mse_test <- mean((y_hat.test - y_test)^2)
  r2_test <- cor(y_hat.test, y_test)^2

  # Return results
  results_StepGWR <- list(Important_vars = selected_variables_step,
                          GWR_y_pred_test = y_hat.test,
                          rrmse=rrmse_test,
                          R_squared = r2_test,
                          mse = mse_test,
                          mae = mae_test)

  return(results_StepGWR)
}
