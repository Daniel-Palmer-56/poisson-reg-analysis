# Daniel Palmer
# Final Project
# STAT 774

# Read in data
psychiatric <- read.csv(file = file.choose(), header = TRUE, na.strings = c("888", "999"))

# Subset questions for response variable (QTOTAL) 
questions <- psychiatric[, 14:82]

# Remove pairs of questions with 0 or negative correlations
corr <- cor(questions)
colnames(corr) <- rownames(corr) <- paste0("Q", 1:69)
index <- function(x){x[x <= 0]}
ind <- apply(corr, 1, FUN = index)
ind2 <- apply(corr, 2, FUN = index)
questions <- questions[, -c(unique(match(ind, ind2)))]
colnames(questions) # Display final set of questions

# Create the response variable (QTOTAL)
QTOTAL <- cbind(rowSums(questions))

# Condense occupation levels into 3 levels
psychiatric$Occupation <- ifelse(psychiatric$Occupation == "1" | psychiatric$Occupation == "2" | psychiatric$Occupation == "3", 1, ifelse(psychiatric$Occupation == "4" | psychiatric$Occupation == "7", 2, 3))
# Add QTOTAL to data
psychiatric <- cbind(psychiatric[, -c(14:82)], QTOTAL)
# Remove level 3 from occupation (administrators and office staff)
psychi <- psychiatric[psychiatric$Occupation != "3", ]
# Remove case with Relative_MI missing
psychi <- psychi[!is.na(psychi$Relative_MI), ]
# Convert appropriate demographics to factors
psychi[, -c(1:2, 5:6, 9, 11, 14)] <- apply(psychi[, -c(1:2, 5:6, 9, 11, 14)], 2, as.factor)

# Create index for sampling train and test data
index <- 1:dim(psychi)[1]
set.seed(772) # Set seed to replicate result
# 90/10 train-test split
train <- sample(index, round(0.9 * length(index)), replace = FALSE)
test <- index[-train]
train <- psychi[train, ]
test <- psychi[test, ]

# Fit the Poisson regression model with all two-way interactions for train data
model <- glm(QTOTAL ~ (Site + Gender + Age + Education + Marital_Status + Race + Years_Employed_MH + Occupation + Hrs_Week_Worked + Setting_Worked + Relative_MI) ^ 2, family = poisson(link = "log"), data = train)
# Perform backward stepwise selection
modelFinal <- step(model)
# Test for goodness-of-fit
deviance(modelFinal) < qchisq(0.99, modelFinal$df.residual)
summary(modelFinal) # Summary output

# Make initial predictions for train and test data
predictionsTrain <- predict(modelFinal, type = "response")
response <- train[, "QTOTAL"]
(trainRMSE <- sqrt(mean((response - predictionsTrain) ^ 2)))
predictionsTest <- predict(modelFinal, newdata = test, type = "response")
pred <- cbind(psychi[names(predictionsTest), ], predictionsTest)
tapply(pred[, 15], pred[, 10], mean) # Predicted for test
tapply(pred[, 14], pred[, 10], mean) # Actual
response <- test[, "QTOTAL"]
(testRMSE <- sqrt(mean((response - predictionsTest) ^ 2)))

# Dispersion
mu <- model.matrix(modelFinal) %*% coef(modelFinal)
response <- train[, "QTOTAL"]
v <- function(x) {x} # v(.)
pearsonResid <- (response - exp(mu)) / sqrt(v(exp(mu)))
sum(pearsonResid ^ 2) / (dim(train)[1] - length(coef(modelFinal)))
sum(resid(modelFinal, type = "pearson") ^ 2) / modelFinal$df.residual

# Show residual plot
plot(resid(modelFinal, type = "pearson") ~ predict(modelFinal, type = "link"), xlab = "Link", ylab = "Pearson Residuals")
abline(h = 0)

# IRWLS
threshold <- 1e-8 # Set the error
eta <- function(x) {log(x)} # Link function (log)
etaInverse <- function(x) {exp(x)} # Inverse link function
v <- function(x) {x ^ 1.08} # v(.)
derivativeEta <- function(x) {1 / x} # Derivative of link function
y <- train[, "QTOTAL"] # Response
designMatrix <- model.matrix(modelFinal) # Design matrix
betaHat <- coef(modelFinal) # Initial estimates for coefficients
muHat <- etaInverse(designMatrix %*% betaHat) # Estimates for mu
WHat <- solve(diag(as.vector(v(muHat) *
                               (derivativeEta(muHat)) ^ 2))) # Weights matrix
zHat <- eta(muHat) + derivativeEta(muHat) * (y - muHat) # Working residual
betaHat <- solve(t(designMatrix) %*%
                   WHat %*% designMatrix) %*% t(designMatrix) %*% WHat %*% zHat # Coef estimates
sets <- list(betaHat) # Create list for each pair of coef estimates

# Solve with IRWLS until convergence
iterations <- 2:100 # Number of iterations
for(i in iterations){
  muHat <- etaInverse(designMatrix %*% betaHat) # Estimates for mu
  WHat <- solve(diag(as.vector(v(muHat) *
                                 (derivativeEta(muHat)) ^ 2))) # Weights matrix
  zHat <- eta(muHat) + derivativeEta(muHat) * (y - muHat) # Working residual
  betaHat <- solve(t(designMatrix) %*%
                     WHat %*% designMatrix) %*% t(designMatrix) %*% WHat %*% zHat # Coef estimates
  sets[[i]] <- betaHat # Add coef estimates to the list
  # Exit and return coefficients when differences < threshold
  if(sqrt(sum((sets[[i]] - sets[[i - 1]]) ^ 2)) < threshold){
    return(print(betaHat))
  }
}

# Summary output
covariance <- solve(t(designMatrix) %*% WHat %*% designMatrix)
rownames(covariance) <- colnames(covariance) <- rownames(betaHat)
se <- sqrt(diag(covariance))
z <- betaHat / se
p <- 2 * pnorm(-abs(betaHat) / se)
summarize <- cbind(betaHat, se, z, p)
colnames(summarize) <- c("Coefficients", "SE", "z", "p-value")
head(summarize)

# Dispersion
mu <- designMatrix %*% betaHat
response <- train[, "QTOTAL"]
pearsonResid <- (response - exp(mu)) / sqrt(v(exp(mu)))
sum(pearsonResid ^ 2) / (dim(train)[1] - length(betaHat))

# Show residual plot
plot(c(mu), c(pearsonResid), xlab = "Link", ylab = "Pearson Residuals")
abline(h = 0)

# Make predictions for train and test data
predictionsTrain2 <- exp(designMatrix %*% betaHat)
response <- train[, "QTOTAL"]
(trainRMSE2 <- sqrt(mean((response - predictionsTrain2) ^ 2)))
designMatrixTest <- model.matrix(glm(formula = QTOTAL ~ Site + Gender + Age + Education + Marital_Status + 
                                       Race + Years_Employed_MH + Occupation + Hrs_Week_Worked + 
                                       Setting_Worked + Relative_MI + Site:Gender + Site:Education + 
                                       Site:Years_Employed_MH + Gender:Years_Employed_MH + Gender:Occupation + 
                                       Gender:Setting_Worked + Age:Years_Employed_MH + Age:Occupation + 
                                       Education:Years_Employed_MH + Education:Occupation + Education:Hrs_Week_Worked + 
                                       Education:Relative_MI + Marital_Status:Years_Employed_MH + 
                                       Marital_Status:Occupation + Marital_Status:Hrs_Week_Worked + 
                                       Years_Employed_MH:Hrs_Week_Worked + Occupation:Hrs_Week_Worked + 
                                       Occupation:Setting_Worked, family = poisson(link = "log"), 
                                     data = test))
predictionsTest2 <- exp(designMatrixTest %*% betaHat)
pred2 <- cbind(psychi[rownames(predictionsTest2), ], predictionsTest2)
tapply(pred[, 15], pred[, 10], mean) # Predicted for initial test
tapply(pred2[, 15], pred2[, 10], mean) # Predicted for adjusted test
tapply(pred2[, 14], pred2[, 10], mean) # Actual
response <- test[, "QTOTAL"]
(testRMSE2 <- sqrt(mean((response - predictionsTest2) ^ 2)))

# Identify significant coefficients
names(p) <- rownames(summarize)
names(p[p <= 0.01])

# Gender by Occupation (nominal-by-nominal)
interact <- data.frame(xtabs(~ Gender + Occupation, data = psychi))
interactModel <- glm(Freq ~ Gender + Occupation, data = interact, family = poisson(link = "log"))
deviance(interactModel) < qchisq(0.99, interactModel$df.residual)
# Evidence against independence

# Gender by Setting_Worked (nominal-by-nominal)
interact2 <- data.frame(xtabs(~ Gender + Setting_Worked, data = psychi))
interactModel2 <- glm(Freq ~ Gender + Setting_Worked, data = interact2, family = poisson(link = "log"))
deviance(interactModel2) < qchisq(0.99, interactModel2$df.residual)
# Evidence against independence

# Occupation by Setting_Worked (nominal-by-nominal)
interact10 <- data.frame(xtabs(~ Occupation + Setting_Worked, data = psychi))
interactModel10 <- glm(Freq ~ Occupation + Setting_Worked, data = interact10, family = poisson(link = "log"))
deviance(interactModel10) < qchisq(0.99, interactModel10$df.residual)
# Evidence against independence

# Education by Years_Employed_MH (ordinal-by-ordinal)
interact3 <- data.frame(xtabs(~ Education + Years_Employed_MH, data = psychi))
interactModel3 <- glm(Freq ~ Education + Years_Employed_MH + I(unclass(Education) * unclass(Years_Employed_MH)), data = interact3, family = poisson(link = "log"))
deviance(interactModel3) < qchisq(0.99, interactModel3$df.residual)
# Evidence against independence

# Age by Years_Employed_MH (ordinal-by-ordinal)
interact8 <- data.frame(xtabs(~ Age + Years_Employed_MH, data = psychi))
interactModel8 <- glm(Freq ~ Age + Years_Employed_MH + I(unclass(Age) * unclass(Years_Employed_MH)), data = interact8, family = poisson(link = "log"))
deviance(interactModel8) < qchisq(0.99, interactModel8$df.residual)
# Evidence against independence

# Occupation by Hrs_Week_Worked (nominal-by-ordinal)
interact4 <- data.frame(xtabs(~ Occupation + Hrs_Week_Worked, data = psychi))
interactModel4 <- glm(Freq ~ Occupation + Hrs_Week_Worked + Occupation:unclass(Hrs_Week_Worked), data = interact4, family = poisson(link = "log"))
deviance(interactModel4) < qchisq(0.99, interactModel4$df.residual)
# Evidence against independence

# Gender by Years_Employed_MH (nominal-by-ordinal)
interact5 <- data.frame(xtabs(~ Gender + Years_Employed_MH, data = psychi))
interactModel5 <- glm(Freq ~ Gender + Years_Employed_MH + Gender:unclass(Years_Employed_MH), data = interact5, family = poisson(link = "log"))
deviance(interactModel5) < qchisq(0.99, interactModel5$df.residual)
# Evidence against independence

# Occupation by Age (nominal-by-ordinal)
interact6 <- data.frame(xtabs(~ Occupation + Age, data = psychi))
interactModel6 <- glm(Freq ~ Occupation + Age + Occupation:unclass(Age), data = interact6, family = poisson(link = "log"))
deviance(interactModel6) < qchisq(0.99, interactModel6$df.residual)
# Evidence against independence

# Site by Years_Employed_MH (nominal-by-ordinal)
interact7 <- data.frame(xtabs(~ Site + Years_Employed_MH, data = psychi))
interactModel7 <- glm(Freq ~ Site + Years_Employed_MH + Site:unclass(Years_Employed_MH), data = interact7, family = poisson(link = "log"))
deviance(interactModel7) < qchisq(0.99, interactModel7$df.residual)
# Evidence against independence

# Relative_MI by Education (nominal-by-ordinal)
interact9 <- data.frame(xtabs(~ Relative_MI + Education, data = psychi))
interactModel9 <- glm(Freq ~ Relative_MI + Education + Relative_MI:unclass(Education), data = interact9, family = poisson(link = "log"))
deviance(interactModel9) < qchisq(0.99, interactModel9$df.residual)
# Evidence against independence

# Site by Education (nominal-by-ordinal)
interact11 <- data.frame(xtabs(~ Site + Education, data = psychi))
interactModel11 <- glm(Freq ~ Site + Education + Site:unclass(Education), data = interact11, family = poisson(link = "log"))
deviance(interactModel11) < qchisq(0.99, interactModel11$df.residual)
# Evidence against independence

# Occupation by Education (nominal-by-ordinal)
interact12 <- data.frame(xtabs(~ Occupation + Education, data = psychi))
interactModel12 <- glm(Freq ~ Occupation + Education + Occupation:unclass(Education), data = interact12, family = poisson(link = "log"))
deviance(interactModel12) < qchisq(0.99, interactModel12$df.residual)
# Evidence against independence