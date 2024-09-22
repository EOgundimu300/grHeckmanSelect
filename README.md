##################################################################################
# Group variable selection in nonignorable missing data using LASSO, SCAD and MCP.
##################################################################################

**Description**

Group variable selection in nonignorable missing data using LASSO, SCAD and MCP.

**Usage**

grHeckSelect(
W,
X,
s,
y,
group_sel,
group_out,
penalty = "grLasso",
nlambda = 100,
lambda,
lambda.min = 0.001,
log.lambda = TRUE,
eps = 1e-04,
max_iter = 10000,
gamma = ifelse(penalty == "grSCAD", 4, 3),
group_multiplier,
init_strat = "MLE",
penalty.factors = NULL
)


**Arguments**

W: A matrix of covariates in selection equation (intercept is not included).

X: A matrix of covariates in outcome equation (intercept is not included).

s: Binary outcome for the selection equation (0/1 and false/true allowed).

y: Continuous outcome for the outcome equation (0 for NA).

group_sel: A vector describing the grouping of the coefficients in the selection equation. It is best if group is a vector of consecutive integers. If there are coefficients to be included in the
model without being penalized, assign them to group 0.

group_out: A vector describing the grouping of the coefficients in the outcome equation. It is best if group is a vector of consecutive integers. If there are coefficients to be included in the model without being penalized, assign them to group 0.

penalty: The penalty to be applied to the model, one of grLasso, grSCAD, or grMCP.

nlambda: The number of lambda values. Default is 100.

lambda: A user supplied sequence of lambda values. Typically, this is left unspecified, and the function automatically computes a grid of lambda values.

lambda.min: The smallest value for lambda, as a fraction of lambda.max. Default is .001.

log.lambda: When TRUE compute the grid values of lambda on log scale (default) or linear scale otherwise.

eps: Convergence threshhold. Default is 1e-4

max_iter: Maximum number of iterations (total across entire path). Default is 10000.

gamma: Tuning parameter of the group MCP/SCAD penalty. Default is 3 for MCP and 4 for SCAD

group_multiplier: A vector of values representing multiplicative factors by which each group's penalty is to be multiplied. The default is the square root of group size.

init_strat: Use MLE to initialize the group descent algorithm otherwise use zeros.

penalty.factors: Allows for the use of weighted versions of grLasso, grSCAD, or grMCP.

Value: class grHeckSelect containing optimal penalized coefficients, lambda values, bic values etc.


**Examples**

*Run examples*

#grHeckSelect(W=W, X=X, s=s, y=y, group_sel, group_out, penalty="grLasso")

#######################################################################################
# Group variable selection in nonignorable missing data using Broken Adaptive Ridge penalty.
###########################################################################################

**Description**

Group variable selection in nonignorable missing data using Broken Adaptive Ridge penalty.
Usage
Arguments
W
A matrix of covariates in selection equation (intercept is not included).
X
A matrix of covariates in outcome equation (intercept is not included).
s
Binary outcome for the selection equation (0/1 and false/true allowed).
y
Continuous outcome for the outcome equation (0 for NA).

group_sel
A vector describing the grouping of the coefficients in the selection equation. It is best if
group is a vector of consecutive integers. If there are coefficients to be included in the
model without being penalized, assign them to group 0.
group_out
A vector describing the grouping of the coefficients in the outcome equation. It is best if
group is a vector of consecutive integers. If there are coefficients to be included in the
model without being penalized, assign them to group 0.
ridge_lambda
A user supplied sequence of lambda values for the initial ridge estimate. Typically, this is
left unspecified, and the function automatically computes a grid of lambda values.
lambda
A user supplied sequence of lambda values for the BAR iteration. Typically, this is left
unspecified, and the function automatically computes a grid of lambda values.
nridge_lambda
The number of lambda values for the ridge. Default is 10.
nlambda
The number of lambda values for the BAR iteration. Default is 50 which leads to 500
combinations in the grid to test.
ridge_lambda.min
The smallest value of ridge_lambda. The default is 0.001.
lambda.min
The smallest value for lambda in BAR iteration. The default is 0.001.
log.ridge_lambda
When TRUE compute the grid values of ridge_lambda on log scale (default) or linear
scale otherwise.
log.lambda
When TRUE compute the grid values of lambda on log scale (default) or linear scale
otherwise.
ridge_eps
The value fo tolerance used in intial ridge estimation.
eps
The tolerance between BAR iteration.
inner_eps
The value of eps for any optimization within each BAR iteration.
max_iter
Maximum number of iterations (total across entire path). Default is 10000.

init_strat
Use MLE to initialize the group descent algorithm otherwise use zeros.
del
threshold_sel The value for which to set parameter values with magnitude less than it to 0.
Default is 0.05.
threshold_out
Same as threshold_sel.
method
Three methods are implemented. Method 2 is the fastest (see details in accompanying
paper).
Value
class grHeckSelect_bar containing optimal penalized coefficients, lambda values, bic values etc.

#Examples
Run examples
#grHeckSelect_bar(W,X,s,y,group_sel, group_out, method=2)

Simulated data with group structure for Heckman model
Description
grHeckman was generated from bivariate normal error with correlation 0.5, and sigma^2 =2.
Usage
Format
A data frame with 1000 rows and 23 variables
yobs
observed outcome variable
ustar
selection indicator
ymiss
underlying outcome variable. It is not needed for the modelling
X1-X20
predictor variables
Details
There are 1000 observations and 20 predictors.
outcome: beta <- c(0.5, 1, 1, 1.5, 1,0.2, 0.2, 0.2, 0.2, 0.2, 0.5, 1, 1.5, 0, 0, 0, 0, 0, 0, 0, 0)
selection: gamma <- c(2.7, 1, 1, 1.5, 1,0.2, 0.2, 0.2, 0.2, 0.2, 0.5, 1, 1.5, 0, 0, 0, 0, 0, 0, 0, 0)
There is no exclusion restriction in the selection equation.
Induce group sctructure:
group_sel <- c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5)
group_out <- c(6,6,6,6,7,7,7,7,8,8,8,8,9,9,9,9,10,10,10,10)

data(grHeckman)
## Not run:
select <- ustar~X1+X2+X3+X4+X5+X6+X7+X8+ X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20
outcome <- yobs~X1+X2+X3+X4+X5+X6+X7+X8+ X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20
data(grHeckman); dd <- grHeckman
mf <- model.frame(select, data = dd)
s <- model.response(mf, "numeric")
W <- model.matrix(select, data = dd)[,-1]
mf2 <- model.frame(outcome, dd)
y <- model.response(mf2, "numeric")
X <- model.matrix(outcome, data = dd)[,-1]
group_sel <- c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5)
group_out <- c(6,6,6,6,7,7,7,7,8,8,8,8,9,9,9,9,10,10,10,10)
#grlasso <- grHeckSelect(W=W, X=X, s=s, y=y, group_sel, group_out, penalty="grLasso")
#coef(grlasso)
#grbar <- grHeckSelect_bar(W,X,s,y,group_sel, group_out, method=2)
## End(Not run)
#coef(grbar)
