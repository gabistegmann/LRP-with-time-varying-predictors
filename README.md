# LRP-with-time-varying-predictors
This is an extension of longRPart2 in order to include time-varying predictions

This method has two stages: in the first one, it fits a longitudinal recursive partitioning 
(as it is done in R's already-existing longRPart2) based on time-invariant predictors, 
and at the second stage it fits a tree to the residuals based on time-varying predictors.

This work is still in progress.
