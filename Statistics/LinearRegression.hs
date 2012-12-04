{-# LANGUAGE BangPatterns #-}

module Statistics.LinearRegression (
    linearRegression,
    linearRegressionRSqr,
    linearRegressionTLS,
    linearRegressionError,
    linearRegressionTLSError,
    EstimationParameters(..),
    estimationQuality,
    converge,
    correl,
    covar,
    ) where

import qualified Data.Vector.Unboxed as U
import Data.Vector.Unboxed ((!))
import Safe (at)
import System.Random
import qualified Statistics.Sample as S
import qualified Statistics.Function as SF
import Data.Function (on)

--- * Simple linear regression

-- | Covariance of two samples
covar :: S.Sample -> S.Sample -> Double
covar xs ys = covar' m1 m2 n xs ys
    where
          !n = fromIntegral $ U.length xs
          !m1 = S.mean xs
          !m2 = S.mean ys
{-# INLINE covar #-}

-- internal function that avoids duplicate calculation of means and lengths where possible
-- Note: trying to make the calculation even more efficient by subtracting m1*m1*n instead of individual subtractions increased errors, probably due to rounding issues.
covar' :: Double -> Double -> Double -> S.Sample -> S.Sample -> Double
covar' m1 m2 n xs ys = U.sum (U.zipWith (*) (U.map (subtract m1) xs) (U.map (subtract m2) ys)) / (n-1)
{-# INLINE covar' #-}

-- | Pearson's product-moment correlation coefficient
correl :: S.Sample -> S.Sample -> Double
correl xs ys = let !c = covar xs ys
                   !sx = S.stdDev xs
                   !sy = S.stdDev ys
               in c / (sx * sy)
{-# INLINE correl #-}

-- | Simple linear regression between 2 samples.
--   Takes two vectors Y={yi} and X={xi} and returns
--   (alpha, beta, r*r) such that Y = alpha + beta*X
--   and where r is the Pearson product-moment correlation
--   coefficient
linearRegressionRSqr :: S.Sample -> S.Sample -> (Double, Double, Double)
linearRegressionRSqr xs ys = (alpha, beta, r2)
    where 
          !c                   = covar' m1 m2 n xs ys
          !r2                  = c*c / (v1*v2)
          !(m1,v1)             = S.meanVarianceUnb xs 
          !(m2,v2)             = S.meanVarianceUnb ys
          !n                   = fromIntegral $ U.length xs
          !beta                = c / v1
          !alpha               = m2 - beta * m1
{-# INLINE linearRegressionRSqr #-}
          
-- | Simple linear regression between 2 samples.
--   Takes two vectors Y={yi} and X={xi} and returns
--   (alpha, beta) such that Y = alpha + beta*X          
linearRegression :: S.Sample -> S.Sample -> (Double, Double)
linearRegression xs ys = (alpha, beta)
    where 
        (alpha, beta, _) = linearRegressionRSqr xs ys
{-# INLINE linearRegression #-}

-- | Total Least Squares (TLS) linear regression.
-- Assumes x-axis values (and not just y-axis values) are random variables and that both variables have similar distributions.
-- interface is the same as linearRegression.
linearRegressionTLS :: S.Sample -> S.Sample -> (Double,Double)
linearRegressionTLS xs ys = (alpha, beta)
    where
          !c                   = covar' m1 m2 n xs ys
          !b                   = (v1 - v2) / c
          !(m1,v1)             = S.meanVarianceUnb xs 
          !(m2,v2)             = S.meanVarianceUnb ys
          !n                   = fromIntegral $ U.length xs
          !betas               = [(-b - sqrt(b^2+4))/2,(-b + sqrt(b^2+4)) /2]
          !beta                = if c > 0 then maximum betas else minimum betas
          !alpha               = m2 - beta * m1
{-# INLINE linearRegressionTLS #-}

---- Robust regression based on "Computing LTS regression for large data sets" by Rousseeuw and Van Driessen 1999.
---- The methods implemented allow using the same scheme both for Least squares (where the errors are only along the vertical axis) and for total least squares, where errors are assumed on both variables.

---- first some helper functions that could probably be made more efficient:

-- get a sub-sample based on indices.
subSample :: [Int] -> S.Sample -> S.Sample
subSample indices original = U.generate (length indices) ((original `U.unsafeIndex`) . (indices `at`))

-- some type definitions to make things clear.
type EstimatedParams = (Double,Double)
type ErrorFunction = (EstimatedParams -> (Double,Double) -> Double)
type Estimator = (S.Sample -> S.Sample -> EstimatedParams)
data EstimationParameters = EstimationParameters { outlierFraction :: Double, estimator :: Estimator,errorFunction :: ErrorFunction }

-- Error functions for the two types of regression currently supported:
linearRegressionError :: ErrorFunction
linearRegressionError (alpha,beta) (x,y) = (y-(beta*x+alpha))^2

linearRegressionTLSError :: ErrorFunction
linearRegressionTLSError (alpha,beta) (x,y) = ey/(1+beta^2)
    where
        ey = linearRegressionError (alpha,beta) (x,y)

-- calculate the size of an expected set with no outliers.
setSize :: EstimationParameters -> S.Sample -> Int
setSize ep xs = round $ (1-(outlierFraction ep)) * (fromIntegral . U.length $ xs)

-- Given an initial estimate of the regression parameters - perform a "concentration" iteration giving at least as good an estimate as the previous one.
concentrationIteration :: EstimationParameters -> S.Sample -> S.Sample -> EstimatedParams -> EstimatedParams
concentrationIteration ep xs ys params = estimator ep good_xs good_ys
    where
        (good_xs,good_ys) = U.unzip good_sample
        set_size = setSize ep xs
        good_sample = U.take set_size . SF.sortBy (compare `on` (errorFunction ep params)) . U.zip xs $ ys

-- Given an estimate of the regression parameters, calculate the "quality" of the estimate by trimming the expected outliers group and summing over the errors of the remaining samples.
estimationQuality :: EstimationParameters -> EstimatedParams -> S.Sample -> S.Sample -> Double
estimationQuality ep params xs ys = U.sum . U.take set_size . SF.sort $ errors
    where
        errors = U.map (errorFunction ep params) . U.zip xs $ ys
        set_size = setSize ep xs

concentrationStep :: EstimatedParams -> S.Sample -> S.Sample -> (EstimatedParams, Double) -> (EstimatedParams, Double)
concentrationStep ep xs ys (prev, prev_err) = (new_estimate, new_err)
    where
        new_estimate = concentrationIteration ep xs ys prev
        new_err = estimationQuality ep new_estimate xs ys

converge :: EstimationParameters -> EstimatedParams -> S.Sample -> S.Sample -> EstimatedParams
converge ep estimate xs ys = fst $ converge' ep (estimate,err) xs ys
    where
        err = estimationQuality ep estimate xs ys

converge' :: EstimationParameters -> (EstimatedParams,Double) -> S.Sample -> S.Sample -> (EstimatedParams,Double)
converge' ep (prev, prev_err) xs ys
    | prev_err == new_err = (prev,prev_err)
    | otherwise = converge' ep (new_estimate, new_err) xs ys
        where
            new_estimate = concentrationIteration ep xs ys prev
            new_err = estimationQuality ep new_estimate xs ys

concentrateNSteps :: EstimationParameters -> EstimatedParams -> S.Sample -> S.Sample -> EstimatedParams
concentrateNSteps ep params xs ys = fst $ concentrateNSteps' ep (params, estimation_quality) xs ys
    where
        estimation_quality = estimationQuality ep params xs ys

{-
concentrateNSteps' :: EstimationParameters -> (EstimatedParams,Double) -> S.Sample -> S.Sample -> (EstimatedParams,Double)
concentrateNSteps' ep (prev, prev_err)

robustFit :: EstimationParameters -> S.Sample -> S.Sample -> EstimatedParams
robustFit ep xs ys
    | n < 2 = error "Robust fit - sample size has to include at least two points"
    | n <= smallSet ep = singleGroupFit ep xs ys
    | otherwise = subgroupFit ep xs ys 
        where
            n = U.length xs

singleGroupFit :: EstimationParameters -> S.Sample -> S.Sample -> EstimatedParams
singleGroupFit ep xs ys
    
    where
        n = U.length xs
        all_pairs = allPairs $ zip (U.toList xs) (U.toList ys)
        initial_sets
            | length all_pairs < (maxSubsetsNum ep) = all_pairs
            | otherwise = fst $ randomSubset (randomGen ep) all_pairs (maxSubsetsNum ep)
-}

allPairs :: [a] -> [(a,a)]
allPairs [] = []
allPairs [x] = []
allPairs [x,y] = [(x,y)]
allPairs (x:xs) = (zip xs . repeat $ x) ++ allPairs xs

randomSubset :: RandomGen g => g -> [a] -> Int -> ([a],g)
randomSubset g xs size
    | size <= 0 = ([],g)
    | otherwise = (head end : remaining_subset, final_g)
        where
            (index,ng) = randomR (0, length xs - 1) g
            (start,end) = splitAt index xs
            (remaining_subset,final_g) = randomSubset ng (start ++ (tail end)) (size-1)

----------------
-- use iterate
-- look for shuffle to generate a random subset
-- use state monad for random-using operations
