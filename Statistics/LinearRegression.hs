{-# LANGUAGE BangPatterns #-}

module Statistics.LinearRegression (
    linearRegression,
    linearRegressionRSqr,
    linearRegressionTLS,
    linearRegressionError,
    linearRegressionTLSError,
    EstimationParameters(..),
    defaultEstimationParameters,
    robustFit,
    EstimatedParams (..),
    estimationQuality,
    converge,
    correl,
    covar,
    ) where

import qualified Data.Vector.Unboxed as U
import Data.Vector.Unboxed ((!))
import Safe (at)
import System.Random
import System.Random.Shuffle (shuffleM)
import Control.Monad.Random.Class
import Control.Monad (liftM)
import Data.Function (on)
import Data.List (minimumBy, sortBy)
import Data.Maybe (fromMaybe)
import qualified Statistics.Sample as S
import qualified Statistics.Function as SF

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
---- The user can also supply his/her own error and regression functions.

---- first some helper functions that could probably be made more efficient:

-- get a sub-sample based on indices.
subSample :: [Int] -> S.Sample -> S.Sample
subSample indices original = U.generate (length indices) ((original `U.unsafeIndex`) . (indices `at`))

-- some type definitions to make things clear.
type EstimatedParams = (Double,Double)
type ErrorFunction = (EstimatedParams -> (Double,Double) -> Double)
type Estimator = (S.Sample -> S.Sample -> EstimatedParams)
data EstimationParameters = EstimationParameters {
    outlierFraction :: Double,  -- what is the maximal fraction of outliers expected in the sample (default 0.25)
    shortIterationSteps :: Int, -- how many concentration steps to take for initial evaluation of a solution (default 3)
    maxSubsetsNum :: Int,       -- what is the maximal number of sampled subsets (pairs of points) to use as starting points for estimation parameters (default 500)
    groupSubsets :: Int,        -- how many candidate-estimations to take from a subgroup of starting samples for complete convergence (default 10)
    mediumSetSize :: Int,       -- what is the maximal size of set that can be analyzed without sub-division (default 600)
    largeSetSize :: Int,        -- what is the maximal size of set that does not require two-step sub-division (default 1500)
    estimator :: Estimator,     -- an estimator function that, given a set of sample points generates a linear approximation of the points (default linearRegression)
    errorFunction :: ErrorFunction -- an error function that, given an estimator and a point, returns the error of the point w.r.t. the estimator (default linearRegressionError)
    }

defaultEstimationParameters = EstimationParameters {
    outlierFraction = 0.75,
    shortIterationSteps = 3,
    maxSubsetsNum = 500,
    groupSubsets = 10,
    mediumSetSize = 600,
    largeSetSize = 1500,
    estimator = linearRegression,
    errorFunction = linearRegressionError
}
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
estimationQuality :: EstimationParameters -> S.Sample -> S.Sample -> EstimatedParams -> Double
estimationQuality ep xs ys params = U.sum . U.take set_size . SF.sort $ errors
    where
        errors = U.map (errorFunction ep params) . U.zip xs $ ys
        set_size = setSize ep xs

concentrationStep :: EstimationParameters -> S.Sample -> S.Sample -> (EstimatedParams, Double) -> (EstimatedParams, Double)
concentrationStep ep xs ys (prev, prev_err) = (new_estimate, new_err)
    where
        new_estimate = concentrationIteration ep xs ys prev
        new_err = estimationQuality ep xs ys new_estimate 

concentration :: EstimationParameters -> S.Sample -> S.Sample -> EstimatedParams -> [(EstimatedParams, Double)]
concentration ep xs ys params = iterate (concentrationStep ep xs ys) (params,err)
    where
        err = estimationQuality ep xs ys params

converge :: EstimationParameters -> S.Sample -> S.Sample -> EstimatedParams -> EstimatedParams
converge ep xs ys = fst . findConvergencePoint . concentration ep xs ys

findConvergencePoint :: Eq a => [(b,a)] -> (b,a)
findConvergencePoint (x:y:ys)
    | snd x == snd y = x
    | otherwise = findConvergencePoint (y:ys)
findConvergencePoint xs = error "Too short a list for conversion (size < 2)"

concentrateNSteps :: EstimationParameters -> S.Sample -> S.Sample -> EstimatedParams -> (EstimatedParams,Double)
concentrateNSteps ep xs ys params = concentration ep xs ys params !! (shortIterationSteps ep)

robustFit :: MonadRandom m => EstimationParameters -> S.Sample -> S.Sample -> m EstimatedParams
robustFit ep xs ys = do
    let n = U.length xs
    if n < 2
        then
            error "cannot fit an input of size < 2"
        else if n == 2
            then return $ lineParams ((U.head xs,U.head ys),(U.last xs,U.last ys))
            else if n < mediumSetSize ep
                then liftM (candidatesToWinner ep xs ys) $ singleGroupFitCandidates ep Nothing xs $ ys
                else if n < largeSetSize ep
                    then largeGroupFit ep xs ys
                    else do
                        (nxs,nys) <- liftM unzip $ randomSubset (zip (U.toList xs) (U.toList ys)) (largeSetSize ep)
                        liftM (candidatesToWinner ep xs ys) $ largeGroupFitCandidates ep (U.fromList nxs) (U.fromList nys)
                
largeGroupFit :: MonadRandom m => EstimationParameters -> S.Sample -> S.Sample -> m EstimatedParams
largeGroupFit ep xs ys = liftM (candidatesToWinner ep xs ys) $ largeGroupFitCandidates ep xs ys

candidatesToWinner :: EstimationParameters -> S.Sample -> S.Sample -> [EstimatedParams] -> EstimatedParams
candidatesToWinner ep xs ys = fst . minimumBy (compare `on` snd) . map (findConvergencePoint . concentration ep xs ys)

largeGroupFitCandidates :: MonadRandom m => EstimationParameters -> S.Sample -> S.Sample -> m [EstimatedParams]
largeGroupFitCandidates ep xs ys = do
    let n = U.length xs
    let sub_groups_num = n `div` (mediumSetSize ep `div` 2)
    let sub_groups_size = n `div` sub_groups_num
    shuffled <- shuffleM $ zip (U.toList xs) (U.toList ys)
    let sub_groups = map (U.unzip . U.fromList) $ splitTo sub_groups_size shuffled
    let sub_groups_candidates = maxSubsetsNum ep `div` sub_groups_num
    candidates_list <- mapM (applyTo $ singleGroupFitCandidates ep (Just sub_groups_candidates)) sub_groups
    let candidates = concat candidates_list
    return . map fst . take (groupSubsets ep) . sortBy (compare `on` snd) . map (findConvergencePoint . concentration ep xs ys) $ candidates

splitTo :: Int -> [a] -> [[a]]
splitTo n = map (take n) . takeWhile (not . null) . iterate (drop n)

applyTo :: (a->b->c) -> (a,b) -> c
applyTo f (x,y) = f x y
    
singleGroupFitCandidates :: MonadRandom m => EstimationParameters -> Maybe Int -> S.Sample -> S.Sample -> m [EstimatedParams]
singleGroupFitCandidates ep m_subsets xs ys = do
    let all_pairs = allPairs $ zip (U.toList xs) (U.toList ys)
    initial_sets <- randomSubset all_pairs $ fromMaybe (maxSubsetsNum ep) m_subsets
    return . map fst . take (groupSubsets ep) . sortBy (compare `on` snd) . map (concentrateNSteps ep xs ys . lineParams) $ initial_sets

lineParams :: ((Double,Double),(Double,Double)) -> EstimatedParams
lineParams ((x1,y1),(x2,y2)) = (alpha,beta)
    where
        beta = (y2-y1)/(x2-x1)
        alpha = y1 - beta*x1

allPairs :: [a] -> [(a,a)]
allPairs [] = []
allPairs [x] = []
allPairs [x,y] = [(x,y)]
allPairs (x:xs) = (zip xs . repeat $ x) ++ allPairs xs

randomSubset :: MonadRandom m => [a] -> Int -> m [a]
randomSubset xs size = liftM (take size) $ shuffleM xs
