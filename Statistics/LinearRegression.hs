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
    EstimatedRelation (..),
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

-- | An estimated linear relation between 2 samples is (alpha,beta) such that Y = alpha + beta*X.
type EstimatedRelation = (Double,Double)

-- | An 'Estimator' is a function that generates an estimated linear regression based on 2 samples. This module provides two estimator functions:
-- 'linearRegression' and 'linearRegressionTLS'
type Estimator = (S.Sample -> S.Sample -> EstimatedRelation)

-- | An 'ErrorFunction' is a function that computes the error of a given point from an estimate. This module provides two error functions correspoinding to the two 'Estimator' functions it defines:
-- * Vertical distance squared via 'linearRegressionError' that should be used with 'linearRegression'
-- 
-- * Total distance squared vie 'linearRegressionTLSError' that should be used with 'linearRegressionTLS'
type ErrorFunction = (EstimatedRelation -> (Double,Double) -> Double)

-- | The robust fit algorithm used has various parameters that can be specified using the 'EstimationParameters' record.
data EstimationParameters = EstimationParameters {
    -- | The maximal fraction of outliers expected in the sample (default 0.25)
    outlierFraction :: Double,
    -- | Number of concentration steps to take for initial evaluation of a solution (default 3)
    shortIterationSteps :: Int,
    -- | Maximal number of sampled subsets (pairs of points) to use as starting points for initial 'EstimationRelation's (default 500)
    maxSubsetsNum :: Int,
    -- | If the initial sample is large, and thus gets subdivided, this is the number of candidate-estimations to take from each subgroup, on which complete convergence will be executed (default 10)
    groupSubsets :: Int,
    -- | Maximal size of sample that can be analyzed without any sub-division (default 600)
    mediumSetSize :: Int,
    -- | Maximal size of set that does not require two-step sub-division (see reference article) (default 1500)
    largeSetSize :: Int,
    -- | Estimator function to use (default linearRegression)
    estimator :: Estimator,
    -- | ErrorFunction to use (default linearRegressionError)
    errorFunction :: ErrorFunction
    }

-- | Default set of parameters to use based on the referenced article.
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

-- | linearRegression error function is the square of the /vertical/ distance of a point from the line.
linearRegressionError :: ErrorFunction
linearRegressionError (alpha,beta) (x,y) = (y-(beta*x+alpha))^2

-- | linearRegressionTLS error function is the square of the /total/ distance of a point from the line.
linearRegressionTLSError :: ErrorFunction
linearRegressionTLSError (alpha,beta) (x,y) = ey/(1+beta^2)
    where
        ey = linearRegressionError (alpha,beta) (x,y)

-- | Helper function to calculate the minimal expected size of uncontaminated data based on the maximal fraction of outliers.
setSize :: EstimationParameters -> S.Sample -> Int
setSize ep xs = round $ (1-outlierFraction ep) * (fromIntegral . U.length $ xs)

-- | Helper function that, given an initial estimated relation, performs a "concentration" step (by generating a new estimate based on a fraction of points laying closest to the previous estimate) giving at least as good an estimate as the previous one.
concentrationIteration :: EstimationParameters -> S.Sample -> S.Sample -> EstimatedRelation -> EstimatedRelation
concentrationIteration ep xs ys params = estimator ep good_xs good_ys
    where
        (good_xs,good_ys) = U.unzip good_sample
        set_size = setSize ep xs
        good_sample = U.take set_size . SF.sortBy (compare `on` errorFunction ep params) . U.zip xs $ ys

-- | Helper function that quantifies the quality of an estimate by summing over the errors of the fraction of the samples that are not considered as outliers (outliers are a predefined fraction of samples that lay the furthest away from an estimate).
estimationQuality :: EstimationParameters -> S.Sample -> S.Sample -> EstimatedRelation -> Double
estimationQuality ep xs ys params = U.sum . U.take set_size . SF.sort $ errors
    where
        errors = U.map (errorFunction ep params) . U.zip xs $ ys
        set_size = setSize ep xs

-- | Helper function that pairs up an iteration and its quality.
----- to be combined with two functions above?
concentrationStep :: EstimationParameters -> S.Sample -> S.Sample -> (EstimatedRelation, Double) -> (EstimatedRelation, Double)
concentrationStep ep xs ys (prev, prev_err) = (new_estimate, new_err)
    where
        new_estimate = concentrationIteration ep xs ys prev
        new_err = estimationQuality ep xs ys new_estimate 

concentration :: EstimationParameters -> S.Sample -> S.Sample -> EstimatedRelation -> [(EstimatedRelation, Double)]
concentration ep xs ys params = iterate (concentrationStep ep xs ys) (params,err)
    where
        err = estimationQuality ep xs ys params

converge :: EstimationParameters -> S.Sample -> S.Sample -> EstimatedRelation -> EstimatedRelation
converge ep xs ys = fst . findConvergencePoint . concentration ep xs ys

findConvergencePoint :: Eq a => [(b,a)] -> (b,a)
findConvergencePoint (x:y:ys)
    | snd x == snd y = x
    | otherwise = findConvergencePoint (y:ys)
findConvergencePoint xs = error "Too short a list for conversion (size < 2)"

concentrateNSteps :: EstimationParameters -> S.Sample -> S.Sample -> EstimatedRelation -> (EstimatedRelation,Double)
concentrateNSteps ep xs ys params = concentration ep xs ys params !! shortIterationSteps ep

robustFit :: MonadRandom m => EstimationParameters -> S.Sample -> S.Sample -> m EstimatedRelation
robustFit ep xs ys = do
    let n = U.length xs
    if n < 2
        then
            error "cannot fit an input of size < 2"
        else if n == 2
            then return $ lineParams ((U.head xs,U.head ys),(U.last xs,U.last ys))
            else if n < mediumSetSize ep
                then liftM (candidatesToWinner ep xs ys) $ singleGroupFitCandidates ep Nothing xs ys
                else if n < largeSetSize ep
                    then largeGroupFit ep xs ys
                    else do
                        (nxs,nys) <- liftM unzip $ randomSubset (zip (U.toList xs) (U.toList ys)) (largeSetSize ep)
                        liftM (candidatesToWinner ep xs ys) $ largeGroupFitCandidates ep (U.fromList nxs) (U.fromList nys)
                
largeGroupFit :: MonadRandom m => EstimationParameters -> S.Sample -> S.Sample -> m EstimatedRelation
largeGroupFit ep xs ys = liftM (candidatesToWinner ep xs ys) $ largeGroupFitCandidates ep xs ys

candidatesToWinner :: EstimationParameters -> S.Sample -> S.Sample -> [EstimatedRelation] -> EstimatedRelation
candidatesToWinner ep xs ys = fst . minimumBy (compare `on` snd) . map (findConvergencePoint . concentration ep xs ys)

largeGroupFitCandidates :: MonadRandom m => EstimationParameters -> S.Sample -> S.Sample -> m [EstimatedRelation]
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
    
singleGroupFitCandidates :: MonadRandom m => EstimationParameters -> Maybe Int -> S.Sample -> S.Sample -> m [EstimatedRelation]
singleGroupFitCandidates ep m_subsets xs ys = do
    let all_pairs = allPairs $ zip (U.toList xs) (U.toList ys)
    initial_sets <- randomSubset all_pairs $ fromMaybe (maxSubsetsNum ep) m_subsets
    return . map fst . take (groupSubsets ep) . sortBy (compare `on` snd) . map (concentrateNSteps ep xs ys . lineParams) $ initial_sets

lineParams :: ((Double,Double),(Double,Double)) -> EstimatedRelation
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
