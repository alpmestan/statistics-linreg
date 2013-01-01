{-# LANGUAGE BangPatterns #-}

module Statistics.LinearRegression (
    -- * Simple linear regression functions
    linearRegression,
    linearRegressionRSqr,
    linearRegressionTLS,
    -- * related functions
    correl,
    covar,
    -- * Robust linear regression
    robustFit,
    nonRandomRobustFit,
    robustFitRSqr,
    -- ** Related types
    EstimationParameters(..),
    ErrorFunction,
    Estimator,
    EstimatedRelation,
    -- ** Provided values
    defaultEstimationParameters,
    linearRegressionError,
    linearRegressionTLSError,
    -- ** Helper functions
    converge,
    -- * References
    -- $references
    ) where

import qualified Data.Vector.Unboxed as U
import Data.Vector.Unboxed ((!))
import Safe (at)
import System.Random
import System.Random.Shuffle (shuffleM)
import Control.Monad.Random.Class
import Control.Monad.Random (evalRand)
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
-- interface is the same as 'linearRegression'.
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

-- | An estimated linear relation between 2 samples is (alpha,beta) such that Y = alpha + beta*X.
type EstimatedRelation = (Double,Double)

-- | An 'Estimator' is a function that generates an estimated linear regression based on 2 samples. This module provides two estimator functions:
-- 'linearRegression' and 'linearRegressionTLS'
type Estimator = (S.Sample -> S.Sample -> EstimatedRelation)

-- | An 'ErrorFunction' is a function that computes the error of a given point from an estimate. This module provides two error functions correspoinding to the two 'Estimator' functions it defines:
-- 
-- * Vertical distance squared via 'linearRegressionError' that should be used with 'linearRegression'
-- 
-- * Total distance squared vie 'linearRegressionTLSError' that should be used with 'linearRegressionTLS'
type ErrorFunction = (EstimatedRelation -> (Double,Double) -> Double)

-- | The robust fit algorithm used has various parameters that can be specified using the 'EstimationParameters' record.
data EstimationParameters = EstimationParameters {
    -- | Maximal fraction of outliers expected in the sample (default 0.25)
    outlierFraction     :: Double,
    -- | Number of concentration steps to take for initial evaluation of a solution (default 3)
    shortIterationSteps :: Int,
    -- | Maximal number of sampled subsets (pairs of points) to use as starting points (default 500)
    maxSubsetsNum       :: Int,
    -- | If the initial sample is large, and thus gets subdivided, this is the number of candidate-estimations to take from each subgroup, on which complete convergence will be executed (default 10)
    groupSubsets        :: Int,
    -- | Maximal size of sample that can be analyzed without any sub-division (default 600)
    mediumSetSize       :: Int,
    -- | Maximal size of sample that does not require two-step sub-division (see reference article) (default 1500)
    largeSetSize        :: Int,
    -- | Estimator function to use (default linearRegression)
    estimator           :: Estimator,
    -- | ErrorFunction to use (default linearRegressionError)
    errorFunction       :: ErrorFunction
    }

-- | Default set of parameters to use (see reference for details).
defaultEstimationParameters = EstimationParameters {
    outlierFraction = 0.25,
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
setSize ep xs = max (n `div` 2 + 1) . round $ (1-outlierFraction ep) * (fromIntegral n)
    where
        n = U.length xs

-- | Helper function that, given an initial estimated relation and the error of the perivous estimation, performs a "concentration" step, generating a new estimate based on a fraction of points laying closest to the previous estimate and estimates the error of the previous estimate based on the same fraction.
-- The result is an estimate that is at least as good as the previous one.
-- The reason the error is calculated for the previous parameters is calculation optimization.
concentrationStep :: EstimationParameters -> S.Sample -> S.Sample -> (EstimatedRelation, Double) -> (EstimatedRelation, Double)
concentrationStep ep xs ys (prev, prev_err) = (new_estimate, new_err)
    where
        set_size = setSize ep xs
        xyerrors = U.map (\p -> (p,errorFunction ep prev p)) $ U.zip xs ys
        (xys,errors) = U.unzip . U.take set_size . SF.sortBy (compare `on` snd) $ xyerrors
        (good_xs,good_ys) = U.unzip xys
        new_estimate = estimator ep good_xs good_ys
        new_err = U.sum errors

-- | Infinite set of consecutive concentration steps.
concentration :: EstimationParameters -> S.Sample -> S.Sample -> EstimatedRelation -> [(EstimatedRelation, Double)]
concentration ep xs ys params = iterate (concentrationStep ep xs ys) (params,-1)

-- | Calculate the optimal (local minimum) estimate based on an initial estimate.
-- The local minimum may not be the global (a.k.a. best) estimate but starting from enough different initial estimates should yield the global optimum eventually.
converge :: EstimationParameters -> S.Sample -> S.Sample -> EstimatedRelation -> EstimatedRelation
converge ep xs ys = fst . findConvergencePoint . concentration ep xs ys

-- | The convergence point is defined as the point the error estimate of which is equal to the next estimate's error.
findConvergencePoint :: Eq a => [(b,a)] -> (b,a)
findConvergencePoint (x:y:ys)
    | snd x == snd y = x
    | otherwise = findConvergencePoint (y:ys)
findConvergencePoint xs = error "Too short a list for conversion (size < 2)"

-- | Many times there is no need for full concentration as bad initial estimates can be discovered after only a few concentration steps.
concentrateNSteps :: EstimationParameters -> S.Sample -> S.Sample -> EstimatedRelation -> (EstimatedRelation,Double)
concentrateNSteps ep xs ys params = concentration ep xs ys params !! shortIterationSteps ep

-- | Finding a robust fit linear estimate between two samples. The procedure requires randomization and is based on the procedure described in the reference.
robustFit :: MonadRandom m => EstimationParameters -> S.Sample -> S.Sample -> m EstimatedRelation
robustFit ep xs ys = do
    let n = U.length xs
-- For optimal performance the exact procedure executed depends on the set size.
    if n < 2
        then
            error "cannot fit an input of size < 2"
        else if n == 2
            then return $ lineParams ((U.head xs,U.head ys),(U.last xs,U.last ys))
            else 
                liftM (candidatesToWinner ep xs ys) $ if n < mediumSetSize ep
                    then
                         singleGroupFitCandidates ep Nothing xs ys
                    else if n < largeSetSize ep
                        then largeGroupFitCandidates ep xs ys
                        else do
                            (nxs,nys) <- liftM unzip $ randomSubset (zip (U.toList xs) (U.toList ys)) (largeSetSize ep)
                            largeGroupFitCandidates ep (U.fromList nxs) (U.fromList nys)

-- | Robust fit yielding also the R-square value of the \"clean\" dataset.
robustFitRSqr :: MonadRandom m => EstimationParameters -> S.Sample -> S.Sample -> m (EstimatedRelation,Double)
robustFitRSqr ep xs ys = do
    er <- robustFit ep xs ys
    let (good_xs,good_ys) = U.unzip . U.take (setSize ep xs) . SF.sortBy (compare `on` errorFunction ep er) $ U.zip xs ys
    return (er,correl good_xs good_ys ^ 2)

-- | A wrapper that executes 'robustFit' using a default random generator (meaning it is only pseudo-random)
nonRandomRobustFit :: EstimationParameters -> S.Sample -> S.Sample -> EstimatedRelation
nonRandomRobustFit ep xs ys = evalRand (robustFit ep xs ys) (mkStdGen 1)

-- | Given a set of initial estimates converge them all and find the optimal one.
candidatesToWinner :: EstimationParameters -> S.Sample -> S.Sample -> [EstimatedRelation] -> EstimatedRelation
candidatesToWinner ep xs ys = fst . minimumBy (compare `on` snd) . map (findConvergencePoint . concentration ep xs ys)

-- | for a large initial sample - subdivide it, then get candidates from each subgroup. Perform full convergence on all the candidates and return the best ones.
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

-- | For a single group (a group that will not be subdivided) pick an initial set of pairs of points, run a few steps on each, then return the most promising candidates.
singleGroupFitCandidates :: MonadRandom m => EstimationParameters -> Maybe Int -> S.Sample -> S.Sample -> m [EstimatedRelation]
singleGroupFitCandidates ep m_subsets xs ys = do
    let all_pairs = allPairs $ zip (U.toList xs) (U.toList ys)
    let return_size = fromMaybe (maxSubsetsNum ep) m_subsets
    initial_sets <- if return_size > length all_pairs
        then return all_pairs
        else randomSubset all_pairs return_size 
    return . map fst . take (groupSubsets ep) . sortBy (compare `on` snd) . map (concentrateNSteps ep xs ys . lineParams) $ initial_sets

-- | Find the line passing between two points. This is the initial estimate to use given two random points.
lineParams :: ((Double,Double),(Double,Double)) -> EstimatedRelation
lineParams ((x1,y1),(x2,y2)) = (alpha,beta)
    where
        beta = (y2-y1)/(x2-x1)
        alpha = y1 - beta*x1

-- | A list of all possible two-element pairs from a list.
allPairs :: [a] -> [(a,a)]
allPairs [] = []
allPairs [x] = []
allPairs [x,y] = [(x,y)]
allPairs (x:xs) = (zip xs . repeat $ x) ++ allPairs xs

-- | Get a random subset of a given size.
randomSubset :: MonadRandom m => [a] -> Int -> m [a]
randomSubset xs size = liftM (take size) $ shuffleM xs

-- | Split a list into sublists of length n.
splitTo :: Int -> [a] -> [[a]]
splitTo n = map (take n) . takeWhile (not . null) . iterate (drop n)

-- | Helper function to adjust parameter handling
applyTo :: (a->b->c) -> (a,b) -> c
applyTo f (x,y) = f x y

-- $references
--
-- * Two Dimensional Euclidean Regression (Stein) <http://www.dspcsp.com/pubs/euclreg.pdf>
--
-- * Computing LTS Regression For Large Data Sets (Rousseeuw and Driessen) <http://agoras.ua.ac.be/abstract/Comlts99.htm>

