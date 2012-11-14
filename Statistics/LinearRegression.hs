{-# LANGUAGE BangPatterns #-}

module Statistics.LinearRegression (
    linearRegressionRSqr,
    linearRegression,
    correl,
    covar,
    linearRegressionTLS
    ) where

import qualified Data.Vector.Unboxed as U
import qualified Statistics.Sample as S

--- * Simple linear regression

-- | Covariance of two samples
covar :: S.Sample -> S.Sample -> Double
covar xs ys = U.sum (U.zipWith (*) (U.map (subtract m1) xs) (U.map (subtract m2) ys)) / (n-1)
    where
          !n = fromIntegral $ U.length xs
          !m1 = S.mean xs
          !m2 = S.mean ys
{-# INLINE covar #-}


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
linearRegressionRSqr xs ys = (alpha, beta, r*r)
    where 
          !c                   = covar xs ys
          !r                   = c / (sx * sy)
          !m1                  = S.mean xs 
          !m2                  = S.mean ys
          !sx                  = S.stdDev xs
          !sy                  = S.stdDev ys
          !n                   = fromIntegral $ U.length xs
          !beta                = r * sy / sx
          !alpha               = m2 - beta * m1
{-# INLINE linearRegressionRSqr #-}
          
-- | Simple linear regression between 2 samples.
--   Takes two vectors Y={yi} and X={xi} and returns
--   (alpha, beta, r*r) such that Y = alpha + beta*X          
linearRegression :: S.Sample -> S.Sample -> (Double, Double)
linearRegression xs ys = (alpha, beta)
    where 
        (alpha, beta, _) = linearRegressionRSqr xs ys
{-# INLINE linearRegression #-}

-- | Total Least Squares (TLS) linear regression assumes x-axis values are also random variables (as opposed to simple linear regression that assumes no errors in the x values).
-- It does assume the distribution of the x values is the same as that of the y values (a.k.a. same error bars/standard deviations).
-- interface is the same as linearRegression.
linearRegressionTLS :: S.Sample -> S.Sample -> (Double, Double, Double,Double)
linearRegressionTLS xs ys = (alpha1, beta1,alpha2, beta2)
    where
          !alpha               = m2 - beta * m1
          !beta                = max (-b - sqrt( b^2-4*a*c)) (-b + sqrt( b^2-4*a*c)) /2*a
          !a                   = covar xs ys
          !b                   = S.variance xs - (S.variance ys)
          !c                   = - covar xs ys
          !m1                  = S.mean xs 
          !m2                  = S.mean ys
          !r                   = (S.variance xs - (S.variance ys)) / (covar xs ys)
          !alpha1               = m2 - beta1 * m1
          !alpha2               = m2 - beta2 * m1
          !beta1                = (-r + sqrt (r^2 + 4)) / 2
          !beta2                = (-r - sqrt (r^2 + 4)) / 2
{-# INLINE linearRegressionTLS #-}
