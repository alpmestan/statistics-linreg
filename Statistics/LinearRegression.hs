{-# LANGUAGE BangPatterns #-}

module Statistics.LinearRegression (linearRegression, correl, covar) where

import qualified Data.Vector.Unboxed as U
import qualified Statistics.Sample as S

--- * Simple linear regression

-- | Covariance of two samples
covar :: S.Sample -> S.Sample -> Double
covar xs ys = U.sum (U.zipWith (*) (U.map f1 xs) (U.map f2 ys)) / (n-1)
    where
          !n = fromIntegral $ U.length xs
          !m1 = S.mean xs
          !m2 = S.mean ys
          f1 = \x -> (x - m1)
          f2 = \x -> (x - m2)
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
          !c                   = U.sum (U.zipWith (*) (U.map (subtract m1) xs) (U.map (subtract m2) ys)) / (n-1)
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
          !c                   = U.sum (U.zipWith (*) (U.map (subtract m1) xs) (U.map (subtract m2) ys)) / (n-1)
          !r                   = c / (sx * sy)
          !m1                  = S.mean xs 
          !m2                  = S.mean ys
          !sx                  = S.stdDev xs
          !sy                  = S.stdDev ys
          !n                   = fromIntegral $ U.length xs
          !beta                = r * sy / sx
          !alpha               = m2 - beta * m1
{-# INLINE linearRegression #-}