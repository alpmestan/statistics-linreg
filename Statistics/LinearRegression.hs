module Statistics.LinearRegression (linearRegression) where

import qualified Data.Vector.Unboxed as U
import Statistics.Sample

--- * Simple linear regression

-- | Covariance
covar :: Sample -> Sample -> (Double, Double, Double)
covar xs ys = (U.sum (U.zipWith (*) (U.map f1 xs) (U.map f2 ys)) / (n-1), m1, m2)
    where
          n = fromIntegral $ U.length xs
          m1 = mean xs
          m2 = mean ys
          f1 = \x -> (x - m1)
          f2 = \x -> (x - m2)

-- | Pearson's product-moment correlation coefficient
correl :: Sample -> Sample -> (Double, Double, Double, Double, Double)
correl xs ys = let (c, m1, m2) = covar xs ys
                   sx = stdDev xs
                   sy = stdDev ys
               in (c / (stdDev xs * stdDev ys), m1, m2, sx, sy)

-- | Simple linear regression between 2 samples.
--   Takes two vectors Y={yi} and X={xi} and returns
--   (alpha, beta) such that Y = alpha + betaX
linearRegression :: Sample -> Sample -> (Double, Double)
linearRegression xs ys = (alpha, beta)
    where 
          (r, m1, m2, sx, sy) = correl xs ys
          beta                = r * sy / sx
          alpha               = m2 - beta * m1
