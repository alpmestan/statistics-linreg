import Statistics.LinearRegression
import Statistics.Sample
import qualified Data.Vector.Unboxed as U
import Control.Monad.Random
import Control.Monad

main = do
--    mapM_ test [1..10]
    test_convergence
    test_robust

    
test k = do  
    let n = 10000000
    let a = k*n + 1
    let b = (k+1)*n
    let xs = U.fromList [a..b]
    let ys = U.map (\x -> x*100 + 2000) xs
    putStrLn "linearRegression:"
    putStrLn . show $ linearRegression xs ys
    putStrLn "linearRegressionTLS:"
    putStrLn . show $ linearRegressionTLS xs ys

test_convergence = do
    let xs = U.fromList [1..10]
    let ys = U.fromList $ [1..5] ++ [0,0] ++ [8..10]
    let iter1 = linearRegression xs ys
    let ep = defaultEstimationParameters
    putStrLn "Initial iteration:"
    putStrLn . show $ iter1
    putStrLn "Successive iteration:"
    putStrLn . show $ converge ep xs ys iter1

test_robust = do
    putStrLn "robustFit test:"
    (simple,non_robust,robust) <- evalRandIO randTest
    putStrLn "linearRegression:"
    putStrLn . show $ simple
    putStrLn "convergedRegression:"
    putStrLn . show $ non_robust
    putStrLn "robustFit:"
    putStrLn . show $ robust

randTest :: MonadRandom m => m (EstimatedParams,EstimatedParams,EstimatedParams)
randTest = do
    first_xs <- liftM (take 800) $ getRandomRs (0.0,100.0)
    first_ys_errs <- getRandomRs (0.0,1.0)
    let first_ys = zipWith (+) first_xs first_ys_errs
    last_xs <- liftM (take 200) $ getRandomRs (100.0,150.0)
    last_ys <- liftM (take 200) $ getRandomRs (25.0,75.0)
    let xs = U.fromList $ first_xs ++ last_xs
    let ys = U.fromList $ first_ys ++ last_ys
    robust <- robustFit defaultEstimationParameters xs ys
    let simple = linearRegression xs ys
    let non_robust = converge defaultEstimationParameters xs ys simple
    return (simple, non_robust, robust)
