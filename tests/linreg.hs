import Statistics.LinearRegression
import Statistics.Sample
import qualified Data.Vector.Unboxed as U

main = do
--    mapM_ test [1..10]
    test_convergence
    
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
    let ep = EstimationParameters { outlierFraction = 0.25, estimator = linearRegression, errorFunction = linearRegressionError}
    putStrLn "Initial iteration:"
    putStrLn . show $ iter1
    putStrLn "Successive iteration:"
    putStrLn . show $ converge ep (iter1, estimationQuality ep iter1 xs ys) xs ys
