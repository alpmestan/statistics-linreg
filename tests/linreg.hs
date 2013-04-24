import Statistics.LinearRegression
import Statistics.Sample
import qualified Data.Vector.Unboxed as U
import Control.Monad.Random
import Control.Monad
import Control.Applicative
import System.Random.MWC
import System.Random.MWC.Distributions
import qualified Data.Packed.Vector as V
import Graphics.Rendering.Plot

main = do
    mapM_ test [1..10]
    test_convergence
    test_robust
    test_variances

    
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

getNormals :: Double -> Double -> Int -> IO [Double]
getNormals mean std n = do
    withSystemRandom . asGenIO $ \gen -> replicateM n (normal mean std gen)
        
testFigure :: U.Vector Double -> U.Vector Double -> (EstimatedRelation, EstimatedRelation, EstimatedRelation) -> Figure ()
testFigure xs ys (simple, non_robust, robust) = do
    let vxs = V.fromList . U.toList $ xs
    let vys = V.fromList . U.toList $ ys
    let dataset = (vxs, [   point vys Cross, line_func simple, line_func non_robust, line_func robust ])

    withTitle . setText $ "linreg test"
    setPlots 1 1
    withPlot (1,1) $ do
        addAxis XAxis (Side Lower) $ do
            setTicks Minor (TickNumber 5)
            withAxisLine $ do
                setLineWidth 1.0 
        addAxis YAxis (Side Lower) $ do
            setTicks Minor (TickNumber 5)
            withAxisLine $ do
                setLineWidth 1.0
        setDataset dataset
        setRangeFromData XAxis Lower Linear
        setRangeFromData YAxis Lower Linear
        setLegend True NorthEast Inside
    where
        line_func (alpha,beta) = line ((\x -> alpha + beta*x) :: Function) (1.0 :: LineWidth)
        
test_robust = do
    putStrLn "generating random dataset for robust fit:"
    first_xs <- getNormals 0.0 10.0 800
    first_ys_errs <- getNormals  0.0 1.0 800
    let first_ys = zipWith (+) first_xs first_ys_errs
    last_xs <- getNormals  50.0 (sqrt 50) 200
    last_ys <- getNormals  0.0 (sqrt 50) 200
    let xs = U.fromList $ first_xs ++ last_xs
    let ys = U.fromList $ first_ys ++ last_ys
    putStrLn "robustFit test results:"
    (simple,non_robust,robust) <- evalRandIO (randTest xs ys)
    putStrLn "linearRegression on dataset:"
    putStrLn . show $ simple
    putStrLn "convergedRegression on dataset:"
    putStrLn . show $ non_robust
    putStrLn "robustFit on dataset:"
    putStrLn . show $ robust
    let filename = "test_robust.png"
    putStrLn $ "Image output is at " ++ filename
    writeFigure PNG filename (800,800) $ testFigure xs ys (simple,non_robust,robust)

randTest :: MonadRandom m => U.Vector Double -> U.Vector Double -> m (EstimatedRelation,EstimatedRelation,EstimatedRelation)
randTest xs ys = do
    robust <- robustFit defaultEstimationParameters xs ys
    let simple = linearRegression xs ys
    let non_robust = converge defaultEstimationParameters xs ys (0.0,0.001) -- simple
    return (simple, non_robust, robust)

test_variances :: IO ()
test_variances = do
    putStrLn "generating random dataset for variance test:"
    let xs = U.fromList [-100..100]
    offsets <- liftM U.fromList $ getNormals 0 10 (U.length xs)
    let ys = U.zipWith (+) xs offsets
    let ab = linearRegression xs ys
    putStrLn $ "estimated fit should be (0,1). It is:" ++ show ab
    let mse = linearRegressionMSE ab xs ys
    putStrLn $ "Calculated MSE of sampled data should be an estimate of 10. it is:" ++ (show . sqrt $ mse)
    let vs = linearRegressionVariances ab xs ys
    putStrLn $ "Calculated variances of the linear fit are estimates of (0.5,1.4777e-4). They are:" ++ show vs
    
