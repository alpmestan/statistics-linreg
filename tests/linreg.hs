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

getNormals :: Double -> Double -> Int -> IO [Double]
getNormals mean var n = do
    withSystemRandom . asGenIO $ \gen -> replicateM n (normal mean var gen)
        
testFigure :: U.Vector Double -> U.Vector Double -> (EstimatedParams, EstimatedParams, EstimatedParams) -> Figure ()
testFigure xs ys (simple, non_robust, robust) = do
    let vxs = V.fromList . U.toList $ xs
    let vys = V.fromList . U.toList $ ys
    let dataset = (vxs, [   point vys Cross, line_func simple, line_func non_robust, line_func robust ])

    withTitle . setText $ "linreg test"
    setPlots 1 1
    withPlot (1,1) $ do
        addAxis XAxis (Side Lower) $ do
            setTicks Minor (Left 5)
            withAxisLine $ do
                setLineWidth 1.0 
        addAxis YAxis (Side Lower) $ do
            setTicks Minor (Left 5)
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
    first_xs <- getNormals 0.0 10.0 8000
    first_ys_errs <- getNormals  0.0 1.0 8000
    let first_ys = zipWith (+) first_xs first_ys_errs
    last_xs <- liftM (map (50+)) $ getNormals  0.0 (sqrt 50) 2000
    last_ys <- getNormals  0.0 (sqrt 50) 2000
    let xs = U.fromList $ first_xs ++ last_xs
    let ys = U.fromList $ first_ys ++ last_ys
    putStrLn "robustFit test:"
    (simple,non_robust,robust) <- evalRandIO (randTest xs ys)
    putStrLn "linearRegression:"
    putStrLn . show $ simple
    putStrLn "convergedRegression:"
    putStrLn . show $ non_robust
    putStrLn "robustFit:"
    putStrLn . show $ robust
    writeFigure PNG ("test.png") (800,800) $ testFigure xs ys (simple,non_robust,robust)

randTest :: MonadRandom m => U.Vector Double -> U.Vector Double -> m (EstimatedParams,EstimatedParams,EstimatedParams)
randTest xs ys = do
    robust <- robustFit defaultEstimationParameters xs ys
    let simple = linearRegression xs ys
    let non_robust = converge defaultEstimationParameters xs ys (0.0,0.001) -- simple
    return (simple, non_robust, robust)
