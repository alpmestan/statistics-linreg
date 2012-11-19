import Statistics.LinearRegression
import Statistics.Sample
import qualified Data.Vector.Unboxed as U

main = do
    mapM_ test [1..10]
    
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
