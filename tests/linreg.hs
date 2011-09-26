import Statistics.LinearRegression
import qualified Data.Vector.Unboxed as U

main = do
    let xs = U.fromList [1..100000]
    let ys = U.map (\x -> x*100 + 2000) xs
    putStrLn . show $ simpleLinReg xs ys