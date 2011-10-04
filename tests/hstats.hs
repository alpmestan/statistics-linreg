import Math.Statistics

main = do
    let xs = [1..10000000]
    let ys = map (\x -> x*100 + 2000) xs
    putStrLn . show $ linreg $ zip xs ys       