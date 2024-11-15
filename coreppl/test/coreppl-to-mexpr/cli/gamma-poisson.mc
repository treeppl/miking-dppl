include "../../cppl-test.mc"
include "../../test.mc"

include "seq.mc"
include "sys.mc"
include "string.mc"
include "common.mc"
include "stats.mc"

-- Utilitaries
recursive
let quantile: Float -> Float -> [(Float, Float)] -> Float = lam inf. lam sup. lam cdf. 
  if eqi (length cdf) 1 then
    0.0
  else
    if gtf (head cdf).0 inf then
      subf (quantile_sup sup (tail cdf)) (head cdf).1  
    else
      quantile inf sup (tail cdf)

let quantile_sup: Float -> [(Float, Float)] -> Float = lam sup. lam cdf.
  if eqi (length cdf) 1 then
    1.0
  else
    if gtf (head cdf).0 sup then
      (head cdf).1
    else
      quantile_sup sup (tail cdf)
end

recursive
let cumul : [Float] -> Float = lam vec.
  if neqi (length vec) 0 then
    addf (head vec) (cumul (tail vec))
  else
    0.0
end

recursive
let distrib: [Float] -> [Float] -> [(Float, Float)] = lam samp. lam weigh. 
  if null samp then [] else 
    concat [((head samp),(head weigh))] (distrib (tail samp) (tail weigh))
end

recursive
let quickSort_distrib : all a. (a -> a -> Float) -> ([(a, a)] -> [(a,a)]) = lam cmp. lam seq.
  if null seq then seq else
    let h = head seq in
    let t = tail seq in
    let lr = partition (lam x. ltf (cmp x.0 h.0) 0.0) t in
    concat (quickSort_distrib cmp lr.0) (cons h (quickSort_distrib cmp lr.1))
end

recursive
let revers_cdf: [(Float, Float)] -> [(Float, Float)] -> [(Float, Float)] = lam distrib. lam vecf.
  let cumul = addf (head distrib).1 (head vecf).1 in
  let vecf = concat [((head distrib).0, cumul)] vecf in
  if neqi (length distrib) 1 then
    revers_cdf (tail distrib) vecf
  else
    vecf
end

-- recursive
-- let obs_cdf: [(Float, Float)] -> Int = lam distrib.
--   print (float2string (head distrib).0);
--   print ",";
--   print (float2string (head distrib).1);
--   print "\n";
--   if neqi (length distrib) 1 then
--     obs_cdf (tail distrib)
--   else
--     print "\n";
--     0
-- end

-- recursive
-- let obs_distr: [Float] -> [Float] -> Int = lam samples. lam wei.
--   print (float2string (head samples));
--   print ",";
--   print (float2string (head wei));
--   print "\n";
--   if neqi (length samples) 1 then
--     obs_distr (tail samples) (tail wei)
--   else
--     print "\n";
--     0
-- end

-- recursive
-- let obs_list: Int -> [Float] -> Int = lam ind. lam list.
--   print (int2string ind);
--   print ",";
--   print (float2string (head list));
--   print "\n";
--   if neqi (length list) 1 then
--     obs_list (addi 1 ind) (tail list)
--   else
--     print "\n";
--     0
-- end

mexpr

let s = 10e-3 in
let eq = eqCoin s in
let t = testCpplMExpr "gamma-poisson.mc" 20000 in
let res = [
  t 0   "-m is-lw --cps none",
  t 0   "-m is-lw --cps partial",
  t 0   "-m is-lw --cps partial --no-early-stop",
  t 0   "-m is-lw --cps full",
  t 0   "-m is-lw --cps full --no-early-stop",
  t 0   "-m smc-bpf --cps partial --resample manual",
  t 0   "-m smc-bpf --cps partial --resample align",
  t 0   "-m smc-bpf --cps partial --resample likelihood",
  t 0   "-m smc-bpf --cps full --resample manual",
  t 0   "-m smc-bpf --cps full --resample align",
  t 0   "-m smc-bpf --cps full --resample likelihood",
  t 0   "-m smc-apf --cps partial --resample manual",
  t 0   "-m smc-apf --cps partial --resample align",
  t 0   "-m smc-apf --cps partial --resample likelihood",
  t 0   "-m smc-apf --cps full --resample manual",
  t 0   "-m smc-apf --cps full --resample align",
  t 0   "-m smc-apf --cps full --resample likelihood",
  t 500 "-m pmcmc-pimh --cps partial",
  t 500 "-m pmcmc-pimh --cps full",
  t 500 "-m mcmc-trace",
  t 500 "-m mcmc-naive",
  t 500 "-m mcmc-lightweight --align --cps none ",
  t 500 "-m mcmc-lightweight --align --cps partial",
  t 500 "-m mcmc-lightweight --align --cps full",
  t 500 "-m mcmc-lightweight",
  t 500 "-m mcmc-lw-dk --align --cps none --kernel --drift 2.0 --mcmc-lw-gprob 0.0",
  t 500 "-m mcmc-lw-dk --align --cps partial --kernel --drift 2.0 --mcmc-lw-gprob 0.0"
] in

recursive
let test : [CpplRes] -> Int = lam resVec. 
  let result = head(resVec) in
  let weights = (map exp result.lweights) in
  -- let test = obs_list weights in
  -- re normalize after truncature of burning
  let acc = cumul weights in
  let distr = distrib (map string2float result.samples) (map (mulf (divf 1.0 acc)) weights)  in
  let sort = quickSort_distrib subf distr in
  let temp = revers_cdf sort [(0.0, 0.0)] in
  let cdf_dist = tail (reverse temp) in

  utest (quantile 100.2653 104.2089 cdf_dist) with 0.95 using eq in
  utest (quantile 100.5786 103.8882 cdf_dist) with 0.9 using eq in
  utest (quantile 100.9405 103.5191 cdf_dist) with 0.8 using eq in
  utest (quantile 101.3801 103.0735 cdf_dist) with 0.6 using eq in
  printLn " EoM";

  if neqi (length (tail resVec)) 0 then
    test (tail resVec)
  else
    print " EoT";
    1
in

let u = test res in
()