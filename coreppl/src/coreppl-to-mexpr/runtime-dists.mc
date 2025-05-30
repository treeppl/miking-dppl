-- Runtime support for first-class distributions in coreppl-to-mexpr compiler.

include "math.mc"
include "seq.mc"
include "map.mc"
include "float.mc"
include "ext/dist-ext.mc"

-- Weiner process
let wienerSample : () -> Float -> Float = lam.
  -- TODO(oerikss, 2024-03-20): We save the observed trace of the process in a
  -- map. This is ofcourse limiting in the sense that we need to keep the
  -- observed traces of all Wiener process in memory. In the futurue we want to
  -- store the trace to disk if it gets too big.

  -- Adapted from:
  -- https://github.com/lazyppl-team/lazyppl/blob/89f162641a27e2db3c8632f9ae32c33b78db5ea2/src/WienerDemo.lhs#L70
  let trace = ref (mapSingleton cmpf 0. 0.) in
  lam t.
    let tr = deref trace in
    mapFindOrElse
      (lam.
        let w =
          -- Brownian bridge: https://en.wikipedia.org/wiki/Brownian_bridge
          switch (mapFindLower t tr, mapFindUpper t tr)
          case (None _, None _) then
            error "impossible, trace should always be initialized with { 0 => 0 }"
          case (Some (t1, a), None _) then gaussianSample a (sqrt (subf t t1))
          case (None _, Some (t2, b)) then gaussianSample b (sqrt (subf t2 t))
          case (Some (t1, a), Some (t2, b)) then
            let mu =
              (addf
                 a
                 (mulf (subf t t1) (divf (subf b a) (subf t2 t1))))
            in
            let sigma =
              sqrt (divf (mulf (subf t2 t) (subf t t1)) (subf t2 t1))
            in
            gaussianSample mu sigma
          end
        in
        modref trace (mapInsert t w tr); w)
      t
      tr

-- Fast version of `wienerSample` that assumes the sampled Wiener realization is
-- evaluated at increasing positive times.
let wienerSampleUnsafe : () -> Float -> Float = lam.
  let state = ref (0., 0.) in
  lam t2.
    match deref state with (w, t1) in
    if eqf t2 t1 then w else
      let w = (gaussianSample w (sqrt (subf t2 t1))) in
      modref state (w, t2); w

-- Base interface
lang RuntimeDistBase
  syn Dist a =

  sem sample : all a. Dist a -> a

  sem expectation : all a. Dist a -> a

  sem logObserve : all a. Dist a -> (a -> Float)
end

-- Elementary distributions
lang RuntimeDistElementary = RuntimeDistBase
  syn Dist a =
  | DistGamma {shape : Float, scale : Float}
  | DistExponential {rate : Float}
  | DistPoisson {lambda : Float}
  | DistBinomial {n : Int, p : Float}
  | DistBernoulli {p : Float}
  | DistBeta {a : Float, b : Float}
  | DistGaussian {mu : Float, sigma : Float}
  | DistMultinomial {n : Int, p : [Float]}
  | DistCategorical {p : [Float]}
  | DistDirichlet {a : [Float]}
  | DistUniform {a : Float, b : Float}
  | DistUniformDiscrete {a : Int, b : Int}
  | DistWiener {cps : Bool, a : ()}
  | DistLomax {scale: Float, shape : Float}
  | DistBetabin {n:Int, a: Float, b: Float}
  | DistNegativeBinomial {n:Int, p: Float}

  sem sample =
  | DistGamma t -> unsafeCoerce (gammaSample t.shape t.scale)
  | DistExponential t -> unsafeCoerce (exponentialSample t.rate)
  | DistPoisson t -> unsafeCoerce (poissonSample t.lambda)
  | DistBinomial t -> unsafeCoerce (binomialSample t.p t.n)
  | DistBernoulli t -> unsafeCoerce (bernoulliSample t.p)
  | DistBeta t -> unsafeCoerce (betaSample t.a t.b)
  | DistGaussian t -> unsafeCoerce (gaussianSample t.mu t.sigma)
  | DistMultinomial t -> unsafeCoerce (multinomialSample t.p t.n)
  | DistCategorical t -> unsafeCoerce (categoricalSample t.p)
  | DistDirichlet t -> unsafeCoerce (dirichletSample t.a)
  | DistUniform t -> unsafeCoerce (uniformContinuousSample t.a t.b)
  | DistUniformDiscrete t -> unsafeCoerce (uniformDiscreteSample t.a t.b)
  | DistWiener {cps = false, a = a} -> unsafeCoerce (wienerSample a)
  | DistWiener {cps = true, a = a} ->
    unsafeCoerce (let w = wienerSample a in lam k. lam x. k (w x))
  | DistLomax t -> unsafeCoerce (lomaxSample t.shape t.scale)
  | DistBetabin t -> unsafeCoerce (betabinSample t.n t.a t.b)
  | DistNegativeBinomial t -> unsafeCoerce (negativeBinomialSample t.n t.p)

  -- Expectation of primitive distributions over real values
  sem expectation =
  | DistGamma t -> unsafeCoerce (mulf t.shape t.scale)
  | DistExponential t -> unsafeCoerce (divf 1. t.rate)
  | DistPoisson t -> unsafeCoerce t.lambda
  | DistBinomial t -> unsafeCoerce (mulf t.p (int2float t.n))
  | DistBernoulli t -> unsafeCoerce t.p
  | DistBeta t -> unsafeCoerce (divf t.a (addf t.a t.b))
  | DistGaussian t -> unsafeCoerce t.mu
  | DistMultinomial t ->
    error "expectation undefined for the multinomial distribution"
  | DistCategorical t ->
    error "expectation undefined for the categorical distribution"
  | DistDirichlet t ->
    error "expectation undefined for the Dirichlet distribution"
  | DistUniform t -> unsafeCoerce (divf (addf t.a t.b) 2.)
  | DistUniformDiscrete t -> unsafeCoerce (divf (int2float (addi t.a t.b)) 2.)
  | DistWiener _ -> error "expectation undefined for the Wiener process"

  sem logObserve =
  | DistGamma t -> unsafeCoerce (gammaLogPdf t.shape t.scale)
  | DistExponential t -> unsafeCoerce (exponentialLogPdf t.rate)
  | DistPoisson t -> unsafeCoerce (poissonLogPmf t.lambda)
  | DistBinomial t -> unsafeCoerce (binomialLogPmf t.p t.n)
  | DistBernoulli t -> unsafeCoerce (bernoulliLogPmf t.p)
  | DistBeta t -> unsafeCoerce (betaLogPdf t.a t.b)
  | DistGaussian t -> unsafeCoerce (gaussianLogPdf t.mu t.sigma)
  | DistMultinomial t ->
    unsafeCoerce (lam o.
      if eqi t.n (foldl1 addi o) then multinomialLogPmf t.p o
      else negf inf)
  | DistCategorical t -> unsafeCoerce (categoricalLogPmf t.p)
  | DistDirichlet t -> unsafeCoerce (dirichletLogPdf t.a)
  | DistUniform t -> unsafeCoerce (uniformContinuousLogPdf t.a t.b)
  | DistUniformDiscrete t -> unsafeCoerce (uniformDiscreteLogPdf t.a t.b)
  | DistWiener _ -> error "logObserve undefined for the Wiener process"
  | DistLomax t -> unsafeCoerce (lomaxLogPdf t.shape t.scale)
  | DistBetabin t -> unsafeCoerce (betabinLogPmf t.n t.a t.b)
  | DistNegativeBinomial t -> unsafeCoerce (negativeBinomialLogPmf t.n t.p)
end

-- Empirical distribution
lang RuntimeDistEmpirical = RuntimeDistBase
  syn Dist a =
  | DistEmpirical {
      cumulativeWeights : [Float],
      logWeights : [Float],
      samples : [a],
      degenerate : Bool,

      -- NOTE(dlunde,2022-10-18): Monte Carlo inference often returns futher
      -- information (not only the samples). This field is used for this.
      extra : EmpiricalExtra
    }

  -- TODO(dlunde,2022-10-18): Would probably be nice to specify this per
  -- inference algorithm (in their respective subfolders)
  syn EmpiricalExtra =
  | EmpNorm { normConst: Float }
  | EmpMCMC { acceptRate: Float }

  sem constructDistEmpiricalHelper =
  | samples ->
    match unzip samples with (logWeights, samples) in
    let extra = EmpNorm {normConst = 0.0} in
    constructDistEmpirical samples logWeights extra

  -- DistEmpirical should always be created via this function
  sem constructDistEmpirical samples logWeights =
  | extra ->

    -- print "INPUT: "; printLn (strJoin "," (map float2string logWeights));

    let samples = toRope samples in
    let logWeights = toRope logWeights in

    -- Compute LSE
    let maxLogWeight = foldl (lam acc. lam lw. if geqf lw acc then lw else acc)
                         (negf inf) logWeights in

    let degenerate = eqf maxLogWeight (negf inf) in

    -- print "MAX: "; printLn (float2string maxLogWeight);
    let lse =
      addf maxLogWeight
        (log (foldl (lam acc. lam lw. addf acc (exp (subf lw maxLogWeight)))
                0. logWeights))
    in
    -- print "LSE: "; printLn (float2string lse);

    -- Normalize logweights (so that the new LSE = 0)
    let logWeights = map (lam lw. subf lw lse) logWeights in
    -- print "LOG-NORMALIZED: "; printLn (strJoin "," (map float2string logWeights));

    -- Computes cumulative (non-log) weights
    let f = lam acc. lam x.
      let acc = addf acc (exp x) in
      (acc, acc)
    in
    match mapAccumL f 0.0 logWeights with (_, cumulativeWeights) in
    -- print "CUMULATIVE NORMALIZED: "; printLn (strJoin "," (map float2string cumulativeWeights));

    -- NOTE(larshum, 2023-10-24): Apply a map to flatten the cumulative
    -- weights, to ensure accessing it (when sampling) is a constant-time
    -- operation.
    map (lam. ()) cumulativeWeights;

    DistEmpirical {
      cumulativeWeights = cumulativeWeights, logWeights = logWeights,
      degenerate = degenerate, samples = samples, extra = extra
    }

  sem empiricalSamples =
  | DistEmpirical t -> (t.samples, t.logWeights)
  | _ -> ([], [])

  sem empiricalDegenerate =
  | DistEmpirical t -> t.degenerate
  | _ -> false

  sem empiricalNormConst =
  | DistEmpirical t ->
    match t.extra with EmpNorm { normConst = normConst } then normConst else nan
  | _ -> nan

  sem empiricalAcceptRate =
  | DistEmpirical t ->
    match t.extra with EmpMCMC { acceptRate = acceptRate }
    then acceptRate else nan
  | _ -> nan

  sem sample =
  | DistEmpirical t ->
    -- NOTE(larshum, 2023-05-03): Sample by choosing a value in range [0, y)
    -- and finding the index of the maximal cumulative weight which is less
    -- than the chosen value. The sample at this index is returned. Note that
    -- we use the final cumulative weight rather than 1 for y, to prevent
    -- rounding errors from skewing the probability of the last sample.
    let x = uniformContinuousSample 0.0 (last t.cumulativeWeights) in

    -- NOTE(larshum, 2023-05-03): The sampled value 'x' cannot be greater than
    -- the last element of the sequence of cumulative weights, so we should
    -- always be able to find an index for which the comparison function
    -- returns a non-negative number.
    let cmp = lam y. if ltf (subf y x) 0.0 then negi 1 else 0 in
    match lowerBoundBinarySearch cmp t.cumulativeWeights with Some idx then
      unsafeCoerce (get t.samples idx)
    else
      error "Sampling from empirical distribution failed"

  sem expectationEmpiricalFloat logWeights =| samples ->
    let weights = map exp logWeights in
    foldl addf 0. (zipWith mulf weights (unsafeCoerce samples))

  sem expectation =
  | DistEmpirical t ->
    -- NOTE(oerikss, 2024-09-13): We assume that samples are floats. The
    -- type-system should reject expectation of distributions over other types.
    unsafeCoerce (expectationEmpiricalFloat t.logWeights t.samples)

  sem logObserve =
  -- TODO(dlunde,2022-10-18): Implement this?
  | DistEmpirical t -> error "Log observe not supported for empirical distribution"

  -- Creates an empirical distribution with a subsample of size n
  -- NOTE(vsenderov, 2024-03-28):
  -- The log-weights are not normalized here. First, normalization in the
  -- returned sub-sample from the program is not strictly desirable.
  -- In addition, as of now, constructDistEmpirical will normalize the weights,
  -- anyway.  For reference, with comments I show how to normalize the weights.
  sem constructDistEmpiricalSubsample n samples logWeights =
  | extra ->
    let d = constructDistEmpirical samples logWeights extra in
    let s = create n (lam. sample d) in
    -- without normalizing the sample weights
    let l = make n 0.0 in
    -- with normalizing
    --let w = negf (log (int2float n)) in
    --let l = make n w in
    constructDistEmpirical s l extra
end

lang RuntimeDist =
  RuntimeDistElementary +
  RuntimeDistEmpirical
end

type Dist a = use RuntimeDist in Dist a

-- We include the below definitions to produce non-mangled functions, which we
-- can refer to in the runtime handler without hard-coding the mangled prefix.
let distEmpiricalSamples : all a. use RuntimeDist in Dist a -> ([a], [Float]) =
  use RuntimeDist in
  empiricalSamples

let distEmpiricalDegenerate : all a. use RuntimeDist in Dist a -> Bool =
  use RuntimeDist in
  empiricalDegenerate

let distEmpiricalNormConst : all a. use RuntimeDist in Dist a -> Float =
  use RuntimeDist in
  empiricalNormConst

let distEmpiricalAcceptRate : all a. use RuntimeDist in Dist a -> Float =
  use RuntimeDist in
  empiricalAcceptRate

let sample : all a. use RuntimeDist in Dist a -> a =
  use RuntimeDist in
  sample

let expectation : all a. use RuntimeDist in Dist a -> a =
  use RuntimeDist in
  expectation

let logObserve : all a. use RuntimeDist in Dist a -> a -> Float =
  use RuntimeDist in
  logObserve

mexpr

let randomTimes = create 10 (lam. gaussianSample 0. 1.) in
let w = wienerSample () in
utest map w randomTimes with map w randomTimes in
utest map w randomTimes with map w randomTimes in
utest map w randomTimes with map w randomTimes in

()
