include "../../../models/diversification-models/crbd-synthetic.mc"
include "../../cppl-test.mc"
include "../../test.mc"

include "seq.mc"
include "sys.mc"
include "string.mc"
include "common.mc"
include "stats.mc"

-- Utilitary

recursive
let fact: Int -> Int = lam num.
  if eqi num 1 then
    1
  else
    muli num (fact (subi num 1))
end

-- These fonction come from

let crbd_lnS: Float -> Float -> Float -> Float -> Float = lam time. lam lambda. lam mu. lam rho.
  let r = subf lambda mu in
  let eRT = exp (mulf (subf 0.0 r) time) in
  let lnNum = log r in
  let lnDenom = log (subf lambda (mulf (subf lambda (divf r rho)) eRT))  in
  subf lnNum lnDenom

let crbd_survivorshipBias: Tree -> Float -> Float -> Float -> Float = lam tree. lam lambda. lam mu. lam rho.
  match tree with Node root in
    mulf -2.0 (crbd_lnS root.age lambda mu rho)

let crbd_lnGhat : Float -> Float -> Float -> Float -> Float = lam time. lam lambda. lam mu. lam rho.
  let r = subf lambda mu in
  let eRT = exp (mulf (subf 0.0 r) time) in
  let f = log (subf lambda (mulf (subf lambda (divf r rho)) eRT)) in
  subf (mulf (subf 0.0 r) time) (mulf 2.0 f)

recursive
let crbd_lnLike: Tree -> Float -> Float -> Float -> Float = lam tree. lam lambda. lam mu. lam rho.
  match tree with Leaf leaf then
    0.0 
  else match tree with Node root in
    let lnLikeLeft = crbd_lnLike root.left lambda mu rho in
    let lnLikeRight = crbd_lnLike root.right lambda mu rho in
    addf (addf lnLikeLeft lnLikeRight) (crbd_lnGhat root.age lambda mu rho)
end

let exactcrbdLikelihood: Tree -> Float -> Float -> Float -> Float = lam tree. lam lambda. lam mu. lam rho.
  -- Compute correction factor from oriented to labelled unoriented trees
  let numLeaves = countLeaves tree in
  let corrFactor = subf (mulf (int2float (subi numLeaves 1)) (log 2.0)) (log (int2float (fact numLeaves))) in

  -- Non-stalked tree, unconditional likelihood
  let ln1 = mulf (int2float (subi numLeaves 2)) (log lambda) in
  let ln2 = mulf (int2float numLeaves) (log rho) in
  match tree with Node root in
    let ln3 = mulf 2.0 (crbd_lnGhat root.age lambda mu rho) in
  let ln4 = crbd_lnLike tree lambda mu rho in
  let ln5 = mulf (int2float (subi 0 numLeaves)) (crbd_lnGhat 0.0 lambda mu rho) in
  -- let ln6 = -2.0*crbd_lnS(tree.age, lambda, mu, rho) in
  let ln6 = 0.0 in
  
  addf (addf (addf (addf (addf (addf corrFactor ln1) ln2) ln3) ln4) ln5) ln6

mexpr

let test = crbd_survivorshipBias tree 0.2 0.1 1.0 in
let test1= exactcrbdLikelihood tree 0.2 0.1 1.0 in
let test2= addf (crbd_survivorshipBias tree 0.2 0.1 1.0) (exactcrbdLikelihood tree 0.2 0.1 1.0) in

print "\n";
print (float2string test);
print "\n";
print (float2string test1);
print "\n";
print (float2string test2);

let en = eqCrbdSynthetic 0.2 2. in
let e = eqCrbdSyntheticMean 0.2 in
let rhs = crbdSyntheticTruth in
let r = resCrbdSynthetic in
let c = cpplResOfDist float2string in

utest r (c 0   (infer (Default {}) model))                                              with rhs using en in

()

