let model: ([Float], [[Float]], [[Float]]) -> () -> Float = lam params. lam.
  let y = assume (Gamma params.0 params.1) in
  iter (lam obs. observe obs (Poisson y)) params.2;
  y