/-
- MExprPPLStaticDelayedANF: normalizes the terms; however, apply special conditions to terms that create the vertices 
- The tests for --static-delay can be run via `make -s -f test-coreppl.mk static-delay`
-/

include "digraph.mc"
include "coreppl.mc"
include "dist.mc"
include "ext/math-ext.mc"
include "mexpr/shallow-patterns.mc"

let debug = false
-- This part is for the ANF transformation for not to ANF distributions with assume and observe but only the parameters
lang MExprPPLStaticDelayedANF = MExprPPL + MExprANFAll
  sem normalize (k:Expr -> Expr) =
  | TmAssume ({ dist = TmDist ({ dist = dist } & td) } & t) ->
    normalizeDist
      (lam dist. k (TmAssume { t with dist = TmDist { td with dist = dist } }))
      dist
  | TmObserve ({ value = value, dist = TmDist ({ dist = dist } & td) } & t) ->
    normalizeName
      (lam value.
        normalizeDist
          (lam dist.
             k (TmObserve {{ t with value = value }
                               with dist = TmDist { td with dist = dist}}))
          dist)
      value
  | TmApp ({lhs=TmApp ({lhs=TmConst ({val=CGet ()}&c),rhs=seq}&a2),rhs=ind}&a1) ->
    normalizeName
      (lam seq.
        normalizeName
          (lam ind.
             k (TmApp {{a1 with lhs= TmApp {{a2 with lhs=TmConst c} with rhs=seq}} with rhs=ind}))
          ind)
      seq
  | TmApp ({lhs=TmApp ({lhs=TmConst ({val=CCreate ()}&c),rhs=rep}&a2),rhs=TmLam l}&a1) ->
    k (TmApp a1)
  | TmApp ({lhs=TmApp ({lhs=TmConst ({val=CIter ()}&c),rhs=TmLam l}&a2),rhs=lst}&a1) ->
    normalizeName
      (lam lst.
             k (TmApp {{a1 with lhs= TmApp {{a2 with lhs=TmConst c} with rhs=
             TmLam {l with body=normalizeTerm l.body} }} with rhs=lst}))
      lst
  | TmApp ({lhs=TmApp ({lhs=TmConst ({val=CIteri ()}&c),rhs=TmLam ({body=TmLam l2}&l1)}&a2),rhs=lst}&a1) ->
    normalizeName
      (lam lst.
             k (TmApp {{a1 with lhs= TmApp {{a2 with lhs=TmConst c} with 
              rhs=TmLam {l1 with body=TmLam {l2 with body=normalizeTerm l2.body}}}} with rhs=lst}))
      lst
end

lang PBNGraph = MExprAst + MExprPPL 

  type Label = Int
  -- m: a mapping from a vertex ident to a corresponding vertex
  type PBN = {
    g:Digraph Vertex Label, -- a graph to keep the dependencies
    m:Map Name Vertex, -- a mapping from vertex id to vertex 
    targets:[Name] -- targets to be sampled at runtime
  }
  -- different types of vertices
  syn Vertex =
  | RandomVarNode {ident:Name,
                    val:Option Expr,
                    state:Ref Int, -- 0:blue (assume), 1:red (stable)
                    dist:Ref Dist,
                    mDist:Ref (Option Dist),
                    listId:Ref (Option Name),
                    plateId:Ref (Option Name)} --if it belongs to a plate (ref to most inner plate)
  | CodeBlockNode {ident:Name,
                    code:Expr,
                    plateId:Ref (Option Name),
                    ret:Bool} --if it belongs to a plate
  | ListNode {ident:Name,
              items:[Vertex],
              dist:Dist,
              create:Option Expr,
              lamParam:Ref (Option Name),
              outParam:Ref (Option Name),
              plateId:Ref (Option Name)}  --if it belongs to a plate, nested lists are not allowed
  | MultiplexerNode {ident:Name,
                      indexId:Name,
                      list:Vertex,
                      mDist:Ref (Option Dist),
                      plateId:Ref (Option Name)} --if it belongs to a plate
  | PlateNode {ident:Name,
               varIds:[Name], -- new variables introduced
               iterlId:Name, -- name of the observations to iterate over
               plateId:Ref (Option Name)}--if it belongs to a plate
  | FoldNode {ident:Name,
               varIds:[Name], -- variables introduced got from plate 
               lamAccId:Name, -- new lam acc id for updated parameters
               --margPIdM:Ref (Map Name Name), -- keep track of the marginalized parameters of the outside values, normally we would acces the params by their distribution but we cannot do that since we assign fold's return value as the marginlaized parameter to the outside dist 
               vToIndex:Ref (Map Name Int), -- parameters tag in the accumulated record
               iterlId:Name, -- name of the observations to iterate over
               plateId:Ref (Option Name),
               retBlockId:Name,
               accId:Name}

  -- to print out the vertex types
  sem v2str: Vertex -> String
  sem v2str =
  | RandomVarNode v -> let id = v.ident in 
                      let plateStr = match deref v.plateId with Some id then join ["(", id.0, ", ",(int2string (sym2hash id.1)), ")"] else "-" in
                       join ["\nRandomVarNode ident: (", id.0, ", ",(int2string (sym2hash id.1)), ")",
                       "\n state:", if eqi (deref v.state) 0 then "initialized" else if eqi (deref v.state) 1 then "marginalized" else "stabilized"
                       ,"\n dist " , (mexprToString (dist_ (deref v.dist)))
                       ,"\n val " , match v.val with Some val then mexprToString val else "-"
                       ,"\nplateId: ", plateStr]
  | CodeBlockNode v -> let id = v.ident in 
                      let plateStr =match deref v.plateId with Some id then join ["(", id.0, ", ",(int2string (sym2hash id.1)), ")"] else "-" in                       let ret = if v.ret then " true" else " false" in
                      join ["\nCodeBlockNode ident: (", id.0, ", ",(int2string (sym2hash id.1)), ")",
                           "\nCode:",expr2str v.code,
                            "\nIsRet:",ret,
                            "\nplateId: ", plateStr,"\n"]
  | ListNode v -> let id = v.ident in
                  let plateStr =match deref v.plateId with Some id then join ["(", id.0, ", ",(int2string (sym2hash id.1)), ")"] else "-" in
                  join ["\nListNode ident: (", id.0, ", ",(int2string (sym2hash id.1)), ")",
                  "\nplateId: ", plateStr,
                  foldl (lam acc. lam v. let v = getId v in join [acc, " ", v.0 ,":",(int2string (sym2hash v.1)),"\t"]) "" v.items]
  | MultiplexerNode v -> let id = v.ident in
                  let plateStr =match deref v.plateId with Some id then join ["(", id.0, ", ",(int2string (sym2hash id.1)), ")"] else "-" in
                        join ["\nMultiplexerNode ident: (", id.0, ", ",(int2string (sym2hash id.1)), ")",
                        "\nplateId: ", plateStr]
  | PlateNode v -> let id = v.ident in
                  let plateStr = match deref v.plateId with Some id then join ["(", id.0, ", ",(int2string (sym2hash id.1)), ")"] else "-" in
                        join ["\nPlateNode ident: (", id.0, ", ",(int2string (sym2hash id.1)), ")",
                              "\nplateId: ", plateStr]
  | FoldNode v -> let id = v.ident in
                  let plateStr = match deref v.plateId with Some id then join ["(", id.0, ", ",(int2string (sym2hash id.1)), ")"] else "-" in

                        join ["\nFoldNode ident: (", id.0, ", ",(int2string (sym2hash id.1)), ")"]

  sem printVertices g =
  | v2str -> iter (lam v. print (join [(v2str v), "\n"])) (digraphVertices g)

  sem cmprVertex: Vertex -> Vertex -> Int
  sem cmprVertex v1 =
  | v2 -> cmprVertexH (v1,v2)

  sem cmprVertexH: (Vertex,Vertex) -> Int
  sem cmprVertexH =
  | (RandomVarNode t1, RandomVarNode t2) -> nameCmp t1.ident t2.ident
  | (CodeBlockNode t1, CodeBlockNode t2) -> nameCmp t1.ident t2.ident
  | (ListNode t1, ListNode t2) -> nameCmp t1.ident t2.ident
  | (MultiplexerNode t1, MultiplexerNode t2) -> nameCmp t1.ident t2.ident
  | (PlateNode t1, PlateNode t2) -> nameCmp t1.ident t2.ident
  | (PlateNode t1, FoldNode t2) -> nameCmp t1.ident t2.ident
  | (FoldNode t1, PlateNode t2) -> nameCmp t1.ident t2.ident
  | (FoldNode t1, FoldNode t2) -> nameCmp t1.ident t2.ident
  | (t1,t2) -> subi (constructorTag t1) (constructorTag t2)

  sem cmprEdge: (Vertex,Vertex,Label) -> (Vertex,Vertex,Label) -> Int
  sem cmprEdge e1 =
  | (v1, v2, _) -> let cmprV1 = cmprVertex e1.0 v1 in
                   if eqi cmprV1 0 then cmprVertex e1.1 v2 else cmprV1

  sem getId: Vertex -> Name
  sem getId =
  | (RandomVarNode {ident=id} | CodeBlockNode {ident=id} | MultiplexerNode {ident=id} |
     ListNode {ident=id} | PlateNode {ident=id} | FoldNode {ident=id})-> id

  sem getPlateId: Vertex -> Option Name
  sem getPlateId =
  | (RandomVarNode {plateId=pid} | CodeBlockNode {plateId=pid} | MultiplexerNode {plateId=pid} |
     ListNode {plateId=pid}| PlateNode {plateId=pid}| FoldNode {plateId=pid})-> deref pid
  
  sem getListId: Vertex -> Option Name
  sem getListId =
  | RandomVarNode l -> deref l.listId
  | _ -> None ()

  sem getDist =
  | RandomVarNode l -> deref l.dist

  -- remove a vertex from PBN
  sem removeVertexPBN: PBN -> Vertex -> PBN
  sem removeVertexPBN pbn =
  | v -> let g = digraphRemoveVertex v pbn.g in
    let m = mapRemove (getId v) pbn.m in
    let pbn = {{pbn with m=m} with g=g} in
    pbn
    
  -- add a vertex to a PBN
  sem addVertexPBN: PBN -> Vertex -> PBN
  sem addVertexPBN pbn =
  | v -> 
  -- here is the error
  let g = digraphAddUpdateVertex v pbn.g in
    let m = mapInsert (getId v) v pbn.m in
    {{pbn with g=g} with m=m}

end

-- a language fragment for posterior and posterior predictive functions
lang ConjugatePrior = CorePPL + MExprAst + MExprPPL + PBNGraph

  -- (d1:likelihood,d2:prior) checks whether d1 and d2 are conjugate
  sem isConjugatePrior: Name -> (Dist,Dist) -> Bool
  sem isConjugatePrior pid =
  | (DBernoulli _, DBeta _) -> true
  | (DGaussian d1, DGaussian _) -> match d1.mu with TmVar v in nameEq v.ident pid
  | (DCategorical _, DDirichlet _) -> true
  | _ -> false

  -- check if two distributions family is equivalent
  sem eqFamilyDist: (Dist,Dist) -> Bool 
  sem eqFamilyDist =
  |(d1, d2) -> eqi (constructorTag d1) (constructorTag d2)

  -- check whether a list consists of rvs with same distribution family
  sem validList: PBN -> [Option Vertex] -> Bool
  sem validList pbn =
  | [Some (RandomVarNode r)] ++ as ->
    match deref r.listId with Some _ then false --already in another list
    else validListH pbn (deref r.dist) as
  | [t] ++ as -> false

  sem validListH: PBN -> Dist -> [Option Vertex] -> Bool
  sem validListH pbn dst =
  | [Some (RandomVarNode r)] ++ as ->
    match deref r.listId with Some _ then false
    else if eqFamilyDist (deref r.dist,dst) then validListH pbn dst as else false
  | [t] ++ as -> false
  | [] -> true

  sem getParams: Dist -> Expr
  sem getParams =
  | DBernoulli d -> utuple_ [d.p]
  | DBeta d -> utuple_ [d.a,d.b]
  | DGaussian d -> utuple_ [d.mu, d.sigma]
  | DCategorical d -> utuple_ [d.p]
  | DDirichlet d -> utuple_ [d.a]
  
  sem changeParams: Name -> Dist -> Dist
  sem changeParams param =
  | DBernoulli d -> DBernoulli {d with p=tupleproj_ 0 (nvar_ param)}
  | DBeta d -> DBeta {{d with a=tupleproj_ 0 (nvar_ param) } with b=tupleproj_ 1 (nvar_ param) }
  | DGaussian d -> DGaussian {{d with mu=tupleproj_ 0 (nvar_ param)} with sigma=tupleproj_ 1 (nvar_ param)}
  | DCategorical d -> DCategorical {d with p=tupleproj_ 0 (nvar_ param)}
  | DDirichlet d -> DDirichlet {d with a=tupleproj_ 0 (nvar_ param)}
  
  -- given the likelihood, the prior and the observartion calculates the posterior
  -- (d1: likelihood, d2: prior)
  sem posterior: Option Expr -> Option (Name,Expr) -> Option Name -> (Dist,Dist) -> (Vertex,Dist,[Name])
  sem posterior obs indices plateId =
  | (DBernoulli d1,DBeta d2) ->
    let val = match obs with Some val then val else never in
    let aName = nameSym "postA" in
    let bName = nameSym "postB" in
    let postAlpha = nulet_ aName (if_ val (addf_ d2.a (float_ 1.)) d2.a) in
    let postBeta = nulet_ bName (if_ val d2.b (addf_ d2.b (float_ 1.))) in
    let code = match indices with Some (mInd, lInd) then
      let eqN =(nameSym "eq") in
       bindall_ [nulet_ eqN ((eqi_ (nvar_ mInd) lInd)),
                nulet_ aName (if_ (nvar_ eqN) (if_ val (addf_ d2.a (float_ 1.)) d2.a) d2.a),
                nulet_ bName (if_ (nvar_ eqN) (if_ val d2.b (addf_ d2.b (float_ 1.))) d2.b)
                ]
      else ((bind_ postAlpha postBeta)) in
    let tName = nameSym "paramR" in
    let rho = CodeBlockNode {ident=tName, code=code, ret=false, plateId=ref plateId} in
    let paramNames = [aName,bName] in
    (rho, DBeta {{d2 with a=(nvar_ aName)} with b=(nvar_ bName)},paramNames)
  | (DGaussian d1, DGaussian d2) ->
    let muName = nameSym "postMu" in
    let sigmaName = nameSym "postSigma" in
    let val = match obs with Some val then val else never in
    let s02 = (mulf_ d2.sigma d2.sigma) in
    let s2 = (mulf_ d1.sigma d1.sigma) in
    let muRHS = addf_ (divf_ d2.mu s02) (divf_ val s2) in
    let muLHS = divf_ (float_ 1.0) (addf_ (divf_ (float_ 1.0) s02) (divf_ (float_ 1.0) s2)) in
    let postMu = nulet_ muName (mulf_ muRHS muLHS) in
    let sigma = divf_ (float_ 1.0) (addf_ (divf_ (float_ 1.0) s02) (divf_ (float_ 1.0) s2)) in
    let postSigma = nulet_ sigmaName (appf1_ (var_ "externalSqrt") sigma) in
    let code = match indices with Some (mInd, lInd) then
      let eqN =(nameSym "eq") in
       bindall_ [nulet_ eqN ((eqi_ (nvar_ mInd) lInd)),
                nulet_ muName ( if_ (nvar_ eqN) (mulf_ muRHS muLHS) d2.mu),
                nulet_ sigmaName (if_ (nvar_ eqN) (appf1_ (var_ "externalSqrt") sigma) d2.sigma)
                ]
      else (bind_ postMu postSigma) in
    let tName = nameSym "paramR" in
    let rho = CodeBlockNode {ident=tName, code=code, ret=false,plateId=ref plateId} in
    let paramNames = [muName,sigmaName] in
    (rho, DGaussian {{d2 with mu= nvar_ muName} with sigma= nvar_ sigmaName},paramNames)
  | (DCategorical d1, DDirichlet d2) ->
    let val = match obs with Some val then val else never in
    let aName = nameSym "postA" in
    let postA = nulet_ aName (mapi_ ( ulam_ "i" (ulam_ "e" (if_ (eqi_ (var_ "i") val) (addf_ (var_ "e") (float_ 1.0)) (var_ "e")))) d2.a) in
    let code = match indices with Some (mInd, lInd) then
      let eqN =(nameSym "eq") in
       bindall_ [nulet_ eqN ((eqi_ (nvar_ mInd) lInd)),
                nulet_ aName (if_ (nvar_ eqN) (mapi_ ( ulam_ "i" (ulam_ "e" (if_ (eqi_ (var_ "i") val) (addf_ (var_ "e") (float_ 1.0)) (var_ "e")))) d2.a) d2.a)  ]
      else postA in
    let tName = nameSym "paramR" in
    let rho = CodeBlockNode {ident=tName, code=code, ret=false,plateId=ref plateId} in
    let paramNames = [aName] in
    (rho, DDirichlet {d2 with a=nvar_ aName},paramNames)
  | _ -> error "posterior:not supported"

  -- input (d1: likelihood, d2: prior)
  -- output (rho:Vertex, q:Expr)
  sem posteriorPredictive: Option Name -> (Dist,Dist) -> Option (Vertex,Dist)
  sem posteriorPredictive plateId =
  | (DBernoulli d1, DBeta d2) ->
    let postP = divf_ d2.a (addf_ d2.a d2.b) in
    let tName = nameSym "param" in
    let pName = nameSym "margP" in
    let letT = nulet_ pName postP in
    let rho = CodeBlockNode {ident=tName, code=letT, ret=false,plateId=ref plateId} in
    Some (rho, DBernoulli {d1 with p=nvar_ pName})

  | (DGaussian d1,DGaussian d2) ->
    let s02 = (mulf_ d2.sigma d2.sigma) in
    let s2 = (mulf_ d1.sigma d1.sigma) in
    let postMu = mulf_ s02 (divf_ d2.mu s02) in
    let postSigma = appf1_ (var_ "externalSqrt") (addf_ s02 s2) in
    let tName = nameSym "param" in
    let mName = nameSym "margMu" in
    let sName = nameSym "margSigma" in
    let letT = bind_ (nulet_ mName postMu) (nulet_ sName postSigma) in
    let rho = CodeBlockNode {ident=tName, code=letT, ret=false,plateId=ref plateId} in
    Some (rho, DGaussian {{d1 with mu=nvar_ mName} with sigma=nvar_ sName})

  | (DCategorical d1,DDirichlet d2) ->
    let sumName = nameSym "sum" in
    let sumai = nulet_ sumName (foldl_ (ulam_ "acc" (ulam_ "i" (addf_ (var_ "acc") (var_ "i")))) (float_ 0.0) (d2.a)) in
    let postP = map_ (ulam_ "ai" (divf_ (var_ "ai") (nvar_ sumName))) d2.a in
    let tName = nameSym "param" in
    let pName = nameSym "margP" in
    let letT = nulet_ pName postP in
    let rho = CodeBlockNode {ident=tName, code=bind_ sumai letT, ret=false,plateId=ref plateId} in
    Some (rho, DCategorical {d1 with p=nvar_ pName})
  | _ -> None ()

  sem posteriorPredictiveL (plateId:Option Name) (indexId:Name) =
  | ([DBernoulli d]++as)&mDists -> 
    let lName = nameSym "" in
    let pName = nameSym "margP" in
    let params = map getParams mDists in
    let lLet = nulet_ lName (seq_ params) in
    let pLet = nulet_ pName (tupleproj_ 0 (get_ (nvar_ lName) (nvar_ indexId))) in
    let code = bindall_ [lLet,pLet] in
    let rho = CodeBlockNode {ident=nameSym "",code=code,ret=false,plateId=ref plateId} in
    Some (rho, DBernoulli {d with p=nvar_ pName})
  | ([DGaussian d]++as)&mDists -> 
    let lName = nameSym "" in
    let muName = nameSym "margMu" in
    let sigmaName = nameSym "margSigma" in
    let params = map getParams mDists in
    let lLet = nulet_ lName (seq_ params) in
    let muLet = nulet_ muName (tupleproj_ 0 (get_ (nvar_ lName) (nvar_ indexId))) in
    let sigmaLet = nulet_ sigmaName (tupleproj_ 1 (get_ (nvar_ lName) (nvar_ indexId))) in
    let code = bindall_ [lLet,muLet,sigmaLet] in
    let rho = CodeBlockNode {ident=nameSym "",code=code,ret=false,plateId=ref plateId} in
    Some (rho, DGaussian {{d with mu=nvar_ muName} with sigma=nvar_ sigmaName})
  | ([DCategorical d]++as)&mDists -> 
    let lName = nameSym "" in
    let pName = nameSym "margP" in
    let params = map getParams mDists in
    let lLet = nulet_ lName (seq_ params) in
    let pLet = nulet_ pName (tupleproj_ 0 (get_ (nvar_ lName) (nvar_ indexId))) in
    let code = bindall_ [lLet,pLet] in
    let rho = CodeBlockNode {ident=nameSym "",code=code,ret=false,plateId=ref plateId} in
    Some (rho, DCategorical {d with p=nvar_ pName})
  | _ -> None ()
end

-- First step of the algorithm (Static PBN constructor): create/construct a PBN from a given program
lang CreatePBN = ConjugatePrior

  -- m: a mapping from a variable name to its corresponding vertex id. Several let bindings can corresspond to a single code block vertex
  type CreateAcc = {
    m:Map Name Name,
    blockIdent:Option Name,
    vertexId:Option Name,
    plateId:Option Name,
    isRet:Bool
  }

  sem emptyCreateAcc: () -> CreateAcc
  sem emptyCreateAcc =
  | _ -> { m=mapEmpty nameCmp, blockIdent=(None ()), plateId=None (), vertexId=None (), isRet=false }

  sem createM : Expr -> PBN
  sem createM =
  | prog -> match createPBN {g=digraphEmpty cmprVertex eqi,targets=[],m=mapEmpty nameCmp} (emptyCreateAcc ()) prog with (pbn,_) in pbn

  -- create edges based on the dependencies of vertex v
  sem createEdges: Vertex -> PBN -> CreateAcc -> Set (Vertex,Vertex,Label) -> Expr -> Set (Vertex,Vertex,Label)
  sem createEdges v pbn cAcc edges =
  | TmVar t ->
    -- find the corresponding vertex ident from the variable ident
    match mapLookup t.ident cAcc.m with Some vertexId then
      let vFrom:Vertex = mapLookupOrElse (lam. error "createEdges:Lookup failed") vertexId pbn.m in
      -- create an edge to the source vertex from the vertex that it depends on
      if digraphEqv pbn.g vFrom v then edges --check if they are in the same codeblock if so no need to create an edge
      else setInsert (vFrom, v, 0) edges
    else edges -- if cannot find id then it must be created with lambda scoping so ignore
  | t -> sfold_Expr_Expr (createEdges v pbn cAcc) edges t

    -- finds the random variable identities within an expression
  sem findTargetRVs: PBN -> Expr -> PBN
  sem findTargetRVs pbn =
  | TmVar t -> 
    match mapLookup t.ident pbn.m with Some v then
      match v with (RandomVarNode _) then {pbn with targets=cons t.ident pbn.targets} else
      -- if a multiplexer is returned, every item becomes target
      match v with (MultiplexerNode v) then
        match v.list with ListNode l in
        {pbn with targets=foldl (lam ids. lam r. cons (getId r) ids) pbn.targets l.items}
      else pbn
    else pbn
  | t -> sfold_Expr_Expr findTargetRVs pbn t

  sem createCodeBlock: PBN -> CreateAcc -> Expr -> (Option Name, Option Name) -> (Vertex,Name)
  sem createCodeBlock pbn cAcc t =
  -- merge with a previously created block with ident 'bid'
  | (Some id, Some bid) -> match mapLookupOrElse (lam. error "createCodeBlock:Lookup failed") bid pbn.m with CodeBlockNode c in
                           let v = CodeBlockNode {c with code=bind_ c.code (nulet_ id t)} in
                           (v,bid)
  | (Some id, None ()) -> let v = CodeBlockNode {ident=id,code=(nulet_ id t),ret=false, plateId=ref cAcc.plateId} in (v,id)
  | _ -> let ident = nameSym "" in
      let isRet = match cAcc.plateId with Some _ then false else cAcc.isRet in
      let v = CodeBlockNode {ident=ident,code=t,ret=isRet, plateId=ref cAcc.plateId} in
      (v,ident)

  -- given vertex, its id for pbn.m and id for cAcc.m2 and expr for env
  sem addVertex: PBN -> CreateAcc -> (Vertex,Name) -> (PBN, CreateAcc)
  sem addVertex pbn cAcc =
  | (v,id2) ->
    let pbn = addVertexPBN pbn v in
    let m = mapInsert id2 (getId v) cAcc.m in
    (pbn, {{cAcc with m=m} with blockIdent=None ()})

  sem replaceInexpr =
  | TmType t -> TmType {t with inexpr=unit_}
  | TmRecLets t -> TmRecLets {t with inexpr=unit_}
  | TmExt t  -> TmExt {t with inexpr=unit_}
  | TmConDef t -> TmConDef {t with inexpr=unit_}

  sem createPBN: PBN -> CreateAcc -> Expr -> (PBN, CreateAcc)
  sem createPBN pbn cAcc =
  | TmLet t ->
    let res = createPBNH pbn {{cAcc with vertexId=(Some t.ident)} with isRet=false} t.body in
    createPBN res.0 res.1 t.inexpr
  -- all other with inexpr
  | (TmRecLets {inexpr=inexpr} | TmType {inexpr=inexpr} | TmExt {inexpr=inexpr} | TmConDef {inexpr=inexpr}) & t -> 
    let res = createPBNH pbn {{cAcc with isRet=false} with vertexId=None ()} (replaceInexpr t) in
    createPBN res.0 res.1 inexpr
  | t -> let res = createPBNH pbn {{cAcc with isRet=true} with vertexId=None ()} t in
         (res.0,res.1)

  sem createPBNH:PBN -> CreateAcc -> Expr -> (PBN, CreateAcc,Option Vertex)
  sem createPBNH pbn cAcc =
  | (TmAssume {dist=TmDist {dist=dist}} | TmObserve {dist=TmDist {dist=dist}}) & t ->
    -- get the ident if it comes from a let expression
    let id = match cAcc.vertexId with Some id then id else nameSym "rv" in
    -- if an observe then get its value
    let val = match t with TmObserve t then Some t.value else None () in
    -- create an initialized (state 0) vertex 
    let v = RandomVarNode {ident = id, val = val, state = ref 0, dist = ref dist, mDist = ref (None ()), plateId=ref cAcc.plateId, listId=ref (None ())} in
    -- add the vertex to the graph and to the context
    match addVertex pbn cAcc (v,id) with (pbn,cAcc) in
    -- if it is an observe, add it to the targets to be conditioned
    let targets = match t with TmObserve _ then cons id pbn.targets else pbn.targets in
    -- create edges to the created random variable node v from the nodes that it depends on
    let edges = setToSeq (createEdges v pbn cAcc (setEmpty cmprEdge) t) in
    let g = digraphMaybeAddEdges edges pbn.g in
    ({{pbn with targets=targets} with g=g},{cAcc with blockIdent=None()} ,Some v)
  | TmVar t -> if cAcc.isRet then createPBNGeneric pbn cAcc (TmVar t) else never -- aliases are removed
  | TmSeq t ->
    -- get the ident if it comes from a let expression
    let id = match cAcc.vertexId with Some id then id else nameSym "seq" in
    -- get the item vertices
    let items = map (lam v. match v with TmVar v in mapLookup v.ident pbn.m) t.tms in
    let v = if validList pbn items then
      let items = map (lam r. match r with Some (RandomVarNode r) in modref r.listId (Some id);RandomVarNode r) items in
      (ListNode {ident=id, items=items,plateId=ref cAcc.plateId,create=None (),dist=getDist (get items 0),lamParam=ref (None ()), outParam=ref (None ())},pbn)
    else -- if items are not valid, then this should be code block so no transformation will be performed
        let res = createCodeBlock pbn cAcc (TmSeq t) (cAcc.vertexId,cAcc.blockIdent) in
        let edges = setToSeq (createEdges res.0 pbn cAcc (setEmpty cmprEdge) (TmSeq t)) in
        let g = digraphMaybeAddEdges edges pbn.g in
        (res.0,{pbn with g=g})
    in
    match addVertex v.1 cAcc (v.0,id) with (pbn,cAcc) in
    let blockIdent = match v.0 with CodeBlockNode c then Some c.ident else None () in
    (pbn,{cAcc with blockIdent=blockIdent} ,Some v.0)
    | TmApp ({lhs=(TmApp ({lhs=TmConst ({val=CCreate()}&c),rhs=rep})&a1),
              rhs=TmLam ({body=TmAssume {dist=TmDist {dist=dist}}}&l)}&a2) ->
    let id = match cAcc.vertexId with Some id then id else nameSym "create" in
    let accH = {{{cAcc with blockIdent=None()} with vertexId=None()} with isRet=false} in
    match createPBNH pbn accH l.body with (pbn, cAcc, Some item) in
    let v = ListNode {ident=id,items=[item],plateId=ref cAcc.plateId,dist=getDist item,lamParam=ref (Some (nameSym "p")), outParam=ref (None ()),create=Some rep} in
    let pbn =addVertexPBN pbn item in
    match addVertex pbn cAcc (v,id) with (pbn,cAcc) in
    let edges = setToSeq (createEdges v pbn cAcc (setEmpty cmprEdge) (TmApp a2)) in
    let g = digraphMaybeAddEdges edges pbn.g in
    ({pbn with g=g},{cAcc with blockIdent=None ()}, Some v)
 | TmApp ({lhs=(TmApp ({lhs=TmConst ({val=CGet ()}&c),rhs=TmVar seq})&t2),rhs=TmVar ind}&a) ->
    let id = match cAcc.vertexId with Some id then id else nameSym "get" in
    let pbn = findTargetRVs pbn (TmVar ind) in
    let v =
      -- if there is no such list node created, create a codeblock
      match mapLookup seq.ident pbn.m with None () then
        (createCodeBlock pbn cAcc (TmApp a) (cAcc.vertexId,cAcc.blockIdent)).0
      else -- there is a list node created which consists of valid items
        let lst = mapLookupOrElse (lam. error "not found") seq.ident pbn.m in
        MultiplexerNode {ident=id,indexId=ind.ident,list=lst,plateId=ref cAcc.plateId,mDist=ref (None ())}
    in
    match addVertex pbn cAcc (v,id) with (pbn,cAcc) in
    let edges = setToSeq (createEdges v pbn cAcc (setEmpty cmprEdge) (TmApp a)) in
    let g = digraphMaybeAddEdges edges pbn.g in
    let blockIdent = match v with CodeBlockNode c then Some c.ident else None () in
    ({pbn with g=g}, {{cAcc with blockIdent=blockIdent} with vertexId=None ()},Some v)
  | TmApp ({lhs=(TmApp ({lhs=TmConst ({val=CIter()}&c),rhs=TmLam l})&a1),rhs=TmVar lst}&a2) ->
    createPlate pbn cAcc [l.ident] l.body lst.ident (TmApp a2)
  | TmApp ({lhs=(TmApp ({lhs=TmConst ({val=CIteri()}&c),rhs=TmLam ({body=TmLam l2}&l1)})&a1),rhs=TmVar lst}&a2) ->
    createPlate pbn cAcc [l1.ident,l2.ident] l2.body lst.ident (TmApp a2)
  | t -> createPBNGeneric pbn cAcc t

  sem createPlate: PBN -> CreateAcc -> [Name] -> Expr -> Name -> Expr -> (PBN,CreateAcc,Option Vertex)
  sem createPlate pbn cAcc idents body iterlId =
  | t ->
    let id = match cAcc.vertexId with Some id then id else never in
    match createPBN pbn {{{{cAcc with blockIdent=None()} with plateId=Some id} with vertexId=None ()} with isRet=false} body with (pbnB,cAccB) in
    let v = PlateNode {ident=id, varIds=idents,iterlId=iterlId, plateId=ref cAcc.plateId} in
    match addVertex pbnB cAcc (v,id) with (pbn,cAcc) in
    let edges = setToSeq (createEdges v pbn cAcc (setEmpty cmprEdge) (nvar_ iterlId)) in
    let g = digraphMaybeAddEdges edges pbn.g in
    ({pbn with g=g},{cAcc with blockIdent=None ()},Some v)

  sem createPBNGeneric: PBN -> CreateAcc -> Expr -> (PBN,CreateAcc,Option Vertex)
  sem createPBNGeneric pbn cAcc =
  | t ->
    let t = match t with TmVar v then
      match mapLookup v.ident pbn.m with Some (PlateNode p) then unit_ else t else t in
    let v = createCodeBlock pbn cAcc t (cAcc.vertexId,cAcc.blockIdent) in
    let id = match cAcc.vertexId with Some id then id else v.1 in
    match addVertex pbn cAcc (v.0,id) with (pbn,cAcc) in
    let edges = setToSeq (createEdges v.0 pbn cAcc (setEmpty cmprEdge) t) in
    let pbn = {pbn with g = digraphMaybeAddEdges edges pbn.g} in
    --let pbn = findTargetRVs pbn t in
    let blockIdent = match v.0 with CodeBlockNode c then Some c.ident else None () in
    (pbn,{cAcc with blockIdent=blockIdent}, Some v.0)

end

-- Third step of the algorithm (Program Reconstructor): reconstructing a probabilistic program from a PBN
lang RecreateProg = PBNGraph + MExprAst + MExprPPL

  sem extractPLVertices m v =
  | Some pid ->
      let vertices = mapLookupOrElse (lam. []) pid m in
      mapInsert pid (snoc vertices v) m
  | _ -> m

  sem recreate: PBN -> Expr
  sem recreate =
  | pbn ->
    let order = digraphTopologicalOrder pbn.g in
    let pbn = modifyGraph pbn in
    let order = digraphTopologicalOrder pbn.g in
    let vRet = filter (lam v.
      match v with CodeBlockNode c then
        let np = match deref c.plateId with Some _ then false else true in
        and c.ret np else false) order in
    let vRet = if eqi (length vRet) 0 then error "recreate:no return" else (get vRet 0) in
    let plateVertices = foldl (lam acc. lam v. extractPLVertices acc v (getPlateId v)) (mapEmpty nameCmp) order in
    let order = filter (lam v. match v with CodeBlockNode c then not c.ret else true) order in
    let order = filter (lam v. match getPlateId v with Some _ then false else true) order in
    let createItems = foldl (lam acc. lam v. match v with ListNode l then match l.create with Some _ then
        match get l.items 0 with item in setInsert (getId item) acc else acc else acc) (setEmpty nameCmp) order in
    let order = filter (lam v. match v with RandomVarNode v then not (setMem v.ident createItems) else true) order in
    recreateVertex plateVertices pbn (snoc order vRet)


  sem modifyGraph: PBN -> PBN
  sem modifyGraph =
  | pbn ->
    match pbn with {g=g,m=m,targets=targets} in
    let lists = filter (lam v. match v with ListNode _ then true else false) (digraphVertices g) in
    let g = foldl (lam g. lam l.
            match l with ListNode r in
            foldl (lam g:Digraph Vertex Label. lam i:Vertex.
                    let edges = digraphEdgesTo i g in
                    digraphMaybeAddEdges (map (lam e. (e.0,l,e.2)) edges) g) g r.items) g lists in
    let g = foldl (lam g. lam v. 
      match v with (PlateNode _ | FoldNode _) then g else
      match getPlateId v with Some pid then let pv = mapLookupOrElse (lam. error "not found") pid pbn.m in
      digraphMaybeAddEdge v pv 0 g else g) g (digraphVertices g) in 
    {pbn with g=g}

  sem recreateCode plateVertices pbn =
  | CodeBlockNode t -> t.code
  | RandomVarNode v -> let body = match v.val with Some val then
    TmObserve {dist=dist_ (deref v.dist), value=val,ty=tyunknown_, info = NoInfo ()}
    else TmAssume {dist=dist_ (deref v.dist), ty=tyunknown_, info = NoInfo ()} in
    nulet_ v.ident body
  | MultiplexerNode m -> nulet_ m.ident (get_ (nvar_ (getId m.list)) (nvar_ m.indexId))
  | ListNode l ->
    match l.create with Some val then 
      match (get l.items 0) with (RandomVarNode v)&rv in
      match deref l.lamParam with Some lamParam in
      match deref l.outParam with Some outparam then
      nulet_ l.ident (map_ (nulam_ lamParam (bind_ (recreateCode plateVertices pbn rv) (nvar_ v.ident))) (nvar_ outparam))
      else 
      nulet_ l.ident (create_ val (nulam_ lamParam (bind_ (recreateCode plateVertices pbn rv) (nvar_ v.ident))))
    else nulet_ l.ident (TmSeq {tms=(map (lam i. nvar_ (getId i)) l.items), ty=tyunknown_,info=NoInfo ()})
  | PlateNode p ->
    let vItems = mapLookupOrElse (lam. error "recreateCode:Lookup failed") p.ident plateVertices in
    let bdyIn = foldl (lam acc. lam v. bind_ acc (recreateCode plateVertices pbn v)) unit_ vItems in
    let body = if eqi (length p.varIds) 1 then (iter_ (nulam_ (get p.varIds 0) bdyIn) (nvar_ p.iterlId))
    else (iteri_ (nulam_ (get p.varIds 0) (nulam_ (get p.varIds 1) bdyIn)) (nvar_ p.iterlId)) in
    nulet_ p.ident body
  | FoldNode f ->
    let vItems = mapLookupOrElse (lam. error "recreateCode:Lookup failed") f.ident plateVertices in
    let vPRet = filter (lam v. match v with CodeBlockNode c then c.ret else false) vItems in   
    let vPRet = get vPRet 0 in
    let vItems = filter (lam v. match v with CodeBlockNode c then not c.ret else true) vItems in
    let bdyIn = foldl (lam acc. lam v. bind_ acc (recreateCode plateVertices pbn v)) unit_ (snoc vItems vPRet) in
    let bdy = foldl (lam acc. lam l. nulam_ l acc) bdyIn (snoc f.varIds f.lamAccId ) in
    nulet_ f.ident (foldl_ bdy (nvar_ f.accId) (nvar_ f.iterlId))
    

  sem recreateVertex: Map Name [Vertex] -> PBN -> [Vertex] -> Expr
  sem recreateVertex plateVertices pbn =
  | [(CodeBlockNode v)&t] ++ as -> let code = (recreateCode plateVertices pbn t) in
    if v.ret then code else bind_ code (recreateVertex plateVertices pbn as)
  | [(RandomVarNode _)&t] ++ as -> bind_ (recreateCode plateVertices pbn t) (recreateVertex plateVertices pbn as)
  | [(MultiplexerNode _)&t] ++ as -> bind_ (recreateCode plateVertices pbn t) (recreateVertex plateVertices pbn as)
  | [(PlateNode _) & t] ++ as -> bind_ (recreateCode plateVertices pbn t) (recreateVertex plateVertices pbn as)
  | [(FoldNode _) & t] ++ as -> bind_ (recreateCode plateVertices pbn t) (recreateVertex plateVertices pbn as)
  | [(ListNode _) & t] ++ as -> bind_ (recreateCode plateVertices pbn t) (recreateVertex plateVertices pbn as)
  | [] -> unit_

end

let modifiedBFS : all v. all l. v -> v -> Digraph v l -> Bool
  = lam source. lam dst. lam g.
  recursive let work = lam fs. lam level. lam dist:Map v Int. lam u.
    if null fs then u else
    match
      foldl (lam acc:([v], Map v Int,Bool). lam f.
        foldl (lam acc:([v], Map v Int,Bool). lam v.
          if mapMem v acc.1 then
            if digraphEqv g dst v then
              (acc.0,acc.1, false)
            else acc
          else (cons v acc.0, mapInsert v level acc.1,acc.2)
        ) acc (digraphSuccessors f g)) ([],dist,u) fs
      with (ns, dist, u) then
        if not u then u
        else
          work ns (addi level 1) dist u
      else never
    in
    work [source] 1 (mapInsert source 1 (mapEmpty (digraphCmpv g))) true

lang TransformPBN = ConjugatePrior

  -- accumulator for the transformation part
  type TAcc =
  {
    accName:Map Name Name,
    plateTBR:Set Name
  }

  sem emptyTAcc: () -> TAcc
  sem emptyTAcc = 
  | _ -> {accName=mapEmpty nameCmp,plateTBR=setEmpty nameCmp}


  sem orderPlates pbn order =
  | [p] ++ as ->
    match p with PlateNode p in
    if any (lam v. match (getPlateId v) with Some pid then nameEq pid p.ident else false) as then 
      snoc (orderPlates pbn order as) p.ident 
    else orderPlates pbn (snoc order p.ident) as 
  | [] -> order

  sem transformPBN: (PBN,TAcc) -> PBN
  sem transformPBN =
  | (pbn,tAcc) -> 
    let plates = filter (lam v. match v with PlateNode _ then true else false) (digraphVertices pbn.g) in
    let plateIds = (map getId plates) in
    let orderedPlates = orderPlates pbn [] plates in
    let plateTargets = map (lam p. 
          let p = mapLookupOrElse (lam. error "Lookup failed") p pbn.m in
          match p with PlateNode p in
          filter (lam v. let v= mapLookupOrElse (lam. error"") v pbn.m in match (getPlateId v) with Some pid then nameEq pid p.ident else false) pbn.targets) orderedPlates in
    match foldl (lam acc. lam targets.
      match acc with (pbn,tAcc) in
      match transformPBNH pbn tAcc targets with (pbn,tAcc) in
      match transformPBNH pbn tAcc (setToSeq tAcc.plateTBR) with (pbn,tAcc) in
      (pbn,tAcc)) (pbn,tAcc) plateTargets  with (pbn,tAcc) in
    let plateIdsSet = setOfSeq nameCmp plateIds in
    let otherTargets = filter (lam t. not (mapMem t plateIdsSet)) pbn.targets in
    (transformPBNH pbn tAcc otherTargets).0
    -- for each target, graft and reorder the target
  sem transformPBNH pbn tAcc =
  | [tId]++as -> 
    let t = mapLookupOrElse (lam. error "lookup failed.") tId pbn.m in
    let graftRes:(PBN,TAcc) = graft pbn tAcc t in
    match graftRes with (pbn,tAcc) in
    let reorderRes:(PBN,TAcc) = reorder pbn tAcc  t in
    match reorderRes with (pbn,tAcc) in
    transformPBNH pbn tAcc as
  | [] -> (pbn,tAcc)

  sem isStabilized: PBN -> Vertex -> Bool
  sem isStabilized pbn =
  | RandomVarNode v -> if eqi (deref v.state) 2 then true else false
  | MultiplexerNode m -> match m.list with ListNode l in
    foldl (lam acc. lam i. or acc (isStabilized pbn i)) false l.items
  
  -- set its state as marginalized
  sem addToMarginalized: Dist -> Vertex -> ()
  sem addToMarginalized q = 
  | RandomVarNode v -> modref v.state 1; modref v.mDist (Some q)

  -- check whether childPID plate is nested in the parentPID plate
  sem isNestedPl: PBN -> Name -> Name -> Bool
  sem isNestedPl pbn parentPID = 
  | childPID -> match mapLookupOrElse (lam. error "isNestedPl:Lookup failed.") childPID pbn.m with 
    (PlateNode {plateId=pid}|FoldNode {plateId=pid}) in
    match deref pid with Some pid then -- if child's one of outer plates is parent plate
      or (nameEq pid parentPID) (isNestedPl pbn parentPID pid) else false

  -- (child, parent)
  sem createMParameter: PBN -> TAcc -> (Vertex,Vertex) -> Option (PBN,TAcc,Vertex,Dist)        
  sem createMParameter pbn tAcc =
  | (RandomVarNode v, RandomVarNode p) & t -> 
    createMParameterH pbn tAcc t (deref v.plateId, deref p.plateId, None ())
  | (RandomVarNode v, MultiplexerNode p) & t -> 
    match p.list with ListNode l in
    createMParameterH pbn tAcc t (deref v.plateId, deref l.plateId, deref p.plateId)

  sem createMParameterH: PBN -> TAcc -> (Vertex,Vertex) -> (Option Name, Option Name, Option Name) -> Option (PBN,TAcc,Vertex,Dist)        
  sem createMParameterH pbn tAcc t =
  | (None (), None (), None ()) -> createMParameterNP pbn tAcc t
  | (Some pid, None (), None ()) -> createMParameterTDP pbn tAcc t 
  | (Some pid, Some pid2, None ()) -> if nameEq pid pid2 then createMParameterNP pbn tAcc t else
    if isNestedPl pbn pid2 pid then
      createMParameterTDP pbn tAcc t else None ()
  | (Some pid, Some pid2, Some pid3) ->
    if nameEq pid pid2 then
      if nameEq pid2 pid3 then (createMParameterNP pbn tAcc t)
      else None ()
    else if nameEq pid2 pid3 then
          if isNestedPl pbn pid3 pid then
              createMParameterTDP pbn tAcc t
          else None ()
          else None ()
  | (Some pid, None (), Some pid2) -> if nameEq pid pid2 then createMParameterTDP pbn tAcc t else 
    if isNestedPl pbn pid2 pid then
      createMParameterTDP pbn tAcc t else None ()
  | _ -> None ()

  sem createMarginalizedListParam pbn = 
  | MultiplexerNode p -> match p.list with ListNode l in
    let paramsId =nameSym "mLstParams" in
    let params = match l.create with Some rep then
        match deref l.outParam with Some outparam then (nvar_ outparam)
        else
          match get l.items 0 with RandomVarNode v in 
          match deref (v.mDist) with Some pMarginalizedDist in create_ rep (ulam_ "" (getParams pMarginalizedDist))
      else
        seq_ (map (lam i. match i with RandomVarNode v in 
            match deref (v.mDist) with Some pMarginalizedDist in getParams pMarginalizedDist
          ) l.items) in
    -- create a list from the marginalized parameters of the list items
    let paramBlock = CodeBlockNode {ident=nameSym "", code=nulet_ paramsId params, ret=false, plateId=l.plateId} in
    let pbn = addVertexPBN pbn paramBlock in
    let pbn = inheritMDependencies pbn paramBlock p.list in
    let pbn = match l.create with Some rep then match deref l.outParam with Some outparam then 
        {pbn with g=digraphMaybeAddEdge (mapLookupOrElse (lam. error "") outparam pbn.m) paramBlock 0 pbn.g} else pbn else pbn in
    -- create a block to select the ith marginalized parameter from the list
    let selectedParamId = nameSym "selectedParam" in
    let pMarginalizedDistParam = nulet_ selectedParamId (get_ (nvar_ paramsId) (nvar_ p.indexId)) in
    let selectedBlock = CodeBlockNode {ident=nameSym "", code=pMarginalizedDistParam, ret=false, plateId=l.plateId} in
    let pbn = addVertexPBN pbn selectedBlock in
    let pbn = {pbn with g=digraphMaybeAddEdge paramBlock selectedBlock 0 pbn.g} in
    let pbn = inheritMDependencies pbn selectedBlock (MultiplexerNode p) in

    -- change the marginalized distribution with the selected one
    match get l.items 0 with RandomVarNode v in 
    match deref v.mDist with Some mdist in
    (pbn, changeParams selectedParamId mdist, selectedBlock)


  sem createMParameterNP: PBN -> TAcc -> (Vertex, Vertex) -> Option (PBN,TAcc,Vertex,Dist)
  sem createMParameterNP pbn tAcc =
  | (RandomVarNode t, RandomVarNode p)&v -> match deref p.mDist with Some pMarginalizedDist in
    match posteriorPredictive (deref t.plateId) (deref t.dist, pMarginalizedDist) with Some (rho,q) then
      let pbn = addVertexPBN pbn rho in
      -- inherit the dependencies
      let pbn = inheritMDependencies pbn rho v.0 in
      let pbn = inheritMDependencies pbn rho v.1 in
      Some (pbn, tAcc,rho, q)
    else None ()
  | (RandomVarNode t, MultiplexerNode p)&v -> match p.list with ListNode l in
    match createMarginalizedListParam pbn v.1 with (pbn, pMarginalizedDist, selectedBlock) in
    modref p.mDist (Some pMarginalizedDist);
    match posteriorPredictive (deref t.plateId) (deref t.dist, pMarginalizedDist) with Some (rho,q) then
      let pbn = addVertexPBN pbn rho in
      -- inherit the dependencies
      let pbn = inheritMDependencies pbn rho v.0 in
      let pbn = inheritMDependencies pbn rho v.1 in
      let pbn = {pbn with g=digraphMaybeAddEdge selectedBlock rho 0 pbn.g} in
      Some (pbn, tAcc,rho, q)
    else None ()

  sem addParam pbn param =
  | CodeBlockNode ({code=TmLet t}&c) -> match t.body with TmRecord r in
    let newRec = record_add (int2string (mapSize r.bindings)) param t.body in
    let newCB = (CodeBlockNode {c with code = nulet_ t.ident newRec}) in
    (addVertexPBN pbn newCB, newCB)
    

  sem createMParameterTDP: PBN -> TAcc -> (Vertex, Vertex) -> Option (PBN,TAcc,Vertex,Dist)
  sem createMParameterTDP pbn tAcc =
  | (RandomVarNode t, RandomVarNode p)&v -> (if debug then print (join ["createMParameterTDP-nolist",v2str v.0,"\n"]) else ());
    let tAcc = {tAcc with plateTBR=setInsert t.ident tAcc.plateTBR} in
    -- Get parent's marginalized distribution parameters
    match (deref p.mDist) with Some pMarginalizedDist in 
    let param = (getParams pMarginalizedDist) in
    let ppId = nameSym "pMargParam" in
    modref p.mDist (Some (changeParams ppId pMarginalizedDist));
    match createMParameterTDPH pbn tAcc param ppId v with (pbn,_) in
    createMParameterNP pbn tAcc v 
  | (RandomVarNode t, MultiplexerNode p)&v ->(if debug then print (join ["createMParameterTDP-muxnode",v2str v.0,"\n"]) else ());
    let tAcc = {tAcc with plateTBR=setInsert t.ident tAcc.plateTBR} in
    match p.list with ListNode l in --get the list input of the multiplexer node
    let ppId = nameSym "pMargParam" in
    let param = match l.create with Some rep then
        match get l.items 0 with RandomVarNode v in 
        match deref (v.mDist) with Some pMarginalizedDist in create_ rep (ulam_ "" (getParams pMarginalizedDist))
      else
        seq_ (map (lam i. match i with RandomVarNode v in 
            match deref (v.mDist) with Some pMarginalizedDist in getParams pMarginalizedDist
          ) l.items) in
    match createMParameterTDPH pbn tAcc param ppId v with (pbn,cb) in
    let selectedParamId = nameSym "selectedParam" in
    let pMarginalizedDistParam = nulet_ selectedParamId (get_ (nvar_ ppId) (nvar_ p.indexId)) in
    let selectedBlock = CodeBlockNode {ident=nameSym "", code=pMarginalizedDistParam, ret=false, plateId=t.plateId} in
    let pbn = addVertexPBN pbn selectedBlock in

    match get l.items 0 with RandomVarNode a in 
    let pbn = {pbn with g=digraphMaybeAddEdge cb selectedBlock 0 pbn.g} in
    match deref a.mDist with Some mdist in
    let pbn = inheritMDependencies pbn selectedBlock (MultiplexerNode p) in
    let pMarginalizedDist = changeParams selectedParamId mdist in
    modref p.mDist (Some pMarginalizedDist);
    match posteriorPredictive (deref t.plateId) (deref t.dist, pMarginalizedDist) with Some (rho,q) then
      let pbn = addVertexPBN pbn rho in
      -- inherit the dependencies
      let pbn = inheritMDependencies pbn rho v.0 in
      let pbn = inheritMDependencies pbn rho v.1 in
    let pbn = {pbn with g=digraphMaybeAddEdge selectedBlock rho 0 pbn.g} in
      Some (pbn, tAcc,rho, q)
    else None ()

  sem propagateThroughPlates pbn initParam innerAccID newParamId v acc =
  | (Some targetPID, ppid)  -> 
    (if debug then print (join ["propagateThroughPlates",v2str v.0,"\n"]) else ());
      match mapLookup targetPID pbn.m with Some container in 
      match (getPlateId container) with Some pid then 
        let res = match ppid with Some ppid then 
          if nameEq pid ppid then Some (createFoldNode pbn (initParam) newParamId innerAccID targetPID v container)
        else None () else None () in match res with Some _ then res
        else let gpContainer = mapLookupOrElse (lam. error "") pid pbn.m in
          let index = match gpContainer with PlateNode _ then 0 else match gpContainer with FoldNode l in mapSize (deref l.vToIndex) in
          let outerAccId = match gpContainer with FoldNode gp then gp.lamAccId else (nameSym "outerParam") in
          let param = tupleproj_ index (nvar_ outerAccId)  in
          match createFoldNode pbn param newParamId innerAccID targetPID v container with (pbn, (FoldNode fl), index) & res in
          propagateThroughPlates pbn initParam innerAccID (Some outerAccId) v (Some res) (deref fl.plateId, None ()) 
      else Some (createFoldNode pbn (initParam) newParamId innerAccID targetPID v container)
  | _ -> acc

  -- TODO: here if already in fold parameters then do not add
  sem createFoldNode pbn param newParamId innerAccID targetPID v =
  | PlateNode pl ->
    (if debug then print (join ["createFoldNode-plate",v2str v.0,"\n"]) else ());
    let accAppId = nameSym "pInitParam" in
    let lamAccId = match newParamId with Some id then id else innerAccID in
    let retBId = nameSym "r" in
    let cbInitParam = CodeBlockNode {ident = accAppId, code = nulet_ accAppId (utuple_ [param]), ret=false, plateId=pl.plateId} in
    let pbn = addVertexPBN pbn cbInitParam in
    let pbn = inheritMDependencies pbn cbInitParam v.1 in
    let index = 0 in
    let vToIndex = mapInsert (getId v.1) index (mapEmpty nameCmp) in
    let f = FoldNode {ident = pl.ident, varIds=pl.varIds, plateId=pl.plateId, lamAccId=lamAccId, accId=accAppId, iterlId=pl.iterlId, retBlockId=retBId, vToIndex=ref vToIndex} in
    let pbn = {pbn with g=digraphRemoveVertex (PlateNode pl) pbn.g} in -- since constructor tag is diff should be removed
    let pbn = addVertexPBN pbn f in
    ({pbn with g=digraphMaybeAddEdge cbInitParam f 0 pbn.g}, f,0)
  | FoldNode fl ->
    match mapLookup (getId v.1) (deref fl.vToIndex) with Some index then (pbn, (FoldNode fl), index)
    else 
      let index = mapSize (deref fl.vToIndex) in
      modref fl.vToIndex (mapInsert (getId v.1) index (deref fl.vToIndex));
      let cbInitParam = mapLookupOrElse (lam. error "") fl.accId pbn.m in
      match addParam pbn param cbInitParam with (pbn,_) in
      (pbn, (FoldNode fl), index)

  sem createMParameterTDPH pbn tAcc param ppId =
  | (RandomVarNode t, (RandomVarNode _ | MultiplexerNode _)) & v  -> 
    (if debug then print (join ["createMParameterTDPH",v2str v.0,"\n"]) else ());
    match (deref t.plateId) with Some tpid in
    let parentPID = match v.1 with MultiplexerNode m then match m.list with ListNode l in deref l.plateId
      else match v.1 with RandomVarNode p in deref p.plateId in
    let innerAccID = match mapLookup tpid pbn.m with Some (FoldNode f) then f.lamAccId else nameSym "p" in 
    match propagateThroughPlates pbn param innerAccID (None ()) v (None ()) (Some tpid, parentPID) with Some (pbn, f, index) in
    let cb = CodeBlockNode {ident = nameSym "", code = nulet_ ppId (tupleproj_ index (nvar_ innerAccID)), ret=false, plateId=t.plateId} in
    let pbn = addVertexPBN pbn cb in
    (pbn,cb)

  sem inheritMDependencies: PBN -> Vertex -> Vertex -> PBN
  sem inheritMDependencies pbn toV =
  | (MultiplexerNode m)&fromV -> 
    let parents = filter (lam v. match v with CodeBlockNode _ then true
                            else match v with RandomVarNode r then true
                            else false) (digraphPredeccessors fromV pbn.g) in
    let g = foldl (lam acc. lam gp. digraphMaybeAddEdge gp toV 0 acc) pbn.g parents in
    {pbn with g=g}
  | ListNode l -> foldl (lam pbn. lam i. inheritMDependencies pbn toV i) pbn l.items
  | fromV ->
    -- get the codeblock parents and stabilized nodes of t
    let parents = filter (lam v. match v with CodeBlockNode _ then true
                            else match v with RandomVarNode r then eqi (deref r.state) 2
                            else false) (digraphPredeccessors fromV pbn.g) in
    let g = foldl (lam acc. lam gp. digraphMaybeAddEdge gp toV 0 acc) pbn.g parents in
    {pbn with g=g}
  
  sem marginalize: PBN -> TAcc -> Vertex -> (PBN,TAcc)
  sem marginalize pbn tAcc =
  | (RandomVarNode v) & t -> (if debug then print (join ["Marginalize ", v2str t, "\n"]) else ());
    -- filter its random variable parents that are not stabilized
    let parents = filter (lam p. match p with RandomVarNode _ | MultiplexerNode _ then 
        not (isStabilized pbn p) else false) (digraphPredeccessors t pbn.g) in
    if null parents then (if debug then print (join ["Marginalize: no parents", "\n"]) else ());
      addToMarginalized (deref v.dist) t; (pbn, tAcc)
    else (if debug then print (join ["Marginalize: has parents", "\n"]) else ());
      let parent = get parents 0 in
      if not (modifiedBFS parent t pbn.g) then
        (if debug then print "Marginalize: can cause cycles reordering the parent\n" else ());
        let res = reorder pbn tAcc parent in
        marginalize res.0 res.1 t
      else match createMParameter pbn tAcc (t, parent) with Some (pbn,tAcc,rho,q) then
          addToMarginalized q t; 
          let g = digraphMaybeAddEdge rho t 0 pbn.g in
          ({pbn with g=g}, tAcc)
        else (if debug then print "Marginalize: no conjugate prior rel\n" else ());
          match reorder pbn tAcc parent with (pbn, tAcc) in
          match marginalize pbn tAcc t with (pbn, tAcc) in
          (pbn, tAcc)
  
  -- the target can be a random variable or a list
  sem graft:PBN -> TAcc -> Vertex -> (PBN, TAcc)
  sem graft pbn tAcc =
  | (RandomVarNode v) & t ->     
    if eqi (deref v.state) 2 then (pbn, tAcc) -- if is stabilized, do nothing
    else (if debug then print (join ["Graft(", v2str t,")\n"]) else ());
      if eqi (deref v.state) 1 then -- if it is marginalized
        (if debug then print "Graft: RV t is already marginalized\n" else ());
        -- get its marginalized random variable child if any
        let child = filter (lam u. match u with RandomVarNode u then eqi (deref u.state) 1 else false) (digraphSuccessors t pbn.g) in
        -- if it does not have a marginalized child, then return the graph
        (if null child then (pbn,tAcc)
        else -- it cannot have more than one marginalized child, throw error
        (if not (eqi (length child) 1) then error "Graft: can only have one marginalized child."
         else -- if it has one marginalized child
          (if debug then print (join ["child node ", (v2str (get child 0)), " to be pruned\n"]) else ());
           -- prune the child so t will become the terminal node on its marginalized path
          pruneD pbn tAcc (get child 0)))
    else -- not marginalized
      (if debug then print "Graft: RV t is not marginalized\n" else ());
      -- get the non-stabilized parents that are either random var or a multiplexer node.
      let parents = filter (lam v. match v with (RandomVarNode _ | MultiplexerNode _) then not (isStabilized pbn v)
                     else false) (digraphPredeccessors t pbn.g) in
      match if null parents then (if debug then print (join ["Graft: t has no parents\n"]) else ());
        marginalize pbn tAcc t
        else let res = switch (get parents 0)
            case RandomVarNode p then (if debug then print (join ["Graft: parent of t is a rv\n"]) else ());
              graft pbn tAcc (RandomVarNode p)
            case MultiplexerNode p then (if debug then print "Graft: t's parent comes from a list\n" else ());
              graft pbn tAcc p.list
            end in marginalize res.0 res.1 t 
      with (pbn, tAcc) in
      -- if t has any child that belongs to a list, prune t; otherwise, the order of the list would not be handled
      if any (lam c. match getListId c with Some _ then true else false) (digraphSuccessors t pbn.g) then pruneD pbn tAcc t else (pbn, tAcc)
  | (ListNode l) & t ->
    let children = digraphSuccessors t pbn.g in
    let res = foldl (lam acc. lam e. graft acc.0 acc.1 e) (pbn,tAcc) l.items in res
    --if gti (length children) 1 then pruneD res.0 res.1 t else res 

  sem pruneD: PBN -> TAcc -> Vertex -> (PBN,TAcc)
  sem pruneD pbn tAcc =
  | (RandomVarNode v) & t -> (if debug then print (join ["Prune(", v2str t,")\n"]) else ());
    if neqi (deref v.state) 1 then error "Prune: t is not marginalized"
    else
      -- get its marginalized child if any
      let children = filter (lam u. match u with RandomVarNode u then eqi (deref u.state) 1 else false) (digraphSuccessors t pbn.g) in
      -- if it does not have a marginalized child then reorder the vertex t.
      (if null children then reorder pbn tAcc t
      else match eqi (length children) 1 with false then error "Prune: t has more than one marginalized child" else
        -- if it has a marginalized child then prune it first.
        match pruneD pbn tAcc (get children 0) with (pbn, tAcc) in
        reorder pbn tAcc t)
  | (ListNode l) & t -> foldl (lam acc. lam e. pruneD acc.0 acc.1 e) (pbn, tAcc) l.items


  sem createRParameter: PBN -> TAcc -> (Vertex,Vertex) -> (PBN,TAcc)
  sem createRParameter pbn tAcc =
  | (RandomVarNode v, RandomVarNode p) & t -> 
    (if debug then print (join ["createRParameter ", v2str t.0, "\n"]) else ());
    createRParameterH pbn tAcc t (deref v.plateId, deref p.plateId, None ())
  | (RandomVarNode v, MultiplexerNode p) & t -> 
    match p.list with ListNode l in
    createRParameterH pbn tAcc t (deref v.plateId, deref l.plateId, deref p.plateId)

  sem createRParameterH: PBN -> TAcc -> (Vertex, Vertex) -> (Option Name, Option Name, Option Name) -> (PBN,TAcc)
  sem createRParameterH pbn tAcc t =
  | (None (), None (), None ()) -> createRParameterNP pbn tAcc (None ()) t
  | (Some pid, None (), None ()) -> createRParameterTDP pbn tAcc t 
  | (Some pid, Some pid2, None ()) -> if nameEq pid pid2 then
      createRParameterNP pbn tAcc (None ()) t else if isNestedPl pbn pid2 pid then
      createRParameterTDP pbn tAcc t else never
  | (Some pid, Some pid2, Some pid3) -> 
    if nameEq pid pid2 then
      if nameEq pid2 pid3 then (createRParameterNP pbn tAcc (None ()) t) else never else never
  | (Some pid, None (), Some pid2) -> if nameEq pid pid2 then 
      match createRParameterTDP pbn tAcc t with (pbn,tAcc) in 
      let pbn = removeVertexPBN pbn t.1 in (pbn, tAcc) else 
      if isNestedPl pbn pid2 pid then
      match createRParameterTDP pbn tAcc t with (pbn,tAcc) in 
      let pbn = removeVertexPBN pbn t.1 in (pbn, tAcc) else  never
  | _ -> never

  sem inheritRDependencies: PBN -> TAcc -> (Vertex, Vertex, Vertex) -> PBN
  sem inheritRDependencies pbn tAcc =
  | (t, p, rho) -> 
    let filterC = (lam v. match v with CodeBlockNode _ then true
                                else match v with RandomVarNode r then neqi (deref r.state) 1
                                else false) in
    let parentsT = filter filterC (digraphPredeccessors t pbn.g) in
    -- get the codeblock parents and stabilized nodes of p
    let parentsP = filter filterC (digraphPredeccessors p pbn.g) in
    -- inherit the dependencies
    let g = foldl (lam acc. lam gp. digraphMaybeAddEdge gp rho 0 acc) pbn.g parentsT in
    let g = foldl (lam acc. lam gp. let g = digraphRemoveEdge gp p 0 acc in digraphMaybeAddEdge gp rho 0 g) g parentsP in
    {pbn with g=g}

  sem createRParameterNP: PBN -> TAcc -> Option (Name,Expr) -> (Vertex, Vertex) -> (PBN,TAcc)
  sem createRParameterNP pbn tAcc indices =
  | (RandomVarNode t, RandomVarNode p)&v ->
    (if debug then print (join ["createRParameterNP ", v2str v.0, "\n"]) else ());
    let obs = match t.val with Some _ then t.val else Some (nvar_ t.ident) in
    match deref p.mDist with Some pMarginalizedDist in
    match posterior obs indices (deref t.plateId) (deref t.dist,pMarginalizedDist) with (rho,q,_) in
    modref p.dist q;
    let pbn = addVertexPBN pbn rho in
    let pbn = inheritRDependencies pbn tAcc (v.0,v.1,rho) in
    let g = digraphMaybeAddEdges [(rho, v.1, 0), (v.0, rho, 0)] pbn.g in
    modref p.mDist (Some q);
    ({pbn with g=g}, tAcc)
  | (RandomVarNode t, MultiplexerNode p)&v -> match p.list with ListNode l in
    match deref p.mDist with (Some pMarginalizedDist) in
    let obs = match t.val with Some _ then t.val else Some (nvar_ t.ident) in
    match posterior obs (None ()) (deref t.plateId) (deref t.dist,pMarginalizedDist) with (rho,q,paramNames) in
    let pbn = addVertexPBN pbn rho in
    let pbn = inheritRDependencies pbn tAcc (v.0,v.1,rho) in
    let pbn = {pbn with g= digraphMaybeAddEdges [(rho, v.1, 0), (v.0, rho, 0)] pbn.g} in
    let reorderedParam  = getParams q in
    let i = nameSym "i" in
    let e = nameSym "e" in
    let params = match l.create with Some rep then
      match deref l.outParam with Some outParam then (nvar_ outParam)
      else match get l.items 0 with RandomVarNode v in 
        match deref (v.mDist) with Some pMarginalizedDist in create_ rep (ulam_ "" (getParams pMarginalizedDist)) 
    else seq_ (map (lam i. match i with RandomVarNode v in match (deref v.mDist) with Some md in getParams md) l.items) in
    let code = mapi_ (nulam_ i 
      (nulam_ e 
        (if_ (eqi_ (nvar_ i) (nvar_ p.indexId)) reorderedParam (nvar_ e)))) (params) in
    let paramName = nameSym "rP" in
    let paramBlock = CodeBlockNode {ident=paramName, code=nulet_ paramName code, ret=false,plateId=t.plateId} in
    let pbn = addVertexPBN pbn paramBlock in
    let edges = [(rho,paramBlock,0)] in
    let pbn = {pbn with g = digraphMaybeAddEdges edges pbn.g} in
    match l.create with Some rep then
      match get l.items 0 with RandomVarNode v in
      match deref l.lamParam with Some id in
      let pbn = match (deref l.outParam) with Some outParam then
        let cb = mapLookupOrElse (lam. error "cannot") outParam pbn.m in
        {pbn with g=digraphMaybeAddEdge cb paramBlock 0 pbn.g}
      else pbn in
      modref l.outParam (Some paramName);
      let q = changeParams id q in
      modref v.dist q;
      modref v.mDist (Some q);
      let edges = [(paramBlock, p.list, 0)] in
      let pbn = {pbn with g = digraphMaybeAddEdges edges pbn.g} in
      (pbn, tAcc)
    else 
      let res = foldl (lam acc. lam e. 
        match e with RandomVarNode v in 
        match acc with (pbn, i) in
        let param = get_ (nvar_ paramName) (int_ i) in
        let paramid = nameSym "param" in
        let cb = CodeBlockNode {ident=nameSym "", code=nulet_ paramid param,ret=false,plateId=l.plateId} in
        let pbn = addVertexPBN pbn cb in
        let edges = [(rho,cb,0),(cb,e,0)] in
        let pbn = {pbn with g=digraphMaybeAddEdges edges pbn.g} in
        let q = changeParams paramid q in
        modref v.dist q;
        modref v.mDist (Some q);
        (pbn,addi i 1)) (pbn, 0) l.items in
      match res with (pbn, _) in
    (pbn, tAcc)

  sem createRParameterTDPH: PBN -> TAcc -> (Name, Vertex) -> Option (Name, Expr) -> (Vertex, Vertex) -> (PBN, TAcc, Name, Vertex, Option Vertex)
  sem createRParameterTDPH pbn tAcc out indices =
  | (RandomVarNode t, RandomVarNode p)&v ->
    match out with (outName, outBlock) in
    let obs = match t.val with Some _ then t.val else Some (nvar_ t.ident) in
    match deref p.mDist with Some pMarginalizedDist in
    match posterior obs indices (deref t.plateId) (deref t.dist, pMarginalizedDist) with (rho, q, paramNames) in
    let pbn = addVertexPBN pbn rho in
    let paramName = nameSym "rP" in
    let paramBlock = CodeBlockNode {ident=nameSym "", code=nulet_ paramName (utuple_ (map nvar_ paramNames)),ret=false,plateId=t.plateId} in
    let pbn = addVertexPBN pbn paramBlock in
    let edges = [(outBlock, v.1, 0), (rho, paramBlock, 0), (v.0, rho, 0)] in
    let pbn = {pbn with g = digraphMaybeAddEdges edges pbn.g} in
    let q = changeParams outName q in
    modref p.dist q;
    modref p.mDist (Some q);
    (pbn, tAcc, paramName, paramBlock, Some rho)
  | (RandomVarNode t, MultiplexerNode p)&v -> match p.list with ListNode l in
    match deref p.mDist with (Some pMarginalizedDist) in
    let obs = match t.val with Some _ then t.val else Some (nvar_ t.ident) in
    match posterior obs (None ()) (deref t.plateId) (deref t.dist,pMarginalizedDist) with (rho,q,paramNames) in
    match deref t.plateId with Some pid in
    match mapLookup pid pbn.m with Some (FoldNode f) in
    match out with (outName, outBlock) in
    let pbn = addVertexPBN pbn rho in
    let pbn = inheritRDependencies pbn tAcc (v.0,v.1,rho) in
    let pbn = {pbn with g= digraphMaybeAddEdges [(rho, v.1, 0), (v.0, rho, 0)] pbn.g} in
    let reorderedParam  = getParams q in
    let i = nameSym "i" in
    let e = nameSym "e" in
    let index = mapLookupOrElse (lam. error "") (getId v.1) (deref f.vToIndex) in
    let params = tupleproj_ index (nvar_ f.lamAccId)  in
    let code = match l.create with Some rep then
      match deref l.outParam with Some outParam then (nvar_ outParam)
      else match get l.items 0 with RandomVarNode v in 
        match deref (v.mDist) with Some pMarginalizedDist in create_ rep (ulam_ "" (getParams pMarginalizedDist)) 
    else seq_ (map (lam i. match i with RandomVarNode v in match (deref v.mDist) with Some md in getParams md) l.items) in
    let code = mapi_ (nulam_ i 
      (nulam_ e 
        (if_ (eqi_ (nvar_ i) (nvar_ p.indexId)) reorderedParam (nvar_ e)))) (params) in
    let paramName = nameSym "rP" in
    let paramBlock = CodeBlockNode {ident=nameSym "", code=nulet_ paramName code, ret=false,plateId=t.plateId} in
    let pbn = addVertexPBN pbn paramBlock in
    let edges = [(rho,paramBlock,0)] in
    let pbn = {pbn with g = digraphMaybeAddEdges edges pbn.g} in
    match l.create with Some rep then
      match get l.items 0 with RandomVarNode v in
      match deref l.lamParam with Some id in
      modref l.outParam (Some outName);
      let q = changeParams id q in
      modref v.dist q;
      modref v.mDist (Some q);
      let edges = [(outBlock, p.list, 0)] in
      let pbn = {pbn with g = digraphMaybeAddEdges edges pbn.g} in
      (pbn, tAcc, paramName, paramBlock, None ())
    else 
      let res = foldl (lam acc. lam i.
        match acc with (pbn, tAcc, edges, cnt) in
        match i with RandomVarNode v in
        let outNameI = nameSym (join ["postParam", int2string cnt]) in
        let outBlockI = CodeBlockNode {ident=nameSym "", code=nulet_ outNameI (get_ (nvar_ outName) (int_ cnt)), ret=false, plateId=l.plateId} in
        let pbn = addVertexPBN pbn outBlockI in
        let newEdges = [(outBlock,outBlockI, 0),(outBlockI, i, 0),((FoldNode f),outBlockI,0)] in
        let q = changeParams outNameI q in
        modref v.dist q; modref v.mDist (Some q);
        (pbn, tAcc, join [newEdges, edges], addi cnt 1)
      ) (pbn, tAcc, [(outBlock, p.list, 0)], 0) l.items in
      match res with (pbn, tAcc, edges, _) in
      let pbn = {pbn with g = digraphMaybeAddEdges edges pbn.g} in
      (pbn, tAcc, paramName, paramBlock, None ())

      -- propagate is not correct when several things are propagated
  sem propagateThroughPlatesR pbn v id =
  | (Some targetPID, None ())  -> 
      match mapLookup targetPID pbn.m with Some (FoldNode f) in
      match deref f.plateId with Some lpid then
        match mapLookup lpid pbn.m with Some (FoldNode gf) in
          let index = mapLookupOrElse (lam. error "") (getId v.1) (deref gf.vToIndex) in
          let param = tupleproj_ index (nvar_ f.ident) in
          let res = match mapLookup gf.retBlockId pbn.m with Some (CodeBlockNode rt) then 
           addParam pbn param (CodeBlockNode rt)
          else let retBlock = CodeBlockNode {ident=gf.retBlockId, code= (nulet_ gf.retBlockId (utuple_ [param])), ret=false,plateId=f.plateId} in
            (addVertexPBN pbn retBlock, retBlock) in 
          match res with (pbn,retBlock) in
          let ret = CodeBlockNode {ident=nameSym "", code=(nvar_ gf.retBlockId), ret=true,plateId=f.plateId} in
          let pbn = addVertexPBN pbn ret in 
          let pbn = {pbn with g=digraphMaybeAddEdge (FoldNode f) retBlock 0 pbn.g} in
          propagateThroughPlatesR pbn v gf.ident (Some lpid, None ())
      else (pbn, id)
  | _ -> (pbn, id)
    
  sem createRParameterTDP: PBN -> TAcc -> (Vertex, Vertex) -> (PBN,TAcc)
  sem createRParameterTDP pbn tAcc =
  | (RandomVarNode t, (RandomVarNode _ | MultiplexerNode _)) & v  -> 
    match deref t.plateId with Some pid in
    match mapLookup pid pbn.m with Some (FoldNode f) in
    let parentPID = match v.1 with MultiplexerNode m then match m.list with ListNode l in deref l.plateId
      else match v.1 with RandomVarNode p in deref p.plateId in
    match propagateThroughPlatesR pbn v f.ident (Some pid, parentPID) with (pbn, outfid) in
    -- create the outblock that gets the result from foldnode
    let outName = nameSym "postParam" in
    let index = mapLookupOrElse (lam. error "index not found") (getId v.1) (deref f.vToIndex) in
    let ppid = match v.1 with MultiplexerNode m then match m.list with ListNode l in l.plateId
      else match v.1 with RandomVarNode p then p.plateId else never in
      -- here f.ident should be the last foldname
    let outBlock = CodeBlockNode {ident=nameSym "", code=nulet_ outName (tupleproj_ index (nvar_ outfid)), ret=false, plateId=ppid} in
    let pbn = addVertexPBN pbn outBlock in
    match createRParameterTDPH pbn tAcc (outName,outBlock) (None ()) v with (pbn, tAcc, paramName, paramBlock) in

    -- here should get the return value as the marginalized
    let res = match mapLookup f.retBlockId pbn.m with Some (CodeBlockNode rt) then
      addParam pbn (nvar_ paramName) (CodeBlockNode rt)
    else
      let retBlock = CodeBlockNode {ident=f.retBlockId, code= (nulet_ f.retBlockId (utuple_ [nvar_ paramName])), ret=false,plateId=t.plateId} in
      (addVertexPBN pbn retBlock, retBlock) in 
    match res with (pbn,retBlock) in
    let ret = CodeBlockNode {ident=nameSym "", code=(nvar_ f.retBlockId), ret=true,plateId=t.plateId} in
    let pbn = addVertexPBN pbn ret in

    match mapLookup outfid pbn.m with Some (FoldNode f) in
    let edges = [(FoldNode f,outBlock,0), (paramBlock,retBlock,0)] in
    let pbn = {pbn with g=digraphMaybeAddEdges edges pbn.g} in
    (pbn, tAcc)


  sem reorder: PBN -> TAcc -> Vertex -> (PBN,TAcc)
  sem reorder pbn tAcc =
  | (RandomVarNode v) & t -> (if debug then print (join ["Reorder ", v2str t, "\n"]) else ());
    -- already stabilized, return the graph
    if eqi (deref v.state) 2 then (pbn,tAcc) else
    let parents = filter (lam v. 
                          match v with (RandomVarNode _ | MultiplexerNode _) then
                            not (isStabilized pbn v)
                          else false) (digraphPredeccessors t pbn.g) in
    if null parents then --if it has no rv parents then directly stabilize the node
      (if debug then print ("Random variable has no parents so directly stabilize") else ());
      -- change its state from 0 to 1 [from assumed to stabilized]
      -- set its distribution as its marginalized distribution
      modref v.state 2;
      match deref v.mDist with Some mDist in
      modref v.dist mDist;
      (pbn, tAcc)
    else 
      let parent = get parents 0 in
      let pbn = {pbn with g=digraphRemoveEdge parent t 0 pbn.g} in
      match createRParameter pbn tAcc (t,parent) with (pbn,tAcc) in
      modref v.state 2;
      match deref v.mDist with Some mDist in
      modref v.dist mDist;
      (pbn, tAcc)
  | MultiplexerNode p -> reorder pbn tAcc p.list
  | ListNode l -> foldl (lam acc. lam i. reorder acc.0 acc.1 i) (pbn, tAcc) l.items
         
end

let modifiedBFS : all v. all l. v -> v -> Digraph v l -> Bool
  = lam source. lam dst. lam g.
  recursive let work = lam fs. lam level. lam dist:Map v Int. lam u.
    if null fs then u else
    match
      foldl (lam acc:([v], Map v Int,Bool). lam f.
        foldl (lam acc:([v], Map v Int,Bool). lam v.
          if mapMem v acc.1 then
            if digraphEqv g dst v then
              (acc.0,acc.1, false)
            else acc
          else (cons v acc.0, mapInsert v level acc.1,acc.2)
        ) acc (digraphSuccessors f g)) ([],dist,u) fs
      with (ns, dist, u) then
        if not u then u
        else
          work ns (addi level 1) dist u
      else never
    in
    work [source] 1 (mapInsert source 1 (mapEmpty (digraphCmpv g))) true

-- the whole algorithm with three steps: Static PBN Constructor + Conjugate Prior Transformer + Program Reconstructor
lang StaticDelay = CreatePBN + TransformPBN + RecreateProg 

  sem removeAlias env =
  | TmLet ({body=TmVar v}&t) ->
    let env = match mapLookup v.ident env with Some id then
    mapInsert t.ident id env else mapInsert t.ident v.ident env in
    removeAlias env t.inexpr
  | TmVar t -> match mapLookup t.ident env with Some id then nvar_ id else TmVar t
  | t -> smap_Expr_Expr (removeAlias env) t

  sem transformModel =
  | prog -> 
    let model = use MExprPPLStaticDelayedANF in (normalizeTerm prog) in
    let model = removeAlias (mapEmpty nameCmp) model in
    let pbn = createM model in
    let pbn = transformPBN ({pbn with targets=(distinct nameEq pbn.targets)},(emptyTAcc ())) in
    recreate pbn

  sem transformLam: Expr -> Expr
  sem transformLam =
  | TmApp ({lhs=TmApp ({lhs=TmConst ({val=CCreate ()}&c),rhs=rep}&a2),rhs=TmLam l}&a1) -> TmApp a1
  | TmApp ({lhs=TmApp ({lhs=TmConst ({val=CIter ()}&c),rhs=TmLam l}&a2),rhs=lst}&a1) -> TmApp a1
  | TmApp ({lhs=TmApp ({lhs=TmConst ({val=CIteri ()}&c),rhs=TmLam ({body=TmLam l2}&l1)}&a2),rhs=lst}&a1) ->
   TmApp a1
  | TmLam l -> let res = TmLam {l with body=transformModel l.body} in
    smap_Expr_Expr transformLam res
  | t -> smap_Expr_Expr transformLam t

  sem transform: Expr -> Expr
  sem transform =
  | prog -> transformModel (transformLam prog)
end


let staticDelay = lam prog. use StaticDelay in
  transform prog