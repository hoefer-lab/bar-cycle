(* in this module, the main module of the package, we make simple tree-based
 * hidden markov model, called BAR (bifurcating autoregressive) model in the
 * article.
 *
 * the BAR model has IPARS.hdim Gaussian hidden variables and allows for a flexible
 * coupling between them upon inheritance. the variables then get linearly
 * combined to give the normal-transformed cell cycle length.
 * *)

open! Containers
open Lacaml.D
open Option.Infix
open Util
open Tree

(* BAR model parametrization using covariances.
 * this is kept generic over the number of hidden variables! *)
module type IPARS = sig
  val hdim : int
  val cmd : mat
  val cmm : mat
  val alpha : vec
  val noise : float
end

(* a typical example *)
module IPars : IPARS = struct
  (*mother daughter hidden covariance *)
  let hdim = 2
  let cmd = Mat.of_array [| [| 0.8;  0.0 |];
                            [| 0.0; -0.4|] |]
  (* mother covariance between hidden vars. *)
  let cmm = Mat.identity 2
  (*the linear combination giving the observation mean*)
  let alpha = Vec.of_array [| 1.; 1. |]
  (*observation noise relative to the intrinsic sigma*)
  let noise = 0.0
end

(* signature for a BAR model parametrization using natural parameters as
 * described below *)
module type PARS = sig
  (** make sampling parameters based on the parametrization
   *
   *
   * x_dtrs = AA.x_mtr + xi_dtrs.
   *
   * x{,i}_dtrs have dimension 2*hdim for both daughters.
   * AA is a stacked version [A;A] of identical single-daughter matrices A.
   * AA == meanScale.
   *
   * xi_dtrs is correlated with covariance Cxi.
   *
   * Cxi usually has the block form
   * [diag(σ^2), γ diag(σ^2); diag(σ^2), γ diag(σ^2)].
   * However, since alpha can scale the outputs, the individual variances are
   * redundant. therefore usually σ==1 are chosen.
   * if xi_dtrs = B . zeta where zeta is independent unit Gaussian, then
   * Cxi = B.B^T. to get a given γ and unit σ, one can choose
   * B = sqrt(1/2 (1 - sqrt(1 - γ^2))) [1, β; β, 1] where
   * β = (1 + sqrt(1 - γ^2))/γ
   *
   * noise is an extra observation noise (not implemented)
   *
   * alpha is the obervarion vector, generating observed times as alpha.x
   *
   * *)

  val hdim : int
  val alpha : vec
  val noise : float
  val meanScale : mat
  val cdGm : mat
end

(* functor for making a natural model parametrization when covariances are given *)
module Pars (Ipars:IPARS) : PARS = struct
  include Ipars
  (* number of hidden variables *)
  let hdim = Mat.dim1 cmd
  let _cmmmd =
    let b = Mat.create (2*hdim) (2*hdim) in
    for i = 1 to hdim do for j = 1 to hdim do
        b.{i,j}           <- cmm.{i,j};
        b.{i,j+hdim}      <- cmd.{i,j};
        b.{i+hdim,j}      <- cmd.{j,i};
        b.{i+hdim,j+hdim} <- cmm.{i,j};
      done; done;
    (* rhomm and rhomd are not independent from each other; the combined block
       matrix [ rhomm, rhomd; rhomd^T, rhomm ] has to be positive definite. *)
    let () =
      let all_evs_positive =
        let b' = lacpy b in
        let all_true f v = Vec.fold (fun x y -> x && (f y)) true v in
        syev b' |> all_true ((>.) 0.) in
      assert all_evs_positive in
    b

  let cmmInv =
    (*inverse as solution to lin eqs. overwrites input! *)
    let cpy = lacpy cmm in
    getri cpy; cpy

  (* compute the covariance of the conditioned daughter variables.
   * formula for drawing from a conditioned Gaussian:
   * If the joint distribution has covariance matrix written in block form as
   * CAA, CAB; CBA, CBB where CAB=CBA^T then we have for the Gaussian
   * distribution of xA conditioned on xB:
   * mean^T = xB^T . CBB^-1 . CBA, or equvalently, mean = CAB . CBB^-1 . xB
   * covariance = CAA - CAB . CBB^-1 . CBA
   * the covariance is always smaller than for unconditioned xA but does not
   * depend on xB. *)
  let cdGm =
    let open Mat in
    gemm cmmInv cmd
    |> gemm ~transa:`T cmd
    |> neg
    |> add cmm
  let meanScale =
    (gemm cmmInv cmd |> Mat.transpose_copy)

  (*finally, the decomp for sampling*)
  (*let cdGmEvs, rotation = U.evs_rot cdGm*)

end


(* example: Exp58 on 5 best-fit parameters taken from the corr fits *)
module Apars_58_5 : PARS = struct

  let hdim = 2

  (*the linear combination giving the observation mean*)
  let alpha = Vec.of_list [ 1.; 1.286 ]
  (*observation noise relative to the intrinsic sigma*)
  let noise = 0.0
  let sigma1, sigma2 = 0.073, 0.850
  (* extra coupling between sisters, in the B mixture matrix *)
  let beta = 0.27

  let meanScale =
    (* the A matrix *)
    let a11 =  0.813 in
    let a22 =  0.475 in
    let a12 = -. 0.384 in
    let a21 = 0.0 in
    let asingle = Mat.of_array [|[|a11; a12|]; [|a21; a22|]|] in
    U.block_matrix [| [|asingle|]; [|asingle|] |]

  (* noise mix between components within a daughter *)
  let dGm_mix = Mat.of_diag
      (Vec.of_list [ sigma1; sigma2])
  (* noise mix between components across daughters *)
  let dGm_xmix = Mat.of_diag
      (Vec.of_list [ beta *. sigma1; beta *. sigma2])
  (* assemble *)
  let cdGm =
    let b = U.block_matrix
      [|[|dGm_mix; dGm_xmix|];
        [|dGm_xmix; dGm_mix|]|] in
    let bbt = gemm ~transb:`T b b in
    bbt
end


(*now let's generate some Gaussian random numbers*)

(*keep the rng outside so it does not get reseeded every time*)
let the_rng = Gsl.Rng.(make MT19937)

module Gauss (Pars:PARS) = struct
  open Pars

  (* some validation *)
  let () =
    assert (Vec.dim alpha = hdim);
    assert (Mat.dim1 meanScale = 2 * hdim);
    assert (Mat.dim2 meanScale =  hdim);
    assert (Mat.dim1 cdGm = 2 * hdim);
    assert (Mat.dim2 cdGm = 2 * hdim)

  let draw_rho_scale cdgm meanscale =
    (* draws Gaussian samples from daughter covar and inheritance matrix *)
    (* the decomp for sampling *)
    let cdGmEvs, rot = U.evs_rot cdgm in
    (* in rare cases we get EVs a bit smaller than 0. Correct for that: *)
    let cdGmEvs = Vec.relu cdGmEvs in
    let sigmas = Vec.sqrt cdGmEvs in
    begin try
        Stat.check_nan_v sigmas
      with Invalid_argument _ ->
        raise @@ Invalid_argument
          (Format.sprintf "cov EVs: %a@." Lacaml.Io.pp_fvec cdGmEvs); end;
    fun mtr_vec ->
      let mean = gemv meanscale mtr_vec in
      U.draw the_rng mean sigmas rot

  let draw_dtrs =
    (* draws Gaussian samples for both daughters in a tall vector *)
    draw_rho_scale cdGm meanScale

  let draw_dtr =
    let cdtr = lacpy ~n:hdim ~m:hdim cdGm in
    let meanscale = lacpy ~n:hdim ~m:hdim meanScale in
    draw_rho_scale cdtr meanscale

  let draw_initial =
    let statEvs, rot =
      let cmm = U.markov_limit_cov meanScale cdGm in
      U.evs_rot cmm in
    let sigmas = Vec.sqrt statEvs in
    let mean = Vec.make0 (Mat.dim2 meanScale) in
    fun () ->
      U.draw the_rng mean sigmas rot

end


(* type of an observed cell cycle time. can be censored -- then the minimum
 * cycle time Min τ is known. otherwise, the value Val τ is known. *)
type tautype = Min of float | Val of float
(*[@@ deriving eq]*)

(* equality *)
let equal_tautype t1 t2 = match t1, t2 with
  | Min f1, Min f2 -> (f1 =. f2)
  | Val f1, Val f2 -> (f1 =. f2)
  | Min _, Val _
  | Val _ , Min _ -> false

(* cell record: observables for a single cell *)
type cell_rec = {tau: tautype; tbirth: float; id: int}

(* pretty printing for cycle times and cells *)
let print_tau = fun fmt tau ->
  let (taustr, tauv) = match tau with
    | Min v -> "τ>", v
    | Val v -> "τ=", v in
  Format.fprintf fmt "%s%.4g" taustr tauv;;
let print_cr = fun fmt {tau; tbirth; id} ->
  let (taustr, tauv) = match tau with
    | Min v -> "τ>", v
    | Val v -> "τ=", v in
  Format.fprintf fmt "%s%.4g tb=%.4g id:%d" taustr tauv tbirth id;;


(* pre-transformation of cycle times to bring them to a desired distribution shape *)
module Transform = struct
  type t = {to_time: float -> float as 'a; from_time: 'a}

  let id_transform = {to_time=Fun.id; from_time=Fun.id}

  (** construct the log transformation from the given mean and var of the
   * logged values; this yields transformed values with mean 0 and unit
   * variance. whether they are truly normal depends on whether the original
   * samples were truly lognormal...*)
  let log_transform (logmean, logvar) =
    let from_time t = (log t -. logmean) /. sqrt logvar in
    let to_time lt = exp (lt *. sqrt logvar +. logmean) in
    {to_time; from_time}

  (** construct the normalizing transform from a given array of (time) values.
   * I.e. the transform t -> tg = C_g(C_e^{-1}(t)) where C_e is the empirical CDF
   * and C_g is the CDF or a standard normal distribution *)
  let normal_transform ?(margins=`Symmetric) ar =
    let sar = Array.sorted Float.compare ar in
    let nar = Stat.normal_transformed_vals ~margins sar in
    let from_time t = U.interpolate_linear t sar nar in
    let to_time nt = U.interpolate_linear nt nar sar in
    {to_time; from_time}

end


(* linear combination of hidden variables to get the (normalized) cycle time *)
module Observe (Pars:PARS) = struct
  let obs vec =
    let proj = dot Pars.alpha vec in
    proj +. U.centered_gauss the_rng Pars.noise

end


(* generate lineage trees with BAR model distribution given by Pars. *)
module Sample (Pars:PARS) = struct

  module Pars = Pars                                               (* public *)
  module Gauss = Gauss (Pars)
  module Observe = Observe (Pars)

  let sample_tree_to ?(child_num=2) ?(synchronize=[]) n x =
    (* height n counts from the current zipper position. synchronize is a quick
       hack to allow one of the two variables to be perfectly correlated
       between sisters *)
    (* only independent or identical children are possible *)
    let tr = T.create x in
    let rec sample n tr = match (n, tr) with
      | 0, _ -> tr                                                   (* done *)
      | _, `Node (x, childlist) when List.length childlist >= child_num ->
        (* all children present; descend *)
        `Node (x, List.map (sample (n-1)) childlist)
      | _, `Node (x, childlist) ->
        (* fill up children by one more *)
        let dx = Gauss.draw_dtr x in
        let () = match (synchronize, childlist) with
          | [], _ | _, [] -> ()
          | synlist, h::_ ->
            match h with `Node (hdx, _) ->
              List.iter (fun xi -> dx.{xi} <- hdx.{xi}) synlist in
        sample n @@ `Node (x, (`Node (dx, [])::childlist)) in
    sample n tr

  let sample_to ?(child_num=2) ?(synchronize=[]) n zp =
    let x = Z.zval zp in
    (* insert back; should be at the same position *)
    let filltr = sample_tree_to ~child_num ~synchronize n x in
    Z.replace filltr zp

  let sample_tree n x =
    let tree = T.create x in
    assert Pars.(Mat.dim1 meanScale / Mat.dim2 meanScale = 2);
    let rec sample n tr = match (n, tr) with
      | 0, _ -> tr
        (* return *)
      | _, `Node (x, (_::_::_ as l)) ->
        (* all children present; descend *)
        `Node (x, List.map (sample (n-1)) l)
      | _, `Node (x, []) ->
        (* fill up with 2 children *)
        let dd = Gauss.draw_dtrs x in
        Stat.check_nan_v dd;
        let dx1, dx2 = U.split2 dd in
        Stat.check_nan_v dx1; Stat.check_nan_v dx2;
        sample n @@ `Node (x, [`Node (dx1, []); `Node (dx2, [])])
      | _ -> raise (Failure "invalid tree geometry encountered") in
    sample n tree

  let sample_joint n zp =
    (* height n counts from the current zipper position. *)
    let start = Z.zval zp in         (* cut out from current zipper position *)
    let filltr = sample_tree n start in
    Z.replace filltr zp    (* insert back; we should be at the same position *)

  let initial_zp () =
    let x = Gauss.draw_initial () in
    Z.zipper (T.create x)

  let corcoeff ar1 ar2 = (Stat.cor 1 ar1 ar2).(0)

  let obs_corr = Z.pair_corr corcoeff Observe.obs

  let x1_corr = Z.pair_corr corcoeff (fun x -> x.{1})
  let x2_corr = Z.pair_corr corcoeff (fun x -> x.{2})
  let x12_cor = Z.pair_xcorr corcoeff
      (fun (x:Vec.t) -> x.{1}) (fun (x:Vec.t) -> x.{2})
  let x21_cor = Z.pair_xcorr corcoeff
      (fun (x:Vec.t) -> x.{2}) (fun (x:Vec.t) -> x.{1})

end


(* cell-specific lineage tree manipulation functions. extends the generic tree
 * functions in module Tree. *)
module Ctree = struct

  let build_cell_tree =
    let open Transform in
    fun ?(transform=id_transform) ?(start_t=0.) obs xtr ->
      let r_id = ref ~-1 in
      let f child_x parent = match parent with
        | {tbirth=tbirth_parent; tau=Val tau_parent; _} ->
          let tbirth = tbirth_parent +. tau_parent in
          let tau = Val (transform.to_time (obs child_x)) in
          let id = incr r_id; !r_id in
          (* attention: child order reversal *)
          Some { tbirth; tau; id }
        | {tau=Min _; _} ->
          None in
      T.fold_filter_tree ~f { tau=Val 0.; tbirth=start_t; id=0 } xtr

  let record_cell_tree ?(max_T=infinity) ?walk_out ctr =
    let disappear_at tbirth = match walk_out with
      | None -> false
      | Some fn -> let prob = fn tbirth in
        (Gsl.Rng.uniform the_rng <. prob) in
    let f child_cr parent_cr = match parent_cr with
      | {tau=Val _; _} ->
        begin match child_cr with
          | {tbirth; _}
            when tbirth >. max_T || disappear_at tbirth ->
            None
          | {tbirth; tau = Val t; _}
          | {tbirth; tau = Min t; _} ->
            if (tbirth +. t <. max_T)
            then Some child_cr
            else Some {child_cr with tau = Min (max_T -. tbirth)} end
      | {tau=Min _; _} ->
        None in
    T.reduce_filter_tree ~f ctr

  let extract_tau_tree =
    let open Transform in
    fun ?(transform=id_transform) (ctr: cell_rec T.t) : tautype T.t ->
      let f = function
        | {tau = Val t; _} -> Val (transform.from_time t)
        | {tau = Min t; _} -> Min (transform.from_time t) in
      T.map f ctr

  let inner_tau =
    fun tau -> match tau with Min _ -> false | Val _ -> true

  let get_inner = function
    | Min _ -> failwith "no survivors?!"
    | Val t -> t

  (* check for cells that have nan as their cycle time. this should happen only
   * at the root! *)
  let prune_survivors_taus tau_zip =
    let inner_cell z =
      let c = Z.zval z in
      if equal_tautype c (Val Float.nan) && Z.depth z > 0
      then failwith "nan at non-root";
      inner_tau c in
    Z.take_while inner_cell tau_zip

  let prune_survivors_cells cell_zip =
    let inner_cell z =
      let c = (Z.zval z).tau in
      if equal_tautype c (Val Float.nan) && Z.depth z > 0
      then failwith "nan at non-root";
      inner_tau c in
    Z.take_while inner_cell cell_zip

  (* erase all extra info and remove survivor cells *)
  let prune_tau_survivors ?(transform=Transform.id_transform) tau_zip =
      prune_survivors_taus tau_zip
      >|= fun zp ->
      Z.fork zp |> Z.map get_inner |> Z.map transform.Transform.from_time

  let prune_survivors ?transform cell_zip =
    let tau_zip = Z.tree cell_zip |> extract_tau_tree |> Z.zipper in
    prune_tau_survivors ?transform tau_zip

  let rev_list project zip =
    let rec fold_ordered acc children =
      let acc' = List.fold_left              (* the immediate children times *)
          (fun ac (`Node (cargo, _)) -> project cargo::ac)
          acc children in
      List.fold_left                                  (* now the descendants *)
        (fun ac (`Node (_, grandch)) -> fold_ordered ac grandch)
        acc' children in
    fold_ordered [] [Z.tree zip]

  (* for already pruned zippers *)
  let tau_vals float_zip = float_zip
    |> rev_list (fun x:float -> x)
    |> List.rev
    |> List.tl     (* matrices are conditioned on first entry; leave it out! *)
    |> Vec.of_list

  (* printing *)
  let prvect = T.print Lacaml.Io.pp_rfvec;;
  let prvecz : Vec.t Z.t T.printer = fun fmt z ->
    let t = Z.tree z in
    Format.fprintf fmt "%a" prvect t
  let prf = T.print U.pf
  let prfz : float Z.t T.printer = fun fmt f ->
    let t = Z.tree f in
    Format.fprintf fmt "%a" prf t
  let prtaut = T.print print_tau
  let prtauz : tautype Z.t T.printer = fun fmt f ->
    let t = Z.tree f in
    Format.fprintf fmt "%a" prtaut t
  let pru = T.print U.pu
  let pruz : unit Z.t T.printer = fun fmt f ->
    let t = Z.tree f in
    Format.fprintf fmt "%a" pru t

  (* printing in a mma compatible nested list format *)
  let rec to_mma print_val tree =
    let pv v = print_val v in
    match tree with
    | `Node (a, []) -> pv a
    | `Node (a, l) ->
      let chs = String.concat ", " (List.map (to_mma print_val) l) in
      Format.sprintf "{%s, {%s}}" (pv a) chs

end


(* correlations between related cells *)

let default_relatives =
  ["md"; "ss"; "cc"; "gmgd"; "cc2"; "an"]

(* string map *)
module SM = Map.Make (String)

let relpairs = SM.of_list
    (List.combine
       default_relatives
       Z.([md_pairs; ss_pairs; cc_pairs; gmgd_pairs; cc2_pairs; an_pairs]))

let make_corrs height ?synchronize ?(rels=default_relatives)
    (module Params:PARS) =
  let zz_top = `Node (Vec.of_list [0.; 0.], []) |> Z.zipper in
  let module Sample = Sample (Params) in

  let rez = Sample.sample_to ?synchronize height zz_top in
  let pairs =
    List.map2
      (fun name pairfun -> (name, pairfun rez))
      rels
      Option.(get_exn @@ sequence_l @@
              List.map (fun k -> SM.get k relpairs) rels)
    |> SM.of_list in
  let ocorrs = SM.map Sample.obs_corr pairs in
  let x1corrs = SM.map Sample.x1_corr pairs in
  let x2corrs = SM.map Sample.x2_corr pairs in
  rez, ocorrs, x1corrs, x2corrs

let rmsdev_corrs cd1 cd2 =
  let sqsum = SM.keys cd1
              |> Sequence.filter_map
                (fun k -> match SM.(get k cd1, get k cd2) with
                   | Some v1, Some v2 -> Some ((v1 -. v2)**2.)
                   | _, _ -> None)
              |> Sequence.fold (+.) 0. in
  sqrt @@ sqsum /. float_of_int (SM.cardinal cd1)


(* the core of the evidence calculation: the module for the evaluation of a
 * likelihood function.
 * for formulas and a prototype see the mathematica notebook recursive-likelihood-new.nb *)
module Likely  = struct
  (*[@@@landmark "auto"]*)

  (* projections *)
  let p1 alpha =
    let proj = Mat.of_col_vecs [|alpha|] in
    Mat.scal (1. /. U.norm alpha) proj; proj

  let p12 alpha =
    let p = p1 alpha in
    U.stack_diag p p

  (* currently not needed *)
  let zero_filled (ml1, ml2, cil1, cil2) =
    match ml1, ml2 with
    | ml1, ml2 when U.is_empty ml2 ->
      (ml1, U.zeros_like ml1, cil1, U.zeros_like cil1,
       Mat.(1, dim1 ml1))
    | ml1, ml2 when U.is_empty ml1 ->
      U.zeros_like ml2, ml2, U.zeros_like cil2, cil2,
      Mat.(dim1 ml2 + 1, 2 * dim1 ml2)
    | _, _ ->
      (ml1, ml2, cil1, cil2,
       Mat.(1, dim1 ml1 + dim1 ml2))

  (* directly calculate small inverses *)
  let stiffness_1 ci12 =
    let d = Mat.dim2 ci12 / 2 in
    let c12 = U.ge_inv ci12 in
    let c1 = lacpy ~m:d ~n:d c12 in
    U.ge_inv c1

  (* rectangular matrices *)
  let tall_c_1_2 c12 =
    let d = Mat.dim2 c12 / 2 in
    let c12_1 = lacpy ~n:d c12 in
    let c12_2 = lacpy ~ac:(d+1) c12 in
    c12_1, c12_2

  (* M^{<} construction *)

  type 'a pair = 'a * 'a

  type children =
    | One of mat pair option
    | Two of mat pair option pair

  let ml (a, alpha, _ci12) =
    let na = U.norm alpha in
    let p = p1 alpha in
    let upper_row = gemm ~alpha:na ~transa:`T p a in
    function
    (* two children *)
    | Two (Some (ml1, _), Some (ml2, _)) ->
      U.block_matrix [| [|upper_row|];
                        [|upper_row|];
                        [|gemm ml1 a|];
                        [|gemm ml2 a|] |]
    | Two (Some (ml, _), None)
    | Two (None, Some (ml, _)) ->
      U.block_matrix [| [|upper_row|];
                        [|upper_row|];
                        [|gemm ml a|] |]
    | Two (None, None) ->
      U.stack_high upper_row upper_row
    (*one child*)
    | One Some (ml, _) ->
      U.stack_high upper_row (gemm ml a)
    | One None ->
      upper_row


  (* C^{<} construction *)

  let trans = `T                                              (* convenience *)

  (* first,  Δ'_{22}^{-1} (in the notation of recursive-likelihood-new.nb) *)

  let delta'_22_inv (_a, _alpha, ci12) =
    (* note: for empty children this makes no sense -> raise exception early *)
    let ci1 = stiffness_1 ci12 in
    let eval ci ml cil =
      let pi = U.cong ~trans ml cil
               |> Mat.add ci
               |> U.ge_inv in                               (* small inverse *)
      let subm = U.cong (gemm cil ml) pi in
      Mat.sub cil subm in
    let invalid = Invalid_argument "delta not defined for empty children" in
    function
    | Two (Some (ml1, cil1), Some (ml2, cil2)) ->
      let ml12 = U.stack_diag ml1 ml2 in
      let cil12 = U.stack_diag cil1 cil2 in
      (eval ci12 ml12 cil12)
    | Two (Some (ml, cil), None)
    | Two (None, Some (ml, cil))
    | One Some (ml, cil) ->
      eval ci1 ml cil
    | Two (None, None)
    | One None ->
      raise invalid


  (* Γ' construction -- the inverse of Δ' *)

  let gamma' (_a, alpha, ci12) =
    (* purely parameter-dependent pre-calculation *)
    let c1 = stiffness_1 ci12 |> U.ge_inv in
    let c12 = ci12 |> U.ge_inv in
    let c12_1, c12_2 = tall_c_1_2 c12 in
    let p1 = p1 alpha in
    let p12 = p12 alpha in
    let d'22i = delta'_22_inv (_a, alpha, ci12) in
    (* core construction with descendants *)
    let eval_no_desc c p =
      let g11 = U.cong ~trans p c |> U.ge_inv in
      g11 in
    let eval c cpart p d'i ml =
      let g11 =
        let cm = gemm cpart ~transb:`T ml in
        let diff = Mat.sub c (U.cong cm d'i) in           (* d'i used later! *)
        U.cong ~trans p diff |> U.ge_inv in
      let auxm = gemm d'i @@ gemm ml @@ gemm ~transa:`T cpart p in
      let g21 = Mat.neg @@ gemm auxm g11 in
      let g12 = Mat.transpose_copy g21 in
      let g22 = Mat.add d'i (U.cong auxm g11) in
      U.block_matrix [|[|g11; g12|];
                       [|g21; g22|]|] in
    function
    | Two (Some (ml1, _), Some (ml2, _)) as two ->
      let d'i12 = d'22i two in
      let ml12 = U.stack_diag ml1 ml2 in
      eval c12 c12 p12 d'i12 ml12
    | Two (Some (ml1, _), None) as first ->
      let d'i = d'22i first in                           (* dimensions as c1 *)
      eval c12 c12_1 p12 d'i ml1
    | Two (None, Some (ml2, _)) as second ->
      let d'i = d'22i second in
      eval c12 c12_2 p12 d'i ml2
    | Two (None, None) ->
      eval_no_desc c12 p12
    | One Some (ml, _) as one ->
      let d'i = d'22i one in
      eval c1 c1 p1 d'i ml
    | One None ->
      eval_no_desc c1 p1


  (* finally, assemble into C^{-1}^{<} *)

  let cil (_a, alpha, ci12) =
    let make_g' = gamma' (_a, alpha, ci12) in
    let nai = 1. /. U.norm alpha in
    fun children ->
      let g' = make_g' children in
      let scale_k = match children with One _ -> 1 | Two _ -> 2 in
      Mat.scal ~m:scale_k nai g';
      Mat.scal ~n:scale_k nai g';
      g'

  let mcil_direct (a, alpha, ci12) =
    let ml = ml (a, alpha, ci12) in
    let cil = cil (a, alpha, ci12) in
    let rec get_mcil z : Mat.t pair option =
      match Gen.to_list (Z.children z) with
      | [] -> None
      | [cz] ->
        let mcil_c = get_mcil cz in                          (* not tail rec *)
        Some (ml @@ One mcil_c, cil @@ One mcil_c)
      | [c1z; c2z] ->
        let mcil_c12 = (get_mcil c1z, get_mcil c2z) in
        Some (ml @@ Two mcil_c12, cil @@ Two mcil_c12)
      | _longer -> raise (Invalid_argument "more than two children!") in
    fun z -> get_mcil z

  (* try to gain some speed through memoizing mcil on the matrices and tree
   * structure. *)

  let skeleton z = Z.map (fun _ -> ()) z

  let mcil = mcil_direct

  (* extract overhanging indices, overhanging cells, internal cells. to be used
   * with U.block_minors
   *
   * convention: overhanging cells come first. the order of tau_vals is the
   * reference for the permutation. (i.e. the reverse order of rev_list) within
   * each of the two sets the ordering is stable. the permutation is given such
   * that
   * U.permute_vector permutation taus ==
   * Vec.of_list (List.rev (Ctree.rev_list Fun.id tau_zip))
   * *)
  let permutation_offs_taus ?(skip_root=true) tau_zip =
    let overhanging, internal =
      tau_zip
      |> Ctree.rev_list Fun.id |> List.rev
      |> (if skip_root then List.tl else Fun.id)  (* if yes, =Ctree.tau_vals *)
      |> List.mapi (fun i tau -> (i+1, tau))        (* index the list from 1 *)
      |> List.partition_map (function
          | (i, Min t) -> `Left (i, t)
          | (i, Val t) -> `Right (i, t)) in
    let offset = List.length overhanging in
    let len = offset + List.length internal in
    let permutation = Array.make len 999 in            (* to be overwritten! *)
    let taus = Vec.create len in
    let blit offs pairs =
      pairs |> List.iteri (fun i (ind, tau) ->
          permutation.(offs+i) <- ind;
          taus.{1+offs+i} <- tau) in
    blit 0 overhanging;
    blit offset internal;
    permutation, offset, taus  (* permutation is 1-based, offset starts at 0 *)

  (* as in: log p = - log Z - 1/2 x.ci.x.
   * here log Z is computed by logZ ci from the stiffness matrix ci. *)
  let log_Z =
    let log2pi = log (acos 0. *. 4.) in
    fun ci ->                                  (* ci is the stiffness matrix *)
      let n = float_of_int (Mat.dim1 ci) in
      (* log det of the covariance *)
      (* this should be correct: the failures happen when det ci ~ 0. *)
      let logdetc = try -. U.log_det ci with
        | Failure s
          when (String.prefix ~pre:"Lacaml.D.potrf: leading minor" s) ->
          +. infinity in
      1./.2. *. (n *. log2pi +. logdetc)

  (* evaluate (-log L) since we want to minimize not maximize *)
  let evaluate_pruned float_zip (a, alpha, ci12) =
    let tv = Ctree.tau_vals float_zip in
    match mcil (a, alpha, ci12) float_zip with
    | None -> None
    | Some (ml, cil) ->
      Some begin fun root_x ->
        let xc = copy root_x in         (* extra safety; not sure if needed. *)
        let dev = Vec.sub tv (gemv ml xc) in
        let quad = 0.5 *. U.vmv dev cil dev in
        quad +. log_Z cil end

  (* we want also a neg. likelihood for the root vector *)
  let evaluate_root (a, _alpha, ci12) =
    (* this is the chain-stationary negative log likelihood for the root
     * vector, not conditioned on anything. in an experiment, the root time is
     * known to be longer than the first observed division time; this requires
     * a bit of work *)
    let cmmi = U.(ci12 |> ge_inv |> markov_limit_cov a |> ge_inv) in
    let lz =
      try
        log_Z cmmi
      with
        Failure s ->
        print_endline s;
        Format.printf "a\n%a\nci12\n%a\ncmmi\n%a\n"
          Lacaml.Io.pp_fmat a
          Lacaml.Io.pp_fmat ci12
          Lacaml.Io.pp_fmat cmmi;
        Format.print_flush ();
        raise Matrix_problem
    in
    fun rootx -> 0.5 *. U.vmv rootx cmmi rootx +. lz

  (* find the root vector that maximizes the joint probability p(τ^<|x)p_c(x)
   * where p_c is the root probability calculated in evaluate_root.
   * formula:
   * x^* = S^{-1}.ML^T.CIL.τ^<
   * where S = ML^T.CIL.ML + CI_c
   * and CI_c is the inv covariance of p_c
   *
   * this neglects the extra info that the root time was higher than the
   * observed root lifetime. deal with it.
   * *)
  let xrootmin tv (ml, cil) (a, _alpha, ci12) =
    let cic = U.(markov_limit_cov a (ge_inv ci12) |> ge_inv) in
    let c = U.cong ~trans ml cil
            |> Mat.add cic
            |> U.ge_inv in
    gemv c @@ gemv ~trans ml @@ gemv cil @@ tv

  let evaluate_pruned_rootmin float_zip (a, alpha, ci12) =
    let tv = Ctree.tau_vals float_zip in
    match mcil (a, alpha, ci12) float_zip with    (* None iff no children *)
    | None -> None
    | Some (ml, cil) ->
      let xroot = xrootmin tv (ml, cil) (a, alpha, ci12) in
      let dev = Vec.sub tv (gemv ml xroot) in
      let quad = 0.5 *. U.vmv dev cil dev in
      let neglogl_root = evaluate_root (a, alpha, ci12) xroot in
      Some (quad +. log_Z cil +. neglogl_root)

  (* evaluate likelihood including integral over the overhanging cells. the
   * root cell is still conditioned on *)

  let evaluate_full_root tau_zip (a, alpha, ci12) =
    match mcil (a, alpha, ci12) tau_zip with
    | None -> None
    | Some (ml, cil) ->
      match permutation_offs_taus tau_zip with
      | _identity, 0, taus ->
        Some begin fun root_x ->
          let dev = Vec.sub taus (gemv ml root_x) in
          let quad = 0.5 *. U.vmv dev cil dev in
          quad +. log_Z cil end
      | permutation, offset, taus ->
        (* ATTENTION: taus are already sorted into overhanging, internal.
         * essentially, taus = permute_vector permute (tauvals tau_zip) *)
        let axx, byy, mxy = U.block_minors cil permutation offset in
        let logZyy = log_Z byy in
        Some begin fun root_x ->
          let mean =
            gemv ml root_x |> U.permute_vector permutation in
          let dev = Vec.sub taus mean in
          let devx = copy ~n:offset dev in
          let devy = copy ~ofsx:(offset+1) dev in
          let quady = 0.5 *. U.vmv devy byy devy in
          let lower = Vec.sub devx (gemv mxy devy) in
          let intxx = -. Integrate_mvnd.log_gauss_integral ~lower axx in
          quady +. logZyy +. intxx end




  (* now combine this with integrating over the root *)

  (* formula. see Stauα in test-root-integration.nb.
   *
   * starting with stiffness cil and mean ml.x where x is the hidden
   * root vector, make an augmented stiffness cilr which corresponds to
   * the joint distribution of τ_0, the time of the root cell, and the times of
   * the descendant cells. this takes as input the assumed stationary
   * distribution of the root cell hidden variables, given by ciroot.
   * this will be taken as the markov chain limit -- a (hopefully small)
   * inaccurracy. *)

  (* _verified_: works as in test-root-integration.nb ! *)

  let ci_prepend_root alpha ci_root (ml, cil) =
    let p1 = p1 alpha in                         (* p1, q1 are both portrait *)
    let q1 = lacpy ~ac:2 (U.fill_basis alpha) in
    let ci_ti_root = Mat.add ci_root (U.cong ~trans ml cil) in
    let c_bar_alpha = U.(cong ~trans q1 ci_ti_root |> ge_inv) in
    let ci_ll =                             (* lower times (i.e. not mother) *)
      Mat.sub cil U.(c_bar_alpha |> cong q1 |> cong ml |> cong cil) in
    let ci_lr = gemm cil @@ gemm ml             (* coupling: lower with root *)
        Mat.(sub (gemm (gemm U.(cong q1 c_bar_alpha) ci_ti_root) p1) p1) in
    let ci_rr = U.cong ~trans p1                                     (* root *)
        (Mat.sub ci_ti_root U.(c_bar_alpha |> cong q1 |> cong ci_ti_root)) in
    let ci_rl = Mat.transpose_copy ci_lr in
    (* combine *)
    let cilr = U.block_matrix [|[|ci_rr; ci_rl|];
                                [|ci_lr; ci_ll|]|] in
    (* rescale to make a stiffness for the actual non-normalized root time *)
    let n_alpha_i = 1. /. U.norm alpha in
    Mat.scal ~n:1 ~ar:1 n_alpha_i cilr;
    Mat.scal ~m:1 ~ac:1 n_alpha_i cilr;
    cilr

  (* validated: this matrix, inverted and with the first row and column
   * removed, is the inverse of ciavg below *)
  let cifull (a, alpha, ci12) (ml, cil) =
    (* using stationary root cell distribution here. *)
    let ci_root = U.(ci12 |> ge_inv |> markov_limit_cov a |> ge_inv) in
    ci_prepend_root alpha ci_root (ml, cil)

  (* evaluate neg log likelihood including integral over the overhanging cells
   * AND integral over the root cell. to do this, we first get the augmented
   * stiffness matrix. *)

  let evaluate_full ?(clip=1e4) tau_zip (a, alpha, ci12) =
    match mcil (a, alpha, ci12) tau_zip with
    | None -> None
    | Some (ml, cil) ->
      (* get full matrix incl root,
       * no conditioning and no ml required in the sequel *)
      let cifull = cifull (a, alpha, ci12) (ml, cil) in
      match permutation_offs_taus ~skip_root:false tau_zip with
      | _identity, 0, taus ->
        let dev = taus in
        let quad = 0.5 *. U.vmv dev cifull dev in
        Some (quad +. log_Z cifull)
      | permute, offset, taus ->              (* now including the root cell *)
        let dev = taus in
        let devx = copy ~n:offset dev in
        let devy = copy ~ofsx:(offset+1) dev in
        let axx, byy, mxy = U.block_minors cifull permute offset in
        let logZyy = log_Z byy in
        let quady = 0.5 *. U.vmv devy byy devy in
        let lower = Vec.sub devx (gemv mxy devy) in      (* mxy contains a - *)
        let intxx = Float.min
            (* call to the FORTRAN multivariate normal integral by genz *)
            clip (-. Integrate_mvnd.log_gauss_integral ~lower axx) in
        Some (quady +. logZyy +. intxx)


(* new alternative: stiffness matrix for all non-root cells where root cell
 * is averaged with the chain limit. this is an alternative to minimizing the
 * root cell latent variables (i.e. fitting the root cell)
 * *)

  (* validated: this is the same as one gets for sampling with a root
   * distribution drawn from the stationary distribution. *)
  let ciavg (a, _alpha, ci12) (ml, cil) =
    (* using stationary root cell distribution to get the averaged stiffness
     * for the descendant cells *)
    let ciroot = U.(ci12 |> ge_inv |> markov_limit_cov a |> ge_inv) in
    (* the woodbury formula for (CL + ML.C∞.ML^t)^{-1} *)
    let inner_ci = U.(ge_inv (Mat.add ciroot (cong ~trans ml cil))) in
    (* cong works since cil is symmetric *)
    Mat.sub cil U.(cong (gemm cil ml) inner_ci)

  let evaluate_pruned_rootavg float_zip (a, alpha, ci12) =
    let tv = Ctree.tau_vals float_zip in
    match mcil (a, alpha, ci12) float_zip with    (* None iff no children *)
    | None -> None
    | Some (ml, cil) ->
      let ciavg = ciavg (a, alpha, ci12) (ml, cil) in
      let dev = tv in
      let quad = 0.5 *. U.vmv dev ciavg dev in
      let logz = log_Z ciavg in
      (* debugging *)
      let res = (quad +. logz) in
      match Float.classify res with
      | FP_normal | FP_subnormal | FP_zero ->
        Some res
      | FP_infinite when (res >. 0.) ->
        (*failwith "see if this gets through";*)
        Some res
      | FP_infinite | FP_nan ->
        failwith
          (Format.sprintf
             "nan or inf log likelihood. quad %f log %f@." quad logz)


  let evaluate_full_avg tau_zip (a, alpha, ci12) =
    match mcil (a, alpha, ci12) tau_zip with
    | None -> None
    | Some (ml, cil) ->
      let ciavg = ciavg (a, alpha, ci12) (ml, cil) in
      match permutation_offs_taus tau_zip with (* this skips the root! *)
      | _identity, 0, taus ->
        let dev = taus in
        let quad = 0.5 *. U.vmv dev ciavg dev in
        Some (quad +. log_Z ciavg)
      | permute, offset, taus when
          (offset = Array.length permute) ->
        let axx = U.permute_matrix permute ciavg in
        let devx = taus in
        let lower = devx in
        let intxx = -. Integrate_mvnd.log_gauss_integral ~lower axx in
        Some intxx
      | permute, offset, taus ->
        (* ATTENTION: taus are already sorted into overhanging, internal.
         * essentially, taus = permute_vector permute (tauvals tau_zip) *)
        let axx, byy, mxy = U.block_minors ciavg permute offset in
        let logZyy = log_Z byy in
        let dev = taus in
        let devx = copy ~n:offset dev in
        let devy = copy ~ofsx:(offset+1) dev in
        let quady = 0.5 *. U.vmv devy byy devy in
        let lower = Vec.sub devx (gemv mxy devy) in
        let intxx = -. Integrate_mvnd.log_gauss_integral ~lower axx in
        Some (quady +. logZyy +. intxx)

end

