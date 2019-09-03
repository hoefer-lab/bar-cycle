(* a utility module for reusable helper functions that depend on Tree and/or
 * Models *)

open! Containers
open Lacaml.D

open Util
open Hmtree
open Tree

module Min = Gsl.Multimin.NoDeriv
open Option.Infix

open Models


(* convenience *)

let opt_get = Option.get_exn

(* iteration over forests *)

(* total time. attention: no auto-rooting of the zipper! *)
let experiment_T zip =
  let open! Sequence in
  Z.to_seq zip
  |> filter_map (fun z ->
      let (r:cell_rec) = Z.zval z in match r.tau with
      | Val tau -> Some (r.tbirth +. tau)   (* should not contribute to max. *)
      | Min tm -> Some (r.tbirth +. tm))
  |> filter (fun f -> not @@ Float.is_nan f)         (* at the root *)
  |> max

(* pick out only the pairs with observed division time: no survivals or deaths.
 * attention: this does not do any filtering for generation etc. *)
let taupairs pairseq =
  let open Sequence in
  pairseq
  |> map (Pair.map_same (fun z -> (Z.zval z).tau))
  |> filter_map (function
      | (Min _, _) | (_, Min _) -> None
      | (Val t1, Val t2) -> if Float.(is_nan t1 || is_nan t2)
        then None else Some (t1, t2));;

let cells_alive_at ?(tolerance=0.1) tt zip =
  Z.leavez zip
  |> Sequence.filter (fun z ->
      let (v:cell_rec) = (Z.zval z) in
      match v.tau with
      | Val _ -> false
      | Min tmin -> v.tbirth +. tmin >=. tt -. tolerance)

(* correlation and n of a pair sequence *)
let cor_n_pair corf taupairseq =
  let ar = taupairseq |> Sequence.to_array in
  let a, b = Array.(map fst ar, map snd ar) in
  (corf 1 a b).(0), Array.length a

(* correlation of a pair sequence *)
let cor_pair corf taupairseq = fst @@ cor_n_pair corf taupairseq

(* coarse-grain bootstrap over families *)
let fam_bootstrap n stat_fun taupairseqlist =
  let seqar = Array.of_list taupairseqlist in
  Stat.seqarray_bs n stat_fun seqar

(* CHANGE: trees without contributions to a particular relatedness pair should
 * be excluded from the start so they do not introduce sparse-sampling noise to
 * the bootstrap interval estimate. *)
let forest_bootstrap_pure n stat_fun pure_pairfun forest =
  forest
  |> List.filter_map (fun z ->
      let p = pure_pairfun z in
      if Sequence.length p = 0 then None else Some p)
  |> fam_bootstrap n stat_fun

(* coarse-grain bootstrap over families of a forest *)
let forest_bootstrap n stat_fun pairfun forest =
  forest_bootstrap_pure n stat_fun (fun z -> taupairs (pairfun z)) forest

(* correlations taking all samples from given forest *)
let frst_full_pure stat_fun pure_pairfun forest =
  let open Sequence in
  forest |> List.to_seq
  |> map pure_pairfun
  |> concat
  |> stat_fun

(* correlations taking all samples from given forest *)
let frst_full stat_fun pairfun forest =
  frst_full_pure stat_fun (fun z -> taupairs (pairfun z)) forest

let cor_n = cor_n_pair Stat.cor
let rcor_n = cor_n_pair Stat.rank_cor
let logcor_n =
  let logcor i a b = Stat.cor i (Array.map log a) (Array.map log b) in
  cor_n_pair logcor
let covar_n = cor_n_pair Stat.covar

(* the correlation and rank correlation *)
let cor = cor_pair Stat.cor
let lcor = cor_pair Stat.log_cor
let rcor = cor_pair Stat.(rank_cor ~tied:true)
let covar = cor_pair Stat.covar

let genfrst_full cor_fun pairfun project12 forest =
  let open Sequence in
  forest
  |> of_list
  |> map pairfun |> concat
  |> map project12
  |> cor_fun

let mean vecs =
  let sum =
    Array.fold (fun v v' -> Vec.add v v') Vec.(make0 (dim vecs.(0))) vecs in
  scal (1. /. float_of_int (Array.length vecs)) sum; sum

(* name clash with function from Util! *)
let covar_array vecs1 vecs2 =
  let m1, m2 = mean vecs1, mean vecs2 in
  let vc1, vc2 =
    Array.(map (Vec.sub m1) vecs1, map (Vec.sub m2) vecs2) in
  let cva = Mat.create (Vec.dim m1) (Vec.dim m2) in
  let l = Array.(make (length vc1) 0.) in
  for i = 1 to Mat.dim1 cva do
    for j = 1 to Mat.dim2 cva do
      for k = 0 to Array.length l - 1 do
        l.(k) <- vc1.(k).{i} *. vc2.(k).{j} done;
      cva.{i,j} <- Stat.mean l done; done;
  cva


(* evaluation of arbitrary correlations *)

exception Pred_corr of
    (Mat.t * Vec.t * Mat.t) * Mat.t * Mat.t * float

(* use the old rule: go along the path from one cell to the other. at each
 * cell, insert the inverse stationary covariance. at each link between cells
 * insert either the mother-daughter, or sister-sister covariance.
 * this has to be the same as the new recursion... *)


let predicted_xcorr a ci12 v1 v2 path =                  (* stationary covar *)
  try
    let c12 = U.ge_inv ci12 in
    let c00 = U.markov_limit_cov a c12 in
    let s00 = U.ge_inv c00 in
    let cmd = U.get_covar_md a c12 in
    let cdm = Mat.transpose_copy cmd in
    let css =
      let aa = gemm (gemm a c00) ~transb:`T a in
      Mat.add aa ~bc:(Mat.dim1 a + 1) c12 in      (* off-diagonal block incl γ *)
    let covar = List.fold_left (fun acc step ->
        let acc' = gemm acc s00 in
        match step with
        | `Up -> gemm acc' cdm
        | `Over -> gemm acc' css
        | `Down -> gemm acc' cmd)
        c00 path in
    let computed_corr =
      let cv = U.vmv v1 covar v2 in
      if cv =. 0.
      then cv
      else
        let var1, var2 = U.(vmv v1 c00 v1, vmv v2 c00 v2) in
        cv /. sqrt var1 /. sqrt var2 in
    (* autocorrelations must be <= 1. only tests when passing the identical
     * vector. *)
    if Equal.physical v1 v2 && computed_corr >. 1. +. Float.epsilon then
      begin let open Format in let open Lacaml.Io in
        printf "computed correlation %f > 1 found\n" computed_corr;
        printf "a=\n%a\nv1, v1=\n%a, %a\nci12=%a\n"
          pp_fmat a pp_fvec v1 pp_fvec v2 pp_fmat ci12;
        printf "covar=\n%a\nc00=\n%a\n"
          pp_fmat covar pp_fmat c00; end;
    computed_corr
  with U.Diverging_recurrence ->
    nan

let predicted_corr (a, alpha, ci12) path =
  predicted_xcorr a ci12 alpha alpha path

let rel_paths = [
  (* convention: from first to second cell in the string. for cc12 go more
   * down than up. *)
  "ss"     , [`Over];
  "md"     , [`Down];
  "gmgd"   , [`Down;`Down];
  "an"     , [`Over;`Down];
  "cc"     , [`Up;`Over;`Down];
  "ggmggd" , [`Down;`Down;`Down];
  "gagn"   , [`Over;`Down;`Down];
  "cc1r"   , [`Up;`Over;`Down;`Down];
  "cc2"    , [`Up;`Up;`Over;`Down;`Down] ]

(* a function to predict the correlations.
 * build_mat may be:
 * * model.build_matrices to use the _unscaled_, (-∞, ∞) parameters; or
 * * model.build_matrices_scaled to use scaled parameters.
 * *)
(* checked with sampling corrs from simulated trees -- correct! *)
let predict_corrs build_mat params relatives =
  let mats = build_mat params in
  let tcor r = predicted_corr mats (List.assoc ~eq:Equal.string r rel_paths) in
  List.map tcor relatives

(* multidim integration over a handful of dimensions *)

(* non-adjustable, generic version *)
let gsl_integrate mc_kind rangevec n f =
  let lo, up = Array.(map fst rangevec, map snd rangevec) in
  Gsl.Monte.integrate mc_kind f ~lo ~up n the_rng


(* minimization!  *)

type it_pars = {mutable step_init_size:float; mutable maxit:int;
                mutable abs_tol:float; mutable succ_steps:int}

let run_min it_pars min_fun pv =
  let tol = it_pars.abs_tol in
  let l = Gslv.length pv in
  let x = Gslv.copy pv in
  let step_size = Gslv.create ~init:it_pars.step_init_size l in
  let nm_min = Min.(make NM_SIMPLEX l) in
  let minim = nm_min min_fun ~x ~step_size in
  let rec runit success_steps it =
    let suc_new =
      let () = Min.iterate minim in
      if Min.test_size minim tol then succ success_steps else 0 in
    match suc_new, it with
    | (_, i) when i >= it_pars.maxit ->
      Format.printf "stopping after %d iterations\n@." it_pars.maxit;
      Error minim
    | (sn, i) when sn >= it_pars.succ_steps ->
      Format.printf "tolerance reached in step %d\n@." i;
      Ok minim
    | (sn, i) ->
      runit sn (succ i) in
  runit 0 0

let get_anyway r = Result.catch ~ok:(fun x -> x) ~err:(fun x -> x) r

(* the following works just as well; above was to test if exceptions are not
 * handled properly.. *)

(*try*)
(*let success_steps = ref 0 in*)
(*for i = 1 to it_pars.maxit do*)
(*Min.iterate minim;*)
(*if Min.test_size minim tol*)
(*then incr success_steps*)
(*else success_steps := 0;*)
(*if !success_steps >= it_pars.succ_steps then raise (Converged i); done;*)
(*Format.printf "stopping after %d iterations\n@." it_pars.maxit;*)
(*Some minim*)
(*with Converged i ->*)
(*Format.printf "tolerance reached in step %d\n@." i;*)
(*Some minim*)

let print_res model fmin res_par =
  print_fmin_vec fmin res_par;
  let (a, alpha, cia) = model.build_matrices res_par in
  print_mats (a, alpha, cia)

let report_res model minimizer =
  let res_par = Gslv.create model.pdim in
  let fmin = Min.minimum ~x:res_par minimizer in
  print_res model fmin res_par;
  fmin, res_par

(* initial conditions *)

let list_product lists =
  let one_more acc l = List.(product cons) l acc in
  List.(fold_left one_more [[]] (rev lists))

let startpars_product bins low up dims =
  let r = Vec.linspace low up bins |> Vec.to_list in
  list_product (List.init dims (fun _ -> r))
  |> List.map Array.of_list
  |> List.map Gslv.of_array

(** latin hypercube sampling *)
let startpars ?(samples=21) ?(low= -3.) ?(up=3.) dims =
  let r = Vec.(linspace low up samples |> to_array) in
  let sampleshuffles =
    Array.(fun _ -> shuffle r; to_list r)
    |> List.init dims in
  let rec cubes acc sl = match List.hd sl with
    | [] -> acc
    | _ ->
      let newparcomb = List.map List.hd sl in
      let newsl = List.map List.tl sl in
      cubes (Gslv.of_array (Array.of_list newparcomb) :: acc) newsl in
  cubes [] sampleshuffles

(* parmap only up to tree size 6 unfortunately *)
type parallelization =
  (* avoid using the exact module names - confusion potential *)
    No | Use_Parmap
  (* Functory works less well; problems with parmap were due to accelerate *)
let cores = ref 3


let parallel_startpars = ref Use_Parmap
let () = parallel_startpars := No

let map_parallel how f l =
  let wrap_f x =
    try Ok (f x) with
    | exn_ -> Result.of_exn exn_ in
  match how with
  | Use_Parmap ->
    let () = Parmap.set_default_ncores !cores in
    let gc_config = Gc.get () in
    (* prevent heap compaction in the children to be safe; this could mess up
     * references to outside memory *)
    let () = Gc.set {gc_config with Gc.max_overhead = 1_000_000} in
    let res = Parmap.(parmap wrap_f @@ L l) in
    let () = Gc.set gc_config in
    res
  | No ->
    List.map wrap_f l

let pick_oks l =
  List.filter_map (function | Error _ -> None | Ok x -> Some x) l

let results ?(keep_stopped=false) it_pars startpars min_fun =
  let resokerr =
    startpars |> map_parallel !parallel_startpars
      begin fun pv ->
        let x = Gslv.copy pv in
        let get_min minimizer =
          let f = (Min.minimum ~x) minimizer in
          (f, x) in
        run_min it_pars min_fun pv
        |> Result.map2 get_min get_min end in
  let oks, lims, exceptions =
    List.fold_left
      (fun (ok, lim, ex) r -> match r with
         | Ok (Ok r)    -> (r::ok, lim, ex)
         | Ok (Error r) -> (ok, r::lim, ex)
         | Error exn_   -> (ok, lim, exn_::ex))
      ([], [], []) resokerr in
  let n, nstop, nex = List.(length resokerr, length lims, length exceptions) in
  Printf.printf "%d stopped iterations (%.0f%%)\n"
    nstop (float_of_int nstop /. float_of_int n *. 100.);
  Printf.printf "%d exceptions (%.0f%%)\n"
    nex (float_of_int nex /. float_of_int n *. 100.);
  if keep_stopped then (oks @ lims, exceptions) else (oks, exceptions)

let comp1f x1 x2 = Float.compare (fst x1) (fst x2)

let make_sl results =
  results |> fst
  |> List.sort comp1f
  |> Array.of_list


(* this is generally a bad idea since evaluation of the function for one tree
 * does not take long. (except maybe for the integration case) *)
let parallelize_forest = ref No

(* this function silently drops exceptions raised in the workers. in this way,
 * parameter sets that lead to exceptions silently get lower results. *)
let forest_eval_lenient ~x eval tau_forest root_x =
  let f z = eval ~x z root_x in
  map_parallel !parallelize_forest f tau_forest
  |> pick_oks                                           (* exeptions dropped *)
  |> List.fold_left
    (fun x -> function None -> x | Some v -> x +. v ) 0.

(* this function raises all exceptions and chokes on empty trees *)
let forest_eval_strict ~x eval tau_forest root_x =
  let f z = eval ~x z root_x in
  List.map f tau_forest
  |> List.map Option.get_exn
  |> List.fold_left (+.) 0.

(* ATTENTION: this always uses a constant root_x for the whole forest! *)
let forest_eval = forest_eval_strict

let find_min_forest it_pars eval tau_forest root_x xa =
  let min_fun = forest_eval eval tau_forest root_x in
  run_min it_pars min_fun xa

(* evaluation with given root hidden vector *)
let eval_build model ~x tau_zip root_x =
  let f =
    model.build_matrices x
    |> Likely.evaluate_pruned tau_zip
    |> Option.get_exn in
  f root_x;;


(* evaluation when the root vector is given as the last entries of the
 * parameter vector *)

(* likelihood of the initial vector added. the root coordinates are now the
 * last two entries of the optimized vector x (ie. x is parameters+root) *)
let eval_build_full build_matrices ranges ~x tau_zip =
  (* convention: build_matrices has to use the lower part of the parameter
   * vector and has to disregard the later part of it *)
  let (a, alpha, ci12) = build_matrices ranges (Gslv.to_array x) in
  let xdim = Vec.dim alpha in
  let rootx =
    Gslv.(subvector ~off:(length x - xdim) ~len:(xdim) x)
    |> U.vec_of_gslvec in
  let lroot = Likely.evaluate_root (a, alpha, ci12) rootx in
  let lcond = ((a, alpha, ci12)
               |> Likely.evaluate_pruned tau_zip
               |> Option.get_exn) @@ rootx in
  lroot +. lcond;;

(* stay with only one header: reading from mma otherwise is messy. *)
let save_to_csv
    ?header ?(to_string=string_of_float)
    filename floats =
  let strings = Array.(map (fun fa -> map to_string fa) floats) in
  let csv = Csv.of_array strings in
  let ch = open_out filename  in
  let csvch = Csv.to_channel ~separator:',' ch in
  let () = match header with
    | None -> ()
    | Some sl -> Csv.output_record csvch sl in
  let () = Csv.output_all csvch csv in
  Csv.close_out csvch

let make_flat_chi2_pars model sl =
  let f (min, parvec) =
    let sc = model.scaled_pars parvec in
    Array.(concat [[|min|]; sc]) in
  Array.map f sl

let make_flat_chi2_pars_corrs model slcorr =
  let f (min, parvec, corrs) =
    let sc = model.scaled_pars parvec in
    Array.(concat [[|min|]; sc; Array.of_list corrs]) in
  Array.map f slcorr


let interpolate frac = Array.map2 (fun v1 v2 -> v1 +. frac *. (v2 -. v1))


(* penalty for wrong ordering on the diagonal. this is a hack to compensate for
 * the fact that the optimizer in gsl can only do rectangular regions. *)

(** enforce increasing diagonal when penalty>1. *)
let a_penalty ~penalty (a, _, _) =
  let f = ref 1. in
  for i = 1 to Mat.dim1 a - 1 do
    if a.{i+1, i+1} <. a.{i, i} then f := penalty; done;
  !f

(* return options: for empty trees we get nothing. *)


(* ATTENTION: ~scp is the scaled, bounded parameter vector,
   ~x is the unbounded one *)

let cell_normalizer zip =
  1. /. float_of_int (Z.length zip)

let eval_build_min_scaled ?(per_cell=false) ?(penalty=2.) model ~scp f_zip =
  let mats = model.build_matrices_scaled scp in
  Likely.evaluate_pruned_rootmin f_zip mats
  >|= ( *.) (a_penalty ~penalty mats)
  >|= ( *.) (if per_cell then cell_normalizer f_zip else 1.)

(* new option: average over the root instead of fitting it *)
let eval_build_avg_scaled ?(per_cell=false) ?(penalty=2.) model ~scp f_zip =
  let mats = model.build_matrices_scaled scp in
  Likely.evaluate_pruned_rootavg f_zip mats
  >|= ( *.) (a_penalty ~penalty mats)
  >|= ( *.) (if per_cell then cell_normalizer f_zip else 1.)

let eval_build_int_scaled ?(per_cell=false) ?(penalty=2.) model ~scp tau_zip =
  let mats = model.build_matrices_scaled scp in
  Likely.evaluate_full tau_zip mats
  >|= ( *.) (a_penalty ~penalty mats)
  >|= ( *.) (if per_cell then cell_normalizer tau_zip else 1.)

let eval_build_full_avg_scaled ?(per_cell=false) ?(penalty=2.) model ~scp tau_zip =
  let mats = model.build_matrices_scaled scp in
  Likely.evaluate_full_avg tau_zip mats
  >|= ( *.) (a_penalty ~penalty mats)
  >|= ( *.) (if per_cell then cell_normalizer tau_zip else 1.)

(* test: evaluate assuming the root cell has root x vector 0,0 *)
let eval_build_0_scaled ?(per_cell=false) ?(penalty=2.) model ~scp f_zip =
  let ((_, alpha, _) as mats) = model.build_matrices_scaled scp in
  let rootx_0 = Vec.(make0 (dim alpha)) in
  Likely.evaluate_pruned f_zip mats
  |> Option.map (fun fn -> fn rootx_0 |> ( *.) (a_penalty ~penalty mats))
  >|= ( *.) (if per_cell then cell_normalizer f_zip else 1.)

(* make the unscaled ones: *)

let _mk_unscaled eval_f =
  fun ?per_cell ?penalty model ~x f_zip ->
    let scp = model.scaled_pars x in
  eval_f ?per_cell ?penalty model ~scp f_zip

let eval_build_min,
    eval_build_avg,
    eval_build_int,
    eval_build_full_avg,
    eval_build_0 =
  _mk_unscaled eval_build_min_scaled,
  _mk_unscaled eval_build_avg_scaled,
  _mk_unscaled eval_build_int_scaled,
  _mk_unscaled eval_build_full_avg_scaled,
  _mk_unscaled eval_build_0_scaled

(*
(* not yet 'modelized' -- integration is not working reliably *)
let eval_build_int_0 build_matrices ranges ~x tau_zip =
  let (a, alpha, ci) = build_matrices ranges (Gslv.to_array x) in
  let xdim = Vec.dim alpha in
  match Likely.evaluate_full_root tau_zip (a, alpha, ci) with
  | None -> raise (Invalid_argument "no mcil possible")
  | Some fn -> fn (Vec.make0 xdim)
   *)

(* experimental data set handling *)


(* map for int tuples to index the experiments *)
module ExpMap = Map.Make (struct
    type t = (int, string) Pair.t
    let compare = Pair.compare Int.compare String.compare
  end)

let exp_merge_both f ma mb =
  let f = fun _ep lr -> match lr with
    | `Left _ | `Right _ -> None
    | `Both (a, b) -> f a b in
  ExpMap.merge_safe ~f ma mb

let exp_filter_map f em =
  ExpMap.mapi f em
  |> ExpMap.filter (fun _ v -> match v with
      | None -> false
      | Some _ -> true)
  |> ExpMap.map Option.get_exn

let list_to_option = function [] -> None | l -> Some l

(* next we need to do a lognormal fit. for that we log the data and then
 * calculate mean and variance *)
let forest_log_mean_vars forest =
  let open Sequence in
  of_list forest
  |> filter_map (Ctree.prune_survivors ?transform:None)
  |> map Z.to_seq
  |> concat
  |> map Z.zval
  |> Stat.log_mean_vars

(* preprocessing pipeline *)
let preprocess_forest ~min_arrival_gen ~max_gen forest =
  let texp = (List.to_seq forest
              |> Sequence.filter_map experiment_T
              |> Sequence.max
              |> Option.get_exn) in
  (* minimum generation numbers that arrive at T *)
  let mingen = forest |> List.map (fun z ->
      cells_alive_at texp z
      |> Sequence.map Z.depth
      |> Sequence.min ?lt:None) in
  (* filter away super low generation arrivals <=3 *)
  let forest' =
    List.map2 begin fun z g ->
      match g with
      | None -> (* non-arriving trees are retained! *)
        Some z
      | Some g when g <= min_arrival_gen ->
        (* arriving low-gen trees indicate slow outliers. cut. *)
        None
      | Some _g -> (* normal arriving trees *)
        Some z end
      forest mingen
    |> List.filter_map (Fun.id) in
  (* filter out overhanging cells over generation gen_lim *)
  let forest'' =
    forest' |> List.filter_map
      (Z.take_while (fun zz -> Z.depth zz <= max_gen)) in
  forest''

let log_transform forest =
  (* log transform the data to be able to fit *)
  Transform.log_transform (forest_log_mean_vars forest)

let log_forest_tau forest =
  List.map (fun z ->
      Z.tree z
      |> Ctree.extract_tau_tree ~transform:(log_transform forest)
      |> Z.zipper)
    forest

let log_forest forest =
  List.filter_map
    (Ctree.prune_survivors ~transform:(log_transform forest))
    forest

let normal_transform ?(margins=`Symmetric) forest =
  (* normal-transform the data to be able to fit
   * without worrying about outliers.
   * margins are for the ECDF: does it start at probability 0, end at 1,
   * or is it symmetrized (default)? *)
  let ar =
      let open Sequence in
      of_list forest
      |> map Z.to_seq
      |> concat
      |> filter_map (fun z -> Ctree.prune_survivors z)
      |> map Z.zval
      |> filter (fun f -> not @@ Float.is_nan f)
      |> to_array in
  Transform.normal_transform ~margins ar

let normal_forest ?margins ?transf_from forest =
  (* if transf_from is given, the transformation is based on it *)
  let transform =
    normal_transform ?margins @@ match transf_from with
    | None -> forest | Some tfrom -> tfrom in
  List.filter_map (Ctree.prune_survivors ~transform) forest

let normal_forest_tau ?margins ?transf_from forest =
  let transform =
    normal_transform ?margins @@ match transf_from with
    | None -> forest | Some tfrom -> tfrom in
  List.map (fun z ->
      z |> Z.tree |> Ctree.extract_tau_tree ~transform |> Z.zipper) forest

(** prune all branches in trees of the forest whose oldest cell has a higher
 * cycle time than the given quantile of the cell cycle time distribution of
 * the whole forest.  *)
let quantile_prune quantile forest =
  let open Sequence in
  let cutoff =
    of_list forest
    |> filter_map (Ctree.prune_survivors ?transform:None)
    |> map Z.to_seq
    |> concat
    |> map Z.zval
    |> filter (fun f -> not @@ Float.is_nan f)
    |> to_array
    |> Stat.quantiles [quantile]
    |> List.hd in
  forest |> List.filter_map (Z.take_while (fun z ->
      match Z.zval z with
      | {tau = Min v; _} | {tau = Val v; _} ->
        (* trick to handle nan case correctly: *)
        not (v >. cutoff)))

(* reusable functions from model_evidence *)

let float_funny flt = match Float.classify flt with
  | FP_normal | FP_zero | FP_subnormal -> false
  | FP_infinite | FP_nan -> true

(** compute the integral of a scaled version of the likeliood. the scaling
 * factor is 10**{exp_offset}. *)
let integrate_evidence' mc_kind n
    ?(exp_offset=0.)
    ?monitor
    ?vegas_state
    pos_log_like
    gaussian_tau_forest
    model =
  let f =
    (* attention: this is a callback passed to the fsl integration routine.
     * exceptions raised in here are questionably propagated. *)
    let pll = pos_log_like model gaussian_tau_forest in     (* allow staging *)
    let counter = ref 0 in
    fun scp ->
      let fscp =
        (* shift in log space to avoid underflow. exp_offset is base 10 *)
        let raw = pll scp in
        let logr = raw +. log 10. *. exp_offset in
        exp logr in
      incr counter;
      let () = match monitor with
        | None -> ()
        | Some mon_fn ->
          let gb = if float_funny fscp then `Bad else `Good in
          (* monitor gets the offset value: *)
          mon_fn !counter gb (Array.copy scp) fscp in
      fscp in
  (* final output is the value including offset - still need to multiply with
   * 1E-exp_offset to get the actual result *)
  let lo, up = model.lo, model.up in
  let open Gsl.Monte in
  match mc_kind, vegas_state with
  | _, None ->
    integrate mc_kind f ~lo ~up n the_rng
  | VEGAS, Some state ->
    integrate_vegas f ~lo ~up n the_rng state
  | (MISER | PLAIN), Some _ ->
    raise (Invalid_argument "Only vegas supports passing a state")


(*exception Offset_wrong of float*)
(** calculate the evidence. as the scaled integral, where the scaling factor is
 * determined by a pre-sampling. the scaling factor is also output. the full
 * unscaled likelihood can be recovered from the output as 10**{-exp_offset} *
 * integral.res, if the float range allows. *)
let evidence'
    ?max_offset
    pos_log_likelihood guess monitor_interval kind n forest model =
  (* init to get an estimate of a good offset. adds 10% to the run time. *)
  let init_offset = guess forest in
  let integral_init = integrate_evidence' kind (n/10) ~exp_offset:init_offset
      pos_log_likelihood forest model in
  let init_res = integral_init.Gsl.Fun.res in
  let exp_offset =
    if float_funny init_res
    then begin
      Format.printf
        "estimated result %f is not finite, using initial guess %f@."
        init_res init_offset;
      init_offset end
    else
      (* with this offset, samples below should be on average close to 1 *)
      let e_o =
        init_offset +.
        Float.round (-. log10 init_res +. log10 (Float.of_int n)) in
      match max_offset with None -> e_o | Some o -> Float.min o e_o in
  (* now the real thing *)
  let samples, values = ref [], ref [] in
  let monitor i gb x fx =
    let is_bad = function `Bad -> true | _ -> false in
    if i mod monitor_interval = 0 || is_bad gb
    then begin
      if is_bad gb then Format.printf "bad sample %f@." fx;
      (samples := x::!samples; values := fx::!values) end in
  let integral = integrate_evidence' kind n
      ~exp_offset ~monitor pos_log_likelihood forest model in
  (integral, exp_offset, !values, !samples)

let args' integrators the_frsts mnms =
  let open List in
  let ex = ExpMap.bindings the_frsts in
  integrators
  |> map (fun i ->
      ex |> map (fun j ->
          mnms |> map (fun k ->
              (i, (j, k)))))
  |> concat |> concat

let integ_string integ =
    (List.assoc ~eq:Pervasives.(=) integ Gsl.Monte.[VEGAS, "vegas"; MISER, "miser"])

(* print to stdout *)
let print_e (exp, description) = Format.sprintf "%d_%s" exp description


(* new addition: write also the shifted results to circumvent floating point
 * reduced range *)

(* write evidences to file *)
let print_header () =
  ["model"; "exp_offset"; "full_result"; "error_estimate"; "samples";
   "scaled_result"; "scaled_error" ]
let print_row = function (mn, (res, offs, vals, _)) ->
  [mn;
   string_of_float offs;
   (* the final, unshifted evidence: *)
   string_of_float (res.Gsl.Fun.res *. 10.**(-.offs));
   string_of_float (res.Gsl.Fun.err *. 10.**(-.offs));
   string_of_int (List.length vals);
   (* the scaled evidence *)
   string_of_float res.Gsl.Fun.res;
   string_of_float res.Gsl.Fun.err;
  ]

(* save the points and values for analysis. *)
(* this now outputs the full non-log, unshifted likelihood *)
let save_points dir dsetname res kind =
  let saver (mn, (_res, offs, vals, sams)) =
    let filename =
      Printf.sprintf
        "%s/samples_values_%s_%s_%s.csv"
        dir dsetname mn kind in
    let vs = List.rev_map2
        (fun vv ss ->
           List.map string_of_float
             ((vv *. 10.**(-.offs)) :: (Array.to_list ss) @ [vv; offs]))
        vals sams in
    Csv.save ~separator:',' filename vs in
  List.iter saver res

(* pair correlations.
 * the root has a nan cycle length; we drop it for the direct ancestors. *)
let thenamedpairs = Z.[
    "ss"     , ss_pairs;
    "md"     , md_pairs';
    "gmgd"   , gmgd_pairs';
    "an"     , an_pairs;
    "cc"     , cc_pairs;
    "ggmggd" , ggmggd_pairs';
    "gagn"   , gagn_pairs;
    "cc1r"   , cc1r_pairs;
    "cc2"    , cc2_pairs; ]

let thepairnames = List.map fst thenamedpairs

(* sanitize: clear nans from float pair sequences. it's better to know where
 * they come from! *)
let sanitize_pair pair z =
  Sequence.filter_map
    (fun (z1, z2) ->
       let v1, v2 = Z.(zval z1, zval z2) in
       if Float.(is_nan v1 || is_nan v2)
       then None
       else Some (v1, v2))
    (pair z)

let model_cors ?(pairnames=thepairnames) model scp =
  let module S = Set.Make(String) in
  if not S.(subset (of_list pairnames) (of_list thepairnames))
  then failwith "unknown relationship encountered";
  let cors = predict_corrs model.build_matrices_scaled scp pairnames in
  Vec.of_list cors

(* a matrix of bootstrap samples of the statistic evaluated for all pair
 * functions. the row index is the sample number. columns are the pair
 * functions. needs a pair function array, which is basically the second column
 * of thenamedpairs. polymorphism trouble -> make that array in the script *)

(* CAUTION: this does not sort out trees that don't contribute to a given
 * pair category! *)
let pair_treewise_bootstrap n stat_fun pairfun_array forest =
  let open Array in
  let seqarrarr =
    forest
    |> List.map (fun z ->
        pairfun_array |> map (fun pf ->
            pf z |> taupairs))
  |> of_list in
  let nforest, nfun = length seqarrarr, length pairfun_array in
  let sample_once () = Sequence.(seqarrarr |> random_array |> take nforest) in
  init n (fun _bs_sample ->
      sample_once ()
      |> Sequence.fold (map2 Sequence.append) (make nfun Sequence.empty)
      |> map stat_fun)
  |> Mat.of_array

(* reject all trees that don't have pairs for all pairfuns! *)
let pair_treewise_bootstrap_filtered n stat_fun pairfun_array forest =
  let open Array in
  let seqarrarr =
    forest
    |> List.map (fun z ->
        map (fun pf -> taupairs (pf z)) pairfun_array)
    |> List.filter (for_all (fun s -> Sequence.length s > 0))
    |> of_list in
  let nforest, nfun = length seqarrarr, length pairfun_array in
  let () = print_int nforest; print_newline () in
  let sample_once () = Sequence.(seqarrarr |> random_array |> take nforest) in
  init n (fun _bs_sample ->
      sample_once ()
      |> Sequence.fold (map2 Sequence.append) (make nfun Sequence.empty)
      |> map stat_fun)
  |> Mat.of_array



(* write evidences to file *)
let print_corr_header () =
  let pn_l = List.map ((^) "l_") thepairnames
  and pn_u = List.map ((^) "u_") thepairnames in
  let corrcols = List.map2 (fun a b -> [a; b]) pn_l pn_u
                 |> List.map2 List.cons thepairnames
                 |> List.concat in
  "model" :: corrcols

let print_corr_row = function (mn, maplulist) ->
  let flatlist =
    maplulist |> List.map (fun (map, (l, u)) -> [map; l; u])
    |> List.concat
    |> List.map string_of_float in
  mn :: flatlist

(* write autocorrelations to file *)
let print_auto_header model =
  let open List in
  let r =
    range model.hdim 1                 (* gotcha: cartesion product reverses *)
    |> map string_of_int in
  let rr =
    cartesian_product [r; r]
    |> map (String.concat "") in
  "tau" :: rr

let print_auto_lines ac_list_list =
  let transpose =
    Mat.(ac_list_list |> of_list |> transpose_copy |> to_list) in
  let strings =
    List.map (fun acl ->
        List.map (fun f -> Format.sprintf "%.3g" f) acl)
      transpose in
  strings

