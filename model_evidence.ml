(*
 * run script to evaluate evidences.
 *
 * this does the following tasks:
 *
 * 1. define a set of models to compare
 * 2. integrate the evidence for them for a given number of datasets
 *
 * raw tree data in csv form have to be present in subdirectories of
 * $HOME/data_trees, see below.
 * output is written to $HOME/tmp/evidence_results/
 *
 * multinormal integration parameters and parameters for the evidence
 * calculation are set in the source code, below.
 *
 * the script can be called as "model_evidence.native batch" in batch mode,
 * which is meant for running it on a cluster.
 *)

open! Containers

open Option.Infix
open Lacaml.D

open Util
open Tree
open Hmtree
open Import_csv
open Test_lib
open Models

module Seq = Sequence
module IM = Integrate_mvnd

let () = Printexc.record_backtrace true;;

(* parallelization. we have no parallel start pars now. *)
let parmap_processes =
  try Sys.getenv "PARMAP_PROCESSES" |> int_of_string
  with Not_found -> 1
let _ = assert (1 <= parmap_processes && 26 >= parmap_processes)

let _ = Format.printf "using %d cores \n" parmap_processes

let () = cores := parmap_processes

let rng_seed = 12n

let batch_mode =
  try String.equal Sys.argv.(1) "batch"
  with Invalid_argument _ -> false
let () = Printf.printf "batch mode: %b\n" batch_mode

(* increase integration sample points -- sometimes needed *)
(* 15_000 is enough for < 0.01% cases where the limit is insufficient *)
let () = Integrate_mvnd.(settings.max_pts_per_dim <- 15_000)

let () = Gsl.Rng.set the_rng (rng_seed)

(* use pre-cut trees or full overlapping trees *)
let (what_trees : [`Precut | `Precut_overhang | `Uncut]) =
  `Precut
  (*`Precut_overhang*)
  (*`Uncut*)

let results_dir =
  if batch_mode then "./"
  else let suffix = match what_trees with
      | `Precut -> "/"
      | `Precut_overhang -> "_precut_overhang/"
      | `Uncut -> "_uncut/" in
    Unix.getenv "HOME" ^ "/tmp/evidence_results" ^ suffix

let () = U.provide_dir results_dir

(* do the integration if possible? *)
type mode = [`Avg | `Min | `Int_avg]
let (eval_mode : mode) =
  `Int_avg
  (*`Avg*)

type treeconf = { dir         : string;
                  tree_suffix : string;
                  integrate   : mode }

(* now decide what to do actually *)
let tconf = match what_trees with
  | `Precut ->
    { dir="precut/"     ; tree_suffix="_Trees_globalgen.csv" ; integrate=`Avg }
  | `Precut_overhang ->
    { dir="precut/"     ; tree_suffix="_Trees_globalgen.csv" ; integrate=eval_mode }
  | `Uncut ->
    { dir="uncut/"      ; tree_suffix="_Trees.csv"           ; integrate=eval_mode }

(* need precut in any case, for determining an unbiased normal transformation *)
let precut_conf = { dir="precut/";
                    tree_suffix="_Trees_globalgen.csv";
                    integrate=`Avg }

(* real or reduced real or simulated or pooled or transf-pooled trees? *)
let (data : [ `Pooled_then_transformed
            | `Transformed_then_pooled
            | `Real
            | `Real_half
            | `Sim ]) =
  `Real
  (*`Pooled_then_transformed*)
  (*`Transformed_then_pooled*)
  (*`ESC*)

let (keys_to_use :
       [ `All
       | `Select of ExpMap.key list]) =
  let elist s = List.map (fun n -> (n, s)) in
  (*`All*)
  (*`Select (elist "on" [58; 60; 62])*)
  (*`Select (elist "off" [58; 60])*)
  (*`Select (elist "rap" [69; 82])*)
  `Select (elist "esc" [1; 2; 3])
  (*`Select [58, "on"]*)
  (*`Select [62, "on"]*)
  (*`Select [62, "on"; 69, "rap"; 82, "rap"]*)
  (*`Select [1, "esc"]*)
  (*`Select [3, "esc"; 58, "on"]*)
  (*`Select [1, "esc"; 2, "esc"; 3, "esc"]*)
  (*`Select [82, "rap"]*)
  (*`Select [586062, "on"]*)
  (*`Select [5862, "on"]*)
  (*`Select [5860, "off"]*)
  (*`Select [6982, "rap"]*)


(* MC integrators *)
let mc_integrators = Gsl.Monte.[
    MISER;
    (*VEGAS*)
  ]

(* outliers *)
let (outlier_treatment :
       [ `Normal | `Quantile of float | `Size of int | `Log_only ]) =
  `Normal
(*`Size 9*)
(*`Quantile 0.98*)

(* likelihood level where to cut off samples taken for the prediction bands *)
let top_level =
 (* for a 1D standard normal: integral from -1.96 to 1.96 gives 0.95. at these
  * points the likelihood is down to 0.146 times the peak value *)
  0.146

(* corresponding to the level, we do the confidence level for the bootstrap *)
let corr_bs_confidence = 0.95
let corr_bs_n = truncate
    (*1e1*)
    1e4

(* evidence integration parameters *)
let monitor_interval = 10
let n_points = truncate
    (*1e3*)
    (*1e4*)
    (*1e5*)
    1e6

(* how many generations to calculate the autocorrelation function *)
let autocorr_gens = 10

(* make some models *)

let adiag_r  = (0., 0.99)               (* < 1 to prevent numerical problems *)
let across_r = (-5., 5.)
let alpha_r  = (0.01, 2.)             (* alpha=0 would open symmetries again *)
let neg_alpha_r = (-. snd alpha_r, -. fst alpha_r)
let chi_r    = (-0.99, 0.99)       (* |.| < 1. is necessary for pos def.ness *)
let gamma_r  = (0., 1.)


(* + *)
let m_1_indep = model_1_indep
    {a11=adiag_r;
     alpha1=alpha_r}

(* + *)
let m_1 = model_1
    {a11=adiag_r;
     alpha1=alpha_r;
     gamma_=gamma_r}

(* + *)
let m_2_uni = model_2_uni
    {a11=adiag_r; a22=adiag_r; a12=across_r;
     alpha1=alpha_r; alpha2=alpha_r;
     gamma_=gamma_r}

(* + *)
let m_2_uni_indep = model_2_uni_indep
    {a11=adiag_r; a22=adiag_r; a12=across_r;
     alpha1=alpha_r; alpha2=alpha_r}

(* + *)
let m_2_indep_indep = model_2_indep_indep
    {a11=adiag_r; a22=adiag_r;
     alpha1=alpha_r; alpha2=alpha_r}

(* + *)
let m_2_indep = model_2_indep
    {a11=adiag_r; a22=adiag_r;
     alpha1=alpha_r; alpha2=alpha_r;
     gamma_=gamma_r}

(* + *)
(* this should be equivalent to m_2_uni *)
let m_2_indep_crossnoise = model_2_indep_crossnoise
    {a11=adiag_r; a22=adiag_r;
     alpha1=alpha_r; alpha2=alpha_r;
     chi_=chi_r; gamma_=gamma_r}

(* + *)
(* this should be equivalent to m_2_uni_indep *)
let m_2_indep_crossnoise_indep = model_2_indep_crossnoise_indep
    {a11=adiag_r; a22=adiag_r;
     alpha1=alpha_r; alpha2=alpha_r;
     chi_=chi_r}

(* + *)
let m_2_uni_same = model_2_uni_same
    {a11=adiag_r; a12=across_r;
     alpha1=alpha_r; alpha2=alpha_r;
     gamma_=gamma_r}

(* + *)
let m_2_uni_only = model_2_uni_only
    {a12=across_r;
     alpha1=alpha_r; alpha2=alpha_r;
     gamma_=gamma_r}

(* + *)
let m_2_uni_a1 = model_2_uni_a1
    {a11=adiag_r; a22=adiag_r; a12=across_r;
     alpha1=alpha_r;
     gamma_=gamma_r}

(* + *)
let m_2_bi = model_2_bi
    {a11=adiag_r; a22=adiag_r; a12=across_r; a21=across_r;
     alpha1=alpha_r; alpha2=alpha_r;
     gamma_=gamma_r}

(* + *)
let m_2_uni_neg_a2 = model_2_uni
    {a11=adiag_r; a22=adiag_r; a12=across_r;
     alpha1=alpha_r; alpha2=neg_alpha_r;
     gamma_=gamma_r}

(* + *)
let m_2_indep_neg_a2 = model_2_indep
    {a11=adiag_r; a22=adiag_r;
     alpha1=alpha_r; alpha2=neg_alpha_r;
     gamma_=gamma_r}

(* + *)
let m_2_indep_posneg = model_2_indep
    {a11=adiag_r; a22=adiag_r;
     alpha1=alpha_r; alpha2=(-. snd alpha_r, snd alpha_r);
     gamma_=gamma_r}

(* + *)
let m_2_uni_a22 = model_2_uni_a22
    {a22=adiag_r; a12=across_r;
     alpha1=alpha_r; alpha2=alpha_r;
     gamma_=gamma_r}

(* synthetic data *)

(* as a point of comparison for Exp58on *)
module Ap_ref : PARS = struct
  let hdim = 2
  let alpha = Vec.of_array [|1.; 1.06|]
  let meanScale =
    let _A = Mat.of_array [|[|0.74 ; -0.5|];
                            [|0.0  ; 0.52|]|] in
    U.block_matrix [|[|_A|];
                     [|_A|]|]
  let cdGm =
    let (sigma1, sigma2) = (1., 1.) in
    let gamma = 0.3 in
    let dia = [sigma1; sigma2] |> Vec.of_list |> Mat.of_diag in
    let odia = [gamma *. sigma1; gamma *. sigma2] |> Vec.of_list |> Mat.of_diag
    in
    U.block_matrix [|[|dia  ; odia|];
                     [|odia ; dia|]|]
  let noise = 0.
end

module Smp_ref = Sample (Ap_ref);;

let pforest =
  let initial_x () = Vec.of_list [0.;0.] in
  let make_pt ~root i =
    let rt = root in
    Z.zipper (T.create rt)
    |> Smp_ref.sample_joint i
    |> Z.map (fun v -> Val (dot Ap_ref.alpha v)) in
  let forest_size = 15 in
  let tree_height = 6 in
  List.init forest_size (fun _i ->
      let root = initial_x () in
      make_pt ~root tree_height)

let sim_forests = ExpMap.(empty |> add (58, "on_sim") pforest)


(* data import. new: import _everything_ here! cut down later. *)

let dir, precut_dir =
  Pair.map_same (fun conf ->
      Sys.getenv "HOME"
      ^ (if not batch_mode
         then  "/tmp/data_trees/"
         else "/data/")
      ^ conf.dir)
    (tconf, precut_conf)

let read_forests dir suffix =
  let f name =
    try exp_forest (dir ^ name ^ suffix)
    with Sys_error _ -> [] in
  let en i = List.assoc ~eq:Int.equal i Import_csv.esc_names in
  let add k f m = match f with
    | [] -> m
    | l -> ExpMap.add k l m in
  let fs = ExpMap.empty
           |> add (62, "on") (f "Exp62on")
           |> add (58, "on") (f "Exp58on")
           |> add (60, "on") (f "Exp60on")
           |> add (58, "off") (f "Exp58off")
           |> add (60, "off") (f "Exp60off")
           |> add (69, "rap") (f "Exp69rap20nM")
           |> add (82, "rap") (f "Exp82rap40nM")
           |> add (1, "esc") (f ("Exp" ^ en 1 ^ "esc"))
           |> add (2, "esc") (f ("Exp" ^ en 2 ^ "esc"))
           |> add (3, "esc") (f ("Exp" ^ en 3 ^ "esc"))
  in                                         (* CUT ONLY IF NOT INTEGRATING! *)
  match tconf.integrate with
  | `Int_avg  -> fs
  | _ -> ExpMap.map (List.filter_map Ctree.prune_survivors_cells) fs


let all_forests = read_forests dir tconf.tree_suffix
let ab = ExpMap.keys all_forests |> Seq.to_list

(* for normal transform in the full-tree case *)
let precut_forests = read_forests precut_dir precut_conf.tree_suffix
let pb = ExpMap.keys precut_forests |> Seq.to_list

(* data pooling. this has to be done here on the level of forests, not in the
 * input csv files -- those rely on global labels which would get mixed up. *)
let pool_data fs =
  let open ExpMap in
  let get k = find k fs in
  empty
  |> add (5862, "on")
    (get (58, "on") @ get (62, "on"))
  |> add (586062, "on")
    (get (58, "on") @ get (60, "on") @ get (62, "on"))
  |> add (5860, "off")
    (get (58, "off") @ get (60, "off"))
  |> add (6982, "rap")
    (get (69, "rap") @ get (82, "rap"))

(* these are pooled before any subtraction of the mean etc *)
let pooled_raw_forests = pool_data all_forests

let pooled_precut_raw_forests = pool_data precut_forests

(* IMPORTANT: reject outliers *)
let preprocess ~precut_forests m =
  let open ExpMap in
  let transformed = begin match outlier_treatment with
    | `Quantile p ->         (* typical value 0.98. this can change results! *)
      m |> map (quantile_prune p) |> map log_forest_tau
    | `Size l ->                                         (* typical value: 9 *)
      m |> map (List.filter (fun z -> Z.length z > l)) |> map log_forest_tau
    | `Normal ->
      m |> mapi (fun k v ->
          (* _always_ use precut for normal transformation. fail if not. *)
          normal_forest_tau ~transf_from:(find k precut_forests) v)
    | `Log_only ->             (* no rejection. not recommended for evidence *)
      m |> map log_forest_tau
  end in
  transformed |> map (List.filter (fun z -> Z.length z > 1))

let transformed_forests = preprocess ~precut_forests all_forests

let pooled_then_transformed_forests = preprocess
    ~precut_forests:pooled_precut_raw_forests pooled_raw_forests

(*let transformed_esc_forests = preprocess esc_data*)

(* alternative pooling: after transformation. i.e. mean etc are subtracted
 * _before_ pooling. this is not the same as data pooling. possibly a bit
 * fishy. *)
let transformed_then_pooled_forests = pool_data transformed_forests

(* downsampling for testing *)
let transformed_forests_smaller =
  let av_length = ExpMap.map List.length transformed_forests
                  |> ExpMap.values
                  |> Seq.map float_of_int
                  |> Seq.to_array
                  |> Stat.mean
                  |> int_of_float in
  transformed_forests
  |> ExpMap.map (List.take (av_length / 2))

(* choose the forests to use *)
let the_forests' = match data with
  | `Pooled_then_transformed -> pooled_then_transformed_forests
  | `Transformed_then_pooled -> transformed_then_pooled_forests
  | `Real                    -> transformed_forests
  | `Real_half               -> transformed_forests_smaller
  | `Sim                     -> sim_forests
  (*| `ESC                     -> transformed_esc_forests*)

(* now finally, cut down to the data we actually want to process *)
let the_forests = let open ExpMap in
  match keys_to_use with
  | `All -> the_forests'
  | `Select l ->
    List.fold_left (fun em k -> add k (find k the_forests') em) empty l


(* make the experimental pair correlations and bootstrap confidence intervals.
 * take stuff from bs_like_data: we don't need a covariance matrix. *)

let theqtile =
  let q = (1. -. corr_bs_confidence) /. 2. in
  [q; 1. -. q]


(* CHANGE: now give the full-data value as center, not the median of the BS. *)
let corbslist frst cor n =
  thenamedpairs
  |> List.map begin fun (_name, pair) ->
    let pure_pair z = Seq.map (Pair.map_same Z.zval) (pair z) in
    let bs =
      forest_bootstrap_pure n cor pure_pair frst
      (* this should not contain nans anymore! *)
      (*|> Seq.of_array*)
      (*|> Seq.filter (fun f -> not (Float.is_nan f))*)
      (*|> Seq.to_array *)
    in
    let value = frst_full_pure cor pure_pair frst in
    let midlohi = match Stat.quantiles theqtile bs with
      | [l; h] -> value, (l, h)
      | _ -> raise (Invalid_argument "need 2 quantiles exactly.") in
    midlohi end

let corfulllist frst cor =
  thenamedpairs
  |> List.map (fun (_name, pair) ->
      let pure_pair z = Seq.map (Pair.map_same Z.zval) (pair z) in
      (*let pure_pair = purify_pair pair in*)
      frst_full_pure cor pure_pair frst)

let to_float_forest taufrst =
  List.filter_map (Ctree.prune_tau_survivors) taufrst

let rank_cor_val_bounds =
  the_forests |> ExpMap.map (fun frst ->
      corbslist (to_float_forest frst) rcor corr_bs_n)

let rank_cor_once =
  the_forests |> ExpMap.map (fun frst ->
      corfulllist (to_float_forest frst) rcor_n)

(* now, evaluate likelihoods *)
(* ranges and true_ranges are set to the same *)

(* choose evaluation *)

let per_cell = false              (* this would not give a correct evidence! *)
let penalty =
  let adiag_order_penalty_factor = 2. in
  adiag_order_penalty_factor

(* now, tautype forests *)
let evaluate model tau_zip scp =
  try match tconf.integrate with
    | `Min ->
      Ctree.prune_tau_survivors tau_zip
      >>= (fun f_zip ->                       (* could be none if tree empty *)
          eval_build_min_scaled ~penalty ~per_cell model ~scp f_zip)
    | `Avg ->
      Ctree.prune_tau_survivors tau_zip
      >>= (fun f_zip ->
          eval_build_avg_scaled ~penalty ~per_cell model ~scp f_zip)
    | `Int_avg ->
      eval_build_full_avg_scaled ~penalty ~per_cell model ~scp tau_zip
  with
  | U.Diverging_recurrence ->
    Some infinity (* accomodate bi models where A eigenvalues are complicated *)
  | IM.Integration_error `Not_enough_points -> Some 1.

let pos_log_like' eval model frst scp =
  frst
  |> List.map (fun zp -> eval model zp scp)
  (* raise on empty trees which give None. within a callback, this just leads
     to a NaN result! *)
  (* debug !! *)
  |> List.filter_map (fun x -> x)
  (*|> List.map opt_get*)
  (* sum of pos natural log likelihoods *)
  |> List.fold_left (-.) 0.

let pos_log_like model frst scp =
  pos_log_like' evaluate model frst scp

let integrate_evidence mc_kind n
    ?(exp_offset=0.)
    ?monitor
    gaussian_tau_forest
    model =
  integrate_evidence' mc_kind n
    ~exp_offset
    ?monitor
    pos_log_like
    gaussian_tau_forest
    model

let integrate_ev_state vegas_state n
    ?(exp_offset=0.)
    ?monitor
    gaussian_tau_forest
    model =
  integrate_evidence' Gsl.Monte.VEGAS n
    ~exp_offset
    ~vegas_state
    ?monitor
    pos_log_like
    gaussian_tau_forest
    model


(* this is nice but we do not need the subdivisions of the parameter space and
 * nonuniform sampling, if we
 * 1. take the MAP from the list of sampled values
 * 2. only use bounds from the set of all samples that are more than x% as
 * likely as the MAP. no way to know what percentage of the posterior mass that
 * is! *)

let guess_offset forest =
  (* trial/error. some of these give NaN especially in vegas *)
  (*3. *. float_of_int (List.length forest)*)
  (* new finer try with actual tree size. it turns out the vegas method is very
     picky. *)
  let totl = List.fold_left (fun ac z -> ac + Z.length z) 0 forest in
  0.5 *. float_of_int totl

(* test! with simulated data... *)
(*let tf = ExpMap.find (58, "on_sim") sim_forests*)
let tf = ExpMap.find (62, "on") the_forests'
let samples = ref []
let values = ref []
let vegas_state ofile dim =
  let open Gsl.Monte in
  let vs = make_vegas_state dim in
  init_vegas vs;
  {(get_vegas_params vs) with
   mode=STRATIFIED;
   verbose=2;
   ostream=Some ofile;}
  |> set_vegas_params vs; vs
let mcres oc =
  let model = m_1 in
  let dim = model.pdim in
  let f () =
    Format.printf "old guess offset %f@." (2. *. float_of_int (List.length tf));
    let exp_offset = guess_offset tf in
    Format.printf "new guess offset %f@." exp_offset;
    let state = vegas_state oc dim in
    let integral = integrate_ev_state state 1_000
        ~exp_offset
        ~monitor:(fun _i _gb x fx ->
            samples := x::!samples; values := fx::!values)
        tf m_1 in
    (integral, exp_offset, state) in
  U.time f ()
let mcr = IO.with_out (results_dir ^ "vegas_info") mcres

let () =
  let i, o, _ = mcr in
  Format.printf "test done; result %f; offset %f@." i.Gsl.Fun.res o


(* ok, now make this for all models *)

(* choice:
   m1_indep: naive.
   m1: minimal extension 1 par
   m2_indep_indep: simplest 2 var. _same would be 1var
   m2_uni: old preferred. 6par
   m2_uni_indepr simplified without sister. 5par
   m2_uni_same: simplified but not trivial. 5par
   m2_indep_crossnoise: should be equiv to m2_uni. 6par
*)

let mnames = [
  m_1                  , "m_1";

  m_2_uni_only         , "m_2_uni_only";
  m_1_indep            , "m_1_indep";
  m_2_indep_indep      , "m_2_indep_indep";
  m_2_uni              , "m_2_uni";
  m_2_uni_indep        , "m_2_uni_indep";
  m_2_uni_same         , "m_2_uni_same";
  m_2_indep_crossnoise , "m_2_indep_crossnoise";
  m_2_uni_a1           , "m_2_uni_a1";
  m_2_bi               , "m_2_bi";
  m_2_indep            , "m_2_indep";
  m_2_uni_neg_a2       , "m_2_uni_neg";
  m_2_indep_neg_a2     , "m_2_indep_neg";
  m_2_indep_posneg     , "m_2_indep_posneg";
  m_2_uni_a22          , "m_2_uni_a22";

]

let models = fst (List.split mnames)
let model_of_name mn =
  List.map Pair.swap mnames
  |> List.assoc ~eq:String.equal mn

let evidence kind n forest model =
  evidence' ~max_offset:600.
    pos_log_like guess_offset monitor_interval kind n forest model





(* a more parallel mapping strategy. unfortunately, it's a hassle to flatten
 * and rebuild... *)

let ncores = !cores

let args = args' mc_integrators the_forests mnames

let map_combinations f =
  if ncores > 1 then
    let l = Parmap.L args in
    Parmap.parmap ~ncores f l
  else
    List.map f args

let print_int_stats () =
  let tot, sm0, not_enough =
    IM.(float_of_int integration_stats.evaluations,
        100. *. float_of_int integration_stats.smaller0,
        100. *. float_of_int integration_stats.not_enough_points) in
  match tconf.integrate with
  | `Int_avg -> Format.printf
                  "total integral evaluations: %0.0f\n\
                   smaller than 0 in %f%%\n\
                   not enough points in %f%%@."
                  tot (sm0/.tot) (not_enough/.tot)
  | _ -> Format.printf
           "no integration requested@."


let calc_integral (integ, ((key, forest), (m, m_name))) =
  Format.printf "expt %d%s, model %s, method %s@."
    (fst key) (snd key) m_name (integ_string integ);
  IM.(
    integration_stats.evaluations <- 0;
    integration_stats.smaller0 <- 0;
    integration_stats.not_enough_points <- 0);
  let n = Gsl.Monte.(match integ with
      | MISER | PLAIN -> n_points
      | VEGAS -> n_points/4) in
  let res = U.time (evidence integ n forest) m in
  print_int_stats ();
  (integ, (key, (m_name, res)))

let ev_res =
  let open List in
  let res = map_combinations calc_integral in
  (* lexical sort. uses polymorphic compare and equality below. ugly but it
   * works at least. *)
  let res = sort Pervasives.compare res in
  (* rebuilding the hierarchical order *)
  let group l =
    let f = function
      | [] -> raise (Invalid_argument "empty sublist")
      | (a, b) :: tl -> (a, b :: map snd tl) in
    l |> group_succ ~eq:(fun (a, _) (a', _) -> Pervasives.(=) a a') |> map f in
  let res_nested =
    res |> group |> map (fun (integ, l) -> (integ, l |> group)) in
  let res_expmap = map (fun (i, ass) -> (i, ExpMap.of_list ass)) res_nested in
  res_expmap

(** get the MAP likelihood and value, and the list of samples whose likelihood
 * is at most down by a factor of "level" *)
let top_v_samples level (_res, _offset, vals, samples) =
  let vs = List.combine vals samples
         |> List.sort Pervasives.compare
         |> List.rev in
  let map_v, map_s = List.hd vs in
  let ok (v, _) = v >. map_v *. level in
  map_v, map_s, List.(take_while ok vs |> map snd)

(* ..factor out the mapping down to the level of results *)
let map_ev f ev_struct =
 ev_struct |> List.map (fun (kind, exp_results) ->
      (kind,
       exp_results |> ExpMap.map (fun models_res ->
                 models_res |> List.map f )))

(* best samples *)
let ev_res_top = map_ev
    (fun (mn, x) ->
       mn, top_v_samples top_level x)
    ev_res

(* returns: model name, (best fit pred, list of acceptable preds) *)
let ev_res_pair_predictions =
  ev_res_top
  |> map_ev (fun (mn, (_v, scp, top_scp)) ->
      let model = model_of_name mn in
      let cors = model_cors model in
      mn, (cors scp, List.map cors top_scp))

(* take the max and min in each correlation over the set of acceptable preds *)
let ev_res_pair_limits =
  ev_res_pair_predictions
  |> map_ev
    begin fun (mn, (cor, cors)) ->
      let d = List.length thepairnames in
      let min' =
        let min_start = Vec.make d infinity in
        List.fold_left (fun v v' -> Vec.min2 v v') min_start cors
      and max' =
        let max_start = Vec.make d (-.infinity) in
        List.fold_left (fun v v' -> Vec.max2 v v') max_start cors in
      Vec.(concat [min'; max']) |> Vec.iter (fun v ->
          if Float.is_nan v then raise (Invalid_argument "found nan"));
      let map' = Vec.to_list cor in
      let lu_list = List.combine (Vec.to_list min') (Vec.to_list max') in
      mn, List.combine map' lu_list end


let _ =
  let print_ev res =
    res |> List.iter (fun (mn, (res, offs, vals, _)) ->
        Format.printf "%s. res: %.15g; offset %f; samples %d\n"
          mn res.Gsl.Fun.res offs (List.length vals)) in
  List.iter
    (fun (integ, em) ->
       Format.printf "\n%s@." (integ_string integ);
       ExpMap.iter (fun k v ->
           Format.printf "%s\n" (print_e k);
           print_ev v) em)
    ev_res

let () =
  let save r f = Csv.save
      ~separator:',' f ((print_header ()) :: List.map print_row r) in
  ev_res |> List.iter begin fun (integ, em) ->
    em |> ExpMap.iter begin fun k v ->
      save v (Format.sprintf "%s/model_evidences_%s_%s.csv"
                results_dir (print_e k) (integ_string integ))
    end end

let () =
  ev_res |> List.iter begin fun (integ, em) ->
    em |> ExpMap.iter begin fun k v ->
      save_points results_dir (print_e k) v (integ_string integ)
    end end

let () =
  let print_exp_corr k =
    let expq = ExpMap.find k rank_cor_val_bounds in
    print_corr_row ("data", expq) in
  let save k r filename = Csv.save
      ~separator:',' filename
      (print_corr_header () ::
       print_exp_corr k ::
       List.map print_corr_row r) in
  ev_res_pair_limits |> List.iter begin fun (integ, em) ->
    em |> ExpMap.iter begin fun k v ->
      save k v (Format.sprintf "%s/model_corr_preds_%s_%s.csv"
                  results_dir (print_e k) (integ_string integ))
    end end


(* we need the autocorrelation functions *)

(* matrices of the best-fits *)
let ev_res_bestfit_mats = map_ev
  (fun (mn, (_v, scp, _scp_list)) ->
    mn, (model_of_name mn).build_matrices_scaled scp)
  ev_res_top

(* write out also the best fit parameter values! *)

let () =
  let fitdir = "best_fits/" in
  let () = U.provide_dir (results_dir ^ fitdir) in
  let save filename (a, alpha, ci12) =
    let write oc =
      Format.(Lacaml.Io.(
          fprintf (of_chan oc)
            "A:\n%a\nalpha\n%a\nC12\n%a\n@."
            pp_fmat a
            pp_fvec alpha
            pp_fmat (U.ge_inv ci12))) in
    IO.with_out filename write in
  ev_res_bestfit_mats |> List.iter begin fun (integ, expmap) ->
    expmap |> ExpMap.iter begin fun k mn_mat_list ->
      mn_mat_list |> List.iter begin fun (mn, mats) ->
        let fn = (Format.sprintf "%s%smats_%s_%s_%s.txt"
                    results_dir fitdir (print_e k) mn (integ_string integ)) in
        save fn mats
      end end end

let path_up i =
  let rec p l = function
    | 0 -> l
    | j -> p (`Up :: l) (j - 1) in
  p [] i

let make_upcorrs (a, alpha, ci12) =
  (* correlations to the ancestors *)
  let d = Vec.dim alpha in
  let paths = List.init autocorr_gens path_up in
  let components = List.range 1 d in
  let uv i = let v = Vec.make0 d in v.{i} <- 1.; v in
  let ac12 =
    List.map (fun i ->
        List.map (fun j ->
            let v1, v2 = uv i, uv j in
            List.map (predicted_xcorr a ci12 v1 v2) paths)
          components)
      components
    |> List.concat in
  let act = List.map (predicted_corr (a, alpha, ci12)) paths in
  act :: ac12

let ev_res_bestfit_autocorrs = map_ev (fun (mn, mats) ->
    mn, make_upcorrs mats)
    ev_res_bestfit_mats

let () =
  let save filename (mn, ac_lines) = Csv.save
      ~separator:',' filename
      ((print_auto_header (model_of_name mn)) :: print_auto_lines ac_lines) in
  ev_res_bestfit_autocorrs
  |> List.iter begin fun (integ, expmap) ->
    expmap |> ExpMap.iter begin fun k mn_acll_list ->
      mn_acll_list |> List.iter begin fun (mn, ac_ll) ->
        save
          (Format.sprintf "%s/autocorrelations_%s_%s_%s.csv"
             results_dir (print_e k) mn (integ_string integ))
          (mn, ac_ll) end end end

