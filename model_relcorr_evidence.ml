(* in this script, we repeat the evidence calculations in model_evidence but
 * this time with an approximate likelihood function that is taken from the
 * correlation coefficients bewteen related cells only.
 *
 * to this end we estimate variances and covariance of all related correlation
 * coeffients from a bootstrap, directly from the data. for this we need a bit
 * of extra machinery. *)


open Containers
open Lacaml.D
open Util
open Tree
open Hmtree
open Import_csv
open Test_lib
open Models

module Seq = Sequence

let batch_mode =
  try String.equal Sys.argv.(1) "batch"
  with Invalid_argument _ -> false
let () = Printf.printf "batch mode: %b\n" batch_mode

(* parallelization. we have no parallel start pars now. *)
let parmap_processes =
  try Sys.getenv "PARMAP_PROCESSES" |> int_of_string
  with Not_found -> 1
let _ = assert (1 <= parmap_processes && 24 >= parmap_processes)

let _ = Format.printf "using %d cores \n" parmap_processes

let () = cores := parmap_processes
let rng_seed =
  Nativeint.of_string Sys.(argv.(Array.length argv - 1))
  |> Option.get_or ~default:1234n

(* convenience for parallel mapping *)
let map_list f args =
  let ncores = !cores in
  if ncores > 1 then
    Parmap.(parmap ~ncores f (L args))
  else
    List.map f args

let () = Gsl.Rng.set the_rng (rng_seed)
let results_dir =
  if batch_mode then "./"
  else "relcor_evidence_results/"
let () =
  let open Unix in
  try
    mkdir results_dir 0o755;
    Printf.printf "made directory\n\n"
  with (Unix_error (EEXIST, _, _)) -> ()


(* simulation parameters *)

let correlation_bootstrap_n = int_of_float
    5e2
    (*1e2*)

let n_points = int_of_float 1e5
let monitor_interval = 10

(* forest_bootstrap from test_lib does it in the wrong order to be able to get
 * covariances. do that again. *)

(* make some models. these are records now! *)

(* ranges *)
let adiag_r  = (0. , 0.99)                     (* prevent numerical problems *)
let across_r = (-5., 5.)
let alpha_r  = (0.01, 2.)                  (* alpha=0 opens symmetries again *)
let neg_alpha_r = (-. snd alpha_r, -. fst alpha_r)
let chi_r    = (-0.99, 0.99)        (* |.| < 1. is necessary for pos defness *)
let gamma_r  = (0. , 1.)


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
    |> Z.map (fun v -> dot Ap_ref.alpha v) in
  let forest_size = 15 in
  let tree_height = 6 in
  List.init forest_size (fun _i ->
      let root = initial_x () in
      make_pt ~root tree_height)

let pforest_rec = List.map (Z.map (fun tau ->
    {tau=Val tau; tbirth=nan; id=0})) pforest


(* data import *)

let dir = if not batch_mode then
    Sys.getenv "HOME" ^ "/tmp/tree_data/precut/"
  else Sys.getenv "HOME" ^ "/data/precut/"

let exp_forests = let open ExpMap in
  let f name = exp_forest (dir ^ name ^ "_Trees_globalgen.csv") in
  let _escf esc_number =
    let name =(List.assoc ~eq:(=) esc_number Import_csv.esc_names) in
    exp_forest (dir ^ "Exp" ^ name ^ "esc_Trees_globalgen.csv") in
  empty
  |> add (62, "on") (f "Exp62on")
  |> add (58, "on") (f "Exp58on")
  (*|> add (58, "on_sim") pforest_rec*)
  |> add (82, "rap") (f "Exp82rap40nM")
  |> add (69, "rap") (f "Exp69rap20nM")
  |> add (60, "off") (f "Exp60off")
  |> add (58, "off") (f "Exp58off")
  (*|> add (1, "esc") (escf 1)*)

  (*|> add (60, "on") (f "Exp60on")*)
  (*|> add (60, "off") (f "Exp60off")*)
  (*|> add (82, "rap") (f "Exp82rap40nM")*)


(* we need also float forests to comply with the interface of model_evidence *)

let the_forests =
  ExpMap.map (List.filter_map Ctree.prune_survivors) exp_forests

(* now make the bs samples *)

(* the full array of relation functions *)
let pairfun_array =
  let thepairs = List.map snd thenamedpairs in
  Array.of_list thepairs

let bs_samples n corf =
  exp_forests |> ExpMap.to_list
  |> map_list (fun (k, frst) ->
      (k, pair_treewise_bootstrap n corf pairfun_array frst))
  |> ExpMap.of_list

(* take only those trees into the bootstrap where each pairfun produces
 * non-emtpy pairs *)
let bs_samples_filtered n pairfuna corf =
  exp_forests |> ExpMap.to_list
  |> map_list (fun (k, frst) ->
      (k, pair_treewise_bootstrap_filtered n corf pairfuna frst))
  |> ExpMap.of_list

(* the full array of relation functions *)
let pairfunsmdc =
  let get a = List.assoc ~eq:Equal.string a thenamedpairs in
  [| get "md"; get "cc"; get "cc2"|]

(* now do some work!
 * compute the rank correlations *)

(* original *)
(*let thesamples = U.time (bs_samples correlation_bootstrap_n) rcor*)

(* this requires each tree to have c2 cousins... *)
let thesamples = U.time (bs_samples_filtered correlation_bootstrap_n pairfunsmdc) rcor

(* this doesn't. still the same number of samples! *)
let make_samples () =
  let s1 = bs_samples_filtered
      correlation_bootstrap_n [|pairfunsmdc.(0); pairfunsmdc.(1)|] rcor in
  let s2 = bs_samples_filtered
      correlation_bootstrap_n [|pairfunsmdc.(0); pairfunsmdc.(2)|] rcor in
  s1, s2
let samplesc1c2 = U.time make_samples ()

let themeans = ExpMap.map Stat.mean_samples thesamples

let thecovs = ExpMap.map Stat.var_covar_samples thesamples

let () = Seq.iter (fun k ->
    Format.printf "key %d%s\n" (fst k) (snd k);
    Format.printf "mean %a\n" pp_vec (ExpMap.find k themeans);
    Format.printf "cov %a\n\n@." pp_mat (ExpMap.find k thecovs);)
    (ExpMap.keys themeans)

(* we also want the ratios of c1, c2 to md,
 * for the prediction verification plot. *)

let mdratio_samples =
  let md mat = Mat.col mat 1 in
  let c1r mat = Vec.div (Mat.col mat 2) (md mat) in
  let c2r mat = Vec.div (Mat.col mat 3) (md mat) in
  let c12r mat = Mat.of_col_vecs_list
      [c1r mat;
       c2r mat] in
  ExpMap.map c12r thesamples

let mdratio_samples =
  let div mat = Vec.div (Mat.col mat 2) (Mat.col mat 1) in
  let s1, s2 = samplesc1c2 in
  let cr1 = ExpMap.map div s1 in
  let cr2 = ExpMap.map div s2 in
  let f key = function
    | `Both (l, r) -> Some (Mat.of_col_vecs_list [l; r])
    | _ -> None in
  ExpMap.merge_safe ~f cr1 cr2

let mdratio_means = ExpMap.map Stat.mean_samples mdratio_samples

let mdratio_covs = ExpMap.map Stat.var_covar_samples mdratio_samples

let mdratio_sdevs =
  ExpMap.map (fun mat -> Mat.copy_diag mat |> Vec.sqrt) mdratio_covs

let mdratio_qts =
  let qsamples mat =
    List.init (Mat.dim2 mat) (fun i ->
        Stat.quantiles
          [0.5; 0.16; 0.84]
          (Mat.col mat (succ i) |> Vec.to_array))
    |> List.concat
    |> Vec.of_list in
  ExpMap.map qsamples mdratio_samples


let () = Seq.iter (fun ((ex, kind) as k) ->
    Format.printf "key %d%s@." ex kind;
    Format.printf "mean ratio %a@." pp_vec (ExpMap.find k mdratio_means);
    Format.printf "m.r. sdev %a@.@." pp_vec (ExpMap.find k mdratio_sdevs);
    Format.printf
      "ratio median lower upper %a@,@." pp_vec (ExpMap.find k mdratio_qts);)
    (ExpMap.keys mdratio_means)


(* export the mdratio data *)
let mdrow key =
  ExpMap.find key mdratio_qts
  |> Vec.to_list
  |> List.map (fun f -> Format.sprintf "%.3f" f)
  |> List.cons (string_of_int (fst key) ^ snd key)

let _ =
  let mdcsv = List.map mdrow (ExpMap.keys mdratio_means |> Seq.to_list)
              |> List.cons ["expt"; "c1/md"; "c1/md lo"; "c1/md hi";
                            "c2/md"; "c2/md lo"; "c2/md hi"]
  in
  Csv.save (results_dir ^ "cousin12_md_rank_corr_ratios_devs.csv") mdcsv


let mdsample_csv () =
  let keys, mats = List.split (ExpMap.bindings mdratio_samples) in
  let bigmat = Mat.of_col_vecs_list
      List.(map (Mat.to_col_vecs_list) mats |> concat) in
  Mat.to_list bigmat |> List.map (fun col ->
      col |> List.map (fun f ->
          Format.sprintf "%.3f" f))
  |> List.cons (List.product (fun (ex, cond) r ->
      r) keys ["c1r"; "c2r"])
  |> List.cons (List.product (fun (ex, cond) r ->
      string_of_int ex ^ cond) keys ["c1r"; "c2r"])

let _ =
  let mdsample = mdsample_csv () in
  Csv.save (results_dir ^ "c1r_c2r_samples.csv") mdsample

(* ok -- scipy confirms *)

(* to reuse the functions from model_evidence, we need evaluation functions
 * that take forests, not experiment keys. to do this we make a hash table with
 * the _stripped_ experimental forests as keys, for memoizing the results above
 * *)

(* allow for a subset of correlations to fit and calculate evidences *)

let frst_means, frst_covs =
  let frst_hash em =
    let t = Hashtbl.create 8 in
    ExpMap.(iter (fun k v -> Hashtbl.add t (find k the_forests) v) em);
    t in
  (frst_hash themeans, frst_hash thecovs)

let frst_mean = Hashtbl.find frst_means
let frst_cov = Hashtbl.find frst_covs


(* next thing: implement a simple Gaussian likelihood with this covariance and
 * mean. it should take a forest as input, calculate all the
 * relative-correlations and plug that in.
 * missing part: for given model parameters, build
 * the mean of the relative correlations (we have this) *)

(*let reduced_pairnames = ["ss"; "md"; "gmgd"; "an"]*)

let reduced_pairnames = ["ss"; "md"; "gmgd"; "an"; "cc"]

let eval_model pairnames mean covar model =
  let idx name =
    1 + (List.find_idx (String.equal name) thepairnames |> opt_get |> fst) in
  let mean = pairnames
             |> List.map (fun name -> mean.{idx name}) |> Vec.of_list in
  let ci =
    let cv = pairnames |> List.map (fun name ->
          (pairnames |> List.map (fun name' ->
               covar.{idx name, idx name'})))
      |> Mat.of_list in
    U.ge_inv cv in
  fun scp ->
    try
      let cors = model_cors ~pairnames model scp in
      let logz = Likely.log_Z ci in
      let dev = Vec.sub cors mean in
      -. logz +. 0.5 *. U.vmv dev ci dev
    with
    (* accomodate bi models where A eigenvalues are complicated and we cannot
     * exclude diverging eigenvalues easily *)
      U.Diverging_recurrence -> infinity

let eval_exp pairnames key model =
  let get m = ExpMap.find key m in
  let (mean, covar) = (get themeans, get thecovs) in
  eval_model pairnames mean covar model


let pos_log_like ?(pairnames=reduced_pairnames) model frst =
  let mean, covar = frst_mean frst, frst_cov frst in
  let ev =  eval_model pairnames mean covar model in
  fun scp ->
    let res = -. ev scp in
    if Float.is_nan res then raise (Invalid_argument "nan!");
    res


(* test! with simulated data... *)
let integrate_evidence ?(pairnames=reduced_pairnames) mc_kind n
    ?(exp_offset=0.)
    ?monitor
    t_forest
    model =
  integrate_evidence' mc_kind n
    ~exp_offset
    ?monitor
    (pos_log_like ~pairnames)
    t_forest
    model

let tf = ExpMap.find (58, "on_sim") the_forests
let samples = ref []
let values = ref []

let mcres =
  let f () =
    let exp_offset = 2. in
    let integral = integrate_evidence
        (*~pairnames:["ss"]*)
        Gsl.Monte.MISER 500
        ~exp_offset
        ~monitor:(fun _i `Good x fx ->  (* TODO check if `Good is correct *)
            samples := x::!samples; values := fx::!values)
        tf m_1 in
    (integral, exp_offset) in
  U.time f ()

(* now for real *)
let mnames = [
  m_1_indep            , "m_1_indep";
  m_1                  , "m_1";
  m_2_indep_indep      , "m_2_indep_indep";
  m_2_uni              , "m_2_uni";
  m_2_uni_indep        , "m_2_uni_indep";
  m_2_uni_same         , "m_2_uni_same";
  m_2_uni_only         , "m_2_uni_only";
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

(* the integration routine from test_like *)
let evidence =
  evidence'
    ~max_offset:100.
    (*(pos_log_like ~pairnames:["ss"])*)
    (pos_log_like ~pairnames:reduced_pairnames)
    (*(pos_log_like ~pairnames:thepairnames)*)
    (fun _frst -> 0.)
    monitor_interval

let args = args' Gsl.Monte.[VEGAS] the_forests mnames

let integ_string = integ_string

(* copied literally *)

let calc_integral (integ, ((key, forest), (m, mn))) =
  Format.printf "expt %d%s, model %s, method %s@."
    (fst key) (snd key) mn
    (integ_string integ);
  let n = Gsl.Monte.(match integ with
      | MISER | PLAIN -> n_points
      | VEGAS -> n_points/4) in
  let res = U.time (evidence integ n forest) m in
  integ, (key, (mn, res))

let ev_res =
  let open List in
  let res = map_list calc_integral args in
  (* lexicographic sort. this works! *)
  let res = sort Pervasives.compare res in
  (* rebuilding the hierarchical order *)
  let group l =
    let f = function
      | [] -> raise (Invalid_argument "empty sublist")
      | (a, b) :: tl -> (a, b :: map snd tl) in
    l |> group_succ ~eq:Equal.(map fst poly) |> map f in
  let res_nested =
    res |> group |> map (fun (integ, l) -> (integ, l |> group)) in
  let res_expmap =
    res_nested |> map (fun (i, ass) -> (i, ExpMap.of_list ass)) in
  res_expmap

let _ =
  let print_ev res =
    res |> List.iter (fun (mn, (res, offs, vals, _)) ->
        Format.printf "%s. res: %5.5g; offset %f; samples %d\n"
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
      save_points
        results_dir (print_e k) v (integ_string integ)
    end end

