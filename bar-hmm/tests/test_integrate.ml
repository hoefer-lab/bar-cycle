(* test the integration over overhanging cells (i.e. cells that have not
 * finished the cycle at the end of the experiment *)

open Containers;;
open Lacaml.D;;
open Util;;
open Tree;;
open Hmtree;;
open Test_lib;;
open Models;;

(* max_overhead 1_000_000 prevents compaction except requested explicitly. *)
(*Gc.(set { (get ()) with*)
          (*verbose = 0x010; max_overhead = 1_000_000 })*)

(* parallelization *)
let parmap_processes =
  try Sys.getenv "PARMAP_PROCESSES" |> int_of_string
  with Not_found -> 1
let _ = assert (1 <= parmap_processes && 8 >= parmap_processes)
let () = Test_lib.parallel_startpars :=
    if parmap_processes > 1 then Use_Parmap else No
let () = Test_lib.cores := parmap_processes

(* how many starting values for the minimizations;
 * parallelization runs over this. *)
let startsamples = 20

let rng_seed =
    Nativeint.of_string Sys.(argv.(Array.length argv - 1))
  |> Option.get_or ~default:1234n
let () = Gsl.Rng.set the_rng (rng_seed)

let batch_mode =
  try String.equal Sys.argv.(1) "batch"
  with Invalid_argument _ -> false

let () = Printf.printf "batch mode: %b\n" batch_mode

let results_dir =
  if batch_mode then "./"
  else "some_meaningful_run_name/"
let () = U.provide_dir results_dir

(* switch to the smaller parameter set, hopefully without redundancy this time,
 * namely a11, a22, a12, α1, α2, β. σ1,σ2 are set to 1. *)

(* take parameters found for the 58pos5 myc on experiment *)
module Ap' = struct
  let hdim = 2
  let noise = 0.
  let alpha = Vec.of_array [|1.;1.06|]
  let meanScale =
    let _A = Mat.of_array [|[|0.52;-0.5|];
                            [|0.0; 0.74|]|] in
    U.block_matrix [|[|_A|]; [|_A|]|]
  let cdGm =
    let (sigma1, sigma2) = (1., 1.) in
    let gamma = 0.3 in
    let dia = [sigma1; sigma2] |> Vec.of_list |> Mat.of_diag in
    let odia = [gamma *. sigma1; gamma *. sigma2] |> Vec.of_list |> Mat.of_diag
    in
    U.block_matrix [|[|dia; odia|]; [|odia; dia|]|]
end

module Smpl' = Sample (Ap');;
module Obs' = Observe (Ap');;

let zvi () = Smpl'.initial_zp ()

(* crude conversion to cell rec *)
let make_tau' vec =
  {tau=Val (dot Ap'.alpha vec); tbirth=nan; id=0}

(* these contain the true parameter ranges *)
let true_ranges = {a11=(-1.,1.);
                   a22=(-1.,1.);
                   a12=(-3.,3.);
                   alpha1=(0.,10.);
                   alpha2=(0.,10.);
                   gamma_=(0.,1.)}


(* old *)
(*{adiag=(-1.,1.); across=(-3.,3.);*)
(*a2=(-0.01,10.);*)
(*sigma=(1e-5,3.); gamma=(0.,1.)}*)

(* do the standard fitting procedure *)
let ranges = true_ranges
(* new parameter building rules *)
let m = model_2_uni true_ranges
let eval = eval_build m
let par_length = 6

let it_pars =
  let step_init_size = 0.5 in
  let maxit, abs_tol, succ_steps = 20, 0.001, 10 in
  {step_init_size; maxit; abs_tol; succ_steps};;

let true_parvec = m.make_parvec Ap'.(meanScale, alpha, U.ge_inv cdGm)

let rebuilt_true_matrices = m.build_matrices true_parvec
let rebuilt_true_parvec = m.scaled_pars true_parvec

let forest_size = 35
let tree_height = 6

(* introduce a mean time to get a cutoff.
   this could be made more realistic by actually exponentiating the time --
   using a log normal distribution effectively *)
let m_tau = 2.
let transform = {Transform.to_time=((+.) m_tau); from_time=((-.) m_tau)}

(* regulate pruning *)
let max_T =  m_tau *. float_of_int tree_height
             *. 1.1

(* forest construction *)
let make_pt ~root i =
  root
  |> Smpl'.sample_tree i
  |> Ctree.build_cell_tree ~transform Obs'.obs
  |> Option.get_exn


let pforest' = List.init forest_size (fun _ ->
    let root = Smpl'.Gauss.draw_initial () in
    (*Vec.fill root 0.;  (* a test using the conditioned eval functions. *)*)
    make_pt ~root tree_height)

(* impose a cutoff time to generate overlapping trees. *)
let pforest_overlapping =
  pforest'
  |> List.map (Ctree.record_cell_tree ~max_T)

(* remove again the mean *)
let pforest_obs =
  pforest_overlapping
  |> List.map (Ctree.extract_tau_tree ~transform)
  |> List.map Z.zipper

(* cut back to a float zipper with values of the internal trees only.
   roughly half the size! *)
let pforest =
  pforest_overlapping
  |> List.filter_map (fun tr ->
      Ctree.prune_survivors ~transform (Z.zipper tr))

let test_sizes =
  List.(map Z.length pforest_obs, map Z.length pforest)

(* now, evaluate likelihoods *)

let per_cell = false              (* this would not give a correct evidence! *)
let penalty =
  let adiag_order_penalty_factor = 2. in
  adiag_order_penalty_factor

(* _dmy argument is there to keep signature the same *)

(* root x optimized; pruned tree *)

(*let eval_rootmin model ~x f_zip _dmy = eval_build_min model ~x f_zip*)
let eval_rootavg model ~x f_zip _dmy = eval_build_avg ~penalty model ~x f_zip

let eval_root0 model ~x f_zip _dmy = eval_build_0 ~penalty model ~x f_zip

(* root x integrated; overhanging integrated, full tree *)

(* integration settings *)

let () = Integrate_mvnd.(
    settings.max_pts_per_dim <- 5000;
    settings.rel_error       <- 1e-2;
    settings.abs_error       <- 1e-6;
    settings.neg_log_clip <- -.1e2)

let eval_int_counter = ref 0
let eval_int model ~x tau_zip _dmy =
  incr eval_int_counter;
  eval_build_int ~penalty model ~x tau_zip

let eval_int_avg model ~x tau_zip _dmy =
  incr eval_int_counter;
  eval_build_full_avg ~penalty model ~x tau_zip

(* done: the int evaluation with averaged root cell produces the same as the
 * non-int evaluation in the case of non truncated forests. *)


(*let eval_int_true ~x tau_zip _dmy =*)
  (*eval_build_int build_matrices_2 true_ranges ~x tau_zip*)


let _ =
  let xtest = Gsl.Vector.create ~init:0. par_length in
  U.time
    (fun () ->
       for _i = 1 to 50 do
         ignore @@ eval_rootavg m ~x:xtest (List.hd pforest) ()
       done) ()

(* with -O3 this is about a factor of 8 slower, in native code. the whole
 * program spends >90% of its time in the mvnd integration routine. *)
let last = !U.last_timed
let _ =
  let xtest = Gsl.Vector.create ~init:0. par_length in
  U.time
    (fun () ->
       for _i = 1 to 10 do
         ignore @@ eval_int_avg m ~x:xtest (List.hd pforest_obs) ()
       done) ()

(* roughly 6X slowdown. not too bad actually.... *)
let _ = 1./. last *. !U.last_timed

let _ = test_sizes;;

let () = parallelize_forest := No

let minimizer =
  U.time
    (find_min_forest
       {it_pars with maxit=1340} (eval_rootavg m) pforest ())
    (Gsl.Vector.create ~init:0. par_length)
  |> get_anyway

let minimizer_int =
  eval_int_counter := 0;
  U.time
    (find_min_forest {it_pars with maxit=1340} (eval_int_avg m) pforest_obs ())
    (*(Gsl.Vector.create ~init:0. par_length)*)
    (Gsl.Vector.copy true_parvec)     (* try to start close to true to see.. *)
  |> get_anyway

(* this still seems not to work at all. something may be wrong with the mvnd
 * wrapper or with the (compilation of the) fortran code itself *)

let fmin, res_par = report_res m minimizer

let fmin_int, res_par_int = report_res m minimizer_int

(* interesting test.
 * the non-integrating branch of cifull already gives
 * something else as best-fit?!?!?!?!? this should not be. what does it
 * actually do? *)


(* test if integration vs truncated evaluation give correlated likelihood
 * over parameters *)

let par_corrs n =
  let rand _ =
    let open Vec in
    let lo, up = (of_array m.lo, of_array m.up) in
    sub up lo
    |> mul (random ~from:0. ~range:1. m.pdim)
    |> add lo
    |> to_array
    (*|> (fun a -> Format.printf "%a@." (Array.pp Float.pp) a; a)*)
    |> m.unconstrained_pars in
  let parsets = List.init n rand in
  let res eval forest = List.map (fun x ->
      try forest_eval (eval m) ~x forest ()
      with Integrate_mvnd.Integration_error _ -> Float.nan)
      parsets in
  res eval_rootavg pforest, res eval_int_avg pforest_obs

let (lla, lli) =
  let la, li = par_corrs 62 in
  Csv.save "llcorrs.csv"
    List.([map string_of_float la; map string_of_float li]);
  (la, li)

(* test if integration vs truncated evaluation give correlated likelihood
 * values over forests *)
let forest_corrs size n =
  let frsts _ =
    let forest_overlapping =
      List.init size (fun _ ->
          make_pt ~root:(Smpl'.Gauss.draw_initial ()) tree_height)
      |> List.map (Ctree.record_cell_tree ~max_T:(max_T/.1.5)) in
    let forest_obs =
      forest_overlapping
      |> List.map (Ctree.extract_tau_tree ~transform)
      |> List.map Z.zipper in
    let forest =
      forest_overlapping
      |> List.filter_map (fun tr ->
          Ctree.prune_survivors ~transform (Z.zipper tr)) in
    forest, forest_obs in
  let res eval1 eval2 forest1 forest2 =
    (*let parvec = Gslv.create ~init:0. m.pdim in*)
    let parvec = Vec.random m.pdim |> Integrate_mvnd.v2c in
    try (forest_eval (eval1 m) ~x:parvec forest1 (),
         forest_eval (eval2 m) ~x:parvec forest2 ())
    with
    | Integrate_mvnd.Integration_error _
    | Invalid_argument _ -> Float.(nan, nan) in
  List.init n (fun _i ->
      let forest, forest_obs = frsts _i in
      (*res eval_rootavg eval_int_avg forest forest_obs)*)
      res eval_rootavg eval_int forest forest_obs)
  |> List.split

let (fla, fli) =
  let la, li = forest_corrs 35 63 in
  Csv.save "flcorrs.csv"
    List.([map string_of_float la; map string_of_float li]);
  (la, li)


(* THIS SEEMS TO WORK NOW! WE CAN DO FULL-TREE BIAS CORRECTIONS! YESS! *)
(* tested for int_avg and int versions... *)

let () = parallelize_forest := No
let _ = !eval_int_counter

(*let () =*)
  (*let true_chi2 =*)
    (*let x = Gsl_vector.of_array true_parvec in*)
    (*forest_eval ~x (eval_rootmin true_ranges) pforest () in*)
  (*Format.printf "true chi2: %.5f@." true_chi2;*)
  (*Format.printf "minimizing chi2: %.5f@." fmin*)

let () =
  let true_chi2_int =
    let x = true_parvec in
    forest_eval ~x (eval_int m) pforest_obs () in
  Format.printf "true chi2: %.5f@." true_chi2_int;
  Format.printf "minimizing chi2: %.5f@." fmin_int

let resf_avg =
  let make () =
    let min_fun = forest_eval (eval_rootavg m) pforest () in
    results
      ~keep_stopped:true
      {it_pars with maxit=560}
      (startpars ~samples:startsamples par_length)
      min_fun in
  U.time make ()
(*1926s with these parameters and with parallel *)

let _ = Format.printf "results best chi squares: \n %a"
    (List.pp Float.pp)
    List.(sort Float.compare (map fst (fst resf_avg)))

let resf =
  eval_int_counter := 0;
  let make () =
    (*let min_fun = forest_eval (eval_int m) pforest_obs () in*)
    let min_fun = forest_eval (eval_int_avg m) pforest_obs () in
    results
      ~keep_stopped:true
      {it_pars with maxit=30}
      (startpars ~samples:startsamples par_length)
      min_fun in
  (* 200s with 10 and 30 *)
  (* 11_000s with 100 and 300 *)
  U.time make ()

let _ = Format.printf "results best chi squares: \n %a"
    (List.pp Float.pp)
    List.(rev @@ sort Float.compare (map fst (fst resf)))

let sl = make_sl resf

(*let true_chi2_forest_roots = forest_eval (eval_rootmin true_ranges) pforest ()*)
    (*~x:(Gsl_vector.of_array true_parvec)*)

let true_chi2_forest_int = forest_eval
    (eval_int m) pforest_obs () ~x:true_parvec

(*let _ =*)
  (*let fmin, res_par = sl.(0) in*)
  (*print_res ~builder:build_matrices ranges fmin res_par;*)
  (*print_res ~builder:build_matrices true_ranges true_chi2_forest_int true_parvec*)

let () =
  let flat = make_flat_chi2_pars m sl in
  let true_string =
    let ss = Array.map string_of_float rebuilt_true_parvec in
    let sss = String.concat "," (Array.to_list ss) in
    ("true parameters: " ^ sss) in
  save_to_csv
    ~header:[Printf.sprintf "time needed: %.2f" !U.last_timed;
             Printf.sprintf
               "true parameters chi2: %.2f" true_chi2_forest_int;
             true_string]
    (results_dir ^ "integrated_forest_best_fits.csv") flat

