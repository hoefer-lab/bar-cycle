(* simulation of trees according to the the growth-progression model. this is a
 * reimplementation of the equivalent code in R, mainly for speed and
 * parallelizability on a cluster. *)


open! Containers

open Lacaml.D

open Util
open Tree
open Hmtree
open Import_csv

module Tl = Test_lib

(* global settings *)
let () = Printexc.record_backtrace true

(* seeding of the rng *)
(*let () = Format.print_string @@ Sys.getenv "PBS_JOBID"*)
let () =
  Sys.getenv_opt "PBS_JOBID"
  |> Option.flat_map (fun s -> String.split_on_char '.' s |> List.head_opt)
  |> Option.flat_map Nativeint.of_string
  |> Option.iter (fun n ->
      Format.printf "Seeding rng to %nd@." n; Gsl.Rng.set the_rng n)

module Global = struct
  let n_bootstrap = 1000
end

type parray =
  [ `One of float array                                  (* one parameter set *)
  | `Arr of (* array of parameter sets,
              array first running par, array second running par *)
      (float array array array * float array * float array)
  | `Flat_list of float array list ]

module Stats = struct
  let get_tau = (fun (_, _, (tau:float)) -> tau)
  let pairnames = Tl.thepairnames

  (* we can keep the roots since these now have a meaningful time. *)
  let thenamedpairs = Z.[
      "ss"     , ss_pairs;
      "md"     , md_pairs;
      "gmgd"   , gmgd_pairs;
      "an"     , an_pairs;
      "cc"     , cc_pairs;
      "ggmggd" , ggmggd_pairs;
      "gagn"   , gagn_pairs;
      "cc1r"   , cc1r_pairs;
      "cc2"    , cc2_pairs; ]

  let make_all_corrs trees =
    let module S = Sequence in
    let pairs pairfun =
      trees
      |> S.of_list
      |> S.map Z.zipper
      |> S.map pairfun
      |> S.concat                        (* concatenating all relative pairs *)
      |> S.map (Pair.map_same (fun z ->
          get_tau (Z.zval z))) in
    thenamedpairs
    |> List.split |> snd
    |> List.map pairs
    |> List.map Tl.rcor


  let growth_fraction mean_progression_time trees =
    let open Sequence in
    let s =
      trees |> of_list |> map (T.to_seq) |> concat
      (* test if tau was longer than progression *)
      |> map (fun (x, _, tau) ->
          let tauprog = exp x *. mean_progression_time in
          if tau >. tauprog then 1 else 0) in
    float_of_int (fold (+) 0 s) /. float_of_int (length s)
  let mean_tau trees =
    let open Sequence in
    let s = trees |> of_list |> map (T.to_seq) |> concat in
    (fold (fun acc (_, _, tau) -> acc +. tau) 0. s)
    /. float_of_int (length s)
  let quartiles_tau trees =
    let open Sequence in
    trees |> of_list |> map (T.to_seq) |> concat
    |> map get_tau
    |> to_array
    |> Stat.quantiles [0.25; 0.5; 0.75]
end

let read_argv argv =
  try
    let a =  Array.(sub argv 1 (length argv - 1)) in
    match Array.length a with
    | 1 ->                            (* input file of flat parameter combos *)
      let sl = Csv.load ~separator:' ' a.(0) in
      `Flat_list List.(map (fun l -> map float_of_string l |> Array.of_list) sl)
    | 10 ->
      `One (Array.map float_of_string a)
    | 14 ->                   (* unfold into an array of full parameter sets *)
      let a = Array.map float_of_string a in
      let aprefix = Array.sub a 0 8 in
      let murange = U.arange a.(8) a.(9) ~step:a.(10)
      and sigrange = U.arange a.(11) a.(12) ~step:a.(13) in
      let amusig =
        Array.map (fun mu ->
            Array.map (fun sigma ->
                Array.concat [aprefix; [|mu; sigma|]])
              sigrange)
          murange in
      `Arr (amusig, murange, sigrange)
    | _  ->
      failwith "unsupported number of arguments"
  with Failure _ ->
    Format.print_string (
      "\n Usage:\n simulate_growth.{byte,native} gen_init gen_tree n_trees" ^
      " threshold_noise growth_rate saturation" ^
      " rho_md rho_ss mean_log_prog_time std_log_prog_time\n\n" ^
      " Alternatively, both last two arguments can be" ^
      " of the form 'start max step\n" ^
      " Finally, a file can be given as the only argument, which contains\n" ^
      " a single parameter combination on each line\n\n" ^
      " For exponential growth set saturation=nan.\n" ^
      " For logistic growth the parameter the law is" ^
      " s' = growth_rate/saturation s (saturation - s)\n");
    exit 1


(*let a = read_argv Sys.argv*)

module type ARGS = sig
  val a : float array
end

module Cl (Args:ARGS) = struct
  open Args
  let gen_init           = int_of_float a.(0)
  let gen_tree           = int_of_float a.(1)
  let n_trees            = int_of_float a.(2)
  let threshold_noise    = a.(3)
  let growth_rate        = a.(4)
  let saturation         = a.(5)
  let md_correlation     = a.(6)
  let ss_correlation     = a.(7)
  let mean_log_prog_time = a.(8)
  let prog_std           = a.(9)
  let growth_law = if Float.is_nan saturation
    then `Exp else `Logistic
  let mean_prog_time = exp mean_log_prog_time
end


module type PG_MODEL = sig
  (* threshold is at 1 by definition. vary the rate only *)
  val mean_prog_time : float
  val threshold_noise : float
  val det_end_size : float -> float -> float
  val rand_growth_time : float -> Gsl.Rng.t -> float
end


module Model (Args:ARGS) : PG_MODEL = struct
  module Cl = Cl (Args)
  let make = function
    | `Exp ->
      (module (struct
        include Cl
        let det_end_size start_size tau =
          start_size *. exp (growth_rate *. tau)
        let rand_growth_time start_size rng =
          let end_size = 1. +. U.centered_gauss rng threshold_noise in
          log (end_size /. start_size) /. growth_rate
      end) : PG_MODEL)
    | `Logistic ->
      (module (struct
        include Cl
        let det_end_size start_size tau =
          let e = exp (growth_rate *. tau) in
          (saturation *. start_size *. e) /.
          (saturation +. start_size *. (e -. 1.))
        let rand_growth_time start_size rng =
          let end_size = 1. +. U.centered_gauss rng threshold_noise in
          let l = log (end_size *. (start_size -. saturation) /.
                       (start_size *. (end_size -. saturation))) in
          l /. growth_rate
      end): PG_MODEL)
  include (val (make Cl.growth_law))
end

module PPars (Args:ARGS) : PARS = struct
  module Cl = Cl (Args) open Cl
  let hdim = 1
  let alpha = Vec.of_list [1.]
  let noise = 0.
  let c =
    let a = md_correlation in
    let s = prog_std in
    s ** 2. *. (1. -. a ** 2.)
  let gamma =
    let a = md_correlation in
    (ss_correlation -. a ** 2.) /. (1. -. a ** 2.)
  let cdGm =
    Mat.of_list [[c ; gamma *. c]; [gamma *. c; c]]
  let meanScale = Mat.of_list [[md_correlation];
                               [md_correlation]]
end


module Sample_Resource' (Args:ARGS) = struct
  module PPars = PPars (Args)
  module PGPars = Model (Args)
  module Cl = Cl (Args)
  (* implementation using fold_filter_tree. do the full tree first only with
     progression, then fold to add growth. *)
  module Prog = Sample (PPars)
  let sample_tree (branching:[`One | `Two]) n =
    let dim = Mat.dim2 PPars.meanScale in
    assert (Mat.(dim1 PPars.meanScale / dim) = 2);
    fun (x_root, s_end_root, tau_root) ->
      (* first construct the pure progression tree *)
      let progression_tree =
        let v = Vec.make 1 x_root in
        match branching with
        | `One -> Prog.sample_tree_to ~child_num:1 n v
        | `Two -> Prog.sample_tree n v in
      let f v_child (_, s, _) =
        let x_child = v_child.{1} in
        let tp = exp x_child *. PGPars.mean_prog_time in
        let s_start = s /. 2. in
        let tg = PGPars.rand_growth_time s_start the_rng in
        let tau = Float.max tg tp in
        if Float.is_nan tg then failwith
            "growth time is nan";
        if Float.is_nan tp then failwith
            (Printf.sprintf "progression time is nan, x_child = %.3f" x_child);
        let s_end = PGPars.det_end_size s_start tau in
        Some (x_child, s_end, tau) in
      (* to get exacly the same semantics as the other sampler we need to start
       * one level lower. a bit ugly -- could revert for production. *)
      (*T.fold_filter_tree ~f (x, s_end, tau) progression_tree*)
      let tchildren = match progression_tree with `Node (_, cl) ->
        cl
        |> List.map (T.fold_filter_tree ~f (x_root, s_end_root, tau_root))
        |> List.map Option.get_exn in
      (`Node ((x_root, s_end_root, tau_root), tchildren): (float * float * float) T.t)
  let init_sample_one () =
    let tree_root =
      sample_tree `One Cl.gen_init (0., 1., Cl.mean_prog_time)
      |> T.to_seq |> Seq.rev |> Seq.head_exn in
    sample_tree `Two Cl.gen_tree tree_root
end


let collect_results parameters =
  let res_table = Hashtbl.create 18 in
  let res_one ba_dim ba_pos a =
    let module S = Sample_Resource' (struct let a = a end) in
    let trees = List.init S.Cl.n_trees (fun _ -> S.init_sample_one ()) in
    let store1 k v =
      let ga = match Hashtbl.get res_table k with
        | None ->
          let ga' = Bigarray.(Genarray.create Float64 C_layout ba_dim) in
          Hashtbl.add res_table k ga'; ga'
        | Some a ->
          a in
      Bigarray.Genarray.set ga ba_pos v in
    store1
      "growth_fraction" (Stats.growth_fraction S.PGPars.mean_prog_time trees);
    store1
      "mean" (Stats.mean_tau trees);
    List.combine ["q1"; "median"; "q3"] (Stats.quartiles_tau trees)
    |> List.iter (Pair.merge store1);
    List.combine Tl.thepairnames (Stats.make_all_corrs trees)
    |> List.iter (Pair.merge store1); in
  let () = match parameters with
    | `One a ->
      res_one [|1|] [|0|] a
    | `Arr (aa, amu, asig) ->
      let nmu = Array.length amu in
      let nsig = Array.length asig in
      for i = 0 to nmu - 1 do for j = 0 to nsig - 1 do
          res_one [|nmu; nsig|] [|i; j|] aa.(i).(j);
        done done in
  res_table


(* collect parameters for storage purposes. *)
let collect_params parameters =
  let open Bigarray in
  let pars = Hashtbl.create 10 in
  let add_int key i =
    let g = (Genarray.create Int16_unsigned C_layout [|1|]) in
    Genarray.set g [|0|] i;
    Hashtbl.add pars key (`Int g) in
  let add_float key f =
    let g = (Genarray.create Float64 C_layout [|1|]) in
    Genarray.set g [|0|] f;
    Hashtbl.add pars key (`Float g) in
  let add_farray key fa =
    let g = (Genarray.create Float64 C_layout [|Array.length fa|]) in
    Array.iteri (fun i f -> Genarray.set g [|i|] f) fa;
    Hashtbl.add pars key (`Float g) in
  let a = match parameters with
    | `One a -> a
    | `Arr (aa, _, _) -> aa.(0).(0) in
  let module Cl =
    Cl (struct let a = a end) in
  add_int "init_generations" Cl.gen_init;
  add_int "tree_generations" Cl.gen_tree;
  add_int "n_trees" Cl.n_trees;
  add_float "growth_rate" Cl.growth_rate;
  add_float "saturation" Cl.saturation;
  add_float "threshold_noise" Cl.threshold_noise;
  add_float "rho_md" Cl.md_correlation;
  add_float "rho_ss" Cl.ss_correlation;
  let _ = match parameters with
    | `One _ ->
      add_float "mean_log_prog_time" Cl.mean_log_prog_time;
      add_float "std_log_prog_time" Cl.prog_std;
    | `Arr (_, amu, asig) ->
      add_farray "mean_log_prog_time" amu;
      add_farray "std_log_prog_time" asig in
  pars


module Export = struct
  let add_ids t =
    let f revpos (x, size, tau) =
      let id = String.concat ""
          (* increase so that we have 1 and 2 and id's do not get mangled *)
          (List.rev_map (fun i -> string_of_int (i + 1)) revpos) in
      (id, x, size, tau) in
    T.map_rev_pos f [] t
  let pair_data t =
    let f m cl acc =
      let c1 = List.nth_opt cl 0 in
      let c2 = List.nth_opt cl 1 in
      (m, c1, c2) :: acc in
    T.fold_mother_children ~f [] t
  let pair_record ((idm, _, _, taum), c1, c2) =
    let fp f = Format.sprintf "%.2f" f in
    let strings = function
      | None -> ("nan", "nan", "nan", "nan")
      | Some (i, x, s, t) ->
        (i, fp x, fp s, fp t) in
    let id1, _x1, _s1, tau1 = strings c1 in
    let id2, _x2, _s2, tau2 = strings c2 in
    [tau1; tau2; fp taum; id1; id2; idm]
  let trees_csv trees =
    let module S = Sequence in
    trees
    |> List.map add_ids
    |> List.map pair_data
    |> S.of_list
    |> S.map (List.map pair_record)
    |> S.mapi (fun i l -> List.map (List.cons (string_of_int i)) l)
    |> S.intersperse [[]]                    (* empty lines between families *)
    |> S.to_list
    |> List.concat
  let write_npz filename pars results =
    let module Nz = Npy.Npz in
    let outfile = Nz.open_out filename in
    Hashtbl.iter (Nz.write outfile) results;
    Hashtbl.iter (fun k -> function
        | `Int gi -> Nz.write outfile k gi
        | `Float gf -> Nz.write outfile k gf)
      pars;
    Nz.close_out outfile
  let external_convert_to_mat filename =
    let matfilename = List.(
        String.split ~by:"." filename
        |> rev |> tl |> cons "mat" |> rev
        |> String.concat ".") in
    let python_command =
      {|python -c "import numpy; import scipy.io; scipy.io.savemat('|} ^
      matfilename ^ {|', numpy.load('|} ^ filename ^ {|'))"|} in
    let status = Sys.command python_command in
    if not (status = 0)
    then failwith "external conversion to MAT format failed"
  let par_order = [
    "init_generations" ;
    "tree_generations" ;
    "n_trees" ;
    "growth_rate" ;
    "saturation" ;
    "threshold_noise" ;
    "rho_md" ;
    "rho_ss" ;
    "mean_log_prog_time";
    "std_log_prog_time"
  ]
  let res_order = [
    "mean" ;
    "q1" ;
    "median" ;
    "q3" ;
    "ss" ;
    "md" ;
    "gmgd" ;
    "an" ;
    "cc" ;
    "ggmggd" ;
    "gagn" ;
    "cc1r" ;
    "cc2" ;
    "growth_fraction" ;
  ]
  let write_csv filename pars results =
    let open Bigarray in
    let get1 ga = Genarray.get ga [|0|] in
    let fixed, running = List.take_drop 8 par_order in
    let fixed_row =
      fixed
      |> List.map (fun k ->
          match Hashtbl.find pars k with
          | `Int ga -> string_of_int (get1 ga)
          | `Float ga -> Format.(sprintf "%.4f" (get1 ga))) in
    let means, stds =
      let f = function
        | `Int _ -> failwith "internal error"
        | `Float a ->
          let a = array1_of_genarray a in
          let d = Array1.dim a in
          Array.init d (fun i -> Format.(sprintf "%.4f" a.{i})) in
      let l = List.map (fun k -> f (Hashtbl.find pars k)) running in
      List.(nth l 0, nth l 1) in
    let make_row_0 () =
      (* case of single-parameter lines *)
      let get ga = Genarray.get ga [|0|] in
      let ms = [means.(0); stds.(0)] in
      let f k = get (Hashtbl.find results k)
              |> Format.sprintf "%.4f" in
      let rs = List.map f res_order in
      fixed_row @ ms @ rs in
    let make_row i j =
      let getij i j ga = Genarray.get ga [|i; j|] in
      let ms = [means.(i); stds.(j)] in
      let f k = getij i j (Hashtbl.find results k)
                |> Format.sprintf "%.4f" in
      let rs = List.map f res_order in
      fixed_row @ ms @ rs in
    let csv = ref [] in
    if Array.length means = 1     (* crude fix for the single-parameter case *)
    then
      csv := make_row_0 () :: !csv
    else
    for i = 0 to Array.length means - 1 do
      for j = 0 to Array.length stds - 1 do
        csv := make_row i j :: !csv
      done done;
    let oc = open_out_gen [Open_append; Open_creat] 0o644 filename in
    let csv_oc = Csv.to_channel oc in
    Csv.output_all csv_oc !csv;
    Csv.close_out csv_oc
end

let matname = "pars_results"
let csvname = "pars_results"

let main_sample = function (`One _ | `Arr _ as aa ) ->
  let results = U.time ~silent:true collect_results aa in
  Format.printf "trees and statistics computed in %.2fs@." !U.last_timed;
  let pars = collect_params aa in
  pars, results

let main_store_mat pars results =
  let npzname =  Sys.getcwd () ^ "/" ^ matname ^ ".npz" in
  let save () =
    Export.write_npz npzname pars results;
    Export.external_convert_to_mat npzname;
    Unix.unlink npzname in
  save ()

let main_store_csv pars results =
  let csvname =  Sys.getcwd () ^ "/" ^ csvname ^ ".csv" in
  let save () =
    Export.write_csv csvname pars results in
  save ()
  (*U.time ~silent:true save ();*)
  (*Format.printf "statistics saved in %.2fs@." !U.last_timed*)


let _ =
  match read_argv Sys.argv with
    (`One _ | `Arr _ as aa ) ->
    let pars, results = main_sample aa in
    main_store_mat pars results;
    main_store_csv pars results
  | `Flat_list l ->
    l |> List.iter (fun p ->
        let pars, results = main_sample (`One p) in
        main_store_csv pars results)


