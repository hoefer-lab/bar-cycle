(* This is a utility script to transform raw cell cycle tree data to
 * a Mathematica-compatible format. The trees are looked for in the directory
 * "$HOME/tmp/data_trees/precut/".
 * Output files are written to "$HOME/tmp/data_trees_mma/"
 * *)

open Containers
open Lacaml.D
open Util
open Tree
open Hmtree
open Import_csv
open Test_lib
open Models

let () = Printexc.record_backtrace true;;

let rng_seed =
    Nativeint.of_string Sys.(argv.(Array.length argv - 1))
  |> Option.get_or ~default:1234n

let () = Gsl.Rng.set the_rng (rng_seed)
let results_dir =
  Filename.concat (Unix.getenv "HOME") "tmp/data_trees_mma/"
let () =
  let open Unix in
  try
    mkdir results_dir 0o755;
    Printf.printf "made directory\n\n"
  with (Unix_error (EEXIST, _, _)) -> ()


(* real or reduced real or simulated or pooled or transf-pooled trees? *)
let (data : [ `Pooled_then_transformed
            | `Transformed_then_pooled
            | `Real
            | `Real_half
            | `Sim
            | `ESC ]) =
  `Real
  (*`Pooled_then_transformed*)
  (*`Transformed_then_pooled*)
  (*`ESC*)

let (keys_to_use :
       [ `All
       | `Select of ExpMap.key list]) =
  `All
  (*`Select [(82, "rap")]*)
  (*`Select [(586062, "on")]*)
  (*`Select [(5862, "on")]*)
  (*`Select [(5860, "off")]*)
  (*`Select [(6982, "rap")]*)
  (*`Select [(62, "on")]*)

(* outliers *)
(*old: `Quantile 0.95*)
let (outlier_treatment :
       [ `Normal | `Quantile of float | `Size of int | `Log_only ]) =
  (*`Size 9*)
  `Normal
  (*`Quantile 0.98*)

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

let sim_forests = ExpMap.(empty |> add (58, "on_sim") pforest)


(* data import. new: import _everything_ here! cut down later. *)

let dir =
  Filename.concat (Unix.getenv "HOME") "tmp/data_trees/precut/"

let all_forests =
  let open ExpMap in
  let f name = exp_forest
      (dir ^ name ^ "_Trees_globalgen.csv") in
  empty
  |> add (62, "on") (f "Exp62on")
  |> add (58, "on") (f "Exp58on")
  |> add (60, "on") (f "Exp60on")
  |> add (58, "off") (f "Exp58off")
  |> add (60, "off") (f "Exp60off")
  |> add (69, "rap") (f "Exp69rap20nM")
  |> add (82, "rap") (f "Exp82rap40nM")
  (*|> add (62, "on_old") (f "Exp62on_old")*)
  |> ExpMap.map (List.filter_map Ctree.prune_survivors_cells)


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

let esc_data =
  let open ExpMap in
  empty |> add (1, "esc")
    (exp_forest ~read_row:read_row_esc_cut (dir ^ "ESCTrees_all_SiblingRows.csv"))

(* these are pooled before any subtraction of the mean etc *)
let pooled_raw_forests = pool_data all_forests


(* reject outliers *)

let preprocess m =
  let transformed = begin match outlier_treatment with
    | `Quantile p ->         (* typical value 0.95. this can change results! *)
      m |> ExpMap.map (quantile_prune p)
      |> ExpMap.map log_forest
    | `Size l ->                                         (* typical value: 9 *)
      m |> ExpMap.map (List.filter (fun z -> Z.length z > l))
      |> ExpMap.map log_forest
    | `Normal ->
      m |> ExpMap.map normal_forest
    | `Log_only -> (* no rejection. not recommended for evidence calculation *)
      m |> ExpMap.map log_forest
  end in
  transformed
  |> ExpMap.map (List.filter (fun z -> Z.length z > 1))

let transformed_forests = all_forests |> preprocess

let pooled_then_transformed_forests = pooled_raw_forests |> preprocess

let transformed_esc_forests = esc_data |> preprocess

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
  | `ESC                     -> transformed_esc_forests

(* now finally, cut down to the data we actually want to process *)
let the_forests = let open ExpMap in
  match keys_to_use with
  | `All -> the_forests'
  | `Select l ->
    List.fold_left (fun em k -> add k (find k the_forests') em) empty l


(* utility *)
let pp_cellrec_tau = function
  | {Hmtree.tau=Min _; _ } -> "nan"
  | {tau=Val f; _} -> Format.sprintf "%a" U.pf f

let pp_cellrec_tb_tau = function
  | {Hmtree.tau=Min _; tbirth} -> "nan"
  | {tau=Val f; tbirth} -> Format.sprintf "{%a, %a}" U.pf tbirth U.pf f

let pp_cellrec_tb = function
  | {Hmtree.tau=Min _; tbirth} -> "nan"
  | {tau=Val f; tbirth} -> Format.sprintf "%a" U.pf tbirth

let pp_cellrec_td = function
  | {Hmtree.tau=Min _; tbirth} -> "nan"
  | {tau=Val f; tbirth} -> Format.sprintf "%a" U.pf (tbirth +. f)

let z_to_mma p z = Ctree.to_mma p (Z.tree z)

let save_forests pp_val prefix emap =
  let f (exp_no, cond) v =
    CCIO.with_out
      (results_dir ^ String.concat "_" [prefix; string_of_int exp_no; cond] ^ ".m")
      (fun oc ->
         let tree_strings = List.map (z_to_mma pp_val) v in
         output_string oc
           ("{" ^ String.concat ", " tree_strings ^ "}")) in
  ExpMap.iter f emap

(* save stuff *)

let () = save_forests pp_cellrec_tau "raw" all_forests

let () = save_forests
    (fun f -> Format.sprintf "%a" U.pf f)
    "processed"
    the_forests

let () = save_forests pp_cellrec_tb_tau "tb_raw" all_forests

(* not so useful: sisters have the same time of birth *)
let () = save_forests pp_cellrec_tb "tb" all_forests

let () = save_forests pp_cellrec_td "td" all_forests


