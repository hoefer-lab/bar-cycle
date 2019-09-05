(* here we collect all model definitions to be compared for evidence. *)

open! Containers
open Lacaml.D

open Util

module Gslv = Gsl.Vector


(* utility: tanh scaling *)

(* scale y from  -∞, ∞  to min, max *)
let scale_par (min, max) y =
  let (s, s') = 0.5 *. (max +. min), 0.5 *. (max -. min) in
  s +. s' *. tanh y

(* and back, checked. *)
let inv_scale_par (min, max) x =
  if x >. max || x <. min
  then raise (Invalid_argument "parameter out of scaling bounds");
  let (s, s') = 0.5 *. (max +. min), 0.5 *. (max -. min) in
  Gsl.Math.atanh ((x -. s) /. s')

(* utility: printing *)

let print_fmin_vec fmin res_par =
  Format.printf "\n min reached: %.3f\n\n parameter vector: @.%a\n@." fmin
    pp_vec (Vec.of_array @@ Gslv.to_array res_par)

let print_mats (a, alpha, cia) =
  Format.printf "parameter matrices:\na@.%a\n@." Lacaml.Io.Toplevel.pp_fmat a;
  Format.printf "\nalpha@.%a@." pp_vec alpha;
  Format.printf "\nca@.%a\n@." pp_mat (U.ge_inv cia)


(* all model parameter sets are records; they each require their own type. this
 * is somewhat verbose and could probably be done with less boilerplate... *)

type m1_indep_r =
      {a11: float * float as 'a;
       alpha1: 'a;}
type m1_r  =
      {a11: float * float as 'a;
       alpha1: 'a;
       gamma_: 'a}
type m2_indep_indep_r =
      {a11: float * float as 'a;
       a22: 'a;
       alpha1: 'a;
       alpha2: 'a;}
type m2_indep_r =
      {a11: float * float as 'a;
       a22: 'a;
       alpha1: 'a;
       alpha2: 'a;
       gamma_: 'a}
type m2_indep_crossnoise_r =
      {a11: float * float as 'a;
       a22: 'a;
       alpha1: 'a;
       alpha2: 'a;
       chi_: 'a;
       gamma_: 'a}
type m2_indep_crossnoise_indep_r =
      {a11: float * float as 'a;
       a22: 'a;
       alpha1: 'a;
       alpha2: 'a;
       chi_: 'a;}
type m2_uni_same_r =
      {a11: float * float as 'a;
       a12: 'a;
       alpha1: 'a;
       alpha2: 'a;
       gamma_: 'a}
type m2_uni_only_r =
      {a12: float * float as 'a;
       alpha1: 'a;
       alpha2: 'a;
       gamma_: 'a}
type m2_uni_a1_r =
      {a11: float * float as 'a;
       a22: 'a;
       a12: 'a;
       alpha1: 'a;
       gamma_: 'a}
type m2_uni_indep_r =
      {a11: float * float as 'a;
       a22: 'a;
       a12: 'a;
       alpha1: 'a;
       alpha2: 'a}
type m2_uni_r =
      {a11: float * float as 'a;
       a22: 'a;
       a12: 'a;
       alpha1: 'a;
       alpha2: 'a;
       gamma_: 'a}
type m2_uni_crossnoise_r =
      {a11: float * float as 'a;
       a22: 'a;
       a12: 'a;
       alpha1: 'a;
       alpha2: 'a;
       chi_: 'a;
       gamma_: 'a}
type m2_uni_a22_r =
      {a22: float * float as 'a;
       a12: 'a;
       alpha1: 'a;
       alpha2: 'a;
       gamma_: 'a}
type m2_bi_indep_r =
      {a11: float * float as 'a;
       a22: 'a;
       a12: 'a;
       a21: 'a;
       alpha1: 'a;
       alpha2: 'a}
type m2_bi_r =
      {a11: float * float as 'a;
       a22: 'a;
       a12: 'a;
       a21: 'a;
       alpha1: 'a;
       alpha2: 'a;
       gamma_: 'a}

type model_ranges =
  | M1_indep of m1_indep_r
  | M1 of m1_r
  | M2_indep_indep of m2_indep_indep_r
  | M2_indep of m2_indep_r
  | M2_indep_crossnoise of m2_indep_crossnoise_r
  | M2_indep_crossnoise_indep of m2_indep_crossnoise_indep_r
  | M2_uni_same of m2_uni_same_r
  | M2_uni_only of m2_uni_only_r
  | M2_uni_a1 of m2_uni_a1_r
  | M2_uni_indep of m2_uni_indep_r
  | M2_uni of m2_uni_r
  | M2_uni_crossnoise of m2_uni_crossnoise_r
  | M2_uni_a22 of m2_uni_a22_r
  | M2_bi_indep of m2_bi_indep_r
  | M2_bi of m2_bi_r


type model = {
  (* scaled-parameter ranges *)
  lo : float array;
  up : float array;
  (* parameter count *)
  pdim : int;
  (* hidden variables count *)
  hdim : int;
  (* map -∞, ∞ par vector to the specified ranges and inverse *)
  scaled_pars : Gslv.vector -> float array;
  unconstrained_pars : float array -> Gslv.vector;
  (* build a matrix from scaled and inverse *)
  build_matrices_scaled : float array -> mat * vec * mat;
  extract_scaled : mat * vec * mat -> float array;
  (* chain the above and inverse *)
  build_matrices : Gslv.vector -> mat * vec * mat;
  make_parvec : mat * vec * mat -> Gslv.vector;
  (* finally the typed ranges *)
  ranges : model_ranges
}


let model_parnames model =
  (match model.ranges with
   | M1_indep _ ->
     ["a11"; "alpha1"]
   | M1 _ ->
     ["a11"; "alpha1"; "gamma"]
   | M2_indep_indep _ ->
     ["a11"; "a22"; "alpha1"; "alpha2"]
   | M2_indep _ ->
     ["a11"; "a22"; "alpha1"; "alpha2"; "gamma"]
   | M2_indep_crossnoise _ ->
     ["a11"; "a22"; "alpha1"; "alpha2"; "chi"; "gamma"]
   | M2_indep_crossnoise_indep _ ->
     ["a11"; "a22"; "alpha1"; "alpha2"; "chi"]
   | M2_uni_same _ ->
     ["a11"; "a12"; "alpha1"; "alpha2" ]
   | M2_uni_only _ ->
     ["a12"; "alpha1"; "alpha2" ]
   | M2_uni_a1 _ ->
     ["a11"; "a22"; "a12"; "alpha1"; "gamma" ]
   | M2_uni_indep _ ->
     ["a11"; "a22"; "a12"; "alpha1"; "alpha2" ]
   | M2_uni _ ->
     ["a11"; "a22"; "a12"; "alpha1"; "alpha2"; "gamma" ]
   | M2_uni_crossnoise _ ->
     ["a11"; "a22"; "a12"; "alpha1"; "alpha2"; "chi"; "gamma" ]
   | M2_uni_a22 _ ->
     ["a22"; "a12"; "alpha1"; "alpha2"; "gamma" ]
   | M2_bi_indep _ ->
     ["a11"; "a22"; "a12"; "a21"; "alpha1"; "alpha2" ]
   | M2_bi _ ->
     ["a11"; "a22"; "a12"; "a21"; "alpha1"; "alpha2"; "gamma" ])



(* evaluation of the likelihood *)

(* generic function to make a model module. the only things we really need to know are
 * how arrange parameters into matrices and get them back *)
let prepare_model hdim build_matrices_scaled extract_scaled ranges =
  let ra =  (match ranges with
      | M1_indep { a11; alpha1 } ->
        [|a11; alpha1|]
      | M1 { a11; alpha1; gamma_ } ->
        [|a11 ;alpha1; gamma_|]
      | M2_indep_indep { a11; a22; alpha1; alpha2 } ->
        [|a11; a22; alpha1; alpha2|]
      | M2_indep { a11; a22; alpha1; alpha2; gamma_ } ->
        [|a11; a22; alpha1; alpha2; gamma_|]
      | M2_indep_crossnoise { a11; a22; alpha1; alpha2; chi_; gamma_ } ->
        [|a11; a22; alpha1; alpha2; chi_; gamma_|]
      | M2_indep_crossnoise_indep { a11; a22; alpha1; alpha2; chi_ } ->
        [|a11; a22; alpha1; alpha2; chi_|]
      | M2_uni_indep { a11; a22; a12; alpha1; alpha2} ->
        [|a11; a22; a12; alpha1; alpha2|]
      | M2_uni_same { a11; a12; alpha1; alpha2; gamma_ } ->
        [|a11; a12; alpha1; alpha2; gamma_|]
      | M2_uni_only { a12; alpha1; alpha2; gamma_ } ->
        [|a12; alpha1; alpha2; gamma_|]
      | M2_uni_a1 { a11; a22; a12; alpha1; gamma_ } ->
        [|a11; a22; a12; alpha1; gamma_|]
      | M2_uni { a11; a22; a12; alpha1; alpha2; gamma_ } ->
        [|a11; a22; a12; alpha1; alpha2; gamma_|]
      | M2_uni_crossnoise { a11; a22; a12; alpha1; alpha2; chi_; gamma_ } ->
        [|a11; a22; a12; alpha1; alpha2; chi_; gamma_|]
      | M2_uni_a22 { a22; a12; alpha1; alpha2; gamma_ } ->
        [| a22; a12; alpha1; alpha2; gamma_|]
      | M2_bi_indep { a11; a22; a12; a21; alpha1; alpha2} ->
        [|a11; a22; a12; a21; alpha1; alpha2|]
      | M2_bi { a11; a22; a12; a21; alpha1; alpha2; gamma_ } ->
        [|a11; a22; a12; a21; alpha1; alpha2; gamma_|]) in
  let hdim = hdim in
  let lo, up = Array.(map fst ra, map snd ra) in
  let pdim = Array.length lo in
  let scaled_pars parvec = Array.map2 scale_par
      ra (Gslv.to_array parvec) in
  let unconstrained_pars sc =
    Array.map2 inv_scale_par ra sc |> Gslv.of_array in
  let build_matrices parvec =
    parvec |> scaled_pars |> build_matrices_scaled in
  let make_parvec (a, alpha, cia) =
    (a, alpha, cia) |> extract_scaled |> unconstrained_pars in
  {lo; up; pdim; hdim;
   scaled_pars; unconstrained_pars;
   build_matrices; build_matrices_scaled;
   make_parvec; extract_scaled;
   ranges}

let check_111 ca =
  assert Float.(equal_precision ~epsilon:(10.*.epsilon) ca.{1,1} 1.)


(* the original model *)

let rec model_2_uni ranges = model_2_uni_standard ranges

and model_2_uni_standard (ranges:m2_uni_r) =
  (* the parvec has the meaning
     [a11, a22, a12, alpha1, alpha2, gamma] *)
  let build_matrices_scaled sc =
    let a = Mat.make0 2 2 in
    a.{1,1} <- sc.(0);
    a.{2,2} <- sc.(1);
    a.{1,2} <- sc.(2);
    let alpha = Vec.of_array [|sc.(3); sc.(4)|] in
    let bd, bo = Mat.(identity 2, identity 2) in
    Mat.scal sc.(5) bo;
    let ca = U.block_matrix [|[|bd; bo|]; [|bo; bd|]|] in
    let cia = U.ge_inv ca in
    (a, alpha, cia) in
  let extract_scaled (a, alpha, cia) =
    let ca = U.ge_inv cia in
    check_111 ca;
    let sc = Array.of_list
        [a.{1,1};
         a.{2,2};
         a.{1,2};
         alpha.{1};
         alpha.{2};
         ca.{1,3};] in
    sc in
  prepare_model 2 build_matrices_scaled extract_scaled (M2_uni ranges)

and model_2_uni_normalized12 (ranges:m2_uni_r) =
  (* the parvec has the meaning
     [a11, a22, a12, alpha1, alpha2, gamma] *)
  let build_matrices_scaled sc =
    let a = Mat.make0 2 2 in
    a.{1,1} <- sc.(0);
    a.{2,2} <- sc.(1);
    (* try fitting with unscaled a12 *)
    let alpha = Vec.of_array [|sc.(3); sc.(4)|] in
    a.{1,2} <- sc.(2) *. alpha.{2} /. alpha.{1} ;
    let bd, bo = Mat.(identity 2, identity 2) in
    Mat.scal sc.(5) bo;
    let ca = U.block_matrix [|[|bd; bo|]; [|bo; bd|]|] in
    let cia = U.ge_inv ca in
    (a, alpha, cia) in
  let extract_scaled (a, alpha, cia) =
    let ca = U.ge_inv cia in
    check_111 ca;
    let sc = Array.of_list
        [a.{1,1};
         a.{2,2};
         (* unscaled a12 see above *)
         a.{1,2} *. alpha.{1} /. alpha.{2};
         alpha.{1};
         alpha.{2};
         ca.{1,3};] in
    sc in
  prepare_model 2 build_matrices_scaled extract_scaled (M2_uni ranges)

(* some more models *)

exception Bad_model of string

(* one-dimensional without sister *)
let model_1_indep ranges =
  (* the parvec has the meaning
     [a11, alpha1] *)
  let build_matrices_scaled sc =
    let a = Mat.make0 1 1 in
    a.{1,1} <- sc.(0);
    let alpha = Vec.of_array [|sc.(1)|] in
    let ca = Mat.identity 2 in
    let cia = U.ge_inv ca in
    (a, alpha, cia) in
  let extract_scaled (a, alpha, _cia) =
    let sc = Array.of_list
        [a.{1,1};
         alpha.{1};] in
    sc in
  prepare_model 1 build_matrices_scaled extract_scaled (M1_indep ranges)

(* one-dimensional with sister *)
let model_1 ranges =
  (* the parvec has the meaning
     [a11, alpha1, gamma] *)
  let build_matrices_scaled sc =
    let a = Mat.make0 1 1 in
    a.{1,1} <- sc.(0);
    let alpha = Vec.of_array [|sc.(1)|] in
    let bd, bo = Mat.(identity 1, identity 1) in
    Mat.scal sc.(2) bo;
    let ca = U.block_matrix [|[|bd; bo|]; [|bo; bd|]|] in
    let cia = U.ge_inv ca in
    (a, alpha, cia) in
  let extract_scaled (a, alpha, cia) =
    let ca = U.ge_inv cia in
    check_111 ca;
    let sc = Array.of_list
        [a.{1,1};
         alpha.{1};
         ca.{1,2};] in
    sc in
  prepare_model 1 build_matrices_scaled extract_scaled (M1 ranges)

(* no sister coupling and no coupled inheritance.*)
let model_2_indep_indep ranges =
   (* the parvec has the meaning
    * [a11, a22, alpha1, alpha2] *)
  let build_matrices_scaled sc =
    let a = Mat.make0 2 2 in
    a.{1,1} <- sc.(0);
    a.{2,2} <- sc.(1);
    let alpha = Vec.of_array [|sc.(2); sc.(3)|] in
    let ca = Mat.identity 4 in
    let cia = U.ge_inv ca in
    (a, alpha, cia) in
  let extract_scaled (a, alpha, cia) =
    let ca = U.ge_inv cia in
    check_111 ca;
    let sc = Array.of_list
        [a.{1,1};
         a.{2,2};
         alpha.{1};
         alpha.{2};] in
    sc in
  prepare_model 2 build_matrices_scaled extract_scaled (M2_indep_indep ranges)

(* no cross inheritance.*)
let model_2_indep ranges =
   (* the parvec has the meaning
    * [a11, a22, alpha1, alpha2, gamma] *)
  let build_matrices_scaled sc =
    let a = Mat.make0 2 2 in
    a.{1,1} <- sc.(0);
    a.{2,2} <- sc.(1);
    let alpha = Vec.of_array [|sc.(2); sc.(3)|] in
    let bd, bo = Mat.(identity 2, identity 2) in
    Mat.scal sc.(4) bo;
    let ca = U.block_matrix [|[|bd; bo|]; [|bo; bd|]|] in
    let cia = U.ge_inv ca in
    (a, alpha, cia) in
  let extract_scaled (a, alpha, cia) =
    let ca = U.ge_inv cia in
    check_111 ca;
    let sc = Array.of_list
        [a.{1,1};
         a.{2,2};
         alpha.{1};
         alpha.{2};
         ca.{1,3};] in
    sc in
  prepare_model 2 build_matrices_scaled extract_scaled
    (M2_indep ranges)

(* no cross inheritance. *)
let model_2_indep_crossnoise ranges =
   (* the parvec has the meaning
    * [a11, a22, alpha1, alpha2, chi, gamma] *)
  let build_matrices_scaled sc =
    let a = Mat.make0 2 2 in
    a.{1,1} <- sc.(0);
    a.{2,2} <- sc.(1);
    let alpha = Vec.of_array [|sc.(2); sc.(3)|] in
    let bd = Mat.identity 2 in
    bd.{1,2} <- sc.(4); bd.{2,1} <- sc.(4);
    let bo = lacpy bd in
    Mat.scal sc.(5) bo;
    let ca = U.block_matrix [|[|bd; bo|]; [|bo; bd|]|] in
    let cia = U.ge_inv ca in
    (a, alpha, cia) in
  let extract_scaled (a, alpha, cia) =
    let ca = U.ge_inv cia in
    let sc = Array.of_list
        [a.{1,1};
         a.{2,2};
         alpha.{1};
         alpha.{2};
         ca.{1,2};
         ca.{1,3} /. ca.{1,1};] in
    sc in
  prepare_model 2 build_matrices_scaled extract_scaled
    (M2_indep_crossnoise ranges)

(* no cross inheritance. *)
let model_2_indep_crossnoise_indep ranges =
  (* the parvec has the meaning
   * [a11, a22, alpha1, alpha2, chi] *)
  let build_matrices_scaled sc =
    let a = Mat.make0 2 2 in
    a.{1,1} <- sc.(0);
    a.{2,2} <- sc.(1);
    let alpha = Vec.of_array [|sc.(2); sc.(3)|] in
    let bd, bo = Mat.(identity 2, make0 2 2) in
    bd.{1,2} <- sc.(4); bd.{2,1} <- sc.(4);
    let ca = U.block_matrix [|[|bd; bo|]; [|bo; bd|]|] in
    let cia = U.ge_inv ca in
    (a, alpha, cia) in
  let extract_scaled (a, alpha, cia) =
    let ca = U.ge_inv cia in
    let sc = Array.of_list
        [a.{1,1};
         a.{2,2};
         alpha.{1};
         alpha.{2};
         ca.{1,2};] in
    sc in
  prepare_model 2 build_matrices_scaled extract_scaled
    (M2_indep_crossnoise_indep ranges)

let model_2_uni_same (ranges:m2_uni_same_r) =
  (* the parvec has the meaning
     [a11, a12, alpha1, alpha2, gamma] *)
  let build_matrices_scaled sc =
    let a = Mat.make0 2 2 in
    a.{1,1} <- sc.(0);
    a.{2,2} <- sc.(0);                                      (* same diagonal *)
    a.{1,2} <- sc.(1);
    let alpha = Vec.of_array [|sc.(2); sc.(3)|] in
    let bd, bo = Mat.(identity 2, identity 2) in
    Mat.scal sc.(4) bo;
    let ca = U.block_matrix [|[|bd; bo|]; [|bo; bd|]|] in
    let cia = U.ge_inv ca in
    (a, alpha, cia) in
  let extract_scaled (a, alpha, cia) =
    let ca = U.ge_inv cia in
    check_111 ca;
    let sc = Array.of_list
        [a.{1,1};
         a.{1,2};
         alpha.{1};
         alpha.{2};
         ca.{1,3};] in
    sc in
  prepare_model 2 build_matrices_scaled extract_scaled (M2_uni_same ranges)

(* only coupling, no persistence of the individual variables. *)
let model_2_uni_only (ranges:m2_uni_only_r) =
  (* the parvec has the meaning
     [a12, alpha1, alpha2, gamma] *)
  let build_matrices_scaled sc =
    let a = Mat.make0 2 2 in
    a.{1,1} <- 0.;
    a.{2,2} <- 0.;                                      (* same diagonal *)
    a.{1,2} <- sc.(0);
    let alpha = Vec.of_array [|sc.(1); sc.(2)|] in
    let bd, bo = Mat.(identity 2, identity 2) in
    Mat.scal sc.(3) bo;
    let ca = U.block_matrix [|[|bd; bo|]; [|bo; bd|]|] in
    let cia = U.ge_inv ca in
    (a, alpha, cia) in
  let extract_scaled (a, alpha, cia) =
    let ca = U.ge_inv cia in
    check_111 ca;
    let sc = Array.of_list
        [a.{1,2};
         alpha.{1};
         alpha.{2};
         ca.{1,3};] in
    sc in
  prepare_model 2 build_matrices_scaled extract_scaled (M2_uni_only ranges)

let model_2_uni_a1 (ranges:m2_uni_a1_r) =
  (* the parvec has the meaning
     [a11, a22, a12, alpha1, gamma] *)
  let build_matrices_scaled sc =
    let a = Mat.make0 2 2 in
    a.{1,1} <- sc.(0);
    a.{2,2} <- sc.(1);
    a.{1,2} <- sc.(2);
    let alpha = Vec.of_array [|sc.(3); 0.|] in
    let bd, bo = Mat.(identity 2, identity 2) in
    Mat.scal sc.(4) bo;
    let ca = U.block_matrix [|[|bd; bo|]; [|bo; bd|]|] in
    let cia = U.ge_inv ca in
    (a, alpha, cia) in
  let extract_scaled (a, alpha, cia) =
    let ca = U.ge_inv cia in
    check_111 ca;
    let sc = Array.of_list
        [a.{1,1};
         a.{2,2};
         a.{1,2};
         alpha.{1};
         ca.{1,3};] in
    sc in
  prepare_model 2 build_matrices_scaled extract_scaled (M2_uni_a1 ranges)

let model_2_uni_indep ranges =
  (* the parvec has the meaning
     [a11, a22, a12, alpha1, alpha2, gamma] *)
  let build_matrices_scaled sc =
    let a = Mat.make0 2 2 in
    a.{1,1} <- sc.(0);
    a.{2,2} <- sc.(1);
    a.{1,2} <- sc.(2);
    let alpha = Vec.of_array [|sc.(3); sc.(4)|] in
    let bd, bo = Mat.(identity 2, make0 2 2) in
    let ca = U.block_matrix [|[|bd; bo|]; [|bo; bd|]|] in
    let cia = U.ge_inv ca in
    (a, alpha, cia) in
  let extract_scaled (a, alpha, _cia) =
    let sc = Array.of_list
        [a.{1,1};
         a.{2,2};
         a.{1,2};
         alpha.{1};
         alpha.{2};] in
    sc in
  prepare_model 2 build_matrices_scaled extract_scaled
    (M2_uni_indep ranges)

let model_2_uni_crossnoise ranges =
  (* the parvec has the meaning
     [a11, a22, a12, alpha1, alpha2, chi, gamma] *)
  let _ = raise (Bad_model "overparametrized") in
  let build_matrices_scaled sc =
    let a = Mat.make0 2 2 in
    a.{1,1} <- sc.(0);
    a.{2,2} <- sc.(1);
    a.{1,2} <- sc.(2);
    let alpha = Vec.of_array [|sc.(3); sc.(4)|] in
    let bd, bo =
      let chi = sc.(5) in
      Mat.(
        of_list [[1.; chi]; [chi; 1.]],
        of_list [[1.; chi]; [chi; 1.]]) in
    Mat.scal sc.(6) bo;
    let ca = U.block_matrix [|[|bd; bo|]; [|bo; bd|]|] in
    let cia = U.ge_inv ca in
    (a, alpha, cia) in
  let extract_scaled (a, alpha, cia) =
    let ca = U.ge_inv cia in
    check_111 ca;
    let sc = Array.of_list
        [a.{1,1};
         a.{2,2};
         a.{1,2};
         alpha.{1};
         alpha.{2};
         ca.{1,2};
         ca.{1,3};] in
    sc in
  prepare_model 2 build_matrices_scaled extract_scaled
    (M2_uni_crossnoise ranges)

let model_2_uni_a22 ranges =
  (* the parvec has the meaning
     [a22, a12, alpha1, alpha2, gamma] *)
  let build_matrices_scaled sc =
    let a = Mat.make0 2 2 in
    a.{2,2} <- sc.(0);
    a.{1,2} <- sc.(1);
    let alpha = Vec.of_array [|sc.(2); sc.(3)|] in
    let bd, bo = Mat.(identity 2, identity 2) in
    Mat.scal sc.(4) bo;
    let ca = U.block_matrix [|[|bd; bo|]; [|bo; bd|]|] in
    let cia = U.ge_inv ca in
    (a, alpha, cia) in
  let extract_scaled (a, alpha, cia) =
    let ca = U.ge_inv cia in
    check_111 ca;
    let sc = Array.of_list
        [a.{2,2};
         a.{1,2};
         alpha.{1};
         alpha.{2};
         ca.{1,3};] in
    sc in
  prepare_model 2 build_matrices_scaled extract_scaled
    (M2_uni_a22 ranges)

let model_2_bi_indep ranges =
  (* the parvec has the meaning
     [a11, a22, a12, a21, alpha1, alpha2, gamma].
     Attention: the A matrix may get eigenvalues that are bigger than 1 in
     absolute value. in the integration: have to return 0... *)
  let _ = raise (Bad_model "overparametrized") in
  let build_matrices_scaled sc =
    let a = Mat.make0 2 2 in
    a.{1,1} <- sc.(0);
    a.{2,2} <- sc.(1);
    a.{1,2} <- sc.(2);
    a.{2,1} <- sc.(3);
    let alpha = Vec.of_array [|sc.(4); sc.(5)|] in
    let bd, bo = Mat.(identity 2, make0 2 2) in
    let ca = U.block_matrix [|[|bd; bo|]; [|bo; bd|]|] in
    let cia = U.ge_inv ca in
    (a, alpha, cia) in
  let extract_scaled (a, alpha, _cia) =
    let sc = Array.of_list
        [a.{1,1};
         a.{2,2};
         a.{1,2};
         a.{2,1};
         alpha.{1};
         alpha.{2};] in
    sc in
  prepare_model 2 build_matrices_scaled extract_scaled
    (M2_bi_indep ranges)

let model_2_bi ranges =
  (* the parvec has the meaning
     [a11, a22, a12, a21, alpha1, alpha2, gamma].
     Attention: the A matrix may get eigenvalues that are bigger than 1 in
     absolute value. in the integration: have to return 0... *)
  (* allow this for evidence calculation test *)
  Format.printf "\n Attention: Model 2_bi is overparameterized!@.";
  (*let _ = raise (Bad_model "overparametrized") in*)
  let build_matrices_scaled sc =
    let a = Mat.make0 2 2 in
    a.{1,1} <- sc.(0);
    a.{2,2} <- sc.(1);
    a.{1,2} <- sc.(2);
    a.{2,1} <- sc.(3);
    let alpha = Vec.of_array [|sc.(4); sc.(5)|] in
    let bd, bo = Mat.(identity 2, identity 2) in
    Mat.scal sc.(6) bo;
    let ca = U.block_matrix [|[|bd; bo|]; [|bo; bd|]|] in
    let cia = U.ge_inv ca in
    (a, alpha, cia) in
  let extract_scaled (a, alpha, cia) =
    let ca = U.ge_inv cia in
    check_111 ca;
    let sc = Array.of_list
        [a.{1,1};
         a.{2,2};
         a.{1,2};
         a.{2,1};
         alpha.{1};
         alpha.{2};
         ca.{1,3};] in
    sc in
  prepare_model 2 build_matrices_scaled extract_scaled
    (M2_bi ranges)

