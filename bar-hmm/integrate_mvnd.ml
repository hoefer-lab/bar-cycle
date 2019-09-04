(* wrapper around the fortan mvndst library from Genz. this must be compiled
 * beforehand as described in the subdirectory multivariate_normal_int and the
 * shared library must be present in $HOME/lib/
 *
 * the correct compilation call is in
 * "multivariate_normal_int/c compat command line"
 * *)


(* this seems to work fine now, after adding the appropriate variable types etc
 * to the FORTRAN source. However, it is slow for calculating the overhang
 * probabilities, since it is an algorithm for general coupled normal variates.
 *)

open Bigarray
open Lacaml.D
open Ctypes
open! Containers
open Util

(*[@@@landmark "auto"]*)
(*let integration = Landmark.register "integration"*)

(* make the integration tolerances dynamic settings *)
type integration_setting = {mutable rel_error: float;
                            mutable abs_error:float;
                            mutable max_pts_per_dim:int;
                            mutable neg_log_clip:float}
let settings =
  {rel_error=1e-2;
   abs_error=1e-6;
   max_pts_per_dim=1000;
   neg_log_clip= -.infinity}

(* get the shared library *)
let os =
  assert (String.equal Sys.os_type "Unix");
  let ic = Unix.open_process_in "uname" in
  let uname = input_line ic in
  let () = close_in ic in
  if String.equal uname "Darwin" then `MacOs else `Linux
let mvnd =
  let filename = Sys.getenv "HOME" ^ "/lib/mvndst." ^
                 (match os with `MacOs -> "dylib" | _ -> "so") in
  Dl.dlopen ~flags:Dl.([RTLD_NOW;]) ~filename


module C_compat = struct
  (* adjustments on the fortran side *)
  let integer = Ctypes.nativeint
  let symbol_name = "mvnormaldist"
  let ap v = bigarray_start array1 v
  let ip k = allocate integer (Nativeint.of_int k)
  let fp f = allocate double f
  let i32_to_int = Nativeint.of_int32
  let make_one_int_array n =
    let a = Array1.create Nativeint c_layout n in
    Array1.fill a 1n; a
  let zero_int = 0n
  let toint = Nativeint.to_int
end

module C_compat2 = struct
  (* alternative adjustments made on the fortran side *)
  let integer = Ctypes.int
  let symbol_name = "mvnormaldist"
  let ap v = bigarray_start array1 v
  let ip k = allocate integer k
  let fp f = allocate double f
  let i32_to_int = Int32.to_int
  let make_one_int_array n =
    let a = Array1.create Int c_layout n in
    Array1.fill a 1; a
  let zero_int = 0
  let toint = Fun.id
end

module Original = struct
  (* original f77 program *)
  let integer = Ctypes.int32_t
  let symbol_name = "mvndst_"
  let ap v = bigarray_start array1 v
  let ip k = allocate integer (Int32.of_int k)
  let fp f = allocate double f
  let i32_to_int = Fun.id
  let make_one_int_array n =
    let a = Array1.create Int32 c_layout n in
    Array1.fill a 1l; a
  let zero_int = 0l
  let toint = Int32.to_int
end


(* choose which! *)

(*open C_compat*)
(*open Original*)
open C_compat2                               (* this version works correctly *)


(* the integration routine which gives the integral in value *)
let integrate =
  Foreign.foreign ~from:mvnd
    symbol_name
    (ptr integer @-> (* N *)
     ptr double  @-> (* LOWER *)
     ptr double  @-> (* UPPER *)
     ptr integer @-> (* INFIN *)
     ptr double  @-> (* CORREL *)
     ptr integer @-> (* MAXPTS *)
     ptr double  @-> (* ABSEPS *)
     ptr double  @-> (* RELEPS *)
     ptr double  @-> (* ERROR *)
     ptr double  @-> (* VALUE *)
     ptr integer @-> (* INFORM *)
     returning void);;

(* c_layout arrays to get the start address *)
let c_make0 n =
  let a = Array1.create float64 c_layout n in
  Array1.fill a 0.;
  a

let v2c fa =
  fa |> genarray_of_array1
  |> (fun fga -> Genarray.change_layout fga c_layout)
  |> array1_of_genarray

(* arrange a symmetric positive matrix in the way the library expects *)
let lower_vec mat =
  let d = Mat.dim1 mat in
  let v = Vec.create ((d*(d-1)) / 2) in
  for i = 1 to d do for j = 1 to i - 1 do
      v.{j + ((i-2)*(i-1))/2} <- mat.{i,j}
    done done;
  v

(* correlation and variances *)
let corr_vars ci =
  (* covariance matrix. po_inv is ok, only lower triangle is used *)
  let c = U.po_inv ci in
  let vars = Mat.copy_diag c in
  let istds = vars |> Vec.sqrt |> Vec.reci in
  Mat.scal_cols c istds;
  Mat.scal_rows istds c;
  lower_vec c, vars

exception Integration_error of
    [ `Not_enough_points
    | `Wrong_dimensionality
    | `Result_nan
    | `Result_inf
    | `Internal ]

type int_stats =
  { mutable evaluations       : int;
    mutable not_enough_points : int;
    mutable smaller0          : int }
let integration_stats =
  {evaluations=0; not_enough_points=0; smaller0=0}

(* attention: this is the POSITIVE log of the integrated probability *)
(* tested - works for 300 dim identity matrix. scaling of variances works as
 * well with the correct sign. *)
let log_gauss_integral
    ?lower ?upper
    ?infin
    ?maxpts
    ?releps ?abseps
    ci =
  let (@||) v default = Option.get_or ~default v in
  let n = Mat.dim1 ci in
  let correl, vars = corr_vars ci in
  (* no jacobian needed: this is in the normalization already *)
  let l, u = lower @|| Vec.make0 n, upper @|| Vec.make0 n in
  if not Vec.(dim l = n && dim u = n) then
    raise (Invalid_argument "need full starting vector");
  let sigmas = Vec.sqrt vars in
  (* now scale them to the unit variance correl used in the integration *)
  let lower, upper = Vec.(div l sigmas, div u sigmas) in
  let infin = infin @|| make_one_int_array n in
  let maxpts = maxpts @|| n * settings.max_pts_per_dim in
  let releps = releps @|| settings.rel_error in
  let abseps = abseps @|| settings.abs_error in
  (* pointers *)
  let _N = ip n in
  let _LOWER, _UPPER, _INFIN = ap (v2c lower), ap (v2c upper), ap infin in
  let _CORREL = ap correl in
  let _MAXPTS, _ABSEPS, _RELEPS = ip maxpts, fp abseps, fp releps in
  let _ERROR,  _VALUE, _INFORM = fp 0., fp 0., ip 0 in
  (*Landmark.enter integration;*)
  let integral = integrate                                         (* do it! *)
      _N
      _LOWER _UPPER _INFIN
      _CORREL
      _MAXPTS _ABSEPS _RELEPS
      _ERROR _VALUE _INFORM;
    begin match toint !@_INFORM with
      | 0 -> ()
      | 1 ->
        integration_stats.not_enough_points <-
          integration_stats.not_enough_points + 1;
        Format.printf
          "not enough evaluation points (%d per dim), error estimate %f@."
          settings.max_pts_per_dim !@_ERROR;
        (* don't raise, continue after recording the incident *)
        (*raise @@ Integration_error `Not_enough_points*)
      | 2  ->
        Format.printf "too low or high dimensionality";
        raise @@ Integration_error `Wrong_dimensionality
      | _  ->
        Format.printf "internal error in gaussian integration routine";
        raise @@ Integration_error `Internal
    end;
    !@_VALUE in
  (*Landmark.exit integration;*)
  (*for some miraculous reason without these checks, get errors from gsl!*)
  if Float.is_nan integral then
    raise @@ Integration_error `Result_nan;
  if not (integral <. Float.max_finite_value) then
    raise @@ Integration_error `Result_inf;
  if integral <. 0. then
    integration_stats.smaller0 <- integration_stats.smaller0 + 1;
  (*raise (Integration_error "< 0.");*)
  let integral = Float.max (exp settings.neg_log_clip) integral in
  (* -inf may be clipped *)
  let res = log integral in
  (*Format.printf "variances %a@." Lacaml.Io.pp_fvec vars;*)
  (*Format.printf "res %.5g@." res;*)
  integration_stats.evaluations <- integration_stats.evaluations + 1;
  res

