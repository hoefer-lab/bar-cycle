(* a script to submit jobs through the PBS batch processing system.
 * like submit_evidence.ml but for simulation runs of the growth-progression
 * model *)

(* use this directly on the server with ocamlrun ocaml submitter.ml *)
(* job directories will be created beneath the current directory! *)

(* submit a bunch of batch jobs for simulating the growth progression model.  *)

#use "topfind";;
#require "unix";;
#require "containers";;
#require "gen";;
open Format
open Containers

let () = Printexc.record_backtrace true;

(* handle parameter set looping {{{*)

type prange = {start:float; stop:float; step:float}

(* repeat so as to not depend on Util.. *)
let range {start; stop; step} =
  assert (step >. 0.);
  if Float.(is_nan start || is_nan stop) then [|nan|]
  else start
       |> Gen.unfold (fun v ->
           if v >. stop +. 1E-5 then None else Some (v, v +. step))
       |> Gen.to_array

let sr start = {start; stop=start; step=1.}

let pp_range fmt ({ start; stop; step } as pr) =
  let r = range pr in
  if Array.length r = 1 then
    fprintf fmt "@[[ %.4f ]@]" start
  else
    fprintf fmt "@[[ %.4f : %.4f : %.4f ]@]" start stop step

let in_range v ar =
  let epsilon = 1e-6 in
  let compare_with_margin =
    let eq = Float.equal_precision ~epsilon in
    fun f f' ->
      if eq f f' then 0
      else if f <. f' then -1
      else 1 in
  match Array.bsearch compare_with_margin v ar with
  | `At _i -> true
  | _ -> false

let not_all_in_some_range arll vl =
  not List.(exists (for_all2 in_range vl) arll)

type full_pars = {
  gen_init:int;
  gen_tree:int;
  n_trees:int;
  threshold_noise:prange;
  growth_rate:prange;
  saturation:float;
  rho_md:prange;
  rho_ss:prange;
  mean_log_prog_time:prange;
  std_log_prog_time:prange;
}

let pp_full_pars fmt
    { gen_init; gen_tree; n_trees; threshold_noise; growth_rate; saturation;
      rho_md; rho_ss; mean_log_prog_time; std_log_prog_time } =
  fprintf fmt
    "@[ gen_init:%d\n gen_tree:%d\n n_trees:%d\n threshold_noise:%a\n growth_rate:%a\n saturation:%.1f\n rho_md:%a\n rho_ss:%a\n mean_log_prog_time:%a\n std_log_prog_time:%a @]"
    gen_init gen_tree n_trees
    pp_range threshold_noise pp_range growth_rate saturation
    pp_range rho_md pp_range rho_ss
    pp_range mean_log_prog_time pp_range std_log_prog_time

let parrays
    {threshold_noise; growth_rate; saturation; rho_md;
     rho_ss; mean_log_prog_time; std_log_prog_time } =
  [threshold_noise; growth_rate; sr saturation; rho_md;
   rho_ss; mean_log_prog_time; std_log_prog_time ]
  |> List.map range

let n_parsets_loops fp =
  let ns = parrays fp
           |> List.map Array.length
           |> List.fold_left (fun ac n -> match ac with
               | [] -> [n]
               | h :: _ -> n * h :: ac) [] in
  List.hd ns, List.nth ns 2

let n_parsets fp = fst @@ n_parsets_loops fp
let n_loops fp = snd @@ n_parsets_loops fp
let n_inner_loop fp = n_parsets fp / n_loops fp

(* transpose *)
(*|> List.(fold_left (Fun.flip (map2 cons)) [[]; []; []; []; []]) *)


(* loop over parameters *)
let args_list_individual
    ?skip
    ({gen_init; gen_tree; n_trees} as fp : full_pars) =
  printf "looping over %d total parameter sets@." (n_parsets fp);
  let skip_arrays = match skip with
    | None -> []
    | Some spl ->
      spl |> List.map begin fun sp ->
        printf "skipping from %d parameter sets@." (n_parsets sp);
        parrays sp end in
  let outer_loop_lists = List.map Array.to_list (parrays fp) in
  let l = List.cartesian_product outer_loop_lists in
  assert (List.length l = n_parsets fp);
  let final_list =
  l          (* don't skip if no range exists which constains all parameters *)
  |> List.filter (not_all_in_some_range skip_arrays)
  |> List.map begin function
    | [threshold_noise; growth_rate; saturation; rho_md; rho_ss; mlp; slp] ->
      List.map string_of_int [gen_init; gen_tree; n_trees] @
      List.map (fun f -> sprintf "%.4f" f)
        [ threshold_noise; growth_rate; saturation; rho_md; rho_ss; mlp; slp]
    | _ -> failwith "incompatible list dimension, need 10 pars per line" end in
  let len = List.length final_list in
  printf "total of %d valid parameter sets@." (len);
  final_list

(* loop over parameters *)
let args_list
    ?skip
    ({gen_init; gen_tree; n_trees;
      mean_log_prog_time=mlp; std_log_prog_time=slp } as fp : full_pars) =
  printf "looping over %d total parameter sets in %d jobs@."
    (n_parsets fp) (n_loops fp);
  let skip_arrays = match skip with
    | None -> []
    | Some spl ->
      spl |> List.map begin fun sp ->
        printf "skipping from %d outer loops@." (n_loops sp);
        List.(take 5 (parrays sp)) end in
  let outer_loop_lists =
    List.(take 5 (map Array.to_list (parrays fp))) in
  let l =
    List.cartesian_product outer_loop_lists in
  assert (List.length l = n_loops fp);
  l          (* don't skip if no range exists which constains all parameters *)
  |> List.filter (not_all_in_some_range skip_arrays)
  |> List.map begin function
    | [threshold_noise; growth_rate; saturation; rho_md; rho_ss] ->
      List.map string_of_int [gen_init; gen_tree; n_trees] @
      List.map (fun f -> sprintf "%.4f" f)
        [ threshold_noise; growth_rate; saturation; rho_md; rho_ss;
          mlp.start; mlp.stop; mlp.step; slp.start; slp.stop; slp.step]
    | _ -> failwith "incompatible list dimension" end


let test_pars : full_pars =
  { gen_init=100;
    gen_tree=7;
    n_trees=20;
    threshold_noise = sr 0.06;
    growth_rate = sr 0.073;
    saturation = nan;
    rho_md = sr 0.88;
    rho_ss = sr 0.88;
    mean_log_prog_time = { start=2.10; stop=2.16; step=0.02 };
    std_log_prog_time = { start=0.30; stop=0.36; step=0.02 }
  }
(*}}}*)

(* Erika input {{{

Als Intervalle habe ich immer genommen:
k in 0.0005 steps  (Wachstumsrate)
sigG1 0.01 (threshold noise)
pmd 0.02 (mother-daughter korrelation)
pss 0.02 (sibling korrelation)
mu 0.02, aber 0.04 ist auch in Ordnung und schneller (Progression mean time)
sig 0.02, aber 0.04 ist auch gut (Progression noise)

Experiments:
ESC110930 - exp 3
pss <- 0.88
sigG1 <- 0.06
k <- 0.073
pmd <- 0.88
mu <-2.14
sig <- 0.32

let () = add_test "3esc" nan 0.88 0.06 0.073 0.88 2.14 0.32

ESC110907 - exp 1
pss <- 0.94
sigG1 <- 0.06
k <- 0.061
pmd <- 0.76
mu <- 2.38
sig <- 0.3

let () = add_test "1esc" nan 0.94 0.06 0.061 0.76 2.38 0.3

ESC110613 - exp2
pss <- 0.98
sigG1 <- 0.06
k <- 0.0695 #0.07 good, 0.069
pmd <- 0.88
mu <-2.18
sig <- 0.36

let () = add_test "2esc" nan 0.98 0.06 0.0695 0.88 2.18 0.36

Exp82rap
pss <- 0.99  (hier habe ich 0.99 auch noch mal mit rein genommen. Später aber eher 0.98 und 0.9999 getestet)
sigG1 <- 0.05
k <- 0.039
pmd <- 0.96
mu <-2.76
sig <- 0.36

let () = add_test "82rap" nan 0.99 0.05 0.039 0.96 2.76 0.36

Exp69rap
pss <- 0.9
sigG1 <- 0.05
k <- 0.03
pmd <- 0.86
mu <- 3.06
sig <- 0.28

let () = add_test "69rap" nan 0.9 0.05 0.03 0.86 3.06 0.28

Exp60on
pss <- 0.64
sigG1 <- 0
k <- 0.0455
pmd <- 0.28
mu <-2.62
sig <- 0.36

let () = add_test "60on" nan 0.64 0.0 0.0455 0.28 2.62 0.36

Exp62on
pss <- 0.64
sigG1 <- 0
k <- 0.042
pmd <- 0.36
mu <-2.68
sig <- 0.36

let () = add_test "62on" nan 0.64 0.0 0.042 0.36 2.68 0.36

Exp58on
pss <- 0.88
sigG1 <- 0.03
k <- 0.04
pmd <- 0.66
mu <- 2.76
sig <- 0.38

let () = add_test "58on" nan 0.88 0.03 0.04 0.66 2.76 0.38

Exp58off
pss <- 0.62
sigG1 <- 0.05
k <- 0.032 #0.032
pmd <- 0.44
mu <-3.14
sig <- 0.2

let () = add_test "58off" nan 0.62 0.05 0.032 0.44 3.14 0.2

Exp60off
pss <- 0.6
sigG1 <- 0.05
k <- 0.033
pmd <- 0.3
mu <-3.14
sig <- 0.24

let () = add_test "60off" nan 0.6 0.05 0.033 0.3 3.14 0.24

}}} *)

(* define parameter sets {{{ *)

let testparsets = Hashtbl.create 12

let add_test k saturation rho_ss threshold_noise growth_rate rho_md
    mean_log_prog_time std_log_prog_time =
  Hashtbl.add testparsets k
    {
      gen_init=100;
      gen_tree=7;
      n_trees=500;
      threshold_noise = sr threshold_noise; growth_rate = sr growth_rate;
      saturation;
      rho_md = sr rho_md; rho_ss = sr rho_ss;
      mean_log_prog_time = sr mean_log_prog_time;
      std_log_prog_time = sr std_log_prog_time }

let () = add_test "3esc" nan 0.88 0.06 0.073 0.88 2.14 0.32
let () = add_test "1esc" nan 0.94 0.06 0.061 0.76 2.38 0.3
let () = add_test "2esc" nan 0.98 0.06 0.0695 0.88 2.18 0.36
let () = add_test "82rap" nan 0.99 0.05 0.039 0.96 2.76 0.36
let () = add_test "69rap" nan 0.9 0.05 0.03 0.86 3.06 0.28
let () = add_test "60on" nan 0.64 0.0 0.0455 0.28 2.62 0.36
let () = add_test "62on" nan 0.64 0.0 0.042 0.36 2.68 0.36
let () = add_test "58on" nan 0.88 0.03 0.04 0.66 2.76 0.38
let () = add_test "58off" nan 0.62 0.05 0.032 0.44 3.14 0.2
let () = add_test "60off" nan 0.6 0.05 0.033 0.3 3.14 0.24

let () = add_test "3esc_s20"  20. 0.88 0.06 0.073 0.88 2.14 0.32
let () = add_test "1esc_s20"  20. 0.94 0.06 0.061 0.76 2.38 0.3
let () = add_test "2esc_s20"  20. 0.98 0.06 0.0695 0.88 2.18 0.36
let () = add_test "82rap_s20" 20. 0.99 0.05 0.039 0.96 2.76 0.36
let () = add_test "69rap_s20" 20. 0.9 0.05 0.03 0.86 3.06 0.28
let () = add_test "60on_s20"  20. 0.64 0.0 0.0455 0.28 2.62 0.36
let () = add_test "62on_s20"  20. 0.64 0.0 0.042 0.36 2.68 0.36
let () = add_test "58on_s20"  20. 0.88 0.03 0.04 0.66 2.76 0.38
let () = add_test "58off_s20" 20. 0.62 0.05 0.032 0.44 3.14 0.2
let () = add_test "60off_s20" 20. 0.6 0.05 0.033 0.3 3.14 0.24


let rangeparset = Hashtbl.create 12

let () =
  if false then
    Hashtbl.add rangeparset "test_58on_range" {
      (Hashtbl.find testparsets "58on") with
      growth_rate={ start=0.039; stop=0.041; step=0.0005 };
      threshold_noise={ start=0.02; stop=0.04; step=0.01 };
      mean_log_prog_time = { start=2.60; stop=2.80; step=0.04 };
      std_log_prog_time = { start=0.24; stop=0.44; step=0.04 };}

let () =
  Hashtbl.add rangeparset "exp_58on_range"
    { gen_init=100;
      gen_tree=7;
      n_trees=500;
      growth_rate={ start=0.039; stop=0.041; step=0.0005 };
      (* test *) (*growth_rate={ start=0.039; stop=0.041; step=0.001 };*)
      threshold_noise={ start=0.02; stop=0.04; step=0.01 };
      saturation=20.;
      rho_md={ start=0.5; stop=0.7; step=0.02 };
      rho_ss={ start=0.8; stop=0.9; step=0.02 };
      mean_log_prog_time={ start=2.64; stop=2.9; step=0.04 };
      std_log_prog_time={ start=0.18; stop=0.60; step=0.04 };
    }

let () =
  Hashtbl.add rangeparset "58on_range"
    { gen_init=100;
      gen_tree=7;
      n_trees=500;
      growth_rate={ start=0.039; stop=0.041; step=0.0005 };
      threshold_noise={ start=0.02; stop=0.04; step=0.01 };
      saturation=20.;
      rho_md={ start=0.5; stop=0.7; step=0.02 };
      rho_ss={ start=0.8; stop=0.9; step=0.02 };
      mean_log_prog_time={ start=2.64; stop=2.9; step=0.04 };
      std_log_prog_time={ start=0.18; stop=0.60; step=0.04 };
    }


let () =
  Hashtbl.add rangeparset "58on_range2"
    { gen_init=100;
      gen_tree=7;
      n_trees=500;
      growth_rate={ start=0.0390; stop=0.0410; step=0.0005 };
      threshold_noise={ start=0.00; stop=0.04; step=0.01 };
      saturation=20.;
      rho_md={ start=0.50; stop=0.70; step=0.02 };
      rho_ss={ start=0.64; stop=0.90; step=0.02 };
      mean_log_prog_time={ start=2.40; stop=2.90; step=0.04 };
      std_log_prog_time={ start=0.18; stop=0.60; step=0.04 };
    }


let () =
  Hashtbl.add rangeparset "58on_range3"
    { gen_init=100;
      gen_tree=7;
      n_trees=500;
      growth_rate={ start=0.0385; stop=0.0415; step=0.0005 };
      threshold_noise={ start=0.00; stop=0.06; step=0.01 };
      saturation=20.;
      rho_md={ start=0.50; stop=0.74; step=0.02 };
      rho_ss={ start=0.64; stop=0.90; step=0.02 };
      mean_log_prog_time={ start=2.40; stop=2.90; step=0.04 };
      std_log_prog_time={ start=0.18; stop=0.64; step=0.04 };
    }

let () =
  Hashtbl.add rangeparset "58on_range4"
    { gen_init=100;
      gen_tree=7;
      n_trees=500;
      growth_rate={ start=0.0402; stop=0.0408; step=0.0001 };
      threshold_noise={ start=0.00; stop=0.06; step=0.01 };
      saturation=20.;
      rho_md={ start=0.50; stop=0.74; step=0.02 };
      rho_ss={ start=0.64; stop=0.90; step=0.02 };
      mean_log_prog_time={ start=2.40; stop=2.90; step=0.04 };
      std_log_prog_time={ start=0.18; stop=0.64; step=0.04 };
    }

let () =
  Hashtbl.add rangeparset "58on_range5"
    { gen_init=100;
      gen_tree=7;
      n_trees=500;
      growth_rate={ start=0.0399; stop=0.0409; step=0.0001 };
      threshold_noise={ start=0.00; stop=0.06; step=0.01 };
      saturation=20.;
      rho_md={ start=0.50; stop=0.74; step=0.02 };
      rho_ss={ start=0.64; stop=0.90; step=0.02 };
      mean_log_prog_time={ start=2.40; stop=2.90; step=0.04 };
      std_log_prog_time={ start=0.18; stop=0.64; step=0.04 };
    }


let () =
  Hashtbl.add rangeparset "58on_range6"
    { gen_init=100;
      gen_tree=7;
      n_trees=500;
      growth_rate={ start=0.0396; stop=0.0398; step=0.0001 };
      threshold_noise={ start=0.00; stop=0.02; step=0.01 };
      saturation=20.;
      rho_md={ start=0.50; stop=0.7; step=0.02 };
      rho_ss={ start=0.7; stop=0.90; step=0.02 };
      mean_log_prog_time={ start=2.40; stop=2.84; step=0.04 };
      std_log_prog_time={ start=0.3; stop=0.7; step=0.04 };
    }

let () =
  Hashtbl.add rangeparset "5860off_range"
    { gen_init=100;
      gen_tree=7;
      n_trees=500;
      growth_rate={ start=0.03; stop=0.035; step=0.0005 };
      threshold_noise={ start=0.0; stop=0.1; step=0.02 };
      saturation=20.;
      rho_md={ start=0.2; stop=0.5; step=0.02 };
      rho_ss={ start=0.5; stop=0.76; step=0.02 };
      mean_log_prog_time={ start=3.0; stop=3.26; step=0.04 };
      std_log_prog_time={ start=0.06; stop=0.36; step=0.04 };
    }

let () =
  Hashtbl.add rangeparset "5860off_range2"
    { gen_init=100;
      gen_tree=7;
      n_trees=500;
      growth_rate={ start=0.03; stop=0.035; step=0.0005 };
      threshold_noise={ start=0.0; stop=0.1; step=0.02 };
      saturation=20.;
      rho_md={ start=0.20; stop=0.70; step=0.02 };
      rho_ss={ start=0.48; stop=0.80; step=0.02 };
      mean_log_prog_time={ start=2.84; stop=3.26; step=0.04 };
      std_log_prog_time={ start=0.06; stop=0.48; step=0.04 };
    }

let () =
  Hashtbl.add rangeparset "5860off_range3"
    { gen_init=100;
      gen_tree=7;
      n_trees=500;
      growth_rate={ start=0.0302; stop=0.0308; step=0.0001 };
      threshold_noise={ start=0.0; stop=0.1; step=0.02 };
      saturation=20.;
      rho_md={ start=0.20; stop=0.70; step=0.02 };
      rho_ss={ start=0.48; stop=0.80; step=0.02 };
      mean_log_prog_time={ start=2.84; stop=3.26; step=0.04 };
      std_log_prog_time={ start=0.06; stop=0.48; step=0.04 };
    }

let () =
  Hashtbl.add rangeparset "62on_range"
    { gen_init=100;
      gen_tree=7;
      n_trees=500;
      growth_rate={ start=0.04; stop=0.044; step=0.0005 };
      threshold_noise={ start=0.0; stop=0.04; step=0.01 };
      saturation=20.;
      rho_md={ start=0.24; stop=0.50; step=0.02 };
      rho_ss={ start=0.54; stop=0.74; step=0.02 };
      mean_log_prog_time={ start=2.54; stop=2.8; step=0.04 };
      std_log_prog_time={ start=0.24; stop=0.60; step=0.04 };
    }

let () =
  Hashtbl.add rangeparset "62on_range2"
    { gen_init=100;
      gen_tree=7;
      n_trees=500;
      growth_rate={ start=0.041; stop=0.043; step=0.0005 };
      threshold_noise={ start=0.0; stop=0.04; step=0.01 };
      saturation=20.;
      rho_md={ start=0.24; stop=0.50; step=0.02 };
      rho_ss={ start=0.54; stop=0.74; step=0.02 };
      mean_log_prog_time={ start=2.54; stop=2.9; step=0.04 };
      std_log_prog_time={ start=0.12; stop=0.60; step=0.04 };
    }

let () =
  Hashtbl.add rangeparset "62on_range3"
    { gen_init=100;
      gen_tree=7;
      n_trees=500;
      growth_rate={ start=0.041; stop=0.043; step=0.0005 };
      threshold_noise={ start=0.0; stop=0.04; step=0.01 };
      saturation=20.;
      rho_md={ start=0.16; stop=0.52; step=0.02 };
      rho_ss={ start=0.40; stop=0.74; step=0.02 };
      mean_log_prog_time={ start=2.34; stop=2.9; step=0.04 };
      std_log_prog_time={ start=0.12; stop=0.60; step=0.04 };
    }

let () =
  Hashtbl.add rangeparset "62on_range4"
    { gen_init=100;
      gen_tree=7;
      n_trees=500;
      growth_rate={ start=0.041; stop=0.043; step=0.0005 };
      threshold_noise={ start=0.0; stop=0.04; step=0.01 };
      saturation=20.;
      rho_md={ start=0.12; stop=0.52; step=0.02 };
      rho_ss={ start=0.40; stop=0.78; step=0.02 };
      mean_log_prog_time={ start=2.34; stop=2.9; step=0.04 };
      std_log_prog_time={ start=0.12; stop=0.64; step=0.04 };
    }

let () =
  Hashtbl.add rangeparset "62on_range5"
    { gen_init=100;
      gen_tree=7;
      n_trees=500;
      growth_rate={ start=0.0417; stop=0.0423; step=0.0001 };
      threshold_noise={ start=0.0; stop=0.04; step=0.01 };
      saturation=20.;
      rho_md={ start=0.12; stop=0.52; step=0.02 };
      rho_ss={ start=0.40; stop=0.78; step=0.02 };
      mean_log_prog_time={ start=2.34; stop=2.9; step=0.04 };
      std_log_prog_time={ start=0.12; stop=0.64; step=0.04 };
    }


let () =
  Hashtbl.add rangeparset "60on_range"
    { gen_init=100;
      gen_tree=7;
      n_trees=500;
      growth_rate={ start=0.044; stop=0.047; step=0.0005 };
      threshold_noise={ start=0.0; stop=0.04; step=0.01 };
      saturation=20.;
      rho_md={ start=0.14; stop=0.42; step=0.02 };
      rho_ss={ start=0.52; stop=0.78; step=0.02 };
      mean_log_prog_time={ start=2.42; stop=2.72; step=0.04 };
      std_log_prog_time={ start=0.08; stop=0.52; step=0.04 };
    }

let () =
  Hashtbl.add rangeparset "60on_range2"
    { gen_init=100;
      gen_tree=7;
      n_trees=500;
      growth_rate={ start=0.0475; stop=0.0485; step=0.0005 };
      threshold_noise={ start=0.0; stop=0.04; step=0.01 };
      saturation=20.;
      rho_md={ start=0.14; stop=0.42; step=0.02 };
      rho_ss={ start=0.52; stop=0.78; step=0.02 };
      mean_log_prog_time={ start=2.42; stop=2.72; step=0.04 };
      std_log_prog_time={ start=0.08; stop=0.52; step=0.04 };
    }

let () =
  Hashtbl.add rangeparset "60on_range3"
    { gen_init=100;
      gen_tree=7;
      n_trees=500;
      growth_rate={ start=0.0475; stop=0.0485; step=0.0005 };
      threshold_noise={ start=0.0; stop=0.04; step=0.01 };
      saturation=20.;
      rho_md={ start=0.06; stop=0.42; step=0.02 };
      rho_ss={ start=0.48; stop=0.78; step=0.02 };
      mean_log_prog_time={ start=2.34; stop=2.78; step=0.04 };
      std_log_prog_time={ start=0.08; stop=0.52; step=0.04 };
    }

let () =
  Hashtbl.add rangeparset "60on_range4"
    (* missing combinations *)
    { gen_init=100;
      gen_tree=7;
      n_trees=500;
      growth_rate={ start=0.0465; stop=0.0470; step=0.0005 };
      threshold_noise={ start=0.0; stop=0.04; step=0.01 };
      saturation=20.;
      rho_md={ start=0.06; stop=0.42; step=0.02 };
      rho_ss={ start=0.48; stop=0.78; step=0.02 };
      mean_log_prog_time={ start=2.34; stop=2.78; step=0.04 };
      std_log_prog_time={ start=0.08; stop=0.52; step=0.04 };
    }

let () =
  Hashtbl.add rangeparset "60on_range5"
    (* missing combinations *)
    { gen_init=100;
      gen_tree=7;
      n_trees=500;
      growth_rate={ start=0.0471; stop=0.0478; step=0.0001 };
      threshold_noise={ start=0.0; stop=0.04; step=0.01 };
      saturation=20.;
      rho_md={ start=0.06; stop=0.42; step=0.02 };
      rho_ss={ start=0.48; stop=0.78; step=0.02 };
      mean_log_prog_time={ start=2.34; stop=2.78; step=0.04 };
      std_log_prog_time={ start=0.08; stop=0.52; step=0.04 };
    }

let () =
  Hashtbl.add rangeparset "69rap_range"
    { gen_init=100;
      gen_tree=7;
      n_trees=500;
      growth_rate={ start=0.028; stop=0.032; step=0.0005 };
      threshold_noise={ start=0.02; stop=0.08; step=0.01 };
      saturation=20.;
      rho_md={ start=0.70; stop=0.96; step=0.02 };
      rho_ss={ start=0.82; stop=0.98; step=0.02 };
      mean_log_prog_time={ start=2.84; stop=3.18; step=0.04 };
      std_log_prog_time={ start=0.10; stop=0.50; step=0.04 };
    }

let () =
  Hashtbl.add rangeparset "69rap_range_inner_expand1"
    { (Hashtbl.find rangeparset "69rap_range") with
      mean_log_prog_time={ start=3.20; stop=3.28; step=0.04 };
    }

let () =
  Hashtbl.add rangeparset "69rap_range_inner_expand2"
    { (Hashtbl.find rangeparset "69rap_range") with
      mean_log_prog_time={ start=2.80; stop=2.80; step=0.04 };
    }


let () =
  Hashtbl.add rangeparset "69rap_range2"
    { gen_init=100;
      gen_tree=7;
      n_trees=500;
      growth_rate={ start=0.028; stop=0.032; step=0.00025 };
      threshold_noise={ start=0.02; stop=0.08; step=0.01 };
      saturation=20.;
      rho_md={ start=0.72; stop=0.96; step=0.02 };
      rho_ss={ start=0.74; stop=0.98; step=0.02 };
      mean_log_prog_time={ start=2.80; stop=3.28; step=0.04 };
      std_log_prog_time={ start=0.10; stop=0.50; step=0.04 };
    }

let () =
  Hashtbl.add rangeparset "69rap_range3"
    { gen_init=100;
      gen_tree=7;
      n_trees=500;
      growth_rate={ start=0.028; stop=0.032; step=0.00025 };
      threshold_noise={ start=0.00; stop=0.10; step=0.01 };
      saturation=20.;
      rho_md={ start=0.72; stop=0.96; step=0.02 };
      rho_ss={ start=0.74; stop=0.98; step=0.02 };
      mean_log_prog_time={ start=2.72; stop=3.28; step=0.04 };
      std_log_prog_time={ start=0.06; stop=0.58; step=0.04 };
    }

let () =
  Hashtbl.add rangeparset "69rap_range4"
    { gen_init=100;
      gen_tree=7;
      n_trees=500;
      growth_rate={ start=0.0302; stop=0.0308; step=0.0001 };
      threshold_noise={ start=0.00; stop=0.10; step=0.01 };
      saturation=20.;
      rho_md={ start=0.72; stop=0.96; step=0.02 };
      rho_ss={ start=0.74; stop=0.98; step=0.02 };
      mean_log_prog_time={ start=2.72; stop=3.28; step=0.04 };
      std_log_prog_time={ start=0.06; stop=0.58; step=0.04 };
    }

let () =
  Hashtbl.add rangeparset "82rap_range"
    { gen_init=100;
      gen_tree=7;
      n_trees=500;
      growth_rate={ start=0.037; stop=0.041; step=0.0005 };
      threshold_noise={ start=0.02; stop=0.09; step=0.01 };
      saturation=20.;
      rho_md={ start=0.88; stop=0.98; step=0.02 };
      rho_ss={ start=0.92; stop=0.98; step=0.02 };
      mean_log_prog_time={ start=2.56; stop=2.9; step=0.04 };
      std_log_prog_time={ start=0.10; stop=0.66; step=0.04 };
    }

let () =
  Hashtbl.add rangeparset "82rap_range2"
    { gen_init=100;
      gen_tree=7;
      n_trees=500;
      growth_rate={ start=0.037; stop=0.0425; step=0.0005 };
      threshold_noise={ start=0.01; stop=0.10; step=0.01 };
      saturation=20.;
      rho_md={ start=0.88; stop=0.98; step=0.02 };
      rho_ss={ start=0.88; stop=0.98; step=0.02 };
      mean_log_prog_time={ start=2.40; stop=2.98; step=0.04 };
      std_log_prog_time={ start=0.06; stop=0.70; step=0.04 };
    }

let () =
  Hashtbl.add rangeparset "82rap_range3"
    { gen_init=100;
      gen_tree=7;
      n_trees=500;
      growth_rate={ start=0.0396; stop=0.0404; step=0.0001 };
      threshold_noise={ start=0.01; stop=0.10; step=0.01 };
      saturation=20.;
      rho_md={ start=0.88; stop=0.98; step=0.02 };
      rho_ss={ start=0.88; stop=0.98; step=0.02 };
      mean_log_prog_time={ start=2.40; stop=2.98; step=0.04 };
      std_log_prog_time={ start=0.06; stop=0.70; step=0.04 };
    }

let () =
  Hashtbl.add rangeparset "82rap_range4"
    { gen_init=100;
      gen_tree=7;
      n_trees=500;
      growth_rate={ start=0.0385; stop=0.0485; step=0.0005 };
      threshold_noise={ start=0.01; stop=0.10; step=0.01 };
      saturation=20.;
      rho_md={ start=0.91; stop=0.99; step=0.01 };
      rho_ss={ start=0.93; stop=0.99; step=0.01 };
      mean_log_prog_time={ start=2.40; stop=2.98; step=0.04 };
      std_log_prog_time={ start=0.06; stop=0.70; step=0.04 };
    }

let () =
  Hashtbl.add rangeparset "esc_range"
    { gen_init=100;
      gen_tree=7;
      n_trees=500;
      growth_rate={ start=0.059; stop=0.075; step=0.0005 };
      threshold_noise={ start=0.03; stop=0.07; step=0.01 };
      saturation=20.;
      rho_md={ start=0.68; stop=0.94; step=0.02 };
      rho_ss={ start=0.80; stop=0.98; step=0.02 };
      mean_log_prog_time={ start=2.00; stop=2.46; step=0.04 };
      std_log_prog_time={ start=0.18; stop=0.60; step=0.04 };
    }

let () =
  Hashtbl.add rangeparset "esc_range2"
    { gen_init=100;
      gen_tree=7;
      n_trees=500;
      growth_rate={ start=0.059; stop=0.075; step=0.0005 };
      threshold_noise={ start=0.01; stop=0.09; step=0.01 };
      saturation=20.;
      rho_md={ start=0.62; stop=0.98; step=0.02 };
      rho_ss={ start=0.74; stop=0.98; step=0.02 };
      mean_log_prog_time={ start=1.88; stop=2.46; step=0.04 };
      std_log_prog_time={ start=0.18; stop=0.68; step=0.04 };
    }

let () =
  Hashtbl.add rangeparset "esc_range3"
    { gen_init=100;
      gen_tree=7;
      n_trees=500;
      growth_rate={ start=0.059; stop=0.075; step=0.0005 };
      threshold_noise={ start=0.01; stop=0.10; step=0.01 };
      saturation=20.;
      rho_md={ start=0.62; stop=0.98; step=0.02 };
      rho_ss={ start=0.74; stop=0.98; step=0.02 };
      mean_log_prog_time={ start=1.80; stop=2.46; step=0.04 };
      std_log_prog_time={ start=0.18; stop=0.68; step=0.04 };
    }

let () =
  Hashtbl.add rangeparset "esc1_range"
    { gen_init=100;
      gen_tree=7;
      n_trees=500;
      growth_rate={ start=0.0618; stop=0.0622; step=0.0001 };
      threshold_noise={ start=0.01; stop=0.10; step=0.01 };
      saturation=20.;
      rho_md={ start=0.62; stop=0.98; step=0.02 };
      rho_ss={ start=0.74; stop=0.98; step=0.02 };
      mean_log_prog_time={ start=1.80; stop=2.46; step=0.04 };
      std_log_prog_time={ start=0.18; stop=0.68; step=0.04 };
    }

(* debug *)
(*threshold_noise={ start=0.01; stop=0.10; step=0.01 };*)
(*rho_md={ start=0.62; stop=0.98; step=0.02 };*)
let () =
  Hashtbl.add rangeparset "esc2_range"
    { gen_init=100;
      gen_tree=7;
      n_trees=500;
      growth_rate={ start=0.0698; stop=0.0702; step=0.0001 };
      threshold_noise={ start=0.01; stop=0.01; step=0.01 };
      saturation=20.;
      rho_md={ start=0.62; stop=0.98; step=0.02 };
      rho_ss={ start=0.74; stop=0.98; step=0.02 };
      mean_log_prog_time={ start=1.80; stop=2.46; step=0.04 };
      std_log_prog_time={ start=0.18; stop=0.68; step=0.04 };
    }

let () =
  Hashtbl.add rangeparset "esc3_range"
    { gen_init=100;
      gen_tree=7;
      n_trees=500;
      growth_rate={ start=0.0733; stop=0.0737; step=0.0001 };
      threshold_noise={ start=0.01; stop=0.10; step=0.01 };
      saturation=20.;
      rho_md={ start=0.62; stop=0.98; step=0.02 };
      rho_ss={ start=0.74; stop=0.98; step=0.02 };
      mean_log_prog_time={ start=1.80; stop=2.46; step=0.04 };
      std_log_prog_time={ start=0.18; stop=0.68; step=0.04 };
    }

(*}}}*)

let code_dir = String.concat "/" [Sys.getenv "HOME"; "bar-cycle";]
let collection_subdir = "data_growth_progression"
let program = "simulate_growth_progression"
let csvname = "pars_results.csv"
let code_type =
  (* byte only for debugging *)
  (*"byte"*)
  "native"

let pars_filename = "pars"
let pars_chunk_per_job = 1_800
let cat_chunk = 100

(* pbs parameters *)
let jobname = program
let processors = 1
let hours = 1
let mem = 4


let shell = ["#!/bin/bash"]
let owd = "$PBS_O_WORKDIR"
let cmd_prefix = "time"


(* compile {{{*)

let compiled_program_to_run =
  (* keep the executable safe from overwriting *)
  sprintf "%d" (Hashtbl.hash (Unix.time ()))

let compile_sequence wd =
  [ "rm -f *.{byte,native}";
    "ocamlbuild -clean";
    sprintf "ocamlbuild -use-ocamlfind %s.%s" program code_type;
    sprintf "mkdir -p %s/%s" wd collection_subdir;
    sprintf "cp %s.%s %s/%s/%s"
      program code_type wd collection_subdir compiled_program_to_run;
    sprintf "git rev-parse HEAD > %s/%s/headrev" wd collection_subdir;]

(* no need to compile for each run:
 * this is handled by the command line arguments *)
let compile_first () =
  let wd = Sys.getcwd () in
  let () = Sys.chdir code_dir in
  let cmd = (String.concat " && " (compile_sequence wd)) in
  let () = print_endline cmd in
  let success = Sys.command cmd in
  let () = Sys.chdir wd in
  success

(* }}} *)

(* submit scripts {{{ *)

let submit_file script =
  let open Unix in
  let scriptfile = open_out "scrp.sh" in
  output_string scriptfile (String.concat "\n" script);
  close_out scriptfile;
  let chan = open_process_in "qsub scrp.sh" in
  let id = input_line chan in
  let () = print_endline id in
  let status = close_process_in chan in
  let status_code = match status with
      WEXITED s | WSIGNALED s | WSTOPPED s -> s in
  Unix.sleepf 0.01;
  (status_code, id)

let submit_pipe script =
  let open Unix in
  let script = String.concat "\n" script in
  let chan = open_process_in
      (* "heredoc" crazyness: *)
      (sprintf "qsub - <<\"EOF\"\n%s\nEOF\n" script) in
  let id = input_line chan in
  let () = print_endline id in
  let status_code = match close_process_in chan with
      WEXITED s | WSIGNALED s | WSTOPPED s -> s in
  Unix.sleepf 0.01;
  (status_code, id)

(* pbs script for a single inner loop given by args *)
let script expname args =
  let jobdir = sprintf "%s/%s/$PBS_JOBID" collection_subdir expname in
  let pbs_directives = [
    sprintf "-N %s" jobname;
    sprintf "-l walltime=%d:00:00" hours;
    sprintf "-l mem=%dgb" mem;
    (* these are the nodes that have the same processor as curry0 *)
    (*sprintf "-l nodes=curry60-0+curry60-1+curry82-0+curry82-1";*)
    (* bigmem excludes the oldest nodes; this should be enough. *)
    sprintf "-l nodes=1:bigmem:ppn=%d" processors;
    sprintf "%s" "-V";                                  (* environment setup *)
    sprintf "-e %s/err" jobdir;
    sprintf "-o %s/out" jobdir;]
    |> List.map ((^) "#PBS ") in
  let dirs = [
    sprintf "cd %s" owd;
    sprintf "mkdir -p %s" jobdir;
    sprintf "cd %s" jobdir;] in
  let pars =
    [ sprintf "echo '%s' >> %s" (String.concat " " args) pars_filename] in
  (* prevent depending on completed jobs if jobs are too fast *)
  let wait = ["sleep 10"] in
  let command =
    [cmd_prefix;
     sprintf "%s/%s/%s" owd collection_subdir compiled_program_to_run]
    @ args
    |> String.concat " "
    |> (fun x -> [x]) in
  List.concat [
    shell;
    pbs_directives;
    dirs;
    pars;
    wait;
    command; ]

(* pbs script for sequential evaluation of a list of arg vectors *)
let seq_script expname (args_list: string list list) =
  let jobdir = sprintf "%s/%s/$PBS_JOBID" collection_subdir expname in
  let email = {|nils.becker@bioquant.uni-heidelberg.de|} in
  let pbs_directives = [
    sprintf "-N %s" jobname;
    sprintf "-l walltime=72:00:00";
    sprintf "-l mem=%dgb" mem;
    (* these are the nodes that have the same processor as curry0: *)
    (*sprintf "-l nodes=curry60-0+curry60-1+curry82-0+curry82-1";*)
    (* bigmem excludes the oldest nodes; this should be enough. *)
    sprintf "-l nodes=1:bigmem:ppn=%d" processors;
    sprintf "%s" "-V";                                  (* environment setup *)
    sprintf "-e %s/err" jobdir;
    sprintf "-o %s/out" jobdir;
    (* debug *)
    sprintf "-m a";                                            (* abort only *)
    sprintf "-M %s" email; ]
    |> List.map ((^) "#PBS ") in
  let dirs = [
    sprintf "cd %s" owd;
    sprintf "mkdir -p %s" jobdir;
    sprintf "cd %s" jobdir;] in
  let pars =
    args_list
    |> List.map (String.concat " ")
    |> List.map (fun p -> sprintf "echo '%s' >> %s" p pars_filename) in
  (* prevent depending on completed jobs if jobs are too fast *)
  let wait = ["sleep 10"] in
  (* parse the input from the pars file to avoid command line length
   * limitations! *)
  let command =
    [sprintf "%s %s/%s/%s %s"
       cmd_prefix
       owd collection_subdir compiled_program_to_run
       pars_filename] in
  List.concat [
    shell;
    pbs_directives;
    dirs;
    pars;
    wait;
    command; ]


let concat_script expname dependids jobids =
  let jobdir id = sprintf "%s/%s/%s" collection_subdir expname id in
  let depend = String.concat ":" (jobids @ dependids) in
  let pbs_directives = [
    sprintf "-N %s%s" jobname "_concat";
    sprintf "-l walltime=00:30:00";
    sprintf "-e %s/%s/concat.err" collection_subdir expname;
    sprintf "-o %s/%s/concat.out" collection_subdir expname;
    sprintf "-W depend=afterok:%s" depend;]
    |> List.map ((^) "#PBS ") in
  let dirs = [ sprintf "cd %s" owd ] in
  let cats = jobids |> List.map (fun id ->
      let j, c = jobdir id, csvname in
      sprintf "cat %s/%s >> %s/../%s" j c j c) in
  let wait = ["sleep 10"] in
  List.concat [
    shell;
    pbs_directives;
    dirs;
    wait;
    cats; ]

(* send email after all are completed *)
let completion_script catids =
  let depend = String.concat ":" catids in
  let email = {|nils.becker@bioquant.uni-heidelberg.de|} in
  let pbs_directives = [
    sprintf "-N %s%s" jobname "_cleanup";
    sprintf "-l walltime=00:01:00";
    sprintf "-e %s/%s/cleanup.err" owd collection_subdir;
    sprintf "-o %s/%s/cleanup.out" owd collection_subdir;
    sprintf "-m ae";                                        (* abort and end *)
    sprintf "-M %s" email;
    sprintf "-W depend=afterany:%s" depend;]
    |> List.map ((^) "#PBS ") in
  let command =                                         (* remove executable *)
    [ sprintf "rm -f %s/%s/%s" owd collection_subdir compiled_program_to_run ;
      sprintf "rm -f %s/scrp.sh" owd ] in
  List.concat [
    shell;
    pbs_directives;
    command; ]

(*}}}*)


(* looping {{{*)
let sublists l = List.sublists_of_len ~last:Option.return l

let loop_exp exp pars skip cids =
    (* chunking to keep the depend lists below the shell limit *)
    let args_list_chunks = args_list ~skip pars
                           (* try chunks of 100 to fit into the pipe script *)
                           |> sublists 100 in
    let last_catid =
      let loop_chunk prev_id al =
        let jobids = al |> List.map (fun args ->
            let stat, id = submit_pipe (script exp args) in
            if stat <> 0 then failwith
                (sprintf "failed submitting work job (%d)" stat);
            id) in
        let catstat, catid =                       (* keep sequential concat *)
          submit_pipe (concat_script exp (Option.to_list prev_id) jobids) in
        if catstat <> 0 then failwith
            (sprintf "failed submitting concat job (%d)" catstat);
        Some catid in
      List.fold_left loop_chunk None args_list_chunks in
    cids @ Option.to_list last_catid

let loop_individual exp pars skip cids =
  let args_list_jobs =
    args_list_individual ~skip pars
    |> sublists pars_chunk_per_job in
  printf "\nsubmitting %d jobs\n@." (List.length args_list_jobs);
  (*printf "continue?\n@.";*)
  (*assert (String.equal_caseless (read_line ()) "y");*)
  let exp_list_id =
    let loop_job prev_ids al =
      (* too long for pipe solution? *)
      let stat, id = submit_file (seq_script exp al) in
      if stat <> 0 then failwith
          (sprintf "failed submitting job (%d)" stat);
      id :: prev_ids in
    List.fold_left loop_job [] args_list_jobs in
  let catids =
    let cat_chunks = sublists cat_chunk exp_list_id in
    List.fold_left (fun cids chunk_ids ->
        (* this could get above 100 dependencies, use submit_file *)
        let new_cstat, new_cid =
          submit_file (concat_script exp cids chunk_ids) in
        if new_cstat <> 0 then failwith
            (sprintf "failed submitting concat job (%d)" new_cstat);
        new_cid :: cids)
      [] cat_chunks in
  catids

(*}}}*)


(* run {{{*)
;;
if compile_first () > 0 then failwith "compilation unsuccessful!"
;;
printf "work dir %s@." (Sys.getcwd ());
printf "collection dir %s/%s@." (Sys.getcwd ()) collection_subdir;
;;
printf "reserving %d cores per job@." processors;
;;
let exp_cat_ids =
  List.fold_left2
    (fun cl exp sko ->
       let pars = Hashtbl.find rangeparset exp in
       let skip = List.map (Hashtbl.find rangeparset) sko in
       loop_individual
         exp pars skip cl)
    []

    (*["58on_range2"] [["58on_range"]] *)
    (*["69rap_range"] [[]]*)
    (*["60on_range2"; ] [["60on_range"];]*)
    (*["62on_range2"] [[]]*)
    (*["esc_range"] [[]]*)
    (*["69rap_range_inner_expand1"; "69rap_range_inner_expand2"] [[]; []]*)
    (*["82rap_range"] [[]]*)
    (*["62on_range3"; "60on_range3"] [["62on_range2"]; ["60on_range2"]]*)
    (*["58on_range3"] [["58on_range2"]]*)
    (*["69rap_range3"] [["69rap_range2"]]*)
    (*["82rap_range2"] [["82rap_range"]]*)
    (*["62on_range4"] [["62on_range3"]]*)
    (*["5860off_range2"] [["5860off_range"]]*)
    (*["esc_range2"] [["esc_range"]]*)
    (*["58on_rangetest"] [["58on_range3"]]*)
    (*["60on_range4"; ] [[]]*)
    (*["69rap_range3"] [[]]*)
    (*["esc_range3"] [["esc_range2"]] *)
    (*["60on_range5"; ] [["60on_range4"; "60on_range3"];]*)
    (*["62on_range5"; ] [["62on_range4"];]*)
    (*["69rap_range4"] [["69rap_range3"]]*)
    (*["82rap_range3"] [["82rap_range2"]]*)
    (*["58on_range4"] [["58on_range3"]]*)
    (*["5860off_range3"] [["5860off_range2"]]*)
    (*["esc1_range"] [[]]*)
    (*["esc2_range"] [[]]*)
    (*["esc3_range"] [[]]*)
    (*["58on_range5"] [["58on_range4"]]*)
    (*["82rap_range4"] [["82rap_range3"]]*)
    ["58on_range6"] [["58on_range5"]]
;;
(* the concat ids are just as many as the experiments in the outermost loop,
 * can submit via pipe. *)
submit_pipe (completion_script exp_cat_ids)

(* vim: set foldmethod=marker: *)
