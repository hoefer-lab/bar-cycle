(* a utility module for building trees from raw csv data *)
open! Containers

open Tree
open Hmtree

(* new version: take into account:
 * survival,
 * walking out
 * cell death *)

(* basic record for a sister pair; corresponds to a row in the csv *)
type pair_rec =
  {fam: int;
   tau1: float; taumin1: float;
   tau2: float; taumin2: float;
   tbirth: float;
   taum: float; taugm: float;
   generation: int;
   id1: int; idm: int; idgm: int; id2: int};;

(* for the trees from erikas experiments *)
let read_row row =
  let get f = function
    | h::t -> f h, t
    | [] -> raise Not_found in
  (* work around csv/export suckiness *)
  let f = get float_of_string in
  let i = get (fun s -> int_of_float (float_of_string s))
  in
  (* the following is a bit silly. should have converted to string array. *)
  let fam, row = i row in
  let tau1, row = f row in
  let row = List.tl row in
  let taumin1, row = f row in
  let tau2, row = f row in
  let row = List.tl row in
  let taumin2, row = f row in
  let tbirth, row = f row in
  let taum, row = f row in
  let taugm, row = f row in
  let generation, row = i row in
  let id1, row = i row in
  let idm, row = i row in
  let idgm, row = i row in
  let id2, _row = i row in
  {fam;
   tau1; taumin1; tau2; taumin2; tbirth;
   taum; taugm;
   generation;
   id1; idm; idgm; id2}

(* for the trees from tim schröder *)
let read_row_esc_cut row =
  let get f = function
    | h::t -> f h, t
    | [] -> raise Not_found in
  (* work around csv/export suckiness *)
  let f = get float_of_string in
  let i = get (fun s -> int_of_float (float_of_string s))
  in
  let fam, row = i row in
  let tau1, row = f row in
  let tau2, row = f row in
  let id1, row = i row in
  let id2, row = i row in
  let tbirth, row = f row in
  let idm, row = i row in
  let generation, _row = i row in
  {fam;
   tau1; taumin1=nan; tau2; taumin2=nan; tbirth;
   taum=nan; taugm=nan;
   generation;
   id1; idm; idgm=0; id2}

(* cell time tree printer *)
let print_tt = T.print print_tau

(* cell zipper printer *)
let print_tz = Z.print print_tt


(* cell tree printer *)
let print_crt = T.print print_cr

(* cell zipper printer *)
let print_crz fmt zp = Format.fprintf fmt "%a" print_crt (Z.tree zp)

(* utility for grouping *)
let groupsort by l =
  let cmp r r' = compare (by r) (by r') in
  let eq r r' = 0 = cmp r r' in
  List.group_succ ~eq (List.sort cmp l)

let build_tree pairlist =
  (* input validation *)
  if (pairlist
      |> List.map (fun r -> r.generation)
      |> List.fold_left min max_int) <> 2
  then None
  else Some begin
      (* worker function *)
      let append_child tau taumin tbirth id zp =
        match Float.(is_nan tau, is_nan taumin) with
        | true, true -> zp                                     (* cell death *)
        | other -> let tr = match other with
            | false, true -> T.create {tau=Val tau; tbirth; id}
            | true, false -> T.create {tau=Min taumin; tbirth; id}
            | _ -> raise @@ Invalid_argument "more than one lifetime?" in
          zp
          |> Z.append_child tr
          |> Z.parent |> Option.get_exn in
      (* root node *)
      let root =
        let top_pair =
          List.find_pred_exn (fun r -> r.generation = 2) pairlist in
        let root_cell =
          {tau=Val top_pair.taum;
           tbirth=top_pair.tbirth -. top_pair.taum;            (* nan anyway *)
           id=top_pair.idm} in
        (*let () = Printf.printf "root: %d\n" root_cell.id in*)     (* debug *)
        Z.zipper @@ T.create root_cell in
      (* now build by generation *)
      let rec grow_leaves zip rl gen =
        let newrecs, later = List.partition (fun r -> r.generation = gen) rl in
        match newrecs with
        | [] -> Z.root zip
        | _ ->
          let zip' = List.fold_left
              begin fun zp prec ->
                zp |> Z.root |> Z.leavez (* all leaves, not only current subtree *)
                |> Sequence.find_pred_exn (fun r ->
                    let cur_id = Z.(zval r).id in
                    cur_id = prec.idm)
                |> append_child prec.tau1 prec.taumin1 prec.tbirth prec.id1
                |> append_child prec.tau2 prec.taumin2 prec.tbirth prec.id2
              end
              zip newrecs in
          grow_leaves zip' later (gen + 1) in
      grow_leaves root pairlist 2 end

(* forest import from erikas csv files *)
let exp_forest ?(read_row=read_row) csv_file =
  let tab = Csv.load csv_file in
  let _cols, vals = List.(hd tab, tl tab) in
  let sspairs = List.map read_row vals in
  (* families grouped together *)
  let families = groupsort (fun r -> r.fam) sspairs in
  let forest = List.filter_map build_tree families in
  forest

(* abbreviations for the ESC experiments *)
let esc_names = [
  1, "110907";
  2, "110613";
  3, "110930"]
