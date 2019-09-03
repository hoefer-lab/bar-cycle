(* this module extends RoseTree by S. Cruanes; it collects extensions to
 * general tree manipulations, without cell-cycle specific ingredients. those
 * are in Ctree. *)

open! Containers
open Option.Infix
open Util

module T = struct
  (*[@@@landmark "auto"]*)
  include RoseTree

  let xval (`Node (a, _)) = a

  let iter f t = to_seq t |> Sequence.iter f

  let length t =  to_seq t |> Sequence.length

  let rec height = function
    | `Node (_, []) -> 0
    | `Node (a, h::t) -> max (height h + 1) (height (`Node (a, t)))

  let rec map f ((`Node (x, cl)):'b t) : 'a t =
    let fx = f x in
    `Node (fx, List.map (map f) cl)

  let rec mapd f d (`Node (x, cl)) =
    let fdx = f d x in
    `Node (fdx, List.map (mapd f (d+1)) cl)

  (* untested... *)
  let mapi f i tr =
    let index = ref i in
    let fidx x =
      let fix = f !index x in
      incr index;
      fix in
    map fidx tr

  (* this is not totally clean, since it uses integer lists as the REVERSED
   * position without referring to Zipper *)
  let rec map_rev_pos
      f revpos ((`Node (x, cl)) : 'b t) : 'a t =
    let fpx = f revpos x in
    `Node (fpx, List.mapi (fun i c ->
        map_rev_pos f (i::revpos) c) cl)

  let map_ordered f ((`Node (x, cl)):'b t) : 'a t =
    let rec aux clist =
      (* first do all children *)
      let cflist = List.map
          (fun (`Node (x, _)) -> f x)
          clist in
      let cclist = List.map
          (fun (`Node (_, cl)) -> cl)
          clist in
      (* now descend recursively *)
      List.map2 (fun fx cl -> (`Node (fx, aux cl)))
        cflist cclist in
    let fx = f x in                        (* this is the first evaluated... *)
    `Node (fx, aux cl)

  (* tested: this is exactly the reverse order of Ctree.rev_list *)
  let mapi_ordered f i tr =
    let index = ref i in
    let fidx x =
      let fix = f !index x in
      incr index;
      fix in
    map_ordered fidx tr

  let rec fold_filter_tree
      ~(f:'a -> 'b -> 'b option)
      (init:'b)
      (t:'a tree) : ('b tree) option =
    (* argument order of f kept consistent with `~f in normal fold, if ugly *)
    match t with `Node (x, l) ->
    match f x init with
    | None -> None
    | Some ac ->
      Some (`Node (ac, List.filter_map (fold_filter_tree ~f ac) l))

  let reduce_filter_tree
      ~(f:'a -> 'a -> 'a option)
      (t:'a tree) : ('a tree) =
    match t with `Node (x, l) ->
      (`Node (x, List.filter_map (fold_filter_tree ~f x) l))

  let rec fold_back
      (f: 'a list -> 'b -> 'a)
      (acc)
      (`Node (a, children)) =
    match children with
    | [] -> f acc a
    | _ ->
      let acc' = List.map (fold_back f acc) children in
      f acc' a
      (* some testing: seems to do the trick of folding backwards up the tree.
       * to get at the local neighbors, it could be useful to define this
       * also for the zipper itself.*)

  let rec fold_mother_children
      ~(f:'a -> 'a list -> 'b -> 'b)
      (init: 'b)
      (`Node (m, dl): 'a tree) : 'b =
    let acc = match dl with
      | [] -> init                                       (* skip leaf nodes! *)
      | l -> f m (List.map xval l) init in
    (List.fold_left (fold_mother_children ~f) acc dl)


    (*| `Node (_, ([] | [_])) -> None*)

  (* sexp serialization *)

  (* utility copied from roseTree *)
  let split_at_length_minus_1 l =
    let rev_list = List.rev l in
    match rev_list with
    | []          -> (l, None)
    | [item]      -> ([], Some item)
    | item::items -> (List.rev items, Some item)

  (* debug *)
  let balance = ref 0

  let to_sexp ?(parens=("(",")")) ?(sep=" ") pp_val tree =
    ignore @@ Format.flush_str_formatter ();
    let fmt = Format.str_formatter in
    let sep () = Format.pp_print_string fmt sep in
    let popen () = incr balance;
      Format.(pp_print_string fmt (fst parens)) in
    let pclose () = decr balance;
      Format.(pp_print_string fmt (snd parens)) in
    let rec print_children children =
      match split_at_length_minus_1 children with
      | (_, None) -> ()
      | (non_last, Some last) ->
        sep ();
        popen ();
        print_non_last_children non_last;
        print_last_child last;
        pclose ();
    and print_child_grandchildren (`Node (ch, gch)) = begin match gch with
      | [] -> pp_val fmt ch;
      | _ ->
        popen ();
        pp_val fmt ch;
        print_children gch;
        pclose () end;
    and print_non_last_children non_last_children =
      List.iter
        (fun n -> print_child_grandchildren n; sep ();)
        non_last_children;
    and print_last_child n = print_child_grandchildren n;
    in
    let print_root (`Node (root_value, root_children)) =
      popen ();
      pp_val fmt root_value;
      print_children root_children;
      pclose ();
    in
    print_root tree;
    Format.flush_str_formatter ()

end


module Z = struct

  include RoseTree.Zipper

  let height zpr = T.height (tree zpr)

  let length zpr = T.length (tree zpr)

  let root_path zpr = parent_fold (fun l zp -> T.xval (tree zp) :: l) [] zpr

  (* how deep down is the zipper position *)
  let depth zpr = parent_fold (fun i _zp -> succ i) 0 zpr

  let zval zpr = T.xval (tree zpr)

  (* this maps the entire tree from the root, then zips to the same position.
   * to map only from the position down, need to fork first *)
  let map f zpr =
    let pos = pos zpr in
    root zpr
    |> tree |> T.map f |> zipper
    |> zip_to (pos:>pos)                          (* this should always work *)
    |> Option.get_exn

  let left_sib_count zp =
    zp |> Gen.unfold
      (fun zp' -> match left_sibling zp' with
         | None -> None
         | Some lzip -> Some ((), lzip))
    |> Gen.length

  let children zp =
    0 |> Gen.unfold
      (fun i -> match nth_child i zp with
         | None -> None
         | Some c -> Some (c, i+1))

  (* relatives as options or sequences. do it in such a way that when iterating
     over a full tree, pairs of relatives are not duplicated *)

  let grandparent zp =
    zp |> parent >>= parent

  let greatgrandparent zp =
    zp |> parent >>= parent >>= parent

  let a_sibling zp =
    match zp |> right_sibling with
    | Some _ as a -> a
    | None -> zp |> left_sibling

  let aunt zp =
    (* both sides since from different generations *)
    zp |> parent >>= a_sibling

  let greataunt zp =
    (* both sides *)
    zp |> grandparent >>= a_sibling

  let left_cousins zp =
    (* one side since same gen *)
    match zp |> parent >>= left_sibling with
    | None -> Gen.empty
    | Some la -> children la

  let left_cousins2 zp =
    (* one side since same generation *)
    match zp |> grandparent >>= left_sibling with
    | None -> Gen.empty
    | Some ga -> ga |> children |> Gen.flat_map children

  let cousins1r zp =
    (* both sides *)
    match zp |> greataunt with
     | None -> Gen.empty
     | Some ga -> ga |> children

  let to_seq zip yield =
    (* pre-order depth-first iterator over the subtree rooted at the zipper
     * position. copied from the one for tree.
     * no hygiene: the root of the yielded elements is still the global root!
     * *)
    let rec iter zp =
      yield zp;
      Gen.iter iter (children zp) in
    iter zip

  let leavez zip =
    (* iterator over the leaves of the subtree rooted at the zipper position.
     * no hygiene. *)
    zip |> to_seq |> Sequence.filter
      (fun zp -> match tree zp with
         | (`Node (_, [])) -> true
         | _ -> false)

  let take_while f zip =
    (* iterate over the subtree of the initial position with to_seq, as long as
     * f gives true on the running zipper pos.
     * if not, delete (prune) the current zipper.
     * restart iteration from the initial position.
     * stop when an iteration runs through; return the last full zipper.
     * no hygiene. in unvisited regions, f may still give false. *)
    let pos = pos zip in
    let rec scan z =
      let s = to_seq z in
      let zf = s |> Sequence.find (fun zz ->
          if f zz then None else Some zz) in
      match zf with
      | None -> Some z
      | Some zz ->
        delete zz               (* if f failed at the root, this gives None! *)
        >>= zip_to (pos:>pos)
        >>= scan in
    scan zip

  (* sequences of relative pairs, which should contain every possible related
     pair in a tree exactly once.
     Hygiene: we exclude pairs whose common ancestor is
     above the input zipper position! to get all pairs, go to the root first*)
  (* no check for nans; these are purely structural! *)
  let make_seq_of_opt rel_opt zip =
    zip
    |> fork                                           (* isolate the subtree *)
    |> to_seq (* this iterates only over the subtree but would allow relations
                 out of the subtree to be found *)
    |> Sequence.filter_map (fun s -> U.pair s (rel_opt s))

  let mm_pairs zip = make_seq_of_opt (fun self -> Some self) zip

  let ss_pairs zip = make_seq_of_opt left_sibling zip

  let md_pairs zip = make_seq_of_opt parent zip

  let an_pairs zip = make_seq_of_opt aunt zip

  let gmgd_pairs zip = make_seq_of_opt grandparent zip

  let ggmggd_pairs zip = make_seq_of_opt greatgrandparent zip

  let gagn_pairs zip = make_seq_of_opt greataunt zip

  let make_seq_of_rel rel_gen zip =
    zip
    |> fork                                  (* isolate the subtree as above *)
    |> to_seq
    |> Sequence.flat_map (fun zp ->
        let cs = rel_gen zp |> Sequence.of_gen in
        Sequence.map (fun c -> zp, c) cs)

  let cc_pairs zip = make_seq_of_rel left_cousins zip

  let cc1r_pairs zip = make_seq_of_rel cousins1r zip

  let cc2_pairs zip = make_seq_of_rel left_cousins2 zip

  (* BEWARE: the root must be kept for all but the direct ancestor pairs! *)
  let drop_root pairs zip =
    Sequence.of_gen (children zip)
      |> Sequence.map pairs
      |> Sequence.concat

  let md_pairs' z = drop_root md_pairs z
  let gmgd_pairs' z = drop_root gmgd_pairs z
  let ggmggd_pairs' z = drop_root ggmggd_pairs z

  (* handling of correlations on sequences of Z.t pairs.
     no check for nans!  *)
  let pair_xcorr cor_fn vec_fun1 vec_fun2
      (pair_sq: ('a t * 'a t) Sequence.t) =
    let fpairs = pair_sq |> Sequence.map (fun (z1, z2) ->
        (vec_fun1 (zval z1), vec_fun2 (zval z2))) in
    let ar1, ar2 =
      let open! Sequence in
      (map fst fpairs |> to_array , map snd fpairs |> to_array) in
    cor_fn ar1 ar2

  let pair_corr cor_fn vec_fun = pair_xcorr cor_fn vec_fun vec_fun

  (* printing *)
  let print print_tree fmt zp = Format.fprintf fmt "%a" print_tree (tree zp);;

end
