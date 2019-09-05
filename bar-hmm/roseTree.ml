(*
copyright (c) 2013-2014, Simon Cruanes, Emmanuel Surleau
all rights reserved.

redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.  redistributions in binary
form must reproduce the above copyright notice, this list of conditions and the
following disclaimer in the documentation and/or other materials provided with
the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*)


type +'a t = [`Node of 'a * 'a t list]

type 'a tree = 'a t

type 'a sequence = ('a -> unit) -> unit
type 'a printer = Format.formatter -> 'a -> unit

let create a = `Node (a, [])

let rec fold ~f init_acc (`Node (value, children)) =
  let acc = f value init_acc in
  List.fold_left (fun acc' child_node -> fold ~f acc' child_node) acc children

let to_seq t yield =
  let rec iter (`Node (value, children)) =
    let () = yield value in
    List.iter iter children in
  iter t


let split_at_length_minus_1 l =
  let rev_list = List.rev l in
  match rev_list with
  | []          -> (l, None)
  | [item]      -> ([], Some item)
  | item::items -> (List.rev items, Some item)

let print pp_val formatter tree =
  let rec print_children children indent_string =
    let non_last_children, maybe_last_child =
      split_at_length_minus_1 children
    in
    print_non_last_children non_last_children indent_string;
    match maybe_last_child with
    | Some last_child -> print_last_child last_child indent_string;
    | None            -> ();
  and print_non_last_children non_last_children indent_string =
    List.iter (fun (`Node (child_value, grandchildren)) ->
        Format.pp_print_string formatter indent_string;
        Format.pp_print_string formatter "├ ";
        pp_val formatter child_value;
        Format.pp_force_newline formatter ();
        let indent_string' = indent_string ^ "│ " in
        print_children grandchildren indent_string'
      ) non_last_children;
  and print_last_child (`Node (last_child_value, last_grandchildren)) indent_string =
    Format.pp_print_string formatter indent_string;
    Format.pp_print_string formatter "└ ";
    pp_val formatter last_child_value;
    Format.pp_force_newline formatter ();
    let indent_string' = indent_string ^ "  " in
    print_children last_grandchildren indent_string'
  in
  let print_root (`Node (root_value, root_children)) =
    pp_val formatter root_value;
    Format.pp_force_newline formatter ();
    print_children root_children ""
  in
  print_root tree;
  Format.pp_print_flush formatter ()

module Zipper = struct

  type 'a parent = {
    left_siblings: ('a tree) list ;
    value: 'a ;
    right_siblings: ('a tree) list ;
  }

  type 'a t = {
    tree: 'a tree ;
    lefts: ('a tree) list ;
    rights: ('a tree) list ;
    parents: ('a parent) list ;
  }

  type abs_pos = [`Abs of int list]
  type rel_pos = [`Rel of int list]
  type pos = [abs_pos | rel_pos]

  let zipper tree = { tree = tree ; lefts = []; rights = []; parents = [] }

  let tree zpr = zpr.tree

  let fork zpr = zipper @@ tree zpr

  let left_sibling zpr =
    match zpr.lefts with
    | []    -> None
    | next_left::farther_lefts ->
      Some { zpr with
        tree = next_left ;
        lefts = farther_lefts;
        rights = zpr.tree::zpr.rights ;
      }

  let right_sibling zpr =
    match zpr.rights with
    | []                  -> None
    | next_right::farther_rights ->
      Some { zpr with
        tree = next_right ;
        lefts = zpr.tree::zpr.lefts ;
        rights = farther_rights ;
      }

  let parent zpr =
    match zpr.parents with
    | []                    -> None
    | { left_siblings ; value ; right_siblings }::other_parents ->
      Some {
        tree = `Node (value, List.rev_append zpr.lefts [zpr.tree] @ zpr.rights) ;
        lefts = left_siblings ;
        rights = right_siblings ;
        parents = other_parents ;
      }

  let rec root zpr =
    let maybe_parent_zpr = parent zpr in
    match maybe_parent_zpr with
    | None                -> zpr
    | Some parent_zpr  -> root parent_zpr

  let parent_seq zpr yield =
    let rec iter zp =
      (* includes the current zipper *)
      let () = yield zp in
      match parent zp with
      | None -> ()
      | Some pzp -> iter pzp in
    iter zpr

  let children_seq zpr yield =
    match tree zpr with
    | `Node (_, []) -> ()
    | `Node (value, left::rest) ->
      let first =
        { tree = left;
          lefts = [];
          rights = rest;
          parents =
            {left_siblings=zpr.lefts; value; right_siblings=zpr.rights}
            ::zpr.parents
        } in
      let rec go_from sibling =
        match right_sibling sibling with
        | None -> ()
        | Some zp ->
          let () = yield zp in go_from zp in
      go_from first

  let nth_child n ({ tree = `Node (value, children) ; _ } as zpr ) =
    let lefts, maybe_child, rev_rights, _counter = List.fold_left (
        fun (lefts, maybe_child, rev_rights, counter) tree ->
          let lefts', maybe_child', rev_rights' =
            match counter with
            | _ when counter == n -> (lefts, Some tree, [])
            | _ when counter < n ->
              (tree::lefts, None, [])
            | _                   ->
              (lefts, maybe_child, tree::rev_rights)
          in
          (lefts', maybe_child', rev_rights', counter+1)
      ) ([], None, [], 0) children
    in
    begin match maybe_child with
      | Some child  ->
        Some {
          tree = child ;
          lefts = lefts;
          rights = List.rev rev_rights ;
          parents = {
            left_siblings = zpr.lefts ;
            value = value ;
            right_siblings = zpr.rights ;
          }::zpr.parents ;
        }
      | None        -> None
    end

  let append_child tree ({ tree = `Node (value, children) ; _ } as zpr ) =
    {
      tree ;
      lefts = children ;
      rights = [] ;
      parents = {
        left_siblings = zpr.lefts ;
        value = value ;
        right_siblings = zpr.rights ;
      }::zpr.parents ;
    }

  let insert_left_sibling tree zpr =
    match zpr.parents with
    | []  -> None
    | _   -> Some { zpr with tree ; rights = zpr.tree::zpr.rights }

  let insert_right_sibling tree zpr =
    match zpr.parents with
    | []  -> None
    | _   -> Some { zpr with tree ; lefts = zpr.tree::zpr.lefts }

  let replace tree zpr =
    { zpr with tree }

  let delete ({ tree = `Node (_value, _children) ; _ } as zpr ) =
    match zpr with
    | { lefts = next_left::farther_lefts ; _  }     ->
      Some { zpr with tree = next_left ; lefts = farther_lefts }
    | { rights = next_right::farther_rights ; _ }  ->
      Some { zpr with tree = next_right ; rights = farther_rights }
    | { parents = { left_siblings ; value ; right_siblings }::other_parents ; _ } ->
      Some {
        tree = `Node (value, List.rev_append zpr.lefts zpr.rights) ;
        lefts = left_siblings ;
        rights = right_siblings ;
        parents = other_parents ;
      }
    | _ -> None

  let at_leaf zpr = match zpr.tree with
    | `Node (_, []) -> true
    | _ -> false

  (* todo: the root is not processed at all! not good! *)
  let parent_fold f start ?(upto) (zpr:'a t) =
    let rec aux (acc, zp) = match parent zp with
      | None -> acc
      | Some _ as pzpo when pzpo = upto -> f acc zp
      | Some pzp -> aux (f acc zp, pzp) in
    aux (start, zpr)

  (* this works now -- should be adapted when parent_fold is fixed *)
  let pos zpr : abs_pos =
    let f ilist zp = List.length zp.lefts :: ilist in
    `Abs (parent_fold f [] zpr)

  (* ditto *)
  let rel_pos_from ancestor zpr : rel_pos =
    let f ilist zp = List.length zp.lefts :: ilist in
    `Rel (parent_fold f [] ~upto:ancestor zpr)

  let rel_path zp1 zp2 : rel_pos =
    if root zp1 <> root zp2 then raise (Invalid_argument "different trees");
    let rec thread (a1, a2) = function
      | ([], [], _) -> List.rev (a2 @ a1)
      | (_h::t, [], _) -> thread (-1::a1, a2) (t, [], true)
      | ([], h::t, _) -> thread (a1, h::a2) ([], t, true)
      | (h1::t1, h2::t2, true) ->
        thread (-1::a1, h2::a2) (t1, t2, true)
      | (h1::t1, h2::t2, false) ->
        if h1 = h2 then thread (a1, a2) (t1, t2, false)
        else thread (-1::a1, h2::a2) (t1, t2, true) in
    let p1, p2 = match (pos zp1, pos zp2) with
      | `Abs l1, `Abs l2 -> l1, l2 in
    `Rel (thread ([], []) (p1, p2, false))

  let pos_up pos : abs_pos option= match pos with
    | `Abs [] -> None
    | `Abs (l:int list) -> Some (`Abs List.(l |> rev |> tl |> rev))

  let normalize_pos pos =
    let rec nlist = function
      | h::h'::t when h >= 0 && h' < 0 -> nlist t
      | h::t -> h :: nlist t
      | [] -> [] in
    match pos with
    | `Rel ll -> Some (`Rel (nlist ll))
    | `Abs ll -> match (nlist ll) with
      | h::_ when h < 0 -> None
      | other -> Some (`Abs other)                  (* should start positive *)

  let combine_pos poslist =
    let contract p1o p2 = match p1o with
      | None -> None
      | Some p1 -> begin match p1, p2 with
          | _, `Abs l2 -> Some (`Abs l2)
          | `Abs l1 , `Rel l2 -> normalize_pos (`Abs (l1 @ l2))
          | `Rel l1 , `Rel l2 -> normalize_pos (`Rel (l1 @ l2)) end in
    match poslist with
    | [] -> None
    | h::t -> List.fold_left contract (Some h) t

  let zip_to =
    let f maybe_zip i = match maybe_zip with
      | None -> None
      | Some zip -> nth_child i zip in
    fun pos zpr -> match pos with
    | `Abs l -> List.fold_left f (Some (root zpr)) l
    | `Rel l -> List.fold_left f (Some zpr) l

end
