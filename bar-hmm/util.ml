(* general utility functions that do not depend on the tree libraries, or on
 * any specifics of cell cycle data. mostly block matrix handling and linear
 * algebra (module U) and statistics (module Stats) *)

open! Containers
open Lacaml.D
open Option.Infix
open Bigarray

exception Matrix_problem

(* not needed anymore: containers 2.0 now has float comparators *)
(*module Float_infix = struct*)
  (*open! Float*)
  (*let (>=.) = (>=)*)
  (*let (<=.) = (<=)*)
  (*let (<.) = (<)*)
  (*let (>.) = (>)*)
  (*let (<>.) = (<>)*)
  (*let (=.) = (=)*)
(*end*)
(*include Float_infix*)

module U = struct
  [@@@landmark "auto"]


  (* utilities and timer *)

  let last_timed = ref nan
  let time ?(silent=false) f x =
    let start = Unix.gettimeofday () in
    let res = f x in
    let stop = Unix.gettimeofday () in
    let () = last_timed := stop -. start in
    let () = if not silent then
        Printf.printf "Execution time: %fs\n%!" !last_timed in
    res

  let pair zp zp' = zp' >|= (fun z -> z, zp)
  let pairs zp gen = Gen.map (fun z -> z, zp) gen

  let vec_of_gslvec gv =
    let l = Gsl.Vector.length gv in
    let res = Vec.create l in
    for i = 1 to l do
      res.{i} <- gv.{i-1}
    done; res

  let gslvec_of_vec v =
    let l = Vec.dim v in
    let res = Gsl.Vector.create l in
    for i = 1 to l do
      res.{i-1} <- v.{i}
    done; res

  let filter_some l = List.filter_map (fun x -> x) l

  let arange ~step ?tol start stop =
    assert (step >. 0.);
    if
      Float.(is_nan start || is_nan stop) then [|nan|]
    else
      let tol = match tol with
        | None -> step /. 2.
        | Some t -> t in
      start
      |> Gen.unfold (fun v ->
          if v >. stop +. tol then None else Some (v, v +. step))
      |> Gen.to_array

  (* generic adjoint representation for 1 - 3 dim symmetric matrices *)

  let get aa (i,j) = Array2.get aa i j

  let addim i = i * (i + 1) / 2

  (* get the coordinates: need only the upper triangle *)
  let coords = function
    | 1 -> [1, (1,1)]
    | 2 -> [1, (1,1); 2, (2,2);
            3, (1,2)]
    | 3 -> [1, (1,1); 2, (2,2); 3, (3,3);
            4, (1,2); 5, (2,3);
            6, (1,3) ]
    | _ -> raise (Invalid_argument "adjoint only in dims 1 .. 3 ")
  let coord d aa i = get aa (List.assoc ~eq:(=) i (coords d))
  (* vector for symmetric matrices. looks only in the upper triangle *)
  let cvec d c = Vec.init (addim d) (coord d c)

  (* make the matrix from the coords *)
  let cmat d vec =
    let out = Mat.make0 d d in
    List.iter
      (fun (i, (j, k)) ->
         out.{j,k} <- vec.{i};
         if j <> k then out.{k,j} <- vec.{i};)
      (coords d); out

  (* adjoint matrix associated with a matrix a; operates on vector
   * representation of symmetric matrices. *)
  let ad d a =
    assert (Mat.dim1 a = d);
    let ad_dim = addim d in
    let img j i =
      let um =
        let v = Vec.make0 ad_dim in
        v.{i} <- 1.;
        cmat d v in
      let ad_image = gemm a @@ gemm um ~transb:`T a in
      coord d ad_image j in
    Mat.init_cols ad_dim ad_dim img


  (* markov tree functions *)

  (* attention: this is the limit for a chain, not a tree at equal times. *)
  exception Diverging_recurrence

  let markov_limit_cov a bbt =
    let d = Mat.dim2 a in
    let a = lacpy ~m:d ~n:d a in
    let bbt = lacpy ~m:d ~n:d bbt in
    (* checking that the abs vals of evs of a are all <1!. for cross-coupling
     * the evs can be complex! *)
    let () =
      let[@landmark] all_evs_less_1 =
        let copy_a = lacpy a in                             (* FAKKING BAGG! *)
        let _, evs_real, evs_im, _ = geev copy_a in
        Vec.(add (sqr evs_real) (sqr evs_im)
             |> fold (fun x y -> x && (y <. 1.)) true) in
      if not all_evs_less_1 then raise Diverging_recurrence in
    (* (1-Ad_a)^{-1} c = c_00; where c = vector of bbt *)
    let c = Mat.of_col_vecs_list [cvec d bbt] in
    gesv Mat.(sub (identity (addim d)) (ad d a)) c;
    Mat.col c 1 |> cmat d

  let get_covar_md a bbt =
    (* mother left = rows, daughter right = cols *)
    let c00 = markov_limit_cov a bbt in               (* this is dim-reduced *)
    let a' = lacpy ~m:(Mat.dim1 c00) a in                 (* cut height of A *)
    gemm c00 ~transb:`T a'


  (* for sampling Gaussians *)

  (*get correlation matrix, normalizing with the stationary variances*)
  let c_to_rho cmm =
    let dv = Mat.copy_diag cmm
             |> Vec.sqrt
             |> Vec.reci
             |> Mat.of_diag in
    function m -> gemm dv (gemm m dv)

  (* cumsum of int array with 0 in front and nothing in the end *)
  let isum ar = Array.fold (+) 0 ar

  (* cumsum of an array with a getter function *)
  let icumf f ar =
    Array.fold
      (fun cum el -> match cum with
         | [] -> [f el]
         | h::_ -> (f el + h)::cum)
      [] ar

  (* leave the last cumulative in - does not hurt *)
  let icum ar =
    let l = icumf (fun (x:int) -> x) ar in
    List.(0 :: (rev l))
    |> Array.of_list

  (* eigenvalues and orthogonal matrix from symmetric matrix *)
  let evs_rot mm =
    let work_cpy = lacpy mm in
    let evs = syev ~vectors:true work_cpy in
    evs, work_cpy

  let centered_gauss rng sigma =
    if sigma =. 0. then 0. else Gsl.Randist.gaussian rng ~sigma

  (* multinormal random variates *)
  let draw rng mean sigmas rot =
    let dgausses = Vec.map (centered_gauss rng) sigmas in
    (* gemv overwrites y! however by accident this did not matter.. *)
    let y, beta = copy mean, 1. in
    let cgausses = gemv ~y ~beta rot dgausses in
    cgausses


  (* matrix and vector manipulation *)


  let is_empty mat = Equal.poly Mat.empty mat

  (* split in 2 *)
  let split2 vec =
    let l = Vec.to_list vec in
    let l1, l2 = List.take_drop (List.length l / 2) l in
    Vec.(of_list l1, of_list l2)

  let block_matrix blocks_array =
    (* blocks are Lacaml.D.mat. 1-based! *)
    let open Array in
    let brows = length blocks_array in
    let bcols = length blocks_array.(0) in
    if not (for_all (fun r -> length r = bcols) blocks_array)
    then raise @@ Invalid_argument "blocks must form a matrix";
    (* now to the inner matrices of Lacaml.D.mat type *)
    let rcums, rdim =
      let l = icumf (fun el -> Mat.dim1 el.(0)) blocks_array in
      Array.of_list (0::List.rev l), List.hd l in
    let ccums, cdim =
      let l = icumf (fun el -> Mat.dim2 el) blocks_array.(0) in
      Array.of_list (0::List.rev l), List.hd l in
    let b = Mat.make0 rdim cdim in
    for r = 0 to brows - 1 do
      for c = 0 to bcols - 1 do
        (* much better: blit here! *)
        ignore @@ lacpy
          ~b ~br:(rcums.(r)+1) ~bc:(ccums.(c)+1)
          blocks_array.(r).(c);
      done; done;
    b

  let norm v = sqrt Vec.(sqr_nrm2 v)

  let zeros_like m =
    Mat.(make0 (dim1 m) (dim2 m))

  let stack_diag m1 m2 =
    let rdim1, cdim1, rdim2, cdim2 =
      Mat.(dim1 m1, dim2 m1, dim1 m2, dim2 m2) in
    let b = Mat.make0 (rdim1 + rdim2) (cdim1 + cdim2) in
    ignore @@ lacpy ~b ~br:1 ~bc:1 m1;
    ignore @@ lacpy ~b ~br:(rdim1 + 1) ~bc:(cdim1 + 1) m2;
    b

  let stack_high m1 m2 =
    block_matrix [| [|m1|]; [|m2|] |]

  (* handle empty matrices *)
  let matmul ?(transa=`N) ?(transb=`N) ?c m1 m2 =
    let[@landmark] emptytest = is_empty m1 && is_empty m2 in
    if emptytest
    then Mat.empty
    else
      let[@landmark] gemres = gemm ~transa ~transb ?c m1 m2 in
      gemres

  (* log det of pos def matrices by cholesky *)
  let log_det m =
    (*if Mat.dim1 m = 0 then*)
    (*0.*)
    (*else*)
    let m' = lacpy m in
    potrf m';
    Mat.copy_diag m'
    |> Vec.log
    |> Vec.sum
    |> ( *.) 2.                      (* since the decomposition is M = U^T.U *)

  let ge_inv m =
    let mi = lacpy m in
    getri mi; mi

  let po_inv m =
    (* for symmetric pos def matrices - only looks at the lower triangle.
     * also, only writes the lower triangle: result is not symmetric. have to
     * use sy-aware functions afterwards ! if not, wrong results! *)
    let up = false in
    let mi = lacpy m in
    potrf ~up mi; potri ~up mi; mi

  (* make an orthogonal basis so that the first vector is the normalized given
   * vector, i.e. v/|v|. basis vectors are in the columns of the output *)
  let fill_basis v =
    let l = Vec.dim v in
    let o = Mat.create l l in
    let _ = v |> Mat.from_col_vec |> lacpy ~b:o in
    let tau = geqrf ~n:1 o in
    let () = orgqr ~tau ~m:l ~n:l ~k:1 o in
    let sign = norm v *. o.{1,1} /. v.{1} in
    Mat.scal ~n:1 sign o;
    o

  (* overall this has helped basically not at all; not a good tradeoff for the
   * bug-danger *)
  (** destroy the inner matrix m' only for equal dimensions, otherwise fall
   * back to the safe version of cong. *)
  let cong_destructive ?(trans=`N) m m' =
    let innerdim = Mat.dim1 m' in
    let transa, transb, outerdim = match trans with
      | `N -> `N, `T, Mat.dim1 m
      | `T -> `T, `N, Mat.dim2 m in
    if innerdim = outerdim
    then
      let m'm = Mat.create innerdim innerdim in
      matmul ~c:m'm m' ~transb m |> ignore;
      let mm'm = m' in
      matmul ~c:mm'm ~transa m m'm |> ignore;
      mm'm
    else
      matmul ~transa m @@ matmul m' ~transb m

  let cong ?(trans=`N) m m' =
    let transa, transb = match trans with
      | `N -> `N, `T
      | `T -> `T, `N in
    matmul ~transa m @@ matmul m' ~transb m

  let vmv v m v' = dot v @@ gemv m v'

  (* functions for completing the square; useful for the integration of
   * overhanging cells. *)
  let ainv_b a b =
    (* b is premultiplied in-place.
     * blas: a is also modified! this is not in the lacaml docs *)
    let a' = lacpy a in
    gesv a' b; ()

  (* checked: this order of the permutation indices seems to correspond with
   * Likely.permutation_offs_taus *)
  let permute_matrix indices a =
    let nx = Array.length indices in
    let xi = Lacaml.Common.create_int32_vec nx in
    Array.iteri (fun i el -> xi.{i+1} <- Int32.of_int el) indices;
    let () = assert (Mat.dim1 a = nx) in
    let a = lacpy a in lapmt a xi;         (* necessary to protect the input *)
    let a = Mat.transpose_copy a in lapmt a xi;           (* other dimension *)
    a

  let permute_vector indices v =                     (* indices are 1-based! *)
    let v' = copy v in
    Array.iteri (fun i ind -> v'.{i+1} <- v.{ind}) indices;
    v'

  let permute_back_vector indices v =                         (* the inverse *)
    let v' = copy v in
    Array.iteri (fun i ind -> v'.{ind} <- v.{i+1}) indices;
    v'

  let block_minors a indices offset =
    (* for a symmetrix matrix m, first permute rows and columns according to
     * indices and offset, to get a block matrix a = [axx, axy; ayx, ayy].
     * finally return (axx, byy, mxy) where
     * byy = ayy - ayx axx^-1 axy and mxy = - axx^{-1} axy.
     * _verified_ in an example! *)
    let nx = offset in
    let a = permute_matrix indices a in
    let axx = lacpy ~m:nx ~n:nx a in         (* need to copy to preserve axx *)
    let mxy =
      let axx', m = lacpy axx, lacpy ~m:nx ~ac:(nx+1) a in
      gesv axx' m;        (* axx^-1 axy. writes stuff into axx', result in m *)
      Mat.scal ~-.1. m;                                         (* swap sign *)
      m in
    let byy =                                        (* ayy - ayx axx^-1 axy *)
      (*if Mat.dim2 mxy = 0 then             (* the case of no internal cells! *)*)
        (*Mat.create 0 0*)
      (*else*)
        gemm ~ar:(nx+1) ~k:nx a mxy
        |> Mat.add ~ar:(nx+1) ~ac:(nx+1) a in
    (axx, byy, mxy)


  (* interpolation *)

  (** linear interpolation between given x values with clipping to the end
   * values of the y array *)
  let interpolate_linear xval xar yar =
    let open Array in
    let l = length xar in
    assert (l = length yar);
    match xval with
    | a when Float.is_nan a -> nan
    | _ -> let i1, i2 = match Array.find_idx ((<.) xval) xar with
        (* all x are smaller than xval *)
        | None -> l - 1, l - 1
        (* all are greater *)
        | Some (0, _) -> 0, 0
        (* within range *)
        | Some (i, _) -> i - 1, i in
      let inter x1 x x2 =
        if x1 =. x2
        then (1., 0.)
        else (x2 -. x) /. (x2 -. x1), (x -. x1) /. (x2 -. x1) in
      let c1, c2 = inter xar.(i1) xval xar.(i2) in
      yar.(i1) *. c1 +. yar.(i2) *. c2


  (* printing *)
  let pf fmt f = Format.fprintf fmt "%.4g" f;;
  let pr_gslv fmt v =
    let a = Gsl.Vector.to_array v in
    Array.pp ~sep:"; " Float.pp fmt a;;
  let pu fmt () = Format.fprintf fmt "*"


  (* file system *)

  let provide_dir dirname =
    let open Unix in
    try
      mkdir dirname 0o755;
      Format.printf "made directory %s@." dirname
    with (Unix_error (EEXIST, _, _)) -> ()

end


module Stat = struct

  let mean ar = Array.(
      fold (+.) 0. ar /. float_of_int (length ar))

  let var ar =
    let m = mean ar in
    let m2 = Array.fold (fun acc a -> acc +. a *. a) 0. ar in
    m2 /. float_of_int (Array.length ar) -. m *. m

  let covar_array ?m1 ?m2 as1 as2 =
    let m1 = Option.get_or ~default:(mean as1) m1 in
    let m2 = Option.get_or ~default:(mean as2) m2 in
    let sum = Array.fold2
        (fun acc e1 e2 -> acc +. (e1-.m1) *. (e2-.m2)) 0. as1 as2 in
    sum /. float_of_int (Array.length as1)

  let cor_array as1 as2 =
    let s1, s2 = sqrt (var as1), sqrt (var as2) in
    let cov = covar_array as1 as2 in
    cov /. s1 /. s2

  let cor_fn stat_on_array n ar1 ar2 =
    (* ar1 shifted forward by 0,...,n-1 *)
    let open Array in
    let l = length ar1 in
    let () = assert (l = length ar2) in
    let cvi i =
      let as1 = sub ar1 i (l-i) in
      let as2 = sub ar2 0 (l-i) in
      stat_on_array as1 as2 in
    init n cvi

  let covar = cor_fn (covar_array ?m1:None ?m2:None)

  let cor = cor_fn cor_array

  let log_cor n ar1 ar2 =
    cor n (Array.map log ar1) (Array.map log ar2)

  (*let ranks ar =*)
    (*(* manual implementation - in the meantime CCArray also got it... *)*)
    (*let open Array in*)
    (*let cf2 a b = Float.compare (snd a) (snd b) in*)
    (*let ci2 a b = Int.compare (snd a) (snd b) in*)
    (*let ar' = mapi (fun i f -> (i, f)) ar in*)
    (*sort cf2 ar';*)
    (*let ar'' = mapi (fun i (j, _f) -> i, j) ar' in*)
    (*sort ci2 ar'';*)
    (*map fst ar''*)

  (* speed: verify first, then compare without nan check! *)

  (*let compare_non_nan f f' =*)
    (*if Float.(is_nan f || is_nan f')*)
    (*then raise (Invalid_argument "comparing with nan");*)
    (*Float.compare f f'*)

  let check_nan ar =
    ar |> Array.iter (fun f ->
        if Float.is_nan f
        then raise (Invalid_argument "encountered nan"))

  let check_nan_v v =
    v |> Vec.iter (fun f ->
        if Float.is_nan f
        then raise @@ Invalid_argument
            (Format.sprintf "encountered nan in %a" Lacaml.Io.pp_fvec v))

  let ranks ar =
    check_nan ar; ar
    |> Array.sort_ranking Float.compare
    |> Array.map float_of_int

  let tied_ranks ar =
    let open Sequence in
    check_nan ar;
    let len = Array.length ar in
    let arr = Array.sort_ranking Float.compare ar in
    (*let ars = Array.sorted Float.compare ar in*)
    (*let ars = Array.map (fun i -> ar.(i)) arr in*)
    (*let eq i j = Float.equal ars.(i) ars.(j) in*)
    let eq i j = Float.equal ar.(arr.(i)) ar.(arr.(j)) in
    let av_seq l =
      let rec aux n s = function
        | h::t -> aux (succ n) (s + h) t
        | [] -> float_of_int s /. float_of_int n, n in
      let av, len = aux 0 0 l in
      repeat av |> take len in
    let ordered_tied =
      0 -- pred len
      |> group_succ_by ~eq
      |> map av_seq
      |> concat |> to_array in
    (*Array.init len (fun i -> ordered_tied.(arr.(i)))*)
    Array.map (Array.get ordered_tied) arr

  let rank_cor ?(tied=true) =
    let ranking = if tied then tied_ranks else ranks in
    cor_fn (fun ar1 ar2 -> cor_array (ranking ar1) (ranking ar2))

  (** get mean and variance of the logs of a sequence of (time) values *)
  let log_mean_vars seq =
    let open Sequence in
    seq
    |> filter_map (fun t ->
        if Float.is_nan t then None else Some (log t))
    |> to_array
    |> fun arr -> (mean arr, var arr)

  let ecdf_vals ?(margins=`Symmetric) ar =
    let n = Array.length ar in
    let fn = float_of_int n in
    let sar = Array.sort_ranking Float.compare ar in
    let cvals =
      Int.range 1 n
      |> Sequence.map (fun i -> float_of_int i /. fn)
      |> begin match margins with
        | `Symmetric -> Sequence.map (fun f -> f -. 1./.2./.fn)
        | `Start0 -> Sequence.map (fun f -> f -. 1./.fn)
        | `End1 -> Fun.id end
      |> Sequence.to_array in
    Array.map (fun i -> cvals.(i)) sar

  let normal_transformed_vals ?(margins=`Symmetric) ar =
    let icdf p = Gsl.Cdf.gaussian_Pinv ~sigma:1. ~p in
    ar |> ecdf_vals ~margins |> Array.map icdf

  (* bootstrapping *)

  let resample_once seed list_ =
    let () = Random.init seed in
    let open Sequence in
    list_ |> Array.of_list
    |> random_array
    |> take (List.length list_)
    |> to_list

  (* here, stat_fun takes a sequence *)
  (* put the original up front for testing etc. *)
  let array_bs n stat_fun ar =
    let l = Array.length ar in
    let samestat = stat_fun (Sequence.of_array ar) in
    let rseq = Sequence.(ar |> random_array |> take l) in
    samestat :: List.init (n-1) (fun _ -> stat_fun rseq)
    |> Array.of_list

  (* over a sequence *)
  let seq_bs n stat_fun taupairseq =
    let ar = taupairseq |> Sequence.to_array in
    array_bs n stat_fun ar

  (* block-wise bootstrap over an array of sequences. This recomputes the
   * stat_fun each time. this looks wasteful but in fact, we need to calculate
   * the statistics over trees always over all trees, never tree-wise. This is
   * because tree-wise statistics do not decorrelate. *)
  let seqarray_bs n stat_fun seqar =
    let l = Array.length seqar in
    let rseq = Sequence.(seqar |> random_array |> take l |> concat) in
    Array.init n (fun _ -> stat_fun rseq)        (* this is lazy evaluation! *)

  let quantiles qlist ar =
    let sar = Array.copy ar in
    Array.sort Float.compare sar;
    let poslist = qlist
                  |> List.map (( *.) (Float.of_int (Array.length ar)))
                  |> List.map Float.to_int in
    List.map (fun p -> sar.(p)) poslist;;


  (* mean without nans *)
  let mean_nan v =
    let l = Vec.dim v in
    let r = ref 0 in
    let nansum = Vec.fold (fun ac f ->
        if Float.is_nan f
        then (incr r; ac)
        else ac +. f) 0. v in
    let norm = float_of_int (l - !r) in
    nansum /. norm

  (* cov without nans *)
  let cov_nan v v' =
    let l = Vec.dim v in
    let cv = Vec.make0 l in
    for i = 1 to l do cv.{i} <- v.{i} *. v'.{i} done;
    mean_nan cv

  (* mean of columns from a sample matrix in which the row
   * index is the sample number and column index is the feature/vector
   * component. drop nans! *)
  let mean_samples mat =
    Vec.init (Mat.dim2 mat) (fun i -> mean_nan (Mat.col mat i))

  (* general variance-covariance matrix of a sample matrix drop nans! *)
  let var_covar_samples mat =
    let d = Mat.dim2 mat in
    let c i j =
      let f k =
        let c = Mat.col mat k in
        let m = mean_nan c in
        Vec.add_const ~-.m c in
      cov_nan (f i) (f j) in
    let res = Mat.create d d in
    for i = 1 to d do
      for j = i to d do
        res.{i,j} <- c i j done;
      for j = 1 to i - 1 do
        res.{i,j} <- res.{j,i} done;
    done; res

end
