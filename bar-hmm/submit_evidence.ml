(* a script to submit jobs through the PBS batch processing system. this is
 * quite specific to the compute cluster used, and is kept here mainly for
 * reference. *)

(* use this directly on the server with ocamlrun ocaml submit_evidence.ml *)
(* job directories will be created beneath the current directory! *)

#use "topfind";;
#require "unix";;
open Printf

(* what do we run where? *)
let program = "model_evidence"

(* run directory *)

let collection_subdir = "fullrun_collect_evidence_size_1e5"

(* new preferred quantile: 98%. *)

let collection_subdir = "esc_collect_evidence_quantile_98_1e5"

let collection_subdir = "fullrun_collect_evidence_quantile_98_1e5"

let collection_subdir = "integrate_evidence_normal_transform_1e5"

let collection_subdir = "esc_collect_evidence_normal_transform_1e5"

let collection_subdir = "on_collect_evidence_normal_transform_1e5"

let collection_subdir = "all_integrate_evidence_normal_transform_1e4"

let collection_subdir = "all_collect_evidence_normal_transform_1e4"

let collection_subdir = "on_integrate_evidence_normal_transform_1e4"

let collection_subdir = "all_collect_evidence_normal_transform_1e6"

let collection_subdir = "all_collect_evidence_normal_transform_1e5"

(* do it for the relative-correlations *)
(*let program = "model_relcorr_evidence"*)
(*let collection_subdir = "collect_relcorr_evidence_1e6"*)

(* parameters *)

(* forking paralellization *)
let parmap_processes = [|1; 5; 7; 13; 18; 26|].(4)
(*  this gets blown up by hyperthreading_factor *)
let hyperthreading_factor = (4, 3)

(* for blas. not clear if this would be actually used. keep at 1 for now *)
let num_threads = 1
(* memory per process.
   for forest_size 15 and tree_height 8, 1 gb is just enough *)
let mem_per_process = 2 (*gb*)

let code_dir = String.concat "/" [Sys.getenv "HOME"; "bar-cycle";]

(* keep the executable safe from overwriting *)
let unique_program_name =
  program ^ "-" ^ string_of_int (Hashtbl.hash (Unix.time ()))

let compiled_program_to_run =
  (* has to be copied to the wd to prevent overwriting by later submission! *)
  sprintf "%s.native" unique_program_name

let cmd_prefix =
  let threads = sprintf "OPENBLAS_NUM_THREADS=%d" num_threads in
  let pps = sprintf "PARMAP_PROCESSES=%d" parmap_processes in
  String.concat " " [threads; pps; "time"]
(* this seeds the rng with the order in the submission, not the job id
 * note that $PBS_JOBID is of the form "nnnnnn.curry0" *)
let args = ["batch";]

(* pbs parameters *)
let jobname = program
(* hyperthreading allowance: 4 HT cores per 3 processes reserved *)
let processors = max 2 ((num_threads * parmap_processes)
                        * fst hyperthreading_factor
                        / snd hyperthreading_factor)
let hours = 72

let shell = ["#!/bin/sh"]

(* pbs script *)
let script compile_id =
  let jobdir = sprintf "%s/$PBS_JOBID" collection_subdir in
  let email = {|user@example.edu|} in
  let pbs_directives = [
    sprintf "-N %s" jobname;
    sprintf "-l walltime=%d:00:00" hours;
    sprintf "-l mem=%dgb" (mem_per_process * parmap_processes);
    (* fat >10 GB per core, normal >2 gb per core *)
    (*sprintf "-l nodes=fat:ppn=%d" processors;*)
    sprintf "-l nodes=1:ppn=%d" processors;
    "-V";                                               (* environment setup *)
    sprintf "-e %s/err" jobdir;
    sprintf "-o %s/out" jobdir;
    sprintf "-m ae";                                        (* abort and end *)
    sprintf "-M %s" email;
    sprintf "-W depend=afterok:%s" compile_id;]
    |> List.map ((^) "#PBS ")
  in
  let dirs = [
    "cd $PBS_O_WORKDIR";
    sprintf "mkdir -p %s" jobdir;
    sprintf "cd %s" jobdir;]
  in
  let wait = ["sleep 3"]
  in
  let reproduce = [
    sprintf "cp %s/%s.ml ." code_dir program;]
  in
  let command =
    [cmd_prefix;
     sprintf "../%s" compiled_program_to_run]
    @ (args)
    |> String.concat " "
    |> (fun x -> [x])
  in
  List.concat [
    shell;
    pbs_directives;
    dirs;
    wait;
    reproduce;
    command; ]


(* setup *)

let compiler_switch = "4.04.2+flambda"
let current_switch =
  let c = Unix.open_process_in "opam switch show" in
  input_line c

let compile_sequence wd =
  [ sprintf "opam switch set %s" compiler_switch;
    sprintf {|eval `opam config env`|};
    sprintf "pushd %s" code_dir;
    "rm -f *.{byte,native}";
    "ocamlbuild -clean";
    sprintf "ocamlbuild -use-ocamlfind %s.native" program;
    sprintf "opam switch set %s" current_switch;
    sprintf {|eval `opam config env`|};
    sprintf "mkdir -p %s/%s" wd collection_subdir;
    sprintf "cp %s.native %s/%s/%s"
      program wd collection_subdir compiled_program_to_run;
    sprintf "git rev-parse HEAD > %s/%s/headrev" wd collection_subdir;]


(* the submitting commands *)

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
  (status_code, id)

let submit_pipe script =
  let open Unix in
  let script = String.concat "\n" script in
  (* "heredoc" crazyness *)
  let chan = open_process_in (sprintf "qsub - <<\"EOF\"\n%s\nEOF\n" script) in
  let id = input_line chan in
  let () = print_endline id in
  let status = close_process_in chan in
  let status_code = match status with
      WEXITED s | WSIGNALED s | WSTOPPED s -> s in
  (status_code, id)

let cleanup () =
  (* not very useful: run at the end of the submitter script! *)
  0

(* old *)
let compile_first () =
  let wd = Sys.getcwd () in
  let () = Sys.chdir code_dir in
  let cmd = (String.concat " && " (compile_sequence wd)) in
  let () = print_endline cmd in
  let success = Sys.command cmd in
  let () = Sys.chdir wd in
  success

let compile_on_node () =
  let wd = Sys.getcwd () in
  let jobdir = sprintf "%s" collection_subdir in
  let script =
    shell @
    ([sprintf "-N compile_evidence";
      sprintf "-v LD_LIBRARY_PATH";
      sprintf "-e %s/compile_err" jobdir;
      sprintf "-o %s/compile_out" jobdir; ]
     |> List.map ((^) "#PBS ")) @
    compile_sequence wd in
  let (comp_submit_status, compile_id) = submit_pipe script in
  comp_submit_status, compile_id

(* compilation on the head node *)
(*;;*)
(*if compile_first () > 0*)
(*then raise (Failure "compilation unsuccessful!")*)


(* run *)
;;
let wd = (Sys.getcwd ())
;;
printf "work dir %s\n" wd;
printf "collection dir %s/%s\n" wd collection_subdir;
;;
let comp_sub_status, compile_id = compile_on_node ()
;;
if comp_sub_status > 0
then failwith "compilation could not be submitted"
;;
printf "reserving %d cores per job\n" processors;
;;
(*printf "%s\n" (String.concat "\n" (script compile_id))*)
(*;;*)
let s = submit_pipe (script compile_id)
;;
let exit_status = (fst s) in
exit exit_status
