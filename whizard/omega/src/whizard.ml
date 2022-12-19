(* whizard.ml --

   Copyright (C) 1999-2022 by

       Wolfgang Kilian <kilian@physik.uni-siegen.de>
       Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
       Juergen Reuter <juergen.reuter@desy.de>
       with contributions from
       Christian Speckner <cnspeckn@googlemail.com>

   WHIZARD is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   WHIZARD is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.  *)

open Printf

module type T =
  sig
    type t
    type amplitude
    val trees : amplitude -> t
    val merge : t -> t
    val write : out_channel -> string -> t -> unit
 
   end

module Make (FM : Fusion.Maker) (P : Momentum.T)
    (PW : Momentum.Whizard with type t = P.t) (M : Model.T) =
  struct
    module F = FM(P)(M)

    type tree = (P.t * F.flavor list) list

    module Poles = Map.Make
        (struct
          type t = int * int 
          let compare (s1, t1) (s2, t2) =
            let c = compare s2 s1 in
            if c <> 0 then
              c
            else
              compare t1 t2
        end)

    let add_tree maps tree trees =
      Poles.add maps
        (try tree :: (Poles.find maps trees) with Not_found -> [tree]) trees

    type t =
        { in1 : F.flavor;
          in2 : F.flavor;
          out : F.flavor list;
          trees : tree list Poles.t }

    type amplitude = F.amplitude

(* \thocwmodulesection{Building Trees} *)

(* A singularity is to be mapped if it is timelike and not the
   overall $s$-channel. *)
    let timelike_map c = P.Scattering.timelike c && not (P.Scattering.s_channel c)

    let count_maps n clist =
      List.fold_left (fun (s, t as cnt) (c, _) ->
        if timelike_map c then
          (succ s, t)
        else if P.Scattering.spacelike c then
          (s, succ t)
        else
          cnt) (0, 0) clist

    let poles_to_whizard n trees poles =
      let tree = List.map (fun wf ->
        (P.Scattering.flip_s_channel_in (F.momentum wf), [F.flavor wf])) poles in
      add_tree (count_maps n tree) tree trees

(* \begin{dubious}
    I must reinstate the [conjugate] eventually!
   \end{dubious} *)

    let trees a =
      match F.externals a with
      | in1 :: in2 :: out ->
          let n = List.length out + 2 in
          { in1 = F.flavor in1;
            in2 = F.flavor in2;
            out = List.map (fun f -> (* [M.conjugate] *) (F.flavor f)) out;
            trees = List.fold_left
              (poles_to_whizard n) Poles.empty (F.poles a) }
      | _ -> invalid_arg "Whizard().trees"

(* \thocwmodulesection{Merging Homomorphic Trees} *)

    module Pole_Map =
      Map.Make (struct type t = P.t list let compare = compare end)
    module Flavor_Set =
      Set.Make (struct type t = F.flavor let compare = compare end)

    let add_flavors flist fset =
      List.fold_right Flavor_Set.add flist fset

    let set_of_flavors flist =
      List.fold_right Flavor_Set.add flist Flavor_Set.empty

    let pack_tree map t =
      let c, f =
        List.split (List.sort (fun (c1, _) (c2, _) ->
          compare (PW.of_momentum c2) (PW.of_momentum c1)) t) in
      let f' = 
        try
          List.map2 add_flavors f (Pole_Map.find c map)
        with
        | Not_found -> List.map set_of_flavors f in
      Pole_Map.add c f' map

    let pack_map trees = List.fold_left pack_tree Pole_Map.empty trees

    let merge_sets clist flist =
      List.map2 (fun c f -> (c, Flavor_Set.elements f)) clist flist

    let unpack_map map =
      Pole_Map.fold (fun c f l -> (merge_sets c f) :: l) map []

(* If a singularity is to be mapped (i.\,e.~if it is timelike and not the
   overall $s$-channel), expand merged particles again: *)
    let unfold1 (c, f) =
      if timelike_map c then
        List.map (fun f' -> (c, [f'])) f
      else
        [(c,f)]

    let unfold_tree tree = Product.list (fun x -> x) (List.map unfold1 tree)

    let unfold trees = ThoList.flatmap unfold_tree trees

    let merge t =
      { t with trees = Poles.map
          (fun t' -> unfold (unpack_map (pack_map t'))) t.trees }

(* \thocwmodulesection{Printing Trees} *)

    let flavors_to_string f =
      String.concat "/" (List.map M.flavor_to_string f)

    let whizard_tree t =
      "tree " ^
      (String.concat " " (List.rev_map (fun (c, _) ->
        (string_of_int (PW.of_momentum c))) t)) ^
      " ! " ^
      (String.concat ", " (List.rev_map (fun (_, f) -> flavors_to_string f) t))

    let whizard_tree_debug t =
      "tree " ^
      (String.concat " " (List.rev_map (fun (c, _) ->
        ("[" ^ (String.concat "+" (List.map string_of_int (P.to_ints c))) ^ "]"))
                            (List.sort (fun (t1,_) (t2,_) ->
                              let c =
                                compare
                                  (List.length (P.to_ints t2))
                                  (List.length (P.to_ints t1)) in
                              if c <> 0 then
                                c
                              else
                                compare t1 t2) t))) ^
      " ! " ^
      (String.concat ", " (List.rev_map (fun (_, f) -> flavors_to_string f) t))

    let format_maps = function 
      | (0, 0) -> "neither mapped timelike nor spacelike poles"
      | (0, 1) -> "no mapped timelike poles, one spacelike pole"
      | (0, n) -> "no mapped timelike poles, " ^
          string_of_int n ^ " spacelike poles"
      | (1, 0) -> "one mapped timelike pole, no spacelike pole"
      | (1, 1) -> "one mapped timelike and spacelike pole each"
      | (1, n) -> "one mapped timelike and " ^
          string_of_int n ^ " spacelike poles"
      | (n, 0) -> string_of_int n ^
          " mapped timelike poles and no spacelike pole"
      | (n, 1) -> string_of_int n ^
          " mapped timelike poles and one spacelike pole"
      | (n, n') -> string_of_int n ^ " mapped timelike and " ^
          string_of_int n' ^ " spacelike poles"

    let format_flavor f =
      match flavors_to_string f with
      | "d" -> "d" | "dbar" -> "D"
      | "u" -> "u" | "ubar" -> "U"
      | "s" -> "s" | "sbar" -> "S"
      | "c" -> "c" | "cbar" -> "C"
      | "b" -> "b" | "bbar" -> "B"
      | "t" -> "t" | "tbar" -> "T"
      | "e-" -> "e1" | "e+" -> "E1"
      | "nue" -> "n1" | "nuebar" -> "N1"
      | "mu-" -> "e2" | "mu+" -> "E2"
      | "numu" -> "n2" | "numubar" -> "N2"
      | "tau-" -> "e3" | "tau+" -> "E3"
      | "nutau" -> "n3" | "nutaubar" -> "N3"
      | "g" -> "G" | "A" -> "A" | "Z" -> "Z"
      | "W+" -> "W+" | "W-" -> "W-"
      | "H" -> "H"
      | s -> s ^ " (not translated)"

    module Mappable =
      Set.Make (struct type t = string let compare = compare end)
    let mappable =
      List.fold_right Mappable.add
        [ "T"; "Z"; "W+"; "W-"; "H" ] Mappable.empty

    let analyze_tree ch t =
      List.iter (fun (c, f) ->
        let f' = format_flavor f
        and c' = PW.of_momentum c in
        if P.Scattering.timelike c then begin
          if P.Scattering.s_channel c then
            fprintf ch "      ! overall s-channel %d %s not mapped\n" c' f'
          else if Mappable.mem f' mappable then
            fprintf ch "      map %d s-channel %s\n" c' f'
          else
            fprintf ch
              "      ! %d s-channel %s can't be mapped by whizard\n"
              c' f'
        end else
          fprintf ch "      ! t-channel %d %s not mapped\n" c' f') t

    let write ch pid t =
      failwith "Whizard.Make().write: incomplete"
(*i
      fprintf ch "! whizard trees by O'Mega\n\n";
      fprintf ch "! %s %s -> %s\n"
        (M.flavor_to_string t.in1) (M.flavor_to_string t.in2) 
        (String.concat " " (List.map M.flavor_to_string t.out));
(*i
      fprintf ch "! %d %d -> %s\n\n"
        (whizard_code1 t.n 1) (whizard_code1 t.n 2)
        (String.concat " " (List.map (fun o ->
          string_of_int (whizard_code1 t.n o)) (ThoList.range 3 t.n)));
i*)
      fprintf ch "process %s\n" pid;
      Poles.iter (fun maps ds ->
        fprintf ch "\n    ! %d times %s:\n"
          (List.length ds) (format_maps maps);
        List.iter (fun d ->
          fprintf ch "\n    grove\n";
          fprintf ch "    %s\n" (whizard_tree d);
          analyze_tree ch d) ds) t.trees;
      fprintf ch "\n"
i*)

  end

(* \thocwmodulesection{Process Dispatcher} *)

let arguments = function
  | [] -> ("", "")
  | args ->
      let arg_list = String.concat ", " (List.map snd args) in
      (arg_list, ", " ^ arg_list)

let import_prefixed ch pid name =
  fprintf ch "    use %s, only: %s_%s => %s !NODEP!\n"
    pid pid name name

let declare_argument ch (arg_type, arg) =
  fprintf ch "    %s, intent(in) :: %s\n" arg_type arg

let call_function ch pid result name args =
  fprintf ch "       case (pr_%s)\n" pid;
  fprintf ch "          %s = %s_%s (%s)\n" result pid name args

let default_function ch result default =
  fprintf ch "       case default\n";
  fprintf ch "          call invalid_process (pid)\n";
  fprintf ch "          %s = %s\n" result default

let call_subroutine ch pid name args =
  fprintf ch "       case (pr_%s)\n" pid;
  fprintf ch "          call %s_%s (%s)\n" pid name args

let default_subroutine ch =
  fprintf ch "       case default\n";
  fprintf ch "          call invalid_process (pid)\n"

let write_interface_subroutine ch wrapper name args processes =
  let arg_list, arg_list' = arguments args in
  fprintf ch "  subroutine %s (pid%s)\n" wrapper arg_list';
  List.iter (fun p -> import_prefixed ch p name) processes;
  List.iter (declare_argument ch) (("character(len=*)", "pid") :: args);
  fprintf ch "    select case (pid)\n";
  List.iter (fun p -> call_subroutine ch p name arg_list) processes;
  default_subroutine ch;
  fprintf ch "    end select\n";
  fprintf ch "  end subroutine %s\n" wrapper

let write_interface_function ch wrapper name
    (result_type, result, default) args processes =
  let arg_list, arg_list' = arguments args in
  fprintf ch "  function %s (pid%s) result (%s)\n" wrapper arg_list' result;
  List.iter (fun p -> import_prefixed ch p name) processes;
  List.iter (declare_argument ch) (("character(len=*)", "pid") :: args);
  fprintf ch "    %s :: %s\n" result_type result;
  fprintf ch "    select case (pid)\n";
  List.iter (fun p -> call_function ch p result name arg_list) processes;
  default_function ch result default;
  fprintf ch "    end select\n";
  fprintf ch "  end function %s\n" wrapper

let write_other_interface_functions ch =
  fprintf ch "  subroutine invalid_process (pid)\n";
  fprintf ch "    character(len=*), intent(in) :: pid\n";
  fprintf ch "    print *, \"PANIC:";
  fprintf ch " process `\"//trim(pid)//\"' not available!\"\n";
  fprintf ch "  end subroutine invalid_process\n";
  fprintf ch "  function n_tot (pid) result (n)\n";
  fprintf ch "    character(len=*), intent(in) :: pid\n";
  fprintf ch "    integer :: n\n";
  fprintf ch "    n = n_in(pid) + n_out(pid)\n";
  fprintf ch "  end function n_tot\n"

let write_other_declarations ch =
  fprintf ch "  public :: n_in, n_out, n_tot, pdg_code\n";
  fprintf ch "  public :: allow_helicities\n";
  fprintf ch "  public :: create, destroy\n";
  fprintf ch "  public :: set_const, sqme\n";
  fprintf ch "  interface create\n";
  fprintf ch "     module procedure process_create\n";
  fprintf ch "  end interface\n";
  fprintf ch "  interface destroy\n";
  fprintf ch "     module procedure process_destroy\n";
  fprintf ch "  end interface\n";
  fprintf ch "  interface set_const\n";
  fprintf ch "     module procedure process_set_const\n";
  fprintf ch "  end interface\n";
  fprintf ch "  interface sqme\n";
  fprintf ch "     module procedure process_sqme\n";
  fprintf ch "  end interface\n"

let write_interface ch names =
  fprintf ch "module process_interface\n";
  fprintf ch "  use kinds, only: default  !NODEP!\n";
  fprintf ch "  use parameters, only: parameter_set\n";
  fprintf ch "  implicit none\n";
  fprintf ch "  private\n";
  List.iter (fun p ->
    fprintf ch
      "  character(len=*), parameter, public :: pr_%s = \"%s\"\n" p p)
    names;
  write_other_declarations ch;
  fprintf ch "contains\n";
  write_interface_function ch "n_in" "n_in" ("integer", "n", "0") [] names;
  write_interface_function ch "n_out" "n_out" ("integer", "n", "0") [] names;
  write_interface_function ch "pdg_code" "pdg_code"
    ("integer", "n", "0") [ "integer", "i" ] names;
  write_interface_function ch "allow_helicities" "allow_helicities"
    ("logical", "yorn", ".false.") [] names;
  write_interface_subroutine ch "process_create" "create" [] names;
  write_interface_subroutine ch "process_destroy" "destroy" [] names;
  write_interface_subroutine ch "process_set_const" "set_const"
    [ "type(parameter_set)", "par"] names;
  write_interface_function ch "process_sqme" "sqme"
    ("real(kind=default)", "sqme", "0")
    [ "real(kind=default), dimension(0:,:)", "p";
      "integer, dimension(:), optional", "h" ] names;
  write_other_interface_functions ch;
  fprintf ch "end module process_interface\n"

(* \thocwmodulesection{Makefile} *)

let write_makefile ch names =
  fprintf ch "KINDS = ../@KINDS@\n";
  fprintf ch "HELAS = ../@HELAS@\n";
  fprintf ch "F90 = @F90@\n";
  fprintf ch "F90FLAGS = @F90FLAGS@\n";
  fprintf ch "F90INCL = -I$(KINDS) -I$(HELAS)\n";
  fprintf ch "F90COMMON = omega_bundle_whizard.f90";
  fprintf ch " file_utils.f90 process_interface.f90\n";
  fprintf ch "include Makefile.processes\n";
  fprintf ch "F90SRC = $(F90COMMON) $(F90PROCESSES)\n";
  fprintf ch "OBJ = $(F90SRC:.f90=.o)\n";
  fprintf ch "MOD = $(F90SRC:.f90=.mod)\n";
  fprintf ch "archive: processes.a\n";
  fprintf ch "processes.a: $(OBJ)\n";
  fprintf ch "\t$(AR) r $@ $(OBJ)\n";
  fprintf ch "\t@RANLIB@ $@\n";
  fprintf ch "clean:\n";
  fprintf ch "\trm -f $(OBJ)\n";
  fprintf ch "realclean:\n";
  fprintf ch "\trm -f processes.a\n";
  fprintf ch "parameters.o: file_utils.o\n";
  fprintf ch "omega_bundle_whizard.o: parameters.o\n";
  fprintf ch "process_interface.o: parameters.o\n";
  fprintf ch "%%.o: %%.f90 $(KINDS)/kinds.f90\n";
  fprintf ch "\t$(F90) $(F90FLAGS) $(F90INCL) -c $<\n"

let write_makefile_processes ch names =
  fprintf ch "F90PROCESSES =";
  List.iter (fun f -> fprintf ch " \\\n  %s.f90" f) names;
  fprintf ch "\n";
  List.iter (fun f ->
    fprintf ch "%s.o: omega_bundle_whizard.o parameters.o\n" f;
    fprintf ch "process_interface.o: %s.o\n" f) names

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
