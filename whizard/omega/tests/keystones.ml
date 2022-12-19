(* keystones.ml --

   Copyright (C) 2019-2022 by

       Wolfgang Kilian <kilian@physik.uni-siegen.de>
       Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
       Juergen Reuter <juergen.reuter@desy.de>

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

type field = Coupling.lorentz * int

type argument =
  | G of int (* complex coupling *)
  | N of int (* negative of complex coupling *)
  | M of int (* real mass (or width) *)
  | P of int (* momentum *)
  | F of field (* field *)
  | V of string (* verbatim *)

type keystone =
  { bra : field;
    name : string;
    args : argument list }

type vertex =
  { tag : string;
    keystones : keystone list }

let order_fields (_, i) (_, j) =
  compare i j

let extract_fields { bra; args } =
  List.sort
    order_fields
    (List.fold_left
       (fun acc arg ->
         match arg with
         | F f -> f :: acc
         | _ -> acc)
       [bra] args)

let check_indices field_list =
  if List.exists
       (fun (n, _) -> n > 1)
       (ThoList.classify (List.map snd field_list)) then
    invalid_arg "check_indices";
  ()

let spin_to_string = function
  | Coupling.Scalar -> "Scalar"
  | Coupling.Spinor -> "Spinor"
  | Coupling.ConjSpinor -> "ConjSpinor"
  | Coupling.Majorana -> "Majorana"
  | Coupling.Vector | Coupling.Massive_Vector -> "Vector"
  | Coupling.Tensor_2 -> "Tensor_2"
  | _ -> failwith "spin_to_string"

let fields_to_string fields =
  "[" ^
    String.concat
      "; " (List.map
              (fun (s, i) -> Printf.sprintf "%s(%d)" (spin_to_string s) i)
              fields) ^ "]"

let check_fields ks_list =
  let fields = List.map extract_fields ks_list in
  if not (ThoList.homogeneous fields) then
    begin
      let spins =
        "[" ^ String.concat "; " (List.map fields_to_string fields) ^ "]" in
      invalid_arg ("check_spins: " ^ spins)
    end;
  check_indices (List.hd fields)

open Format_Fortran

let spin_type = function
  | Coupling.Scalar -> "complex(kind=default)"
  | Coupling.Spinor -> "type(spinor)"
  | Coupling.ConjSpinor -> "type(conjspinor)"
  | Coupling.Majorana -> "type(bispinor)"
  | Coupling.Vector | Coupling.Massive_Vector -> "type(vector)"
  | Coupling.Tensor_2 -> "type(tensor)"
  | _ -> failwith "spin_type"

let type_arg = function
  | G _ -> Some "complex(kind=default)"
  | N _ -> Some "complex(kind=default)"
  | M _ -> Some "real(kind=default)"
  | P _ -> Some "type(momentum)"
  | F (s, _) -> Some (spin_type s)
  | V _ -> None

let spin_mnemonic = function
  | Coupling.Scalar -> "phi"
  | Coupling.Spinor -> "psi"
  | Coupling.ConjSpinor -> "psibar"
  | Coupling.Majorana -> "chi"
  | Coupling.Maj_Ghost -> "???"
  | Coupling.Vector -> "a"
  | Coupling.Massive_Vector -> "v"
  | Coupling.Tensor_2 -> "h"
  | _ -> failwith "spin_mnemonic"

let format_coupling i =
  Printf.sprintf "g%d" i

let format_negative_coupling i =
  Printf.sprintf "(-g%d)" i

let format_momentum i =
  Printf.sprintf "p%d" i

let format_mass i =
  Printf.sprintf "m%d" i

let format_field (s, i) =
  Printf.sprintf "%s%d" (spin_mnemonic s) i

let format_declaration = function
  | G i -> format_coupling i
  | N i -> format_coupling i
  | M i -> format_mass i
  | P i -> format_momentum i
  | F f -> format_field f
  | V s -> s

let format_arg = function
  | G i -> format_coupling i
  | N i -> format_negative_coupling i
  | M i -> format_mass i
  | P i -> format_momentum i
  | F f -> format_field f
  | V s -> s

let fusion_to_fortran ff name args =
  let printf fmt = fprintf ff fmt in
  match args with
  | [] -> invalid_arg "fusion_to_fortran"
  | arg1 :: arg2n ->
     printf "%s (%s" name (format_arg arg1);
     List.iter (fun arg -> printf ",@ %s" (format_arg arg)) arg2n;
     printf ")"

(* \begin{dubious}
     The ordering here works for Dirac spinors, but fails for
     Majorana spinors, leading to a sign ambiguity in this test
     \ldots
   \end{dubious} *)
let keystone_to_fortran ff (ksv, { bra; name; args }) =
  let printf fmt = fprintf ff fmt
  and nl = pp_newline ff in
  printf "      @[<2>%s =@ " ksv;
  begin match bra with
  | Coupling.Spinor, _ ->
     fusion_to_fortran ff name args;
     printf "@ * %s" (format_field bra)
  | Coupling.Majorana, _ ->
     begin match args with
     | _ :: F (Coupling.Majorana, _) :: _ ->
        fusion_to_fortran ff name args;
        printf "@ * %s" (format_field bra)
     | _ ->
        printf "%s@ * " (format_field bra);
        fusion_to_fortran ff name args
     end
  | _, _ -> 
     printf "%s@ * " (format_field bra);
     fusion_to_fortran ff name args
  end;
  printf "@]"; nl()

let keystones_to_subroutine ff { tag; keystones } =
  check_fields keystones;
  let printf fmt = fprintf ff fmt
  and nl = pp_newline ff in
  printf "  @[<4>subroutine@ testks_%s@ (repetitions," tag;
  printf "@ passed,@ threshold,@ quiet,@ abs_threshold)@]"; nl ();
  printf "    integer, intent(in) :: repetitions"; nl ();
  printf "    logical, intent(inout) :: passed"; nl ();
  printf "    logical, intent(in), optional :: quiet"; nl ();
  printf "    @[<2>real(kind=default),@ intent(in),@ optional ::";
  printf "@ threshold,@ abs_threshold@]"; nl ();
  printf "    integer :: i"; nl ();
  let ks1 = List.hd keystones in
  let all_momenta =
    List.map
      (fun i -> P i)
      (ThoList.range 0 (List.length (extract_fields ks1) - 1)) in
  let variables =
    ThoList.uniq (List.sort compare (F (ks1.bra) :: ks1.args @ all_momenta)) in
  List.iter
    (fun a ->
      match type_arg a with
      | None -> ()
      | Some t -> printf "    @[<2>%s :: %s@]" t (format_declaration a); nl ())
    variables;
  let ks_list =
    List.map
      (fun (n, ks) -> (Printf.sprintf "ks%d" n, ks))
      (ThoList.enumerate 0 keystones) in
  begin match ks_list with
  | [] -> failwith "keystones_to_fortran"
  | (ksv1, _) :: ks2n ->
     printf "    @[<2>complex(kind=default) ::@ %s" ksv1;
     List.iter (fun (ksv, _) -> printf ",@ %s" ksv) ks2n;
     printf "@]"; nl ()
  end;
  printf "    do i = 1, repetitions"; nl ();
  List.iter
    (fun a ->
      match a with
      | P 0 -> () (* this will be determined by momentum conservation! *)
      | V _ -> ()
      | a ->
         printf "      @[<2>call@ make_random@ (%s)@]" (format_arg a); nl ())
    variables;
  begin match all_momenta with
  | [] -> failwith "keystones_to_fortran"
  | p1 :: p2n ->
     printf "      @[<2>%s =" (format_arg p1);
     List.iter (fun p -> printf "@ - %s" (format_arg p)) p2n;
     printf "@]"; nl ()
  end;
  List.iter (keystone_to_fortran ff) ks_list;
  begin match ks_list with
  | [] -> failwith "keystones_to_fortran"
  | (ksv1, ks1) :: ks2n ->
     List.iter
       (fun (ksv, ks) ->
         printf "      @[<8>call@ expect@ (%s,@ %s," ksv ksv1;
         printf "@ '%s: %s <> %s'," tag ks.name ks1.name;
         printf "@ passed,@ threshold, quiet, abs_threshold)@]";
         nl ())
       ks2n
  end;
  printf "    end do"; nl ();
  printf "  @[<2>end@ subroutine@ testks_%s@]" tag; nl ()

let keystones_to_fortran
      ff ?(reps=1000) ?(threshold=0.85)
      ?(program="keystones_test") ?(omega_module="omega95")
      ?(modules=[]) vertices =
  let printf fmt = fprintf ff fmt
  and nl = pp_newline ff in
  printf "program %s" program; nl ();
  List.iter
    (fun m -> printf "  use %s" m; nl ())
    ("kinds" :: "constants" :: omega_module ::
       "omega_testtools" :: "keystones_tools" :: modules);
  printf "  implicit none"; nl ();
  printf "  logical :: passed"; nl ();
  printf "  logical, parameter :: quiet = .false."; nl ();
  printf "  integer, parameter :: reps = %d" reps; nl ();
  printf "  real(kind=default), parameter :: threshold = %f" threshold; nl ();
  printf "  real(kind=default), parameter :: abs_threshold = 1E-17"; nl ();
  printf "  integer, dimension(8) :: date_time"; nl ();
  printf "  integer :: rsize"; nl ();
  printf "  call date_and_time (values = date_time)"; nl ();
  printf "  call random_seed (size = rsize)"; nl ();
  printf "  @[<8>call random_seed@ (put = spread (product (date_time),";
  printf "@ dim = 1,@ ncopies = rsize))@]"; nl ();
  printf "  passed = .true."; nl ();
  List.iter
    (fun v ->
      printf "  @[<8>call testks_%s@ (reps,@ passed," v.tag;
      printf "@ threshold, quiet, abs_threshold)@]"; nl ())
    vertices;
  printf "  if (passed) then"; nl ();
  printf "    stop 0"; nl ();
  printf "  else"; nl ();
  printf "    stop 1"; nl ();
  printf "  end if"; nl ();
  printf "contains"; nl ();
  List.iter (keystones_to_subroutine ff) vertices;
  printf "end program %s" program; nl ()

let generate ?reps ?threshold ?program ?omega_module ?modules vertices =
  let my_name = Sys.argv.(0) in
  let verbose = ref false
  and cat = ref false
  and usage = "usage: " ^ my_name ^ " ..." in
  Arg.parse
    (Arg.align 
       [ ("-cat", Arg.Set cat, " print test snippets");
	 ("-v", Arg.Set verbose, " be more verbose");
	 ("-verbose", Arg.Set verbose, " be more verbose") ])
    (fun s -> raise (Arg.Bad s))
    usage;
  if !cat then
    keystones_to_fortran
      std_formatter ?reps ?threshold ?program ?omega_module ?modules vertices

type ufo_vertex =
  { v_tag : string;
    v_spins : Coupling.lorentz array;
    v_tensor : UFO_Lorentz.t;
    v_flines : Coupling.fermion_lines }

type ufo_propagator =
  { p_tag : string;
    p_omega : string;
    p_spins : Coupling.lorentz * Coupling.lorentz;
    p_propagator : UFO.Propagator.t }

let transpose p =
  { p_tag = p.p_tag;
    p_omega = p.p_omega;
    p_spins = (snd p.p_spins, fst p.p_spins);
    p_propagator = UFO.Propagator.transpose p.p_propagator }

let equivalent_tensors ?(fermion_lines=[]) v_spins alternatives =
  List.map
    (fun (v_tag, tensor) ->
      let v_tensor =
        UFO_Lorentz.parse
          (Array.to_list v_spins)
          (UFOx.Lorentz.of_string tensor) in
      { v_tag; v_spins; v_tensor; v_flines = fermion_lines })
    alternatives

module P = Permutation.Default

let permute_spins p s = P.array p s

(* We must permute only the free indices, of course.
   Note that we apply the \emph{inverse} permutation to
   the indices in order to match the permutation of the
   particles/spins. *)

(* The following is copied from [UFO.Lorentz.permute].
   We can't simply call it, because the types [UFO.Lorentz.t]
   and [vertex] differ.
   This should be changed to make sure that we're also
   testing [UFO.Lorentz.permute], but note
   that only the ["_p" ^ permutation] naming convention
   and the simultaneous exchange of indices in Lorentz structures
   and fermion lines is relevant for applications! *)
let permute_structure n p (l, f) =
  let permuted = P.array (P.inverse p) (Array.init n succ) in
  let permute_index i =
    if i > 0 then
      UFOx.Index.map_position (fun pos -> permuted.(pred pos)) i
    else
      i in
  (UFO_Lorentz.map_indices permute_index l,
   UFO_Lorentz.map_fermion_lines permute_index f)

let permute_vertex n v p =
  let v_tensor, v_flines = permute_structure n p (v.v_tensor, v.v_flines) in
  { v_tag = v.v_tag ^ "_p" ^ P.to_string (P.inverse p);
    v_spins = permute_spins p v.v_spins;
    v_tensor;
    v_flines }

let vertex_permutations v =
  let n = Array.length v.v_spins in
  List.map (permute_vertex n v) (P.cyclic n)

(* The following is mostly copied from
   [UFO.Lorentz.all_charge_conjugates].
   We can't simply call it, because the types [UFO.Lorentz.t]
   and [vertex] differ.
   This should be changed to make sure that we're also
   testing [UFO.Lorentz.all_charge_conjugates], but note
   that only [UFO_Lorentz.charge_conjugate] and the ["_c%x%x"]
   naming convention is relevant for applications!
   Note also that we're \emph{only} charge conjugating
   fermion lines involving Majoranas. *)
let charge_conjugate1 v (bra, ket as fermion_line) =
  { v_tag = v.v_tag ^ Printf.sprintf "_c%x%x" bra ket;
    v_spins = v.v_spins;
    v_tensor = UFO_Lorentz.charge_conjugate fermion_line v.v_tensor;
    v_flines = v.v_flines }

let charge_conjugate l fermion_lines =
  List.fold_left charge_conjugate1 l fermion_lines

let is_majorana = function
  | Coupling.Majorana | Coupling.Vectorspinor | Coupling.Maj_Ghost -> true
  | _ -> false

let is_majorana_fline v_spins (bra, ket) =
  is_majorana v_spins.(pred bra) || is_majorana v_spins.(pred ket)

(*i
let all_charge_conjugates l =
  List.map
    (charge_conjugate l)
    (ThoList.power (List.filter (is_majorana_fline l.v_spins) l.v_flines))
i*)

let required_charge_conjugates l =
  let saturated_fermion_lines =
    List.filter (fun (bra, ket) -> bra != 1 && ket != 1) l.v_flines in
  List.map
    (charge_conjugate l)
    (ThoList.power
       (List.filter (is_majorana_fline l.v_spins) saturated_fermion_lines))

let keystones_of_ufo_vertex { v_tag; v_spins } =
  { tag = v_tag;
    keystones =
      let fields = Array.mapi (fun i s -> (s, i)) v_spins in
      let n = Array.length fields in
      List.map
        (fun p ->
          let permuted = P.array p fields in
          match Array.to_list permuted with
          | [] -> invalid_arg "keystones_of_ufo_vertex"
          | bra :: args ->
             { bra;
               name = v_tag ^ "_p" ^ P.to_string (P.inverse p);
               args =
                 G (0) ::
                   (ThoList.flatmap (fun (s, i) -> [ F (s, i); P (i) ]) args) })
        (P.cyclic n) }

let keystones_of_propagator { p_tag; p_omega; p_spins } =
  let s0, s1 = p_spins in
  let keystone omega name =
    match omega, s1, name with
    | _, (Coupling.Scalar|Coupling.Tensor_2), _
    | false, _, _ ->
       { bra = (s0, 0);
         name;
         args = [P (1); M (0); M (1); F (s1, 1) ] }
    | _, Coupling.Vector, "pr_gauge" ->
       { bra = (s0, 0);
         name;
         args = [P (1); V ("42.0_default"); F (s1, 1) ] }
    | _, Coupling.Vector, "pr_rxi" ->
       { bra = (s0, 0);
         name;
         args = [P (1); M (0); M (1); V ("42.0_default"); F (s1, 1) ] }
    | _, Coupling.Vector, _ ->
       { bra = (s0, 0);
         name;
         args = [P (1); F (s1, 1) ] }
    | true, _, _ ->
       { bra = (s0, 0);
         name;
         args = [P (1); M (0); M (1); V (".false."); F (s1, 1) ] } in
  { tag = p_tag;
    keystones = [keystone false ("pr_U_" ^ p_tag); keystone true p_omega] }

let merge (ufo_list, omegalib) =
  match ufo_list with
  | [] -> omegalib
  | ufo1 :: _ ->
     { tag = ufo1.v_tag;
       keystones =
         (omegalib.keystones
          @ ThoList.flatmap
              (fun ufo -> (keystones_of_ufo_vertex ufo).keystones)
              ufo_list) }

let fusions ff ?(omega_module="omega95") module_name vertices propagators =
  let printf fmt = fprintf ff fmt
  and nl () = pp_newline ff () in
  printf "module %s" module_name; nl ();
  printf "  use kinds"; nl ();
  printf "  use %s" omega_module; nl ();
  printf "  implicit none"; nl ();
  printf "  private"; nl ();
  let permuted_vertices =
    ThoList.flatmap
      required_charge_conjugates
      (ThoList.flatmap vertex_permutations vertices) in
  List.iter
    (fun v -> printf "  public :: %s" v.v_tag; nl ())
    permuted_vertices;
  List.iter
    (fun p -> printf "  public :: pr_U_%s" p.p_tag; nl ())
    propagators;
  UFO_targets.Fortran.eps4_g4_g44_decl std_formatter ();
  UFO_targets.Fortran.eps4_g4_g44_init std_formatter ();
  printf "contains"; nl ();
  UFO_targets.Fortran.inner_product_functions std_formatter ();
  List.iter
    (fun v ->
      printf "  ! %s" (String.make 68 '='); nl ();
      printf "  ! %s" (UFO_Lorentz.to_string v.v_tensor); nl ();
      UFO_targets.Fortran.lorentz std_formatter v.v_tag v.v_spins v.v_tensor)
    permuted_vertices;
  List.iter
    (fun p ->
      UFO_targets.Fortran.propagator
        std_formatter p.p_tag
        "parameters" p.p_propagator.UFO.Propagator.variables p.p_spins
        p.p_propagator.UFO.Propagator.numerator
        p.p_propagator.UFO.Propagator.denominator)
    propagators;
  printf "end module %s" module_name; nl ()

let generate_ufo ?program ?omega_module ?reps ?threshold
      ?(only_fusions=[]) module_name vertices propagators =
  fusions
    ?omega_module std_formatter module_name
    (only_fusions @ ThoList.flatmap fst vertices) propagators;
  generate
    ?reps ?threshold ?program ?omega_module ~modules:[module_name]
    (List.map merge vertices @ List.map keystones_of_propagator propagators)

(* \begin{dubious}
     placeholder:
   \end{dubious} *)

let generate_ufo_bispinors ?program ?omega_module ?reps ?threshold
      ?(only_fusions=[]) module_name vertices propagators =
  fusions
    ?omega_module std_formatter module_name
    (only_fusions @ ThoList.flatmap fst vertices) propagators;
  generate
    ?reps ?threshold ?program ?omega_module ~modules:[module_name]
    (List.map merge vertices @ List.map keystones_of_propagator propagators)

