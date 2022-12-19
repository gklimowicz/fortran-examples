(* UFO_targets.ml --

   Copyright (C) 1999-2017 by

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

(* \thocwmodulesection{Generating Code for UFO Lorentz Structures} *)

(* \begin{dubious}
     O'Caml before 4.02 had a module typing bug that forced us to put
     these definitions outside of [Lorentz_Fusion].  Since then, they
     might have appeared in more places.  Investigate, if it is
     worthwhile to encapsulate them again.
   \end{dubious} *)

module Q = Algebra.Q
module QC = Algebra.QC

module type T =
  sig

    (* [lorentz formatter name spins v]
       writes a representation of the Lorentz structure [v] of
       particles with the Lorentz representations [spins] as a
       (Fortran) function [name] to [formatter]. *)
    val lorentz :
      Format_Fortran.formatter -> string ->
      Coupling.lorentz array -> UFO_Lorentz.t -> unit

    val propagator :
      Format_Fortran.formatter -> string -> string -> string list ->
      Coupling.lorentz * Coupling.lorentz ->
      UFO_Lorentz.t -> UFO_Lorentz.t -> unit

    val fusion_name :
      string -> Permutation.Default.t -> Coupling.fermion_lines -> string

    val fuse :
      Algebra.QC.t -> string ->
      Coupling.lorentzn -> Coupling.fermion_lines ->
      string -> string list -> string list -> Coupling.fusen -> unit

    val eps4_g4_g44_decl : Format_Fortran.formatter -> unit -> unit
    val eps4_g4_g44_init : Format_Fortran.formatter -> unit -> unit
    val inner_product_functions : Format_Fortran.formatter -> unit -> unit

    module type Test =
      sig
        val suite : OUnit.test
      end

    module Test : Test

  end

module Fortran : T =
  struct

    open Format_Fortran

    let pp_divide ?(indent=0) ff () =
      fprintf ff "%*s! %s" indent "" (String.make (70 - indent) '-');
      pp_newline ff ()

    let conjugate = function
      | Coupling.Spinor -> Coupling.ConjSpinor
      | Coupling.ConjSpinor -> Coupling.Spinor
      | r -> r

    let spin_mnemonic = function
      | Coupling.Scalar -> "phi"
      | Coupling.Spinor -> "psi"
      | Coupling.ConjSpinor -> "psibar"
      | Coupling.Majorana -> "chi"
      | Coupling.Maj_Ghost ->
         invalid_arg "UFO_targets: Maj_Ghost"
      | Coupling.Vector -> "a"
      | Coupling.Massive_Vector -> "v"
      | Coupling.Vectorspinor -> "grav" (* itino *)
      | Coupling.Tensor_1 ->
         invalid_arg "UFO_targets: Tensor_1"
      | Coupling.Tensor_2 -> "h"
      | Coupling.BRS l ->
         invalid_arg "UFO_targets: BRS"

    let fortran_type = function
      | Coupling.Scalar -> "complex(kind=default)"
      | Coupling.Spinor -> "type(spinor)"
      | Coupling.ConjSpinor -> "type(conjspinor)"
      | Coupling.Majorana -> "type(bispinor)"
      | Coupling.Maj_Ghost ->
         invalid_arg "UFO_targets: Maj_Ghost"
      | Coupling.Vector -> "type(vector)"
      | Coupling.Massive_Vector -> "type(vector)"
      | Coupling.Vectorspinor -> "type(vectorspinor)"
      | Coupling.Tensor_1 ->
         invalid_arg "UFO_targets: Tensor_1"
      | Coupling.Tensor_2 -> "type(tensor)"
      | Coupling.BRS l ->
         invalid_arg "UFO_targets: BRS"

    (* The \texttt{omegalib} separates time from space.  Maybe
       not a good idea after all.  Mend it locally \ldots *)
    type wf =
      { pos : int;
        spin : Coupling.lorentz;
        name : string;
        local_array : string option;
        momentum : string;
        momentum_array : string;
        fortran_type : string }

    let wf_table spins =
      Array.mapi
        (fun i s ->
          let spin =
            if i = 0 then
              conjugate s
            else
              s in
          let pos = succ i in
          let i = string_of_int pos in
          let name = spin_mnemonic s ^ i in
          let local_array =
            begin match spin with
            | Coupling.Vector | Coupling.Massive_Vector -> Some (name ^ "a")
            | _ -> None
            end in
          { pos;
            spin;
            name;
            local_array;
            momentum = "k" ^ i;
            momentum_array = "p" ^ i;
            fortran_type = fortran_type spin } )
        spins

    module L = UFO_Lorentz

    (* Format rational ([Q.t]) and complex rational ([QC.t])
       numbers as fortran values. *)
    let format_rational q =
      if Q.is_integer q then
        string_of_int (Q.to_integer q)
      else
        let n, d = Q.to_ratio q in
        Printf.sprintf "%d.0_default/%d" n d

    let format_complex_rational cq =
      let real = QC.real cq
      and imag = QC.imag cq in
      if Q.is_null imag then
        begin
          if Q.is_negative real then
            "(" ^ format_rational real ^ ")"
          else
            format_rational real
        end
      else if Q.is_integer real && Q.is_integer imag then
        Printf.sprintf "(%d,%d)" (Q.to_integer real) (Q.to_integer imag)
      else
        Printf.sprintf
          "cmplx(%s,%s,kind=default)"
          (format_rational real) (format_rational imag)

    (* Optimize the representation if used as a prefactor of
       a summand in a sum. *)
    let format_rational_factor q =
      if Q.is_unit q then
        "+ "
      else if Q.is_unit (Q.neg q) then
        "- "
      else if Q.is_negative q then
        "- " ^ format_rational (Q.neg q) ^ "*"
      else
        "+ " ^ format_rational q ^ "*"

    let format_complex_rational_factor cq =
      let real = QC.real cq
      and imag = QC.imag cq in
      if Q.is_null imag then
        begin
          if Q.is_unit real then
            "+ "
          else if Q.is_unit (Q.neg real) then
            "- "
          else if Q.is_negative real then
            "- " ^ format_rational (Q.neg real) ^ "*"
          else
            "+ " ^ format_rational real ^ "*"
        end
      else if Q.is_integer real && Q.is_integer imag then
        Printf.sprintf "+ (%d,%d)*" (Q.to_integer real) (Q.to_integer imag)
      else
        Printf.sprintf
          "+ cmplx(%s,%s,kind=default)*"
          (format_rational real) (format_rational imag)

    (* Append a formatted list of indices to [name]. *)
    let append_indices name = function
      | [] -> name
      | indices ->
         name ^ "(" ^ String.concat "," (List.map string_of_int indices) ^ ")"

    (* Dirac string variables and their names. *)
    type dsv =
      | Ket of int
      | Bra of int
      | Braket of int

    let dsv_name = function
      | Ket n -> Printf.sprintf "ket%02d" n
      | Bra n -> Printf.sprintf "bra%02d" n
      | Braket n -> Printf.sprintf "bkt%02d" n

    let dirac_dimension dsv indices =
      let tail ilist =
        String.concat "," (List.map (fun _ -> "0:3") ilist) ^ ")" in
      match dsv, indices with
      | Braket _, [] -> ""
      | (Ket _ | Bra _), [] -> ", dimension(1:4)"
      | Braket _, indices -> ", dimension(" ^ tail indices
      | (Ket _ | Bra _), indices -> ", dimension(1:4," ^ tail indices

    (* Write Fortran code to [decl] and [eval]: apply the Dirac matrix
       [gamma] with complex rational entries to the spinor [ket] from
       the left. [ket] must be the name of a scalar variable and cannot
       be an array element.  The result is stored in [dsv_name (Ket n)]
       which can have additional [indices].  Return [Ket n] for further
       processing. *)
    let dirac_ket_to_fortran_decl ff n indices =
      let printf fmt = fprintf ff fmt
      and nl = pp_newline ff in
      let dsv = Ket n in
      printf
        "    @[<2>complex(kind=default)%s ::@ %s@]"
        (dirac_dimension dsv indices) (dsv_name dsv);
      nl ()

    let dirac_ket_to_fortran_eval ff n indices gamma ket =
      let printf fmt = fprintf ff fmt
      and nl = pp_newline ff in
      let dsv = Ket n in
      for i = 0 to 3 do
        let name = append_indices (dsv_name dsv) (succ i :: indices) in
        printf "    @[<%d>%s = 0" (String.length name + 4) name;
        for j = 0 to 3 do
          if not (QC.is_null gamma.(i).(j)) then
            printf
              "@ %s%s%%a(%d)"
              (format_complex_rational_factor gamma.(i).(j))
              ket.name (succ j)
        done;
        printf "@]";
        nl ()
      done;
      dsv

    (* The same as [dirac_ket_to_fortran], but apply the Dirac matrix
       [gamma] to [bra] from the right and return [Bra n]. *)
    let dirac_bra_to_fortran_decl ff n indices =
      let printf fmt = fprintf ff fmt
      and nl = pp_newline ff in
      let dsv = Bra n in
      printf
        "    @[<2>complex(kind=default)%s ::@ %s@]"
        (dirac_dimension dsv indices) (dsv_name dsv);
      nl ()

    let dirac_bra_to_fortran_eval ff n indices bra gamma =
      let printf fmt = fprintf ff fmt
      and nl = pp_newline ff in
      let dsv = Bra n in
      for j = 0 to 3 do
        let name = append_indices (dsv_name dsv) (succ j :: indices) in
        printf "    @[<%d>%s = 0" (String.length name + 4) name;
        for i = 0 to 3 do
          if not (QC.is_null gamma.(i).(j)) then
            printf
              "@ %s%s%%a(%d)"
              (format_complex_rational_factor gamma.(i).(j))
              bra.name (succ i)
        done;
        printf "@]";
        nl ()
      done;
      dsv

    (* More of the same, but evaluating a spinor sandwich and
       returning [Braket n]. *)
    let dirac_braket_to_fortran_decl ff n indices =
      let printf fmt = fprintf ff fmt
      and nl = pp_newline ff in
      let dsv = Braket n in
      printf
        "    @[<2>complex(kind=default)%s ::@ %s@]"
        (dirac_dimension dsv indices) (dsv_name dsv);
      nl ()

    let dirac_braket_to_fortran_eval ff n indices bra gamma ket =
      let printf fmt = fprintf ff fmt
      and nl = pp_newline ff in
      let dsv = Braket n in
      let name = append_indices (dsv_name dsv) indices in
      printf "    @[<%d>%s = 0" (String.length name + 4) name;
      for i = 0 to 3 do
        for j = 0 to 3 do
          if not (QC.is_null gamma.(i).(j)) then
            printf
              "@ %s%s%%a(%d)*%s%%a(%d)"
              (format_complex_rational_factor gamma.(i).(j))
              bra.name (succ i) ket.name (succ j)
        done
      done;
      printf "@]";
      nl ();
      dsv

    (* Choose among the previous functions according to the position
       of [bra] and [ket] among the wavefunctions.  If any is in the
       first position evaluate the spinor expression with the
       corresponding spinor removed, otherwise evaluate the
       spinir sandwich. *)
    let dirac_bra_or_ket_to_fortran_decl ff n indices bra ket =
      if bra = 1 then
        dirac_ket_to_fortran_decl ff n indices
      else if ket = 1 then
        dirac_bra_to_fortran_decl ff n indices
      else
        dirac_braket_to_fortran_decl ff n indices

    let dirac_bra_or_ket_to_fortran_eval ff n indices wfs bra gamma ket =
      if bra = 1 then
        dirac_ket_to_fortran_eval ff n indices gamma wfs.(pred ket)
      else if ket = 1 then
        dirac_bra_to_fortran_eval ff n indices wfs.(pred bra) gamma
      else
        dirac_braket_to_fortran_eval
          ff n indices wfs.(pred bra) gamma wfs.(pred ket)

    (* UFO summation indices are negative integers.  Derive a valid Fortran
       variable name. *)
    let prefix_summation = "mu"
    let prefix_polarization = "nu"
    let index_spinor = "alpha"
    let index_tensor = "nu"

    let index_variable mu =
      if mu < 0 then
        Printf.sprintf "%s%d" prefix_summation (- mu)
      else if mu == 0 then
       prefix_polarization
      else
        Printf.sprintf "%s%d" prefix_polarization mu

    let format_indices indices =
      String.concat "," (List.map index_variable indices)

    module IntPM =
      Partial.Make (struct type t = int let compare = compare end)

    type tensor =
      | DS of dsv
      | V of string
      | T of UFOx.Lorentz_Atom.vector
      | S of UFOx.Lorentz_Atom.scalar
      | Inv of UFOx.Lorentz_Atom.scalar

    (* Transform the Dirac strings if we have Majorana
       fermions involved, in order to implement the algorithm
       from JRR's thesis. NB:
       The following is for reference only, to better understand what JRR
       was doing\ldots *)

    (* If the vertex is (suppressing the Lorentz indices of~$\phi_2$ and~$\Gamma$)
       \begin{equation}
       \label{eq:FVF-Vertex}
         \bar\psi \Gamma\phi \psi
            = \Gamma_{\alpha\beta} \bar\psi_{\alpha} \phi \psi_{\beta}
       \end{equation}
       (cf.~[Coupling.FBF] in the hardcoded O'Mega models),
       then this is the version implemented by [fuse] below. *)

    let tho_print_dirac_current f c wf1 wf2 fusion =
      match fusion with
      | [1; 3] -> printf "%s_ff(%s,%s,%s)" f c wf1 wf2 (* $\Gamma_{\alpha\beta} \bar\psi_{1,\alpha} \psi_{2,\beta}$ *)
      | [3; 1] -> printf "%s_ff(%s,%s,%s)" f c wf2 wf1 (* $\Gamma_{\alpha\beta} \bar\psi_{1,\alpha} \psi_{2,\beta}$ *)
      | [2; 3] -> printf "f_%sf(%s,%s,%s)" f c wf1 wf2 (* $\Gamma_{\alpha\beta} \phi_1 \psi_{2,\beta}$ *)
      | [3; 2] -> printf "f_%sf(%s,%s,%s)" f c wf2 wf1 (* $\Gamma_{\alpha\beta} \phi_1 \psi_{2,\beta}$ *)
      | [1; 2] -> printf "f_f%s(%s,%s,%s)" f c wf1 wf2 (* $\Gamma_{\alpha\beta} \bar\psi_{1,\alpha} \phi_2$ *)
      | [2; 1] -> printf "f_f%s(%s,%s,%s)" f c wf2 wf1 (* $\Gamma_{\alpha\beta} \bar\psi_{1,\alpha} \phi_2$ *)
      | _ -> ()

    (* The corresponding UFO [fuse] exchanges the arguments in the case
       of two fermions.  This is the natural choice for cyclic permutations. *)

    let tho_print_FBF_current f c wf1 wf2 fusion =
      match fusion with
      | [3; 1] -> printf "f%sf_p120(%s,%s,%s)" f c wf1 wf2 (* $\Gamma_{\alpha\beta} \psi_{1,\beta} \bar\psi_{2,\alpha}$ *)
      | [1; 3] -> printf "f%sf_p120(%s,%s,%s)" f c wf2 wf1 (* $\Gamma_{\alpha\beta} \psi_{1,\beta} \bar\psi_{2,\alpha}$ *)
      | [2; 3] -> printf "f%sf_p012(%s,%s,%s)" f c wf1 wf2 (* $\Gamma_{\alpha\beta} \phi_1 \psi_{2,\beta}$ *)
      | [3; 2] -> printf "f%sf_p012(%s,%s,%s)" f c wf2 wf1 (* $\Gamma_{\alpha\beta} \phi_1 \psi_{2,\beta}$ *)
      | [1; 2] -> printf "f%sf_p201(%s,%s,%s)" f c wf1 wf2 (* $\Gamma_{\alpha\beta} \bar\psi_{1,\alpha} \phi_2$ *)
      | [2; 1] -> printf "f%sf_p201(%s,%s,%s)" f c wf2 wf1 (* $\Gamma_{\alpha\beta} \bar\psi_{1,\alpha} \phi_2$ *)
      | _ -> ()

    (* This is how JRR implemented
       (see subsection~\ref{sec:dirac-matrices-jrr}) the Dirac matrices
       that don't change sign under $C\Gamma^T C^{-1} = \Gamma$,
       i.\,e.~$\mathbf{1}$, $\gamma_5$ and~$\gamma_5\gamma_\mu$
       (see [Targets.Fortran_Majorana_Fermions.print_fermion_current])
       \begin{itemize}
         \item In the case of two fermions, the second wave
           function [wf2] is always put into the second slot,
           as described in JRR's thesis.
           \label{pg:JRR-Fusions}
         \item In the case of a boson and a fermion, there is no
           need for both ["f_%sf"] and ["f_f%s"], since the
           latter can be obtained by exchanging arguments.
       \end{itemize} *)

    let jrr_print_majorana_current_S_P_A f c wf1 wf2 fusion =
      match fusion with
      | [1; 3] -> printf "%s_ff(%s,%s,%s)" f c wf1 wf2 (*
        $(C\Gamma)_{\alpha\beta} \bar\psi_{1,\alpha} \psi_{2,\beta} \cong
         C\Gamma $ *)
      | [3; 1] -> printf "%s_ff(%s,%s,%s)" f c wf1 wf2 (*
        $(C\Gamma)_{\alpha\beta} \psi_{1,\alpha} \bar\psi_{2,\beta} \cong
         C\Gamma = C\,C\Gamma^T C^{-1} $ *)
      | [2; 3] -> printf "f_%sf(%s,%s,%s)" f c wf1 wf2 (*
        $\Gamma_{\alpha\beta} \phi_1 \psi_{2,\beta} \cong
         \Gamma $ *)
      | [3; 2] -> printf "f_%sf(%s,%s,%s)" f c wf2 wf1 (*
        $\Gamma_{\alpha\beta} \phi_1 \psi_{2,\beta} \cong
         \Gamma $ *)
      | [1; 2] -> printf "f_%sf(%s,%s,%s)" f c wf2 wf1 (*
        $\Gamma_{\alpha\beta} \phi_1 \bar\psi_{2,\beta} \cong
         \Gamma = C\Gamma^T C^{-1} $ *)
      | [2; 1] -> printf "f_%sf(%s,%s,%s)" f c wf1 wf2 (*
        $\Gamma_{\alpha\beta} \phi_1 \bar\psi_{2,\beta} \cong
         \Gamma = C\Gamma^T C^{-1} $ *)
      | _ -> ()

    (* This is how JRR implemented the Dirac matrices
       that do change sign under $C\Gamma^T C^{-1} = - \Gamma$,
       i.\,e.~$\gamma_\mu$ and~$\sigma_{\mu\nu}$
       (see [Targets.Fortran_Majorana_Fermions.print_fermion_current_vector]). *)

    let jrr_print_majorana_current_V f c wf1 wf2 fusion =
      match fusion with
      | [1; 3] -> printf "%s_ff( %s,%s,%s)" f c wf1 wf2 (*
        $ (C\Gamma)_{\alpha\beta} \bar\psi_{1,\alpha} \psi_{2,\beta} \cong
          C\Gamma $ *)
      | [3; 1] -> printf "%s_ff(-%s,%s,%s)" f c wf1 wf2 (*
        $-(C\Gamma)_{\alpha\beta} \psi_{1,\alpha} \bar\psi_{2,\beta}  \cong
         -C\Gamma = C\,C\Gamma^T C^{-1} $ *)
      | [2; 3] -> printf "f_%sf( %s,%s,%s)" f c wf1 wf2 (*
        $ \Gamma_{\alpha\beta} \phi_1 \psi_{2,\beta} \cong
          \Gamma $ *)
      | [3; 2] -> printf "f_%sf( %s,%s,%s)" f c wf2 wf1 (*
        $ \Gamma_{\alpha\beta} \phi_1 \psi_{2,\beta} \cong
          \Gamma $ *)
      | [1; 2] -> printf "f_%sf(-%s,%s,%s)" f c wf2 wf1 (*
        $-\Gamma_{\alpha\beta} \phi_1 \bar\psi_{2,\beta} \cong
         -\Gamma = C\Gamma^T C^{-1} $ *)
      | [2; 1] -> printf "f_%sf(-%s,%s,%s)" f c wf1 wf2 (*
        $-\Gamma_{\alpha\beta} \phi_1 \bar\psi_{2,\beta} \cong
         -\Gamma = C\Gamma^T C^{-1} $ *)
      | _ -> ()

    (* These two can be unified, if the \texttt{\_c} functions
       implement~$\Gamma'=C\Gamma^T C^{-1}$, but we \emph{must}
       make sure that the multiplication with~$C$ from the left
       happens \emph{after} the transformation~$\Gamma\to\Gamma'$. *)
    let jrr_print_majorana_current f c wf1 wf2 fusion =
      match fusion with
      | [1; 3] -> printf "%s_ff  (%s,%s,%s)" f c wf1 wf2 (*
        $ (C\Gamma)_{\alpha\beta} \bar\psi_{1,\alpha} \psi_{2,\beta} \cong
          C\Gamma $ *)
      | [3; 1] -> printf "%s_ff_c(%s,%s,%s)" f c wf1 wf2 (*
        $(C\Gamma')_{\alpha\beta} \psi_{1,\alpha} \bar\psi_{2,\beta} \cong
         C\Gamma' = C\,C\Gamma^T C^{-1} $ *)
      | [2; 3] -> printf "f_%sf  (%s,%s,%s)" f c wf1 wf2 (*
        $ \Gamma_{\alpha\beta} \phi_1 \psi_{2,\beta} \cong
          \Gamma $ *)
      | [3; 2] -> printf "f_%sf  (%s,%s,%s)" f c wf2 wf1 (*
        $ \Gamma_{\alpha\beta} \phi_1 \psi_{2,\beta} \cong
          \Gamma $ *)
      | [1; 2] -> printf "f_%sf_c(%s,%s,%s)" f c wf2 wf1 (*
        $\Gamma'_{\alpha\beta} \phi_1 \bar\psi_{2,\beta} \cong
         \Gamma' = C\Gamma^T C^{-1} $ *)
      | [2; 1] -> printf "f_%sf_c(%s,%s,%s)" f c wf1 wf2 (*
        $\Gamma'_{\alpha\beta} \phi_1 \bar\psi_{2,\beta} \cong
         \Gamma' = C\Gamma^T C^{-1} $ *)
      | _ -> ()

    (* Since we may assume~$C^{-1}=-C=C^T$, this can be rewritten
       if the \texttt{\_c} functions implement
       \begin{equation}
          \Gamma^{\prime\,T}
            = \left(C\Gamma^T C^{-1}\right)^T
            = \left(C^{-1}\right)^T \Gamma C^T
            = C \Gamma C^{-1} 
       \end{equation}
       instead. *)

    let jrr_print_majorana_current_transposing f c wf1 wf2 fusion =
      match fusion with
      | [1; 3] -> printf "%s_ff  (%s,%s,%s)" f c wf1 wf2 (*
        $ (C\Gamma)_{\alpha\beta} \bar\psi_{1,\alpha} \psi_{2,\beta} \cong
          C\Gamma $ *)
      | [3; 1] -> printf "%s_ff_c(%s,%s,%s)" f c wf2 wf1 (*
        $(C\Gamma')^T_{\alpha\beta} \bar\psi_{1,\alpha} \psi_{2,\beta}  \cong
         (C\Gamma')^T = - C\Gamma $ *)
      | [2; 3] -> printf "f_%sf  (%s,%s,%s)" f c wf1 wf2 (*
        $ \Gamma_{\alpha\beta} \phi_1 \psi_{2,\beta} \cong
          \Gamma $ *)
      | [3; 2] -> printf "f_%sf  (%s,%s,%s)" f c wf2 wf1 (*
        $ \Gamma_{\alpha\beta} \phi_1 \psi_{2,\beta} \cong
          \Gamma $ *)
      | [1; 2] -> printf "f_f%s_c(%s,%s,%s)" f c wf1 wf2 (*
        $\Gamma^{\prime\,T}_{\alpha\beta} \bar\psi_{1,\alpha} \phi_2 \cong
         \Gamma^{\prime\,T} = C\Gamma C^{-1}$ *)
      | [2; 1] -> printf "f_f%s_c(%s,%s,%s)" f c wf2 wf1 (*
        $\Gamma^{\prime\,T}_{\alpha\beta} \bar\psi_{1,\alpha} \phi_2 \cong
         \Gamma^{\prime\,T} = C\Gamma C^{-1} $ *)
      | _ -> ()

    (* where we have used
       \begin{equation}
         (C\Gamma')^T = \Gamma^{\prime,T}C^T
           = C\Gamma C^{-1} C^T = C\Gamma C^{-1} (-C) = - C\Gamma\,.
       \end{equation} *)

    (* This puts the arguments in the same slots as [tho_print_dirac_current]
       above and can be implemented by [fuse], iff we inject the proper
       transformations in [dennerize] below.
       We notice that we do \emph{not} need the conjugated version for
       all combinations, but only for the case of two fermions.
       In the two cases of one column spinor~$\psi$, only the original
       version appears and in the two cases of one row spinor~$\bar\psi$,
       only the conjugated version appears. *)

    (* Before we continue, we must however generalize from the
       assumption~\eqref{eq:FVF-Vertex} that the fields in the
       vertex are always ordered as in~[Coupling.FBF].  First,
       even in this case the slots of the fermions must be exchanged
       to accomodate the cyclic permutations. Therefore we exchange the
       arguments of the [[1; 3]] and [[3; 1]] fusions. *)

    let jrr_print_majorana_FBF f c wf1 wf2 fusion =
      match fusion with (* [fline = (3, 1)] *)
      | [3; 1] -> printf "f%sf_p120_c(%s,%s,%s)" f c wf1 wf2 (*
        $(C\Gamma')^T_{\alpha\beta}
          \psi_{1,\beta} \bar\psi_{2,\alpha}   \cong
         (C\Gamma')^T = - C\Gamma $ *)
      | [1; 3] -> printf "f%sf_p120  (%s,%s,%s)" f c wf2 wf1 (*
        $ (C\Gamma)_{\alpha\beta} \psi_{1,\beta} \bar\psi_{2,\alpha} \cong
          C\Gamma $ *)
      | [2; 3] -> printf "f%sf_p012  (%s,%s,%s)" f c wf1 wf2 (*
        $ \Gamma_{\alpha\beta} \phi_1 \psi_{2,\beta} \cong
          \Gamma $ *)
      | [3; 2] -> printf "f%sf_p012  (%s,%s,%s)" f c wf2 wf1 (*
        $ \Gamma_{\alpha\beta} \phi_1 \psi_{2,\beta} \cong
          \Gamma $ *)
      | [1; 2] -> printf "f%sf_p201  (%s,%s,%s)" f c wf1 wf2 (*
        $\Gamma^{\prime\,T}_{\alpha\beta} \bar\psi_{1,\alpha} \phi_2 \cong
         \Gamma^{\prime\,T} = C\Gamma C^{-1}$ *)
      | [2; 1] -> printf "f%sf_p201  (%s,%s,%s)" f c wf2 wf1 (*
        $\Gamma^{\prime\,T}_{\alpha\beta} \bar\psi_{1,\alpha} \phi_2 \cong
         \Gamma^{\prime\,T} = C\Gamma C^{-1} $ *)
      | _ -> ()

    (* The other two permutations: *)

    let jrr_print_majorana_FFB f c wf1 wf2 fusion =
      match fusion with (* [fline = (1, 2)] *)
      | [3; 1] -> printf "ff%s_p120  (%s,%s,%s)" f c wf1 wf2 (*
        $ \Gamma_{\alpha\beta} \phi_1 \psi_{2,\beta} \cong
          \Gamma $ *)
      | [1; 3] -> printf "ff%s_p120  (%s,%s,%s)" f c wf2 wf1 (*
        $ \Gamma_{\alpha\beta} \phi_1 \psi_{2,\beta} \cong
          \Gamma $ *)
      | [2; 3] -> printf "ff%s_p012  (%s,%s,%s)" f c wf1 wf2 (*
        $\Gamma^{\prime\,T}_{\alpha\beta} \bar\psi_{1,\alpha} \phi_2 \cong
         \Gamma^{\prime\,T} = C\Gamma C^{-1}$ *)
      | [3; 2] -> printf "ff%s_p012  (%s,%s,%s)" f c wf2 wf1 (*
        $\Gamma^{\prime\,T}_{\alpha\beta} \bar\psi_{1,\alpha} \phi_2 \cong
         \Gamma^{\prime\,T} = C\Gamma C^{-1} $ *)
      | [1; 2] -> printf "ff%s_p201  (%s,%s,%s)" f c wf1 wf2 (*
        $ (C\Gamma)_{\alpha\beta} \psi_{1,\beta} \bar\psi_{2,\alpha} \cong
          C\Gamma $ *)
      | [2; 1] -> printf "ff%s_p201_c(%s,%s,%s)" f c wf2 wf1 (*
        $(C\Gamma')^T_{\alpha\beta}
           \psi_{1,\beta} \bar\psi_{2,\alpha} \cong
         (C\Gamma')^T = - C\Gamma $ *)
      | _ -> ()

    let jrr_print_majorana_BFF f c wf1 wf2 fusion =
      match fusion with (* [fline = (2, 3)] *)
      | [3; 1] -> printf "%sff_p120  (%s,%s,%s)" f c wf1 wf2 (*
        $\Gamma^{\prime\,T}_{\alpha\beta} \bar\psi_{1,\alpha} \phi_2 \cong
         \Gamma^{\prime\,T} = C\Gamma C^{-1} $ *)
      | [1; 3] -> printf "%sff_p120  (%s,%s,%s)" f c wf2 wf1 (*
        $\Gamma^{\prime\,T}_{\alpha\beta} \bar\psi_{1,\alpha} \phi_2 \cong
         \Gamma^{\prime\,T} = C\Gamma C^{-1}$ *)
      | [2; 3] -> printf "%sff_p012  (%s,%s,%s)" f c wf1 wf2 (*
        $ (C\Gamma)_{\alpha\beta} \psi_{1,\beta} \bar\psi_{2,\alpha} \cong
          C\Gamma $ *)
      | [3; 2] -> printf "%sff_p012_c(%s,%s,%s)" f c wf2 wf1 (*
        $(C\Gamma')^T_{\alpha\beta} \psi_{1,\beta} \bar\psi_{2,\alpha} \cong
         (C\Gamma')^T = - C\Gamma $ *)
      | [1; 2] -> printf "%sff_p201  (%s,%s,%s)" f c wf1 wf2 (*
        $ \Gamma_{\alpha\beta} \phi_1 \psi_{2,\beta} \cong
          \Gamma $ *)
      | [2; 1] -> printf "%sff_p201  (%s,%s,%s)" f c wf2 wf1 (*
        $ \Gamma_{\alpha\beta} \phi_1 \psi_{2,\beta} \cong
          \Gamma $ *)
      | _ -> ()

    (* In the model, the necessary
       information is provided as [Coupling.fermion_lines], encoded as
       [(right,left)] in the usual direction of the lines.
       E.\,g.~the case of~\eqref{eq:FVF-Vertex} is~[(3,1)].
       Equivalent information is available
       as~[(ket, bra)] in [UFO_Lorentz.dirac_string]. *)

    let is_majorana = function
      | Coupling.Majorana | Coupling.Vectorspinor | Coupling.Maj_Ghost -> true
      | _ -> false

    let is_dirac = function
      | Coupling.Spinor | Coupling.ConjSpinor -> true
      | _ -> false

    let dennerize ~eval wfs atom =
      let printf fmt = fprintf eval fmt
      and nl = pp_newline eval in
      if is_majorana wfs.(pred atom.L.bra).spin ||
           is_majorana wfs.(pred atom.L.ket).spin then
        if atom.L.bra = 1 then
          (* Fusing one or more bosons with a ket like fermion:
             $\chi \leftarrow \Gamma\chi$. *)
          (* Don't do anything,
             as per subsection~\ref{sec:dirac-matrices-jrr}. *)
          atom
        else if atom.L.ket = 1 then
          (* We fuse one or more bosons with a bra like fermion:
             $\bar\chi \leftarrow \bar\chi\Gamma$. *)
          (* $\Gamma\to C \Gamma C^{-1}$. *)
          begin
            let atom = L.conjugate atom in
            printf "    ! conjugated for Majorana"; nl ();
            printf "    ! %s" (L.dirac_string_to_string atom); nl ();
            atom
          end
        else if not atom.L.conjugated then
          (* We fuse zero or more bosons with a sandwich of fermions.
             $\phi \leftarrow \bar\chi\gamma\chi$.*)
          (* Multiply by~$C$ from the left,
             as per subsection~\ref{sec:dirac-matrices-jrr}. *)
          begin
            let atom = L.cc_times atom in
            printf "    ! multiplied by CC for Majorana"; nl ();
            printf "    ! %s" (L.dirac_string_to_string atom); nl ();
            atom
          end
        else
          (* Transposed: multiply by~$-C$ from the left. *)
          begin
            let atom = L.minus (L.cc_times atom) in
            printf "    ! multiplied by -CC for Majorana"; nl ();
            printf "    ! %s" (L.dirac_string_to_string atom); nl ();
            atom
          end
      else
        atom

    (* Write the [i]th Dirac string [ds] as Fortran code to [eval], including
       a shorthand representation as a comment.  Return [ds] with
       [ds.L.atom] replaced by the dirac string variable,
       i,\,e.~[DS dsv] annotated with the internal and external indices.
       In addition write the declaration to [decl].  *)
    let dirac_string_to_fortran ~decl ~eval i wfs ds =
      let printf fmt = fprintf eval fmt
      and nl = pp_newline eval in
      let bra = ds.L.atom.L.bra
      and ket = ds.L.atom.L.ket in
      pp_divide ~indent:4 eval ();
      printf "    ! %s" (L.dirac_string_to_string ds.L.atom); nl ();
      let atom = dennerize ~eval wfs ds.L.atom in
      begin match ds.L.indices with
      | [] ->
         let gamma = L.dirac_string_to_matrix (fun _ -> 0) atom in
         dirac_bra_or_ket_to_fortran_decl decl i [] bra ket;
         let dsv =
           dirac_bra_or_ket_to_fortran_eval eval i [] wfs bra gamma ket in
         L.map_atom (fun _ -> DS dsv) ds
      | indices ->
         dirac_bra_or_ket_to_fortran_decl decl i indices bra ket;
         let combinations = Product.power (List.length indices) [0; 1; 2; 3] in
         let dsv =
           List.map
             (fun combination ->
               let substitution = IntPM.of_lists indices combination in
               let substitute = IntPM.apply substitution in
               let indices = List.map substitute indices in
               let gamma = L.dirac_string_to_matrix substitute atom in
               dirac_bra_or_ket_to_fortran_eval eval i indices wfs bra gamma ket)
             combinations in
         begin match ThoList.uniq (List.sort compare dsv) with
         | [dsv] -> L.map_atom (fun _ -> DS dsv) ds
         | _ -> failwith "dirac_string_to_fortran: impossible"
         end
      end

    (* Write the Dirac strings in the list [ds_list] as Fortran code to
       [eval], including shorthand representations as comments.
       Return the list of variables and corresponding indices to
       be contracted. *)
    let dirac_strings_to_fortran ~decl ~eval wfs last ds_list =
      List.fold_left
        (fun (i, acc) ds ->
          let i = succ i in
          (i, dirac_string_to_fortran ~decl ~eval i wfs ds :: acc))
        (last, []) ds_list

    (* Perform a nested sum of terms, as printed by [print_term]
       (which takes the number of spaces to indent as only argument)
       of the cartesian product of [indices] running from 0 to 3. *)
    let nested_sums ~decl ~eval initial_indent indices print_term =
      let rec nested_sums' indent = function
        | [] -> print_term indent
        | index :: indices ->
           let var = index_variable index in
           fprintf eval "%*s@[<2>do %s = 0, 3@]" indent "" var;
           pp_newline eval ();
           nested_sums' (indent + 2) indices; pp_newline eval ();
           fprintf eval "%*s@[<2>end do@]" indent "" in
      nested_sums' (initial_indent + 2) indices

    (* Polarization indices also need to be summed over, but they
       appear only once. *)
    let indices_of_contractions contractions =
      let index_pairs, polarizations =
        L.classify_indices
          (ThoList.flatmap (fun ds -> ds.L.indices) contractions) in
      try
        ThoList.pairs index_pairs @ ThoList.uniq (List.sort compare polarizations)
      with
      | Invalid_argument s ->
         invalid_arg
           ("indices_of_contractions: " ^
              ThoList.to_string string_of_int index_pairs)

(*i   Printf.eprintf
        "indices_of_contractions: %s / %s\n"
        (ThoList.to_string string_of_int index_pairs)
        (ThoList.to_string string_of_int polarizations);
i*)

    let format_dsv dsv indices =
      match dsv, indices with
      | Braket _, [] -> dsv_name dsv
      | Braket _, ilist ->
         Printf.sprintf "%s(%s)" (dsv_name dsv) (format_indices indices)
      | (Bra _ | Ket _), [] ->
         Printf.sprintf "%s(%s)" (dsv_name dsv) index_spinor
      | (Bra _ | Ket _), ilist ->
         Printf.sprintf
           "%s(%s,%s)" (dsv_name dsv) index_spinor (format_indices indices)

    let denominator_name = "denom_"
    let mass_name = "m_"
    let width_name = "w_"

    let format_tensor t =
      let indices = t.L.indices in
      match t.L.atom with
      | DS dsv -> format_dsv dsv indices
      | V vector -> Printf.sprintf "%s(%s)" vector (format_indices indices)
      | T UFOx.Lorentz_Atom.P (mu, n) ->
         Printf.sprintf "p%d(%s)" n (index_variable mu)
      | T UFOx.Lorentz_Atom.Epsilon (mu1, mu2, mu3, mu4) ->
         Printf.sprintf "eps4_(%s)" (format_indices [mu1; mu2; mu3; mu4])
      | T UFOx.Lorentz_Atom.Metric (mu1, mu2) ->
         if mu1 > 0 && mu2 > 0 then
           Printf.sprintf "g44_(%s)" (format_indices [mu1; mu2])
         else
           failwith "format_tensor: compress_metrics has failed!"
      | S (UFOx.Lorentz_Atom.Mass _) -> mass_name
      | S (UFOx.Lorentz_Atom.Width _) -> width_name
      | S (UFOx.Lorentz_Atom.P2 i) -> Printf.sprintf "g2_(p%d)" i
      | S (UFOx.Lorentz_Atom.P12 (i, j)) -> Printf.sprintf "g12_(p%d,p%d)" i j
      | Inv (UFOx.Lorentz_Atom.Mass _) -> "1/" ^ mass_name
      | Inv (UFOx.Lorentz_Atom.Width _) -> "1/" ^ width_name
      | Inv (UFOx.Lorentz_Atom.P2 i) -> Printf.sprintf "1/g2_(p%d)" i
      | Inv (UFOx.Lorentz_Atom.P12 (i, j)) ->
         Printf.sprintf "1/g12_(p%d,p%d)" i j
      | S (UFOx.Lorentz_Atom.Variable s) -> s
      | Inv (UFOx.Lorentz_Atom.Variable s) -> "1/" ^ s
      | S (UFOx.Lorentz_Atom.Coeff c) -> UFOx.Value.to_string c
      | Inv (UFOx.Lorentz_Atom.Coeff c) -> "1/(" ^ UFOx.Value.to_string c ^ ")"

    let rec multiply_tensors ~decl ~eval = function
      | [] -> fprintf eval "1";
      | [t] -> fprintf eval "%s" (format_tensor t)
      | t :: tensors ->
         fprintf eval "%s@,*" (format_tensor t);
         multiply_tensors ~decl ~eval tensors

    let pseudo_wfs_for_denominator =
      Array.init
        2
        (fun i ->
          let ii = string_of_int i in
          { pos = i;
            spin = Coupling.Scalar;
            name = denominator_name;
            local_array = None;
            momentum = "k" ^ ii;
            momentum_array = "p" ^ ii;
            fortran_type = fortran_type Coupling.Scalar })

    let contract_indices ~decl ~eval indent wf_indices wfs (fusion, contractees) =
      let printf fmt = fprintf eval fmt
      and nl = pp_newline eval in
      let sum_var =
        begin match wf_indices with
        | [] -> wfs.(0).name
        | ilist ->
           let indices = String.concat "," ilist in
           begin match wfs.(0).local_array with
           | None ->
              let component =
                begin match wfs.(0).spin with
                | Coupling.Spinor | Coupling.ConjSpinor | Coupling.Majorana -> "a"
                | Coupling.Tensor_2 -> "t"
                | Coupling.Vector | Coupling.Massive_Vector ->
                   failwith "contract_indices: expected local_array for vectors"
                | _ -> failwith "contract_indices: unexpected spin"
                end in
              Printf.sprintf "%s%%%s(%s)" wfs.(0).name component indices
           | Some a -> Printf.sprintf "%s(%s)" a indices
           end
        end in
      let indices =
        List.filter
          (fun i -> UFOx.Index.position i <> 1)
          (indices_of_contractions contractees) in
      nested_sums
        ~decl ~eval
        indent indices
        (fun indent ->
          printf "%*s@[<2>%s = %s" indent "" sum_var sum_var;
          printf "@ %s" (format_complex_rational_factor fusion.L.coeff);
          List.iter (fun i -> printf "@,g4_(%s)*" (index_variable i)) indices;
          printf "@,(";
          multiply_tensors ~decl ~eval contractees;
          printf ")";
          begin match fusion.L.denominator with
          | [] -> ()
          | d -> printf " / %s" denominator_name
          end;
          printf "@]");
      printf "@]";
      nl ()

    let scalar_expression1 ~decl ~eval fusion =
      let printf fmt = fprintf eval fmt in
      match fusion.L.dirac, fusion.L.vector with
      | [], [] ->
         let scalars =
           List.map (fun t -> { L.atom = S t; L.indices = [] }) fusion.L.scalar
         and inverses =
           List.map (fun t -> { L.atom = Inv t; L.indices = [] }) fusion.L.inverse in
         let contractees = scalars @ inverses in
         printf "@ %s" (format_complex_rational_factor fusion.L.coeff);
         multiply_tensors ~decl ~eval contractees
      | _, [] ->
         invalid_arg
           "UFO_targets.Fortran.scalar_expression1: unexpected spinor indices"
      | [], _ ->
         invalid_arg
           "UFO_targets.Fortran.scalar_expression1: unexpected vector indices"
      | _, _ ->
         invalid_arg
           "UFO_targets.Fortran.scalar_expression1: unexpected indices"

    let scalar_expression ~decl ~eval indent name fusions =
      let printf fmt = fprintf eval fmt
      and nl = pp_newline eval in
      let sum_var = name in
      printf "%*s@[<2>%s =" indent "" sum_var;
      List.iter (scalar_expression1 ~decl ~eval) fusions;
      printf "@]";
      nl ()

    let local_vector_copies ~decl ~eval wfs =
      begin match wfs.(0).local_array with
      | None -> ()
      | Some a ->
         fprintf
           decl "    @[<2>complex(kind=default),@ dimension(0:3) ::@ %s@]" a;
         pp_newline decl ()
      end;
      let n = Array.length wfs in
      for i = 1 to n - 1 do
        match wfs.(i).local_array with
        | None -> ()
        | Some a ->
           fprintf
             decl "    @[<2>complex(kind=default),@ dimension(0:3) ::@ %s@]" a;
           pp_newline decl ();
           fprintf eval "    @[<2>%s(0) = %s%%t@]" a wfs.(i).name;
           pp_newline eval ();
           fprintf eval "    @[<2>%s(1:3) = %s%%x@]" a wfs.(i).name;
           pp_newline eval ()
      done

    let return_vector ff wfs =
      let printf fmt = fprintf ff fmt
      and nl = pp_newline ff in
      match wfs.(0).local_array with
      | None -> ()
      | Some a ->
         pp_divide ~indent:4 ff ();
         printf "    @[<2>%s%%t = %s(0)@]" wfs.(0).name a; nl ();
         printf "    @[<2>%s%%x = %s(1:3)@]" wfs.(0).name a; nl ()

    let multiply_coupling_and_scalars ff g_opt wfs =
      let printf fmt = fprintf ff fmt
      and nl = pp_newline ff in
      pp_divide ~indent:4 ff ();
      let g =
        match g_opt with
        | None -> ""
        | Some g -> g ^ "*" in
      let wfs0name =
        match wfs.(0).local_array with
        | None -> wfs.(0).name
        | Some a -> a in
      printf "    @[<2>%s = %s%s" wfs0name g wfs0name;
      for i = 1 to Array.length wfs - 1 do
        match wfs.(i).spin with
        | Coupling.Scalar -> printf "@,*%s" wfs.(i).name
        | _ -> ()
      done;
      printf "@]"; nl ()

    let local_momentum_copies ~decl ~eval wfs =
      let n = Array.length wfs in
      fprintf
        decl "    @[<2>real(kind=default),@ dimension(0:3) ::@ %s"
        wfs.(0).momentum_array;
      for i = 1 to n - 1 do
        fprintf decl ",@ %s" wfs.(i).momentum_array;
        fprintf
          eval "    @[<2>%s(0) = %s%%t@]"
          wfs.(i).momentum_array wfs.(i).momentum;
        pp_newline eval ();
        fprintf
          eval "    @[<2>%s(1:3) = %s%%x@]"
          wfs.(i).momentum_array wfs.(i).momentum;
        pp_newline eval ()
      done;
      fprintf eval "    @[<2>%s =" wfs.(0).momentum_array;
      for i = 1 to n - 1 do
        fprintf eval "@ - %s" wfs.(i).momentum_array
      done;
      fprintf decl "@]";
      pp_newline decl ();
      fprintf eval "@]";
      pp_newline eval ()

    let contractees_of_fusion
          ~decl ~eval wfs (max_dsv, indices_seen, contractees) fusion =
      let max_dsv', dirac_strings =
        dirac_strings_to_fortran ~decl ~eval wfs max_dsv fusion.L.dirac
      and vectors =
        List.fold_left
          (fun acc wf ->
            match wf.spin, wf.local_array with
            | Coupling.Tensor_2, None ->
               { L.atom =
                   V (Printf.sprintf "%s%d%%t" (spin_mnemonic wf.spin) wf.pos);
                 L.indices = [UFOx.Index.pack wf.pos 1;
                              UFOx.Index.pack wf.pos 2] } :: acc
            | _, None -> acc
            | _, Some a -> { L.atom = V a; L.indices = [wf.pos] } :: acc)
          [] (List.tl (Array.to_list wfs))
      and tensors =
        List.map (L.map_atom (fun t -> T t)) fusion.L.vector
      and scalars =
        List.map (fun t -> { L.atom = S t; L.indices = [] }) fusion.L.scalar
      and inverses =
        List.map (fun t -> { L.atom = Inv t; L.indices = [] }) fusion.L.inverse in
      let contractees' = dirac_strings @ vectors @ tensors @ scalars @ inverses in
      let indices_seen' =
        Sets.Int.of_list (indices_of_contractions contractees') in
      (max_dsv',
       Sets.Int.union indices_seen indices_seen',
       (fusion, contractees') :: contractees)

    let local_name wf =
      match wf.local_array with
      | Some a -> a
      | None ->
         match wf.spin with
         | Coupling.Spinor | Coupling.ConjSpinor | Coupling.Majorana ->
            wf.name ^ "%a"
         | Coupling.Scalar -> wf.name
         | Coupling.Tensor_2 -> wf.name ^ "%t"
         | Coupling.Vector | Coupling.Massive_Vector ->
            failwith "UFO_targets.Fortran.local_name: unexpected spin 1"
         | _ ->
            failwith "UFO_targets.Fortran.local_name: unhandled spin"

    let external_wf_loop ~decl ~eval ~indent wfs (fusion, _ as contractees) =
      pp_divide ~indent eval ();
      fprintf eval "%*s! %s" indent "" (L.to_string [fusion]); pp_newline eval ();
      pp_divide ~indent eval ();
      begin match fusion.L.denominator with
      | [] -> ()
      | denominator ->
         scalar_expression ~decl ~eval 4 denominator_name denominator
      end;
      match wfs.(0).spin with
      | Coupling.Scalar ->
         contract_indices ~decl ~eval 2 [] wfs contractees
      | Coupling.Spinor | Coupling.ConjSpinor | Coupling.Majorana ->
         let idx = index_spinor in
         fprintf eval "%*s@[<2>do %s = 1, 4@]" indent "" idx; pp_newline eval ();
         contract_indices ~decl ~eval 4 [idx] wfs contractees;
         fprintf eval "%*send do@]" indent ""; pp_newline eval ()
      | Coupling.Vector | Coupling.Massive_Vector ->
         let idx = index_variable 1 in
         fprintf eval "%*s@[<2>do %s = 0, 3@]" indent "" idx; pp_newline eval ();
         contract_indices ~decl ~eval 4 [idx] wfs contractees;
         fprintf eval "%*send do@]" indent ""; pp_newline eval ()
      | Coupling.Tensor_2 ->
         let idx1 = index_variable (UFOx.Index.pack 1 1)
         and idx2 = index_variable (UFOx.Index.pack 1 2) in
         fprintf eval "%*s@[<2>do %s = 0, 3@]" indent "" idx1;
         pp_newline eval ();
         fprintf eval "%*s@[<2>do %s = 0, 3@]" (indent + 2) "" idx2;
         pp_newline eval ();
         contract_indices ~decl ~eval 6 [idx1; idx2] wfs contractees;
         fprintf eval "%*send do@]" (indent + 2) ""; pp_newline eval ();
         fprintf eval "%*send do@]" indent ""; pp_newline eval ()
      | Coupling.Vectorspinor ->
         failwith "external_wf_loop: Vectorspinor not supported yet!"
      | Coupling.Maj_Ghost ->
         failwith "external_wf_loop: unexpected Maj_Ghost"
      | Coupling.Tensor_1 ->
         failwith "external_wf_loop: unexpected Tensor_1"
      | Coupling.BRS _ ->
         failwith "external_wf_loop: unexpected BRS"

    let fusions_to_fortran ~decl ~eval wfs ?(denominator=[]) ?coupling fusions =
      local_vector_copies ~decl ~eval wfs;
      local_momentum_copies ~decl ~eval wfs;
      begin match denominator with
      | [] -> ()
      | _ ->
         fprintf decl "    @[<2>complex(kind=default) :: %s@]" denominator_name;
         pp_newline decl ()
      end;
      let max_dsv, indices_used, contractions =
        List.fold_left
          (contractees_of_fusion ~decl ~eval wfs)
          (0, Sets.Int.empty, [])
          fusions in
      Sets.Int.iter
        (fun index ->
          fprintf decl "    @[<2>integer ::@ %s@]" (index_variable index);
          pp_newline decl ())
        indices_used;
      begin match wfs.(0).spin with
      | Coupling.Spinor | Coupling.ConjSpinor | Coupling.Majorana ->
         fprintf decl "    @[<2>integer ::@ %s@]" index_spinor;
         pp_newline decl ()
      | _ -> ()
      end;
      pp_divide ~indent:4 eval ();
      let wfs0name = local_name wfs.(0) in
      fprintf eval "    %s = 0" wfs0name;
      pp_newline eval ();
      List.iter (external_wf_loop ~decl ~eval ~indent:4 wfs) contractions;
      multiply_coupling_and_scalars eval coupling wfs;
      begin match denominator with
      | [] -> ()
      | denominator ->
         pp_divide ~indent:4 eval ();
         fprintf eval "%*s! %s" 4 "" (L.to_string denominator);
         pp_newline eval ();
         scalar_expression ~decl ~eval 4 denominator_name denominator;
         fprintf eval
           "    @[<2>%s =@ %s / %s@]" wfs0name wfs0name denominator_name;
         pp_newline eval ()
      end;
      return_vector eval wfs

    (* TODO: eventually, we should include the momentum among
       the arguments only if required.  But this can wait for
       another day. *)
    let lorentz ff name spins lorentz =
      let printf fmt = fprintf ff fmt
      and nl = pp_newline ff in
      let wfs = wf_table spins in
      let n = Array.length wfs in
      printf "  @[<4>pure function %s@ (g,@ " name;
      for i = 1 to n - 2 do
        printf "%s,@ %s,@ " wfs.(i).name wfs.(i).momentum
      done;
      printf "%s,@ %s" wfs.(n - 1).name wfs.(n - 1).momentum;
      printf ")@ result (%s)@]" wfs.(0).name; nl ();
      printf "    @[<2>%s ::@ %s@]" wfs.(0).fortran_type wfs.(0).name; nl();
      printf "    @[<2>complex(kind=default),@ intent(in) ::@ g@]"; nl();
      for i = 1 to n - 1 do
        printf
          "    @[<2>%s, intent(in) :: %s@]"
          wfs.(i).fortran_type wfs.(i).name; nl();
      done;
      printf "    @[<2>type(momentum), intent(in) ::@ %s" wfs.(1).momentum;
      for i = 2 to n - 1 do
        printf ",@ %s" wfs.(i).momentum
      done;
      printf "@]";
      nl ();
      let width = 80 in (* get this from the default formatter instead! *)
      let decl_buf = Buffer.create 1024
      and eval_buf = Buffer.create 1024 in
      let decl = formatter_of_buffer ~width decl_buf
      and eval = formatter_of_buffer ~width eval_buf in
      fusions_to_fortran ~decl ~eval ~coupling:"g" wfs lorentz;
      pp_flush decl ();
      pp_flush eval ();
      pp_divide ~indent:4 ff ();
(*i   printf "    ! %s" (L.to_string lorentz); nl ();
      pp_divide ~indent:4 ff (); i*)
      printf "%s" (Buffer.contents decl_buf);
      pp_divide ~indent:4 ff ();
      printf "    if (g == 0) then"; nl ();
      printf "      call set_zero (%s)" wfs.(0).name; nl ();
      printf "      return"; nl ();
      printf "    end if"; nl ();
      pp_divide ~indent:4 ff ();
      printf "%s" (Buffer.contents eval_buf);
      printf "  end function %s@]" name; nl ();
      Buffer.reset decl_buf;
      Buffer.reset eval_buf;
      ()

    let use_variables ff parameter_module variables =
      let printf fmt = fprintf ff fmt
      and nl = pp_newline ff in
      match variables with
      | [] -> ()
      | v :: v_list ->
         printf "    @[<2>use %s, only: %s" parameter_module v;
         List.iter (fun s -> printf ", %s" s) v_list;
         printf "@]"; nl ()

    let propagator ff name parameter_module variables
          (bra_spin, ket_spin) numerator denominator =
      let printf fmt = fprintf ff fmt
      and nl = pp_newline ff in
      let width = 80 in (* get this from the default formatter instead! *)
      let wf_name = spin_mnemonic ket_spin
      and wf_type = fortran_type ket_spin in
      let wfs = wf_table [| ket_spin; ket_spin |] in
      printf
        "  @[<4>pure function pr_U_%s@ (k2, %s, %s, %s2)"
        name mass_name width_name wf_name;
      printf " result (%s1)@]" wf_name; nl ();
      use_variables ff parameter_module variables;
      printf "    %s :: %s1" wf_type wf_name; nl ();
      printf "    type(momentum), intent(in) :: k2"; nl ();
      printf
        "    real(kind=default), intent(in) :: %s, %s"
        mass_name width_name; nl ();
      printf "    %s, intent(in) :: %s2" wf_type wf_name; nl ();
      let decl_buf = Buffer.create 1024
      and eval_buf = Buffer.create 1024 in
      let decl = formatter_of_buffer ~width decl_buf
      and eval = formatter_of_buffer ~width eval_buf in
      fusions_to_fortran ~decl ~eval wfs ~denominator numerator;
      pp_flush decl ();
      pp_flush eval ();
      pp_divide ~indent:4 ff ();
      printf "%s" (Buffer.contents decl_buf);
      pp_divide ~indent:4 ff ();
      printf "%s" (Buffer.contents eval_buf);
      printf "  end function pr_U_%s@]" name; nl ();
      Buffer.reset decl_buf;
      Buffer.reset eval_buf;
      ()

    let scale_coupling c g =
      if c = 1 then
        g
      else if c = -1 then
        "-" ^ g
      else
        Printf.sprintf "%d*%s" c g

    let scale_coupling z g =
      format_complex_rational_factor z ^ g

    (* As a prototypical example consider the vertex
       \begin{subequations}
       \label{eq:cyclic-UFO-fusions}
       \begin{equation}
         \bar\psi\fmslash{A}\psi =
            \tr\left(\psi\otimes\bar\psi\fmslash{A}\right)
       \end{equation}
       encoded as \texttt{FFV} in the SM UFO file.  This example
       is useful, because all three fields have different type
       and we can use the Fortran compiler to check our
       implementation.

       In this case we need to generate the following function
       calls with the arguments in the following order
       \begin{center}
         \begin{tabular}{lcl}
           \texttt{F12}:&$\psi_1\bar\psi_2\to A$&
              \texttt{FFV\_p201(g,psi1,p1,psibar2,p2)} \\
           \texttt{F21}:&$\bar\psi_1\psi_2\to A$&
              \texttt{FFV\_p201(g,psi2,p2,psibar1,p1)} \\
           \texttt{F23}:&$\bar\psi_1 A_2 \to \bar\psi$&
              \texttt{FFV\_p012(g,psibar1,p1,A2,p2)} \\
           \texttt{F32}:&$A_1\bar\psi_2 \to \bar\psi$&
              \texttt{FFV\_p012(g,psibar2,p2,A1,p1)} \\
           \texttt{F31}:&$A_1\psi_2\to \psi$&
              \texttt{FFV\_p120(g,A1,p1,psi2,p2)} \\
           \texttt{F13}:&$\psi_1A_2\to \psi$&
              \texttt{FFV\_p120(g,A2,p2,psi1,p1)}
         \end{tabular}
       \end{center} *)

    (* Fortunately, all Fermi signs have been taken
       care of by [Fusions] and we can concentrate on
       injecting the wave functions into the correct slots. *)

    (* The other possible cases are
       \begin{equation}
         \bar\psi\fmslash{A}\psi
       \end{equation}
       which would be encoded as \texttt{FVF} in a UFO file
       \begin{center}
         \begin{tabular}{lcl}
           \texttt{F12}:&$\bar\psi_1 A_2 \to \bar\psi$&
              \texttt{FVF\_p201(g,psibar1,p1,A2,p2)} \\
           \texttt{F21}:&$A_1\bar\psi_2 \to \bar\psi$&
              \texttt{FVF\_p201(g,psibar2,p2,A1,p1)} \\
           \texttt{F23}:&$A_1\psi_2\to \psi$&
              \texttt{FVF\_p012(g,A1,p1,psi2,p2)} \\
           \texttt{F32}:&$\psi_1A_2\to \psi$&
              \texttt{FVF\_p012(g,A2,p2,psi1,p1)} \\
           \texttt{F31}:&$\psi_1\bar\psi_2\to A$&
              \texttt{FVF\_p120(g,psi1,p1,psibar2,p2)} \\
           \texttt{F13}:&$\bar\psi_1\psi_2\to A$&
              \texttt{FVF\_p120(g,psi2,p2,psibar1,p1)}
         \end{tabular}
       \end{center}
       and
       \begin{equation}
         \bar\psi\fmslash{A}\psi =
            \tr\left(\fmslash{A}\psi\otimes\bar\psi\right)\,,
       \end{equation}
       corresponding to \texttt{VFF}
       \begin{center}
         \begin{tabular}{lcl}
           \texttt{F12}:&$A_1\psi_2\to \psi$&
              \texttt{VFF\_p201(g,A1,p1,psi2,p2)} \\
           \texttt{F21}:&$\psi_1A_2\to \psi$&
              \texttt{VFF\_p201(g,A2,p2,psi1,p1)} \\
           \texttt{F23}:&$\psi_1\bar\psi_2\to A$&
              \texttt{VFF\_p012(g,psi1,p1,psibar2,p2)} \\
           \texttt{F32}:&$\bar\psi_1\psi_2\to A$&
              \texttt{VFF\_p012(g,psi2,p2,psibar1,p1)} \\
           \texttt{F31}:&$\bar\psi_1 A_2 \to \bar\psi$&
              \texttt{VFF\_p120(g,psibar1,p1,A2,p2)} \\
           \texttt{F13}:&$A_1\bar\psi_2 \to \bar\psi$&
              \texttt{VFF\_p120(g,psibar2,p2,A1,p1)}
         \end{tabular}
       \end{center}
       \end{subequations} *)

    (* \begin{dubious}
         Once the Majorana code generation is fully debugged,
         we should replace the lists by reverted lists everywhere
         in order to become a bit more efficient.
       \end{dubious} *)

    module P = Permutation.Default

    let factor_cyclic f12__n =
      let f12__, fn = ThoList.split_last f12__n in
      let cyclic = ThoList.cycle_until fn (List.sort compare f12__n) in
      (P.of_list (List.map pred cyclic),
       P.of_lists (List.tl cyclic) f12__)

    let ccs_to_string ccs =
      String.concat "" (List.map (fun (f, i) -> Printf.sprintf "_c%x%x" i f) ccs)

    let fusion_name v perm ccs =
      Printf.sprintf "%s_p%s%s" v (P.to_string perm) (ccs_to_string ccs)

    let fuse_dirac c v s fl g wfs ps fusion =
      let g = scale_coupling c g
      and cyclic, factor = factor_cyclic fusion in
      let wfs_ps = List.map2 (fun wf p -> (wf, p)) wfs ps in
      let args = P.list (P.inverse factor) wfs_ps in
      let args_string =
        String.concat "," (List.map (fun (wf, p) -> wf ^ "," ^ p) args) in
      printf "%s(%s,%s)" (fusion_name v cyclic []) g args_string

    (* We need to look at the permuted fermion lines in order to
       decide wether to apply charge conjugations.  *)

    (* It is not enough to look at the cyclic permutation used
       to move the fields into the correct arguments of
       the fusions \ldots *)
    let map_indices perm unit =
      let pmap = IntPM.of_lists unit (P.list perm unit) in
      IntPM.apply pmap

    (* \ldots{} we also need to inspect the full permutation of
       the fields. *)
    let map_indices2 perm unit =
      let pmap =
        IntPM.of_lists unit (1 :: P.list (P.inverse perm) (List.tl unit)) in
      IntPM.apply pmap

    (* This is a more direct implementation of the composition
       of [map_indices2] and [map_indices], that is used in the
       unit tests. *)
    let map_indices_raw fusion =
      let unit = ThoList.range 1 (List.length fusion) in
      let f12__, fn = ThoList.split_last fusion in
      let fusion = fn :: f12__ in
      let map_index = IntPM.of_lists fusion unit in
      IntPM.apply map_index

    (* Map the fermion line indices in [fl] according to [map_index]. *)
    let map_fermion_lines map_index fl =
      List.map (fun (i, f) -> (map_index i, map_index f)) fl

    (* Map the fermion line indices in [fl] according to [map_index],
       but keep a copy of the original. *)
    let map_fermion_lines2 map_index fl =
      List.map (fun (i, f) -> ((i, f), (map_index i, map_index f))) fl

    let permute_fermion_lines cyclic unit fl =
      map_fermion_lines (map_indices cyclic unit) fl

    let permute_fermion_lines2 cyclic factor unit fl =
      map_fermion_lines2
        (map_indices2 factor unit)
        (map_fermion_lines (map_indices cyclic unit) fl)

    (* \begin{dubious}
         TODO: this needs more more work for the fully
         general case with 4-fermion operators involving Majoranas.
       \end{dubious} *)
    let charge_conjugations fl2 =
      ThoList.filtermap
        (fun ((i, f), (i', f')) ->
          match (i, f), (i', f') with
          | (1, 2), _ | (2, 1), _ -> Some (f, i) (* $\chi^T\Gamma'$ *)
          | _, (2, 3) -> Some (f, i)             (* $\chi^T(C\Gamma')\chi$ *)
          | _ -> None)
        fl2

(*i
    let fuse_majorana c v s fl g wfs ps fusion =
      let g = scale_coupling c g
      and cyclic, factor = factor_cyclic fusion in
      let wfs_ps = List.map2 (fun wf p -> (wf, p)) wfs ps in
      let wfs_ps_string =
        String.concat "," (List.map (fun (wf, p) -> wf ^ "," ^ p) wfs_ps) in
      let args = P.list (P.inverse factor) wfs_ps in
      let args_string =
        String.concat "," (List.map (fun (wf, p) -> wf ^ "," ^ p) args) in
      let f12__, fn = ThoList.split_last fusion in
      Printf.eprintf
        "fusion : %d < %s\n" fn (ThoList.to_string string_of_int f12__);
      Printf.eprintf "cyclic : %s\n" (P.to_string cyclic);
      Printf.eprintf "factor : %s\n" (P.to_string factor);
      let unit = ThoList.range 1 (List.length fusion) in
      Printf.eprintf "permutation     : %s -> %s\n"
        (ThoList.to_string string_of_int unit)
        (ThoList.to_string
           string_of_int (List.map (map_indices cyclic unit) unit));
      Printf.eprintf "permutation raw : %s -> %s\n"
        (ThoList.to_string string_of_int unit)
        (ThoList.to_string
           string_of_int (List.map (map_indices_raw fusion) unit));
      Printf.eprintf "fermion lines : %s\n"
        (ThoList.to_string (fun (i, f) -> Printf.sprintf "%d>%d" i f) fl);
      let fl2 = permute_fermion_lines2 cyclic factor unit fl in
      let fl = permute_fermion_lines cyclic unit fl in
      Printf.eprintf "permuted      : %s\n"
        (ThoList.to_string (fun (i, f) -> Printf.sprintf "%d>%d" i f) fl);
      Printf.eprintf "arguments : %s\n" wfs_ps_string;
      Printf.eprintf "permuted  : %s\n" args_string;
      Printf.eprintf
        ">> %s(%s,%s)\n"
        (fusion_name v cyclic (charge_conjugations fl2)) g args_string;
      printf "%s(%s,%s)" (fusion_name v cyclic (charge_conjugations fl2)) g args_string
i*)

    let charge_conjugations fl2 =
      ThoList.filtermap
        (fun ((i, f), (i', f')) ->
          match (i, f), (i', f') with
          | _, (2, 3) -> Some (f, i)
          | _ -> None)
        fl2

    let fuse_majorana c v s fl g wfs ps fusion =
      let g = scale_coupling c g
      and cyclic, factor = factor_cyclic fusion in
      let wfs_ps = List.map2 (fun wf p -> (wf, p)) wfs ps in
      let args = P.list (P.inverse factor) wfs_ps in
      let args_string =
        String.concat "," (List.map (fun (wf, p) -> wf ^ "," ^ p) args) in
      let unit = ThoList.range 1 (List.length fusion) in
      let ccs =
        charge_conjugations (permute_fermion_lines2 cyclic factor unit fl) in
      printf "%s(%s,%s)" (fusion_name v cyclic ccs) g args_string

    let fuse c v s fl g wfs ps fusion =
      if List.exists is_majorana s then
        fuse_majorana c v s fl g wfs ps fusion
      else
        fuse_dirac c v s fl g wfs ps fusion

    let eps4_g4_g44_decl ff () =
      let printf fmt = fprintf ff fmt
      and nl = pp_newline ff in
      printf "  @[<2>integer,@ dimension(0:3)";
      printf ",@ save,@ private ::@ g4_@]"; nl ();
      printf "  @[<2>integer,@ dimension(0:3,0:3)";
      printf ",@ save,@ private ::@ g44_@]"; nl ();
      printf "  @[<2>integer,@ dimension(0:3,0:3,0:3,0:3)";
      printf ",@ save,@ private ::@ eps4_@]"; nl ()

    let eps4_g4_g44_init ff () =
      let printf fmt = fprintf ff fmt
      and nl = pp_newline ff in
      printf "  @[<2>data g4_@            /@  1, -1, -1, -1 /@]"; nl ();
      printf "  @[<2>data g44_(0,:)@      /@  1,  0,  0,  0 /@]"; nl ();
      printf "  @[<2>data g44_(1,:)@      /@  0, -1,  0,  0 /@]"; nl ();
      printf "  @[<2>data g44_(2,:)@      /@  0,  0, -1,  0 /@]"; nl ();
      printf "  @[<2>data g44_(3,:)@      /@  0,  0,  0, -1 /@]"; nl ();
      for mu1 = 0 to 3 do
        for mu2 = 0 to 3 do
          for mu3 = 0 to 3 do
            printf "  @[<2>data eps4_(%d,%d,%d,:)@ /@ " mu1 mu2 mu3;
            for mu4 = 0 to 3 do
              if mu4 <> 0 then
                printf ",@ ";
              let mus = [mu1; mu2; mu3; mu4] in
              if List.sort compare mus = [0; 1; 2; 3] then
                printf "%2d" (Combinatorics.sign mus)
              else
                printf "%2d" 0;
            done;
            printf " /@]";
            nl ()
          done
        done
      done

    let inner_product_functions ff () =
      let printf fmt = fprintf ff fmt
      and nl = pp_newline ff in
      printf "  pure function g2_ (p) result (p2)"; nl();
      printf "    real(kind=default), dimension(0:3), intent(in) :: p"; nl();
      printf "    real(kind=default) :: p2"; nl();
      printf "    p2 = p(0)*p(0) - p(1)*p(1) - p(2)*p(2) - p(3)*p(3)"; nl();
      printf "  end function g2_"; nl();
      printf "  pure function g12_ (p1, p2) result (p12)"; nl();
      printf "    real(kind=default), dimension(0:3), intent(in) :: p1, p2"; nl();
      printf "    real(kind=default) :: p12"; nl();
      printf "    p12 = p1(0)*p2(0) - p1(1)*p2(1) - p1(2)*p2(2) - p1(3)*p2(3)"; nl();
      printf "  end function g12_"; nl()

    module type Test =
      sig
        val suite : OUnit.test
      end

    module Test : Test =
      struct

        open OUnit

        let assert_mappings fusion =
          let unit = ThoList.range 1 (List.length fusion) in
          let cyclic, factor = factor_cyclic fusion in
          let raw = map_indices_raw fusion
          and map1 = map_indices cyclic unit
          and map2 = map_indices2 factor unit in
          let map i = map2 (map1 i) in
          assert_equal ~printer:(ThoList.to_string string_of_int)
            (List.map raw unit) (List.map map unit)

        let suite_mappings =
          "mappings" >:::

            [ "1<-2" >::
                (fun () ->
                  List.iter assert_mappings (Combinatorics.permute [1;2;3]));

              "1<-3" >::
                (fun () ->
                  List.iter assert_mappings (Combinatorics.permute [1;2;3;4])) ]

        let suite =
          "UFO_targets" >:::
            [suite_mappings]

      end
  end
    
