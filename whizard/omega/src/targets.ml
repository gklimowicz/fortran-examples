(* targets.ml --

   Copyright (C) 1999-2022 by

       Wolfgang Kilian <kilian@physik.uni-siegen.de>
       Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
       Juergen Reuter <juergen.reuter@desy.de>
       with contributions from
       Christian Speckner <cnspeckn@googlemail.com>
       Fabian Bach <fabian.bach@t-online.de> (only parts of this file)
       Marco Sekulla <marco.sekulla@kit.edu> (only parts of this file)
       Bijan Chokoufe Nejad <bijan.chokoufe@desy.de> (only parts of this file)
       So Young Shim <soyoung.shim@desy.de>

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

module Dummy (F : Fusion.Maker) (P : Momentum.T) (M : Model.T) =
  struct
    type amplitudes = Fusion.Multi(F)(P)(M).amplitudes
    type diagnostic = All | Arguments | Momenta | Gauge
    let options = Options.empty
    let amplitudes_to_channel _ _ _ = failwith "Targets.Dummy"
    let parameters_to_channel _ = failwith "Targets.Dummy"
  end

(* \thocwmodulesection{O'Mega Virtual Machine with \texttt{Fortran\;90/95}} *)

(* \thocwmodulesubsection{Preliminaries} *)

module VM (Fusion_Maker : Fusion.Maker) (P : Momentum.T) (M : Model.T) =
  struct

    open Coupling
    open Format

    module CM = Colorize.It(M)
    module F = Fusion_Maker(P)(M)
    module CF = Fusion.Multi(Fusion_Maker)(P)(M)
    module CFlow = Color.Flow
    type amplitudes = CF.amplitudes

(* Options. *)
    type diagnostic = All | Arguments | Momenta | Gauge

    let wrapper_module = ref "ovm_wrapper"
    let parameter_module_external = ref "some_external_module_with_model_info"
    let bytecode_file = ref "bytecode.hbc"
    let md5sum = ref None
    let openmp = ref false
    let kind = ref "default"
    let whizard = ref false

    let options = Options.create
      [ "wrapper_module", Arg.String (fun s -> wrapper_module := s),
          "name of wrapper module";
        "bytecode_file", Arg.String (fun s -> bytecode_file := s),
          "bytecode file to be used in wrapper";
        "parameter_module_external", Arg.String (fun s ->
                                     parameter_module_external := s),
          "external parameter module to be used in wrapper";
        "md5sum", Arg.String (fun s -> md5sum := Some s),
          "transfer MD5 checksum in wrapper";
        "whizard", Arg.Set whizard, "include WHIZARD interface in wrapper";
        "openmp", Arg.Set openmp,
          "activate parallel computation of amplitude with OpenMP"]

(* Integers encode the opcodes (operation codes). *)
    let ovm_ADD_MOMENTA = 1
    let ovm_CALC_BRAKET = 2

    let ovm_LOAD_SCALAR = 10
    let ovm_LOAD_SPINOR_INC = 11
    let ovm_LOAD_SPINOR_OUT = 12
    let ovm_LOAD_CONJSPINOR_INC = 13
    let ovm_LOAD_CONJSPINOR_OUT = 14
    let ovm_LOAD_MAJORANA_INC = 15
    let ovm_LOAD_MAJORANA_OUT = 16
    let ovm_LOAD_VECTOR_INC = 17
    let ovm_LOAD_VECTOR_OUT = 18
    let ovm_LOAD_VECTORSPINOR_INC = 19
    let ovm_LOAD_VECTORSPINOR_OUT = 20
    let ovm_LOAD_TENSOR2_INC = 21
    let ovm_LOAD_TENSOR2_OUT = 22
    let ovm_LOAD_BRS_SCALAR = 30
    let ovm_LOAD_BRS_SPINOR_INC = 31
    let ovm_LOAD_BRS_SPINOR_OUT = 32
    let ovm_LOAD_BRS_CONJSPINOR_INC = 33
    let ovm_LOAD_BRS_CONJSPINOR_OUT = 34
    let ovm_LOAD_BRS_VECTOR_INC = 37
    let ovm_LOAD_BRS_VECTOR_OUT = 38
    let ovm_LOAD_MAJORANA_GHOST_INC = 23
    let ovm_LOAD_MAJORANA_GHOST_OUT = 24
    let ovm_LOAD_BRS_MAJORANA_INC = 35
    let ovm_LOAD_BRS_MAJORANA_OUT = 36

    let ovm_PROPAGATE_SCALAR = 51
    let ovm_PROPAGATE_COL_SCALAR = 52
    let ovm_PROPAGATE_GHOST = 53
    let ovm_PROPAGATE_SPINOR = 54
    let ovm_PROPAGATE_CONJSPINOR = 55
    let ovm_PROPAGATE_MAJORANA = 56
    let ovm_PROPAGATE_COL_MAJORANA = 57
    let ovm_PROPAGATE_UNITARITY = 58
    let ovm_PROPAGATE_COL_UNITARITY = 59
    let ovm_PROPAGATE_FEYNMAN = 60
    let ovm_PROPAGATE_COL_FEYNMAN = 61
    let ovm_PROPAGATE_VECTORSPINOR = 62
    let ovm_PROPAGATE_TENSOR2 = 63

(* \begin{dubious}
    [ovm_PROPAGATE_NONE] has to be split up to different types to work
    in conjunction with color MC \dots
   \end{dubious} *)
    let ovm_PROPAGATE_NONE = 64

    let ovm_FUSE_V_FF = -1
    let ovm_FUSE_F_VF = -2
    let ovm_FUSE_F_FV = -3
    let ovm_FUSE_VA_FF = -4
    let ovm_FUSE_F_VAF = -5
    let ovm_FUSE_F_FVA = -6
    let ovm_FUSE_VA2_FF = -7
    let ovm_FUSE_F_VA2F = -8
    let ovm_FUSE_F_FVA2 = -9
    let ovm_FUSE_A_FF = -10
    let ovm_FUSE_F_AF = -11
    let ovm_FUSE_F_FA = -12
    let ovm_FUSE_VL_FF = -13
    let ovm_FUSE_F_VLF = -14
    let ovm_FUSE_F_FVL = -15
    let ovm_FUSE_VR_FF = -16
    let ovm_FUSE_F_VRF = -17
    let ovm_FUSE_F_FVR = -18
    let ovm_FUSE_VLR_FF = -19
    let ovm_FUSE_F_VLRF = -20
    let ovm_FUSE_F_FVLR = -21
    let ovm_FUSE_SP_FF = -22
    let ovm_FUSE_F_SPF = -23
    let ovm_FUSE_F_FSP = -24
    let ovm_FUSE_S_FF = -25
    let ovm_FUSE_F_SF = -26
    let ovm_FUSE_F_FS = -27
    let ovm_FUSE_P_FF = -28
    let ovm_FUSE_F_PF = -29
    let ovm_FUSE_F_FP = -30
    let ovm_FUSE_SL_FF = -31
    let ovm_FUSE_F_SLF = -32
    let ovm_FUSE_F_FSL = -33
    let ovm_FUSE_SR_FF = -34
    let ovm_FUSE_F_SRF = -35
    let ovm_FUSE_F_FSR = -36
    let ovm_FUSE_SLR_FF = -37
    let ovm_FUSE_F_SLRF = -38
    let ovm_FUSE_F_FSLR = -39

    let ovm_FUSE_G_GG = -40
    let ovm_FUSE_V_SS = -41
    let ovm_FUSE_S_VV = -42
    let ovm_FUSE_S_VS = -43
    let ovm_FUSE_V_SV = -44
    let ovm_FUSE_S_SS = -45
    let ovm_FUSE_S_SVV = -46
    let ovm_FUSE_V_SSV = -47
    let ovm_FUSE_S_SSS = -48
    let ovm_FUSE_V_VVV = -49

    let ovm_FUSE_S_G2 = -50
    let ovm_FUSE_G_SG = -51
    let ovm_FUSE_G_GS = -52
    let ovm_FUSE_S_G2_SKEW = -53
    let ovm_FUSE_G_SG_SKEW = -54
    let ovm_FUSE_G_GS_SKEW = -55

    let inst_length = 8

(* Some helper functions. *)
    let printi ~lhs:l ~rhs1:r1 ?coupl:(cp = 0) ?coeff:(co = 0)
               ?rhs2:(r2 = 0) ?rhs3:(r3 = 0) ?rhs4:(r4 = 0) code =
      printf "@\n%d %d %d %d %d %d %d %d" code cp co l r1 r2 r3 r4

    let nl () = printf "@\n"

    let print_int_lst lst = nl (); lst |> List.iter (printf "%d   ")

    let print_str_lst lst = nl (); lst |> List.iter (printf "%s ")

    let break () = printi ~lhs:0 ~rhs1:0 0

(* Copied from below. Needed for header. *)
(* \begin{dubious}
     Could be fused with [lorentz_ordering].
   \end{dubious} *)
    type declarations =
      { scalars : F.wf list;
        spinors : F.wf list;
        conjspinors : F.wf list;
        realspinors : F.wf list;
        ghostspinors : F.wf list;
        vectorspinors : F.wf list;
        vectors : F.wf list;
        ward_vectors : F.wf list;
        massive_vectors : F.wf list;
        tensors_1 : F.wf list;
        tensors_2 : F.wf list;
        brs_scalars : F.wf list;
        brs_spinors : F.wf list;
        brs_conjspinors : F.wf list;
        brs_realspinors : F.wf list;
        brs_vectorspinors : F.wf list;
        brs_vectors : F.wf list;
        brs_massive_vectors : F.wf list }

    let rec classify_wfs' acc = function
      | [] -> acc
      | wf :: rest ->
          classify_wfs'
            (match CM.lorentz (F.flavor wf) with
            | Scalar -> {acc with scalars = wf :: acc.scalars}
            | Spinor -> {acc with spinors = wf :: acc.spinors}
            | ConjSpinor -> {acc with conjspinors = wf :: acc.conjspinors}
            | Majorana -> {acc with realspinors = wf :: acc.realspinors}
            | Maj_Ghost -> {acc with ghostspinors = wf :: acc.ghostspinors}
            | Vectorspinor ->
                {acc with vectorspinors = wf :: acc.vectorspinors}
            | Vector -> {acc with vectors = wf :: acc.vectors}
            | Massive_Vector ->
                {acc with massive_vectors = wf :: acc.massive_vectors}
            | Tensor_1 -> {acc with tensors_1 = wf :: acc.tensors_1}
            | Tensor_2 -> {acc with tensors_2 = wf :: acc.tensors_2}
            | BRS Scalar -> {acc with brs_scalars = wf :: acc.brs_scalars}
            | BRS Spinor -> {acc with brs_spinors = wf :: acc.brs_spinors}
            | BRS ConjSpinor -> {acc with brs_conjspinors =
                                 wf :: acc.brs_conjspinors}
            | BRS Majorana -> {acc with brs_realspinors =
                               wf :: acc.brs_realspinors}
            | BRS Vectorspinor -> {acc with brs_vectorspinors =
                                   wf :: acc.brs_vectorspinors}
            | BRS Vector -> {acc with brs_vectors = wf :: acc.brs_vectors}
            | BRS Massive_Vector -> {acc with brs_massive_vectors =
                                     wf :: acc.brs_massive_vectors}
            | BRS _ -> invalid_arg "Targets.classify_wfs': not needed here")
            rest

    let classify_wfs wfs = classify_wfs'
      { scalars = [];
        spinors = [];
        conjspinors = [];
        realspinors = [];
        ghostspinors = [];
        vectorspinors = [];
        vectors = [];
        ward_vectors = [];
        massive_vectors = [];
        tensors_1 = [];
        tensors_2 = [];
        brs_scalars = [];
        brs_spinors = [];
        brs_conjspinors = [];
        brs_realspinors = [];
        brs_vectorspinors = [];
        brs_vectors = [];
        brs_massive_vectors = [] } wfs

(* \thocwmodulesubsection{Sets and maps} *)

(* The OVM identifies all objects via integers. Therefore, we need maps
   which assign the abstract object a unique ID. *)

(* I want [int list]s with less elements to come first. Used in conjunction
   with the int list representation of momenta, this will set the outer
   particles at first position and allows the OVM to set them without further
   instructions. *)

(* \begin{dubious}
      Using the Momentum module might give better performance than integer lists?
   \end{dubious} *)
    let rec int_lst_compare (e1 : int list) (e2 : int list) =
      match e1,e2 with
      | [], []  -> 0
      | _, [] -> +1
      | [], _ -> -1
      | [_;_], [_] -> +1
      | [_], [_;_] -> -1
      | hd1 :: tl1, hd2 :: tl2 ->
          let c = compare hd1 hd2 in
          if (c != 0 && List.length tl1 = List.length tl2) then
            c
          else
            int_lst_compare tl1 tl2

(* We need a canonical ordering for the different types
   of wfs. Copied, and slightly modified to order [wf]s, from
   \texttt{fusion.ml}. *)

    let lorentz_ordering wf =
      match CM.lorentz (F.flavor wf) with
      | Scalar -> 0
      | Spinor -> 1
      | ConjSpinor -> 2
      | Majorana -> 3
      | Vector -> 4
      | Massive_Vector -> 5
      | Tensor_2 -> 6
      | Tensor_1 -> 7
      | Vectorspinor -> 8
      | BRS Scalar -> 9
      | BRS Spinor -> 10
      | BRS ConjSpinor -> 11
      | BRS Majorana -> 12
      | BRS Vector -> 13
      | BRS Massive_Vector -> 14
      | BRS Tensor_2 -> 15
      | BRS Tensor_1 -> 16
      | BRS Vectorspinor -> 17
      | Maj_Ghost -> invalid_arg "lorentz_ordering: not implemented"
      | BRS _ -> invalid_arg "lorentz_ordering: not needed"

    let wf_compare (wf1, mult1) (wf2, mult2) =
      let c1 = compare (lorentz_ordering wf1) (lorentz_ordering wf2) in
      if c1 <> 0 then
        c1
      else
        let c2 = compare wf1 wf2 in
        if c2 <> 0 then
          c2
        else
          compare mult1 mult2

    let amp_compare amp1 amp2 =
      let cflow a = CM.flow (F.incoming a) (F.outgoing a) in
      let c1 = compare (cflow amp1) (cflow amp2) in
      if c1 <> 0 then
        c1
      else
        let process_sans_color a =
          (List.map CM.flavor_sans_color (F.incoming a),
           List.map CM.flavor_sans_color (F.outgoing a)) in
        compare (process_sans_color amp1) (process_sans_color amp2)

    let level_compare (f1, amp1) (f2, amp2) =
      let p1 = F.momentum_list (F.lhs f1)
      and p2 = F.momentum_list (F.lhs f2) in
      let c1 = int_lst_compare p1 p2 in
      if c1 <> 0 then
        c1
      else
        let c2 = compare f1 f2 in
        if c2 <> 0 then
          c2
        else
          amp_compare amp1 amp2

    module ISet = Set.Make (struct type t = int list
                            let compare = int_lst_compare end)

    module WFSet = Set.Make (struct type t = CF.wf * int
                             let compare = wf_compare end)

    module CSet = Set.Make (struct type t = CM.constant
                            let compare = compare end)

    module FSet = Set.Make (struct type t = F.fusion * F.amplitude
                            let compare = level_compare end)

(* \begin{dubious}
     It might be preferable to use a [PMap] which maps mom to int, instead of
     this way. More standard functions like [mem] could be used. Also, [get_ID]
     would be faster, $\mathcal{O}(\log N)$ instead of $\mathcal{O}(N)$, and
     simpler.  For 8 gluons: N=127 momenta. Minor performance issue.
   \end{dubious} *)

    module IMap = Map.Make (struct type t = int let compare = compare end)

(* For [wf]s it is crucial for the performance to use a different type of
   [Map]s. *)

    module WFMap = Map.Make (struct type t = CF.wf * int
                             let compare = wf_compare end)

    type lookups = { pmap : int list IMap.t;
                     wfmap : int WFMap.t;
                     cmap : CM.constant IMap.t * CM.constant IMap.t;
                     amap : F.amplitude IMap.t;
                     n_wfs : int list;
                     amplitudes : CF.amplitudes;
                     dict : F.amplitude -> F.wf -> int }

    let largest_key imap =
      if (IMap.is_empty imap) then
        failwith "largest_key: Map is empty!"
      else
        fst (IMap.max_binding imap)

(* OCaml's [compare] from pervasives cannot compare functional types, e.g.
   for type [amplitude], if no specific equality function is given ("equal:
   functional value"). Therefore, we allow to specify the ordering. *)

    let get_ID' comp map elt : int =
      let smallmap = IMap.filter (fun _ x -> (comp x elt) = 0 ) map in
      if IMap.is_empty smallmap then
        raise Not_found
      else
        fst (IMap.min_binding smallmap)

(* \begin{dubious}
     Trying to curry [map] here leads to type errors of the
     polymorphic function [get_ID]?
   \end{dubious} *)

    let get_ID map = match map with
      | map -> get_ID' compare map

    let get_const_ID map x = match map with
      | (map1, map2) -> try get_ID' compare map1 x with
                       _ -> try get_ID' compare map2 x with
                       _ -> failwith "Impossible"

(* Creating an integer map of a list with an optional argument that
   indicates where the map should start counting. *)

    let map_of_list ?start:(st=1) lst =
      let g (ind, map) wf = (succ ind, IMap.add ind wf map) in
      lst |> List.fold_left g (st, IMap.empty) |> snd

    let wf_map_of_list ?start:(st=1) lst =
      let g (ind, map) wf = (succ ind, WFMap.add wf ind map) in
      lst |> List.fold_left g (st, WFMap.empty) |> snd

(* \thocwmodulesubsection{Header} *)

(* \begin{dubious}
     [Bijan:]
     It would be nice to save the creation date as comment. However, the Unix
     module doesn't seem to be loaded on default.
   \end{dubious} *)

    let version =
      String.concat " " [Config.version; Config.status; Config.date]
    let model_name =
      let basename = Filename.basename Sys.executable_name in
      try
        Filename.chop_extension basename
      with
      | _ -> basename


    let print_description cmdline  =
      printf "Model %s\n" model_name;
      printf "OVM %s\n" version;
      printf "@\nBytecode file generated automatically by O'Mega for OVM";
      printf "@\nDo not delete any lines. You called O'Mega with";
      printf "@\n  %s" cmdline;
      (*i
      let t = Unix.localtime (Unix.time() ) in
        printf "@\n on %5d %5d %5d" (succ t.Unix.tm_mon) t.Unix.tm_mday
               t.Unix.tm_year;
       i*)
      printf "@\n"

    let num_classified_wfs wfs =
      let wfs' = classify_wfs wfs in
      List.map List.length
        [ wfs'.scalars @ wfs'.brs_scalars;
          wfs'.spinors @ wfs'.brs_spinors;
          wfs'.conjspinors @ wfs'.brs_conjspinors;
          wfs'.realspinors @ wfs'.brs_realspinors @ wfs'.ghostspinors;
          wfs'.vectors @ wfs'.massive_vectors @ wfs'.brs_vectors
            @ wfs'.brs_massive_vectors @ wfs'.ward_vectors;
          wfs'.tensors_2;
          wfs'.tensors_1;
          wfs'.vectorspinors ]

    let description_classified_wfs =
      [ "N_scalars";
        "N_spinors";
        "N_conjspinors";
        "N_bispinors";
        "N_vectors";
        "N_tensors_2";
        "N_tensors_1";
        "N_vectorspinors" ]

    let num_particles_in amp =
      match CF.flavors amp with
      | [] -> 0
      | (fin, _) :: _ -> List.length fin

    let num_particles_out amp =
      match CF.flavors amp with
      | [] -> 0
      | (_, fout) :: _ -> List.length fout

    let num_particles amp =
      match CF.flavors amp with
      | [] -> 0
      | (fin, fout) :: _ -> List.length fin + List.length fout

    let num_color_indices_default = 2 (* Standard model and non-color-exotica *)

    let num_color_indices amp =
      try CFlow.rank (List.hd (CF.color_flows amp)) with
      _ -> num_color_indices_default

    let num_color_factors amp =
      let table = CF.color_factors amp in
      let n_cflow = Array.length table
      and n_cfactors = ref 0 in
      for c1 = 0 to pred n_cflow do
        for c2 = 0 to pred n_cflow do
          if c1 <= c2 then begin
            match table.(c1).(c2) with
            | [] -> ()
            | _ -> incr n_cfactors
          end
        done
      done;
      !n_cfactors

    let num_helicities amp = amp |> CF.helicities |> List.length

    let num_flavors amp = amp |> CF.flavors |> List.length

    let num_ks amp = amp |> CF.processes |> List.length

    let num_color_flows amp = amp |> CF.color_flows |> List.length

(* Use [fst] since [WFSet.t = F.wf * int]. *)
    let num_wfs wfset = wfset |> WFSet.elements |> List.map fst
                              |> num_classified_wfs

(* [largest_key] gives the number of momenta if applied to [pmap]. *)

    let num_lst lookups wfset =
      [ largest_key lookups.pmap;
        num_particles lookups.amplitudes;
        num_particles_in lookups.amplitudes;
        num_particles_out lookups.amplitudes;
        num_ks lookups.amplitudes;
        num_helicities lookups.amplitudes;
        num_color_flows lookups.amplitudes;
        num_color_indices lookups.amplitudes;
        num_flavors lookups.amplitudes;
        num_color_factors lookups.amplitudes ] @ num_wfs wfset

    let description_lst =
      [ "N_momenta";
        "N_particles";
        "N_prt_in";
        "N_prt_out";
        "N_amplitudes";
        "N_helicities";
        "N_col_flows";
        "N_col_indices";
        "N_flavors";
        "N_col_factors" ] @ description_classified_wfs

    let print_header' numbers =
      let chopped_num_lst = ThoList.chopn inst_length numbers
      and chopped_desc_lst = ThoList.chopn inst_length description_lst
      and printer a b = print_str_lst a; print_int_lst b in
      List.iter2 printer chopped_desc_lst chopped_num_lst

    let print_header lookups wfset = print_header' (num_lst lookups wfset)

    let print_zero_header () =
      let rec zero_list' j =
        if j < 1 then []
        else 0 :: zero_list' (j - 1) in
      let zero_list i = zero_list' (i + 1) in
      description_lst |> List.length |> zero_list |> print_header'

(* \thocwmodulesubsection{Tables} *)

    let print_spin_table' tuples =
      match tuples with
      | [] -> ()
      | _ -> tuples |> List.iter ( fun (tuple1, tuple2) ->
          tuple1 @ tuple2 |> List.map (Printf.sprintf "%d ")
                          |> String.concat "" |> printf "@\n%s" )

    let print_spin_table amplitudes =
      printf "@\nSpin states table";
      print_spin_table' @@ CF.helicities amplitudes

    let print_flavor_table tuples =
      match tuples with
      | [] -> ()
      | _ -> List.iter ( fun tuple -> tuple
                        |> List.map (fun f -> Printf.sprintf "%d " @@ M.pdg f)
                        |> String.concat "" |> printf "@\n%s"
                       ) tuples

    let print_flavor_tables amplitudes =
      printf "@\nFlavor states table";
      print_flavor_table @@ List.map (fun (fin, fout) -> fin @ fout)
                         @@ CF.flavors amplitudes

    let print_color_flows_table' tuple =
        match CFlow.to_lists tuple with
        | [] -> ()
        | cfs -> printf "@\n%s" @@ String.concat "" @@ List.map
                  ( fun cf -> cf |> List.map (Printf.sprintf "%d ")
                                 |> String.concat ""
                  ) cfs

    let print_color_flows_table tuples =
      match tuples with
      | [] -> ()
      | _ -> List.iter print_color_flows_table' tuples

    let print_ghost_flags_table tuples =
      match tuples with
      | [] -> ()
      | _ ->
        List.iter (fun tuple ->
        match CFlow.ghost_flags tuple with
            | [] -> ()
            | gfs -> printf "@\n"; List.iter (fun gf -> printf "%s "
              (if gf then "1" else "0") ) gfs
        ) tuples

    let format_power
      { CFlow.num = num; CFlow.den = den; CFlow.power = pwr } =
      match num, den, pwr with
      | _, 0, _ -> invalid_arg "targets.format_power: zero denominator"
      | n, d, p -> [n; d; p]

    let format_powers = function
      | [] -> [0]
      | powers -> List.flatten (List.map format_power powers)

    (*i
    (* We go through the array line by line and collect all colorfactors which
     * are nonzero because their corresponding color flows match.
     * With the gained intset, we would be able to print only the necessary
     * coefficients of the symmetric matrix and indicate from where the OVM
     * can copy the rest. However, this approach gets really slow for many
     * gluons and we can save at most 3 numbers per line.*)

    let print_color_factor_table_funct table =
      let n_cflow = Array.length table in
      let (intset, _, _ ) =
        let rec fold_array (set, cf1, cf2) =
          if cf1 > pred n_cflow then (set, 0, 0)
          else
              let returnset =
              match table.(cf1).(cf2) with
                  | [] -> set
                  | cf ->
                      ISet.add ([succ cf1; succ cf2] @ (format_powers cf)) set
              in
              if cf2 < pred n_cflow then
                fold_array (returnset, cf1, succ cf2) else
                fold_array (returnset, succ cf1, 0)
        in
        fold_array (ISet.empty, 0, 0)
      in
      let map = map_of_list (ISet.elements intset) in
      List.iter (fun x -> printf "@\n"; let xth = List.nth x in
      if (xth 0 <= xth 1) then List.iter (printf "%d ") x
      else printf "%d %d" 0 (get_ID map x))
        (ISet.elements intset)

    let print_color_factor_table_old table =
      let n_cflow = Array.length table in
      let (intlsts, _, _ ) =
        let rec fold_array (lsts, cf1, cf2) =
          if cf1 > pred n_cflow then (lsts, 0, 0)
          else
              let returnlsts =
              match table.(cf1).(cf2) with
                  | [] -> lsts
                  | cf -> ([succ cf1; succ cf2] @ (format_powers cf)) :: lsts
              in
              if cf2 < pred n_cflow then
                fold_array (returnlsts, cf1, succ cf2) else
                fold_array (returnlsts, succ cf1, 0)
        in
        fold_array ([], 0, 0)
      in
      let intlsts = List.rev intlsts in
      List.iter (fun x -> printf "@\n"; List.iter (printf "%d ") x ) intlsts
      i*)

(* Straightforward iteration gives a great speedup compared to the fancier
   approach which only collects nonzero colorfactors. *)

    let print_color_factor_table table =
      let n_cflow = Array.length table in
      if n_cflow > 0 then begin
        for c1 = 0 to pred n_cflow do
          for c2 = 0 to pred n_cflow do
            if c1 <= c2 then begin
              match table.(c1).(c2) with
              | [] -> ()
              | cf -> printf "@\n"; List.iter (printf "%9d")
                ([succ c1; succ c2] @ (format_powers cf));
            end
          done
        done
      end

    let option_to_binary = function
      | Some _ -> "1"
      | None -> "0"

    let print_flavor_color_table n_flv n_cflow table =
      if n_flv > 0 then begin
        for c = 0 to pred n_cflow do
          printf "@\n";
          for f = 0 to pred n_flv do
            printf "%s " (option_to_binary table.(f).(c))
          done;
        done;
      end

    let print_color_tables amplitudes =
      let cflows =  CF.color_flows amplitudes
      and cfactors = CF.color_factors amplitudes in
      printf "@\nColor flows table: [ (i, j) (k, l) -> (m, n) ...]";
      print_color_flows_table cflows;
      printf "@\nColor ghost flags table:";
      print_ghost_flags_table cflows;
      printf "@\nColor factors table: [ i, j: num den power], %s"
        "i, j are indexed color flows";
      print_color_factor_table cfactors;
      printf "@\nFlavor color combination is allowed:";
      print_flavor_color_table (num_flavors amplitudes) (List.length
        (CF.color_flows amplitudes)) (CF.process_table amplitudes)

(* \thocwmodulesubsection{Momenta} *)

(* Add the momenta of a WFSet to a Iset. For now, we are throwing away the
   information to which amplitude the momentum belongs. This could be optimized
   for random color flow computations. *)

    let momenta_set wfset =
      let get_mom wf = wf |> fst |> F.momentum_list in
      let momenta = List.map get_mom (WFSet.elements wfset) in
      momenta |> List.fold_left (fun set x -> set |> ISet.add x) ISet.empty

    let chop_in_3 lst =
      let ceil_div i j = if (i mod j = 0) then i/j else i/j + 1 in
      ThoList.chopn (ceil_div (List.length lst) 3) lst

(* Assign momenta via instruction code. External momenta [[_]] are already
   set by the OVM. To avoid unnecessary look-ups of IDs we seperate two cases.
   If we have more, we split up in two or three parts. *)

    let add_mom p pmap =
      let print_mom lhs rhs1 rhs2 rhs3 = if (rhs1!= 0) then
        printi ~lhs:lhs ~rhs1:rhs1 ~rhs2:rhs2 ~rhs3:rhs3 ovm_ADD_MOMENTA in
      let get_p_ID = get_ID pmap in
      match p with
      | [] | [_] -> print_mom 0 0 0 0
      | [rhs1;rhs2] -> print_mom (get_p_ID [rhs1;rhs2]) rhs1 rhs2 0
      | [rhs1;rhs2;rhs3] -> print_mom (get_p_ID [rhs1;rhs2;rhs3]) rhs1 rhs2 rhs3
      | more ->
          let ids = List.map get_p_ID (chop_in_3 more) in
          if (List.length ids = 3) then
            print_mom (get_p_ID more) (List.nth ids 0) (List.nth ids 1)
              (List.nth ids 2)
          else
            print_mom (get_p_ID more) (List.nth ids 0) (List.nth ids 1) 0

(* Hand through the current level and print level seperators if necessary. *)

    let add_all_mom lookups pset =
      let add_all' level p =
        let level' = List.length p in
        if (level' > level && level' > 3) then break ();
        add_mom p lookups.pmap; level'
      in
      ignore (pset |> ISet.elements |> List.fold_left add_all' 1)

(* Expand a set of momenta to contain all needed momenta for the computation
   in the OVM. For this, we create a list of sets which contains the chopped
   momenta and unify them afterwards. If the set has become larger, we
   expand again. *)

    let rec expand_pset p =
      let momlst = ISet.elements p in
      let pset_of lst = List.fold_left (fun s x -> ISet.add x s) ISet.empty
        lst in
      let sets = List.map (fun x -> pset_of (chop_in_3 x) ) momlst in
      let bigset = List.fold_left ISet.union ISet.empty sets in
      let biggerset = ISet.union bigset p in
      if (List.length momlst < List.length (ISet.elements biggerset) ) then
        expand_pset biggerset
      else
        biggerset

    let mom_ID pmap wf = get_ID pmap (F.momentum_list wf)

(* \thocwmodulesubsection{Wavefunctions and externals} *)

(* [mult_wf] is needed because the [wf] with same combination of flavor and
   momentum can have different dependencies and content. *)

    let mult_wf dict amplitude wf =
      try
        wf, dict amplitude wf
      with
        | Not_found -> wf, 0

(* Build the union of all [wf]s of all amplitudes and a map of the amplitudes. *)

    let wfset_amps amplitudes =
      let amap = amplitudes |> CF.processes |> List.sort amp_compare
                            |> map_of_list
      and dict = CF.dictionary amplitudes in
      let wfset_amp amp =
        let f = mult_wf dict amp in
        let lst = List.map f ((F.externals amp) @ (F.variables amp)) in
        lst |> List.fold_left (fun s x -> WFSet.add x s) WFSet.empty in
      let list_of_sets = amplitudes |> CF.processes |> List.map wfset_amp in
        List.fold_left WFSet.union WFSet.empty list_of_sets, amap

(* To obtain the Fortran index, we substract the number of precedent wave
   functions. *)

    let lorentz_ordering_reduced wf =
      match CM.lorentz (F.flavor wf) with
      | Scalar | BRS Scalar -> 0
      | Spinor | BRS Spinor -> 1
      | ConjSpinor | BRS ConjSpinor -> 2
      | Majorana | BRS Majorana -> 3
      | Vector | BRS Vector | Massive_Vector | BRS Massive_Vector -> 4
      | Tensor_2 | BRS Tensor_2 -> 5
      | Tensor_1 | BRS Tensor_1 -> 6
      | Vectorspinor | BRS Vectorspinor -> 7
      | Maj_Ghost -> invalid_arg "lorentz_ordering: not implemented"
      | BRS _ -> invalid_arg "lorentz_ordering: not needed"

    let wf_index wfmap num_lst (wf, i) =
      let wf_ID = WFMap.find (wf, i) wfmap
      and sum lst = List.fold_left (fun x y -> x+y) 0 lst in
        wf_ID - sum (ThoList.hdn (lorentz_ordering_reduced wf) num_lst)

    let print_ext lookups amp_ID inc (wf, i) =
      let mom = (F.momentum_list wf) in
      let outer_index = if List.length mom = 1 then List.hd mom else
        failwith "targets.print_ext: called with non-external particle"
      and f = F.flavor wf in
      let pdg = CM.pdg f
      and wf_code =
        match CM.lorentz f with
        | Scalar -> ovm_LOAD_SCALAR
        | BRS Scalar -> ovm_LOAD_BRS_SCALAR
        | Spinor ->
            if inc then ovm_LOAD_SPINOR_INC
            else ovm_LOAD_SPINOR_OUT
        | BRS Spinor ->
            if inc then ovm_LOAD_BRS_SPINOR_INC
            else ovm_LOAD_BRS_SPINOR_OUT
        | ConjSpinor ->
            if inc then ovm_LOAD_CONJSPINOR_INC
            else ovm_LOAD_CONJSPINOR_OUT
        | BRS ConjSpinor ->
            if inc then ovm_LOAD_BRS_CONJSPINOR_INC
            else ovm_LOAD_BRS_CONJSPINOR_OUT
        | Vector | Massive_Vector ->
            if inc then ovm_LOAD_VECTOR_INC
            else ovm_LOAD_VECTOR_OUT
        | BRS Vector | BRS Massive_Vector ->
            if inc then ovm_LOAD_BRS_VECTOR_INC
            else ovm_LOAD_BRS_VECTOR_OUT
        | Tensor_2 ->
            if inc then ovm_LOAD_TENSOR2_INC
            else ovm_LOAD_TENSOR2_OUT
        | Vectorspinor | BRS Vectorspinor ->
            if inc then ovm_LOAD_VECTORSPINOR_INC
            else ovm_LOAD_VECTORSPINOR_OUT
        | Majorana ->
            if inc then ovm_LOAD_MAJORANA_INC
            else ovm_LOAD_MAJORANA_OUT
        | BRS Majorana ->
            if inc then ovm_LOAD_BRS_MAJORANA_INC
            else ovm_LOAD_BRS_MAJORANA_OUT
        | Maj_Ghost ->
            if inc then ovm_LOAD_MAJORANA_GHOST_INC
            else ovm_LOAD_MAJORANA_GHOST_OUT
        | Tensor_1 ->
            invalid_arg "targets.print_ext: Tensor_1 only internal"
        | BRS _ ->
            failwith "targets.print_ext: Not implemented"
      and wf_ind = wf_index lookups.wfmap lookups.n_wfs (wf, i)
      in
        printi wf_code ~lhs:wf_ind ~coupl:(abs(pdg)) ~rhs1:outer_index ~rhs4:amp_ID

    let print_ext_amp lookups amplitude =
      let incoming = (List.map (fun _ -> true) (F.incoming amplitude) @
                      List.map (fun _ -> false) (F.outgoing amplitude))
      and amp_ID = get_ID' amp_compare lookups.amap amplitude in
      let wf_tpl wf = mult_wf lookups.dict amplitude wf in
      let print_ext_wf inc wf = wf |> wf_tpl |> print_ext lookups amp_ID inc in
        List.iter2 print_ext_wf incoming (F.externals amplitude)

    let print_externals lookups seen_wfs amplitude =
      let externals =
        List.combine
          (F.externals amplitude)
          (List.map (fun _ -> true) (F.incoming amplitude) @
           List.map (fun _ -> false) (F.outgoing amplitude)) in
      List.fold_left (fun seen (wf, incoming) ->
        let amp_ID = get_ID' amp_compare lookups.amap amplitude in
        let wf_tpl = mult_wf lookups.dict amplitude wf in
        if not (WFSet.mem wf_tpl seen) then begin
          wf_tpl |> print_ext lookups amp_ID incoming
        end;
        WFSet.add wf_tpl seen) seen_wfs externals

(* [print_externals] and [print_ext_amp] do in principle the same thing but
   [print_externals] filters out dublicate external wave functions. Even with
   [print_externals] the same (numerically) external wave function will be
   loaded if it belongs to a different color flow, just as in the native Fortran
   code.  For color MC, [print_ext_amp] has to be used (redundant instructions
   but only one flow is computed) and the filtering of duplicate fusions has to
   be disabled. *)

    let print_ext_amps lookups =
      let print_external_amp s x = print_externals lookups s x in
      ignore (
        List.fold_left print_external_amp WFSet.empty
          (CF.processes lookups.amplitudes)
        )
      (*i
      List.iter (print_ext_amp lookups) (CF.processes lookups.amplitudes)
      i*)

(* \thocwmodulesubsection{Currents} *)

(* Parallelization issues: All fusions have to be completed before the
   propagation takes place. Preferably each fusion and propagation is done
   by one thread.  Solution: All fusions are subinstructions, i.e. if
   they are read by the main loop they are skipped. If a propagation
   occurs, all fusions have to be computed first. The additional control
   bit is the sign of the first int of an instruction. *)

    (*i TODO: (bcn 2014-07-21) Majorana support will come some day maybe i*)
    let print_fermion_current code_a code_b code_c coeff lhs c wf1 wf2 fusion =
      let printc code r1 r2 = printi code ~lhs:lhs ~coupl:c ~coeff:coeff
        ~rhs1:r1 ~rhs2:r2 in
      match fusion with
      | F13 -> printc code_a wf1 wf2
      | F31 -> printc code_a wf2 wf1
      | F23 -> printc code_b wf1 wf2
      | F32 -> printc code_b wf2 wf1
      | F12 -> printc code_c wf1 wf2
      | F21 -> printc code_c wf2 wf1

      let ferm_print_current = function
        | coeff, Psibar, V, Psi -> print_fermion_current
          ovm_FUSE_V_FF ovm_FUSE_F_VF ovm_FUSE_F_FV coeff
        | coeff, Psibar, VA, Psi -> print_fermion_current
          ovm_FUSE_VA_FF ovm_FUSE_F_VAF ovm_FUSE_F_FVA coeff
        | coeff, Psibar, VA2, Psi -> print_fermion_current
          ovm_FUSE_VA2_FF ovm_FUSE_F_VA2F ovm_FUSE_F_FVA2 coeff
        | coeff, Psibar, A, Psi -> print_fermion_current
          ovm_FUSE_A_FF ovm_FUSE_F_AF ovm_FUSE_F_FA coeff
        | coeff, Psibar, VL, Psi -> print_fermion_current
          ovm_FUSE_VL_FF ovm_FUSE_F_VLF ovm_FUSE_F_FVL coeff
        | coeff, Psibar, VR, Psi -> print_fermion_current
          ovm_FUSE_VR_FF ovm_FUSE_F_VRF ovm_FUSE_F_FVR coeff
        | coeff, Psibar, VLR, Psi -> print_fermion_current
          ovm_FUSE_VLR_FF ovm_FUSE_F_VLRF ovm_FUSE_F_FVLR coeff
        | coeff, Psibar, SP, Psi -> print_fermion_current
          ovm_FUSE_SP_FF ovm_FUSE_F_SPF ovm_FUSE_F_FSP coeff
        | coeff, Psibar, S, Psi -> print_fermion_current
          ovm_FUSE_S_FF ovm_FUSE_F_SF ovm_FUSE_F_FS coeff
        | coeff, Psibar, P, Psi -> print_fermion_current
          ovm_FUSE_P_FF ovm_FUSE_F_PF ovm_FUSE_F_FP coeff
        | coeff, Psibar, SL, Psi -> print_fermion_current
          ovm_FUSE_SL_FF ovm_FUSE_F_SLF ovm_FUSE_F_FSL coeff
        | coeff, Psibar, SR, Psi -> print_fermion_current
          ovm_FUSE_SR_FF ovm_FUSE_F_SRF ovm_FUSE_F_FSR coeff
        | coeff, Psibar, SLR, Psi -> print_fermion_current
          ovm_FUSE_SLR_FF ovm_FUSE_F_SLRF ovm_FUSE_F_FSLR coeff
        | _, Psibar, _, Psi -> invalid_arg
          "Targets.Fortran.VM: no superpotential here"
        | _, Chibar, _, _ | _, _, _, Chi -> invalid_arg
          "Targets.Fortran.VM: Majorana spinors not handled"
        | _, Gravbar, _, _ | _, _, _, Grav -> invalid_arg
          "Targets.Fortran.VM: Gravitinos not handled"

    let children2 rhs =
      match F.children rhs with
      | [wf1; wf2] -> (wf1, wf2)
      | _ -> failwith "Targets.children2: can't happen"

    let children3 rhs =
      match F.children rhs with
      | [wf1; wf2; wf3] -> (wf1, wf2, wf3)
      | _ -> invalid_arg "Targets.children3: can't happen"

    let print_vector4 c lhs wf1 wf2 wf3 fusion (coeff, contraction) =
      let printc r1 r2 r3 = printi ovm_FUSE_V_VVV ~lhs:lhs ~coupl:c
        ~coeff:coeff ~rhs1:r1 ~rhs2:r2 ~rhs3:r3 in
      match contraction, fusion with
      | C_12_34, (F341|F431|F342|F432|F123|F213|F124|F214)
      | C_13_42, (F241|F421|F243|F423|F132|F312|F134|F314)
      | C_14_23, (F231|F321|F234|F324|F142|F412|F143|F413) ->
          printc wf1 wf2 wf3
      | C_12_34, (F134|F143|F234|F243|F312|F321|F412|F421)
      | C_13_42, (F124|F142|F324|F342|F213|F231|F413|F431)
      | C_14_23, (F123|F132|F423|F432|F214|F241|F314|F341) ->
          printc wf2 wf3 wf1
      | C_12_34, (F314|F413|F324|F423|F132|F231|F142|F241)
      | C_13_42, (F214|F412|F234|F432|F123|F321|F143|F341)
      | C_14_23, (F213|F312|F243|F342|F124|F421|F134|F431) ->
          printc wf1 wf3 wf2

    let print_current lookups lhs amplitude rhs =
      let f = mult_wf lookups.dict amplitude in
      match F.coupling rhs with
      | V3 (vertex, fusion, constant) ->
          let ch1, ch2 = children2 rhs in
          let wf1 = wf_index lookups.wfmap lookups.n_wfs (f ch1)
          and wf2 = wf_index lookups.wfmap lookups.n_wfs (f ch2)
          and p1 = mom_ID lookups.pmap ch1
          and p2 = mom_ID lookups.pmap ch2
          and const_ID = get_const_ID lookups.cmap constant in
          let c = if (F.sign rhs) < 0 then - const_ID else const_ID in
          begin match vertex with
          | FBF (coeff, fb, b, f) ->
              begin match coeff, fb, b, f with
              | _, Psibar, VLRM, Psi | _, Psibar, SPM, Psi
              | _, Psibar, TVA, Psi | _, Psibar, TVAM, Psi
              | _, Psibar, TLR, Psi | _, Psibar, TLRM, Psi
              | _, Psibar, TRL, Psi | _, Psibar, TRLM, Psi -> failwith
       "print_current: V3: Momentum dependent fermion couplings not implemented"
              | _, _, _, _ ->
                  ferm_print_current (coeff, fb, b, f) lhs c wf1 wf2 fusion
              end
          | PBP (_, _, _, _) ->
              failwith "print_current: V3: PBP not implemented"
          | BBB (_, _, _, _) ->
              failwith "print_current: V3: BBB not implemented"
          | GBG (_, _, _, _) ->
              failwith "print_current: V3: GBG not implemented"

          | Gauge_Gauge_Gauge coeff ->
              let printc r1 r2 r3 r4 = printi ovm_FUSE_G_GG
                ~lhs:lhs ~coupl:c ~coeff:coeff ~rhs1:r1 ~rhs2:r2 ~rhs3:r3
                ~rhs4:r4 in
              begin match fusion with
              | (F23|F31|F12) -> printc wf1 p1 wf2 p2
              | (F32|F13|F21) -> printc wf2 p2 wf1 p1
              end

          | I_Gauge_Gauge_Gauge _ ->
              failwith "print_current: I_Gauge_Gauge_Gauge: not implemented"

          | Scalar_Vector_Vector coeff ->
              let printc code r1 r2 = printi code
                ~lhs:lhs ~coupl:c ~coeff:coeff ~rhs1:r1 ~rhs2:r2 in
              begin match fusion with
              | (F23|F32) -> printc ovm_FUSE_S_VV wf1 wf2
              | (F12|F13) -> printc ovm_FUSE_V_SV wf1 wf2
              | (F21|F31) -> printc ovm_FUSE_V_SV wf2 wf1
              end

          | Scalar_Scalar_Scalar coeff ->
              printi ovm_FUSE_S_SS ~lhs:lhs ~coupl:c ~coeff:coeff ~rhs1:wf1 ~rhs2:wf2

          | Vector_Scalar_Scalar coeff ->
              let printc code ?flip:(f = 1) r1 r2 r3 r4 = printi code
                ~lhs:lhs ~coupl:(c*f) ~coeff:coeff ~rhs1:r1 ~rhs2:r2 ~rhs3:r3
                ~rhs4:r4 in
              begin match fusion with
              | F23 -> printc ovm_FUSE_V_SS wf1 p1 wf2 p2
              | F32 -> printc ovm_FUSE_V_SS wf2 p2 wf1 p1
              | F12 -> printc ovm_FUSE_S_VS wf1 p1 wf2 p2
              | F21 -> printc ovm_FUSE_S_VS wf2 p2 wf1 p1
              | F13 -> printc ovm_FUSE_S_VS wf1 p1 wf2 p2 ~flip:(-1)
              | F31 -> printc ovm_FUSE_S_VS wf2 p2 wf1 p1 ~flip:(-1)
              end

          | Aux_Vector_Vector _ ->
              failwith "print_current: V3: not implemented"

          | Aux_Scalar_Scalar _ ->
              failwith "print_current: V3: not implemented"

          | Aux_Scalar_Vector _ ->
              failwith "print_current: V3: not implemented"

          | Graviton_Scalar_Scalar _ ->
              failwith "print_current: V3: not implemented"

          | Graviton_Vector_Vector _ ->
              failwith "print_current: V3: not implemented"

          | Graviton_Spinor_Spinor _ ->
              failwith "print_current: V3: not implemented"

          | Dim4_Vector_Vector_Vector_T _ ->
              failwith "print_current: V3: not implemented"

          | Dim4_Vector_Vector_Vector_L _ ->
              failwith "print_current: V3: not implemented"

          | Dim6_Gauge_Gauge_Gauge _ ->
              failwith "print_current: V3: not implemented"

          | Dim4_Vector_Vector_Vector_T5 _ ->
              failwith "print_current: V3: not implemented"

          | Dim4_Vector_Vector_Vector_L5 _ ->
              failwith "print_current: V3: not implemented"

          | Dim6_Gauge_Gauge_Gauge_5 _ ->
              failwith "print_current: V3: not implemented"

          | Aux_DScalar_DScalar _ ->
              failwith "print_current: V3: not implemented"

          | Aux_Vector_DScalar _ ->
              failwith "print_current: V3: not implemented"

          | Dim5_Scalar_Gauge2 coeff ->
              let printc code r1 r2 r3 r4 =  printi code
                ~lhs:lhs ~coupl:c ~coeff:coeff ~rhs1:r1 ~rhs2:r2 ~rhs3:r3
                ~rhs4:r4  in
              begin match fusion with
              | (F23|F32) -> printc ovm_FUSE_S_G2 wf1 p1 wf2 p2
              | (F12|F13) -> printc ovm_FUSE_G_SG wf1 p1 wf2 p2
              | (F21|F31) -> printc ovm_FUSE_G_GS wf2 p2 wf1 p1
              end

          | Dim5_Scalar_Gauge2_Skew coeff ->
              let printc code ?flip:(f = 1) r1 r2 r3 r4 =  printi code
                ~lhs:lhs ~coupl:(c*f) ~coeff:coeff ~rhs1:r1 ~rhs2:r2 ~rhs3:r3
                ~rhs4:r4  in
              begin match fusion with
              | (F23|F32) -> printc ovm_FUSE_S_G2_SKEW wf1 p1 wf2 p2
              | (F12|F13) -> printc ovm_FUSE_G_SG_SKEW wf1 p1 wf2 p2
              | (F21|F31) -> printc ovm_FUSE_G_GS_SKEW wf2 p1 wf1 p2 ~flip:(-1)
              end

          | Dim5_Scalar_Vector_Vector_T _ ->
              failwith "print_current: V3: not implemented"

          | Dim5_Scalar_Vector_Vector_U _ ->
              failwith "print_current: V3: not implemented"

          | Dim5_Scalar_Scalar2 _ ->
              failwith "print_current: V3: not implemented"

          | Dim6_Vector_Vector_Vector_T _ ->
              failwith "print_current: V3: not implemented"

          | Tensor_2_Vector_Vector _ ->
              failwith "print_current: V3: not implemented"

          | Tensor_2_Scalar_Scalar _ ->
              failwith "print_current: V3: not implemented"

          | Dim5_Tensor_2_Vector_Vector_1 _ ->
              failwith "print_current: V3: not implemented"

          | Dim5_Tensor_2_Vector_Vector_2 _ ->
              failwith "print_current: V3: not implemented"

          | Dim7_Tensor_2_Vector_Vector_T _ ->
              failwith "print_current: V3: not implemented"

          | Dim5_Scalar_Vector_Vector_TU _ ->
              failwith "print_current: V3: not implemented"

          | Scalar_Vector_Vector_t _ ->
              failwith "print_current: V3: not implemented"

          | Tensor_2_Vector_Vector_cf _ ->
              failwith "print_current: V3: not implemented"

          | Tensor_2_Scalar_Scalar_cf _ ->
              failwith "print_current: V3: not implemented"

          | Tensor_2_Vector_Vector_1 _ ->
              failwith "print_current: V3: not implemented"

          | Tensor_2_Vector_Vector_t _ ->
              failwith "print_current: V3: not implemented"

          | TensorVector_Vector_Vector _ ->
              failwith "print_current: V3: not implemented"

          | TensorVector_Vector_Vector_cf _ ->
              failwith "print_current: V3: not implemented"

          | TensorVector_Scalar_Scalar _ ->
              failwith "print_current: V3: not implemented"

          | TensorVector_Scalar_Scalar_cf _ ->
              failwith "print_current: V3: not implemented"

          | TensorScalar_Vector_Vector _ ->
              failwith "print_current: V3: not implemented"

          | TensorScalar_Vector_Vector_cf _ ->
              failwith "print_current: V3: not implemented"

          | TensorScalar_Scalar_Scalar _ ->
              failwith "print_current: V3: not implemented"

          | TensorScalar_Scalar_Scalar_cf _ ->
              failwith "print_current: V3: not implemented"

          | Dim6_Scalar_Vector_Vector_D _ ->
              failwith "print_current: V3: not implemented"

          | Dim6_Scalar_Vector_Vector_DP _ ->
              failwith "print_current: V3: not implemented"

          | Dim6_HAZ_D _ ->
              failwith "print_current: V3: not implemented"

          | Dim6_HAZ_DP _ ->
              failwith "print_current: V3: not implemented"

          | Dim6_HHH _ ->
              failwith "print_current: V3: not implemented"  

          | Dim6_Gauge_Gauge_Gauge_i _ ->
              failwith "print_current: V3: not implemented"  

          | Gauge_Gauge_Gauge_i _ ->
              failwith "print_current: V3: not implemented"

          | Dim6_GGG _ ->
              failwith "print_current: V3: not implemented"	

          | Dim6_AWW_DP _ ->
              failwith "print_current: V3: not implemented"

          | Dim6_AWW_DW _ ->
              failwith "print_current: V3: not implemented"

          | Dim6_WWZ_DPWDW _ ->
              failwith "print_current: V3: not implemented"
 
          | Dim6_WWZ_DW _ ->
              failwith "print_current: V3: not implemented"
 
          | Dim6_WWZ_D _ ->
              failwith "print_current: V3: not implemented"

          | Aux_Gauge_Gauge _ ->
              failwith "print_current: V3 (Aux_Gauge_Gauge): not implemented"

   end

(* Flip the sign in [c] to account for the~$\mathrm{i}^2$ relative to diagrams
   with only cubic couplings. *)
      | V4 (vertex, fusion, constant) ->
          let ch1, ch2, ch3 = children3 rhs in
          let wf1 = wf_index lookups.wfmap lookups.n_wfs (f ch1)
          and wf2 = wf_index lookups.wfmap lookups.n_wfs (f ch2)
          and wf3 = wf_index lookups.wfmap lookups.n_wfs (f ch3)
          (*i
          (*and p1 = mom_ID lookups.pmap ch1*)
          (*and p2 = mom_ID lookups.pmap ch2*)
          (*and p3 = mom_ID lookups.pmap ch2*)
          i*)
          and const_ID = get_const_ID lookups.cmap constant in
          let c =
            if (F.sign rhs) < 0 then const_ID else - const_ID in
          begin match vertex with
          | Scalar4 coeff ->
              printi ovm_FUSE_S_SSS ~lhs:lhs ~coupl:c ~coeff:coeff ~rhs1:wf1
                ~rhs2:wf2 ~rhs3:wf3
          | Scalar2_Vector2 coeff ->
              let printc code r1 r2 r3 = printi code
                ~lhs:lhs ~coupl:c ~coeff:coeff ~rhs1:r1 ~rhs2:r2 ~rhs3:r3 in
              begin match fusion with
              | F134 | F143 | F234 | F243 ->
                  printc ovm_FUSE_S_SVV wf1 wf2 wf3
              | F314 | F413 | F324 | F423 ->
                  printc ovm_FUSE_S_SVV wf2 wf1 wf3
              | F341 | F431 | F342 | F432 ->
                  printc ovm_FUSE_S_SVV wf3 wf1 wf2
              | F312 | F321 | F412 | F421 ->
                  printc ovm_FUSE_V_SSV wf2 wf3 wf1
              | F231 | F132 | F241 | F142 ->
                  printc ovm_FUSE_V_SSV wf1 wf3 wf2
              | F123 | F213 | F124 | F214 ->
                  printc ovm_FUSE_V_SSV wf1 wf2 wf3
              end

          | Vector4 contractions ->
              List.iter (print_vector4 c lhs wf1 wf2 wf3 fusion) contractions

          | Vector4_K_Matrix_tho _
          | Vector4_K_Matrix_jr _
          | Vector4_K_Matrix_cf_t0 _          
          | Vector4_K_Matrix_cf_t1 _
          | Vector4_K_Matrix_cf_t2 _
          | Vector4_K_Matrix_cf_t_rsi _
          | Vector4_K_Matrix_cf_m0 _
          | Vector4_K_Matrix_cf_m1 _
          | Vector4_K_Matrix_cf_m7 _          
          | DScalar2_Vector2_K_Matrix_ms _
          | DScalar2_Vector2_m_0_K_Matrix_cf _
          | DScalar2_Vector2_m_1_K_Matrix_cf _
          | DScalar2_Vector2_m_7_K_Matrix_cf _
          | DScalar4_K_Matrix_ms _ ->
              failwith "print_current: V4: K_Matrix not implemented"
          | Dim8_Scalar2_Vector2_1 _ 
          | Dim8_Scalar2_Vector2_2 _
          | Dim8_Scalar2_Vector2_m_0 _
          | Dim8_Scalar2_Vector2_m_1 _
          | Dim8_Scalar2_Vector2_m_7 _
          | Dim8_Scalar4 _ ->
              failwith "print_current: V4: not implemented"
          | Dim8_Vector4_t_0 _ ->
              failwith "print_current: V4: not implemented"
          | Dim8_Vector4_t_1 _ ->
              failwith "print_current: V4: not implemented"              
          | Dim8_Vector4_t_2 _ ->
              failwith "print_current: V4: not implemented"
          | Dim8_Vector4_m_0 _ ->
              failwith "print_current: V4: not implemented"
          | Dim8_Vector4_m_1 _ ->
              failwith "print_current: V4: not implemented"
          | Dim8_Vector4_m_7 _ ->
              failwith "print_current: V4: not implemented"    
          | GBBG _ ->
              failwith "print_current: V4: GBBG not implemented"
          | DScalar4 _
          | DScalar2_Vector2 _ ->
              failwith "print_current: V4: DScalars not implemented"
          | Dim6_H4_P2 _ ->  
              failwith "print_current: V4: not implemented"
          | Dim6_AHWW_DPB _ ->
              failwith "print_current: V4: not implemented"
          | Dim6_AHWW_DPW _ ->
              failwith "print_current: V4: not implemented"
          | Dim6_AHWW_DW _ ->
              failwith "print_current: V4: not implemented"
          | Dim6_Vector4_DW _ ->
              failwith "print_current: V4: not implemented"
          | Dim6_Vector4_W _ ->
              failwith "print_current: V4: not implemented"
          | Dim6_Scalar2_Vector2_D _ ->
              failwith "print_current: V4: not implemented"
          | Dim6_Scalar2_Vector2_DP _ ->
              failwith "print_current: V4: not implemented"
          | Dim6_HWWZ_DW _ ->
              failwith "print_current: V4: not implemented"
          | Dim6_HWWZ_DPB _ ->
              failwith "print_current: V4: not implemented"
          | Dim6_HWWZ_DDPW _ ->
              failwith "print_current: V4: not implemented"
          | Dim6_HWWZ_DPW _ ->
              failwith "print_current: V4: not implemented"
          | Dim6_AHHZ_D _ ->
              failwith "print_current: V4: not implemented"
          | Dim6_AHHZ_DP _ ->
              failwith "print_current: V4: not implemented"
          | Dim6_AHHZ_PB _ ->
              failwith "print_current: V4: not implemented"
          | Dim6_Scalar2_Vector2_PB _ ->           
              failwith "print_current: V4: not implemented"
          | Dim6_HHZZ_T _ ->   
              failwith "print_current: V4: not implemented"

          end

      | Vn (_, _, _) -> invalid_arg "Targets.print_current: n-ary fusion."

(* \thocwmodulesubsection{Fusions} *)

    let print_fusion lookups lhs_momID fusion amplitude =
      if F.on_shell amplitude (F.lhs fusion) then
        failwith "print_fusion: on_shell projectors not implemented!";
      if F.is_gauss amplitude (F.lhs fusion) then
        failwith "print_fusion: gauss amplitudes not implemented!";
      let lhs_wf = mult_wf lookups.dict amplitude (F.lhs fusion) in
      let lhs_wfID = wf_index lookups.wfmap lookups.n_wfs lhs_wf in
      let f = F.flavor (F.lhs fusion) in
      let pdg = CM.pdg f in
      let w =
        begin match CM.width f with
        | Vanishing | Fudged -> 0
        | Constant -> 1
        | Timelike -> 2
        | Complex_Mass -> 3
        | Running -> 4
        | Custom _ -> failwith "Targets.VM: custom width not available"
        end
      in
      let propagate code = printi code ~lhs:lhs_wfID ~rhs1:lhs_momID
        ~coupl:(abs(pdg)) ~coeff:w ~rhs4:(get_ID' amp_compare lookups.amap amplitude)
      in
      begin match CM.propagator f with
      | Prop_Scalar ->
          propagate ovm_PROPAGATE_SCALAR
      | Prop_Col_Scalar ->
          propagate ovm_PROPAGATE_COL_SCALAR
      | Prop_Ghost ->
          propagate ovm_PROPAGATE_GHOST
      | Prop_Spinor ->
          propagate ovm_PROPAGATE_SPINOR
      | Prop_ConjSpinor ->
          propagate ovm_PROPAGATE_CONJSPINOR
      | Prop_Majorana ->
          propagate ovm_PROPAGATE_MAJORANA
      | Prop_Col_Majorana ->
          propagate ovm_PROPAGATE_COL_MAJORANA
      | Prop_Unitarity ->
          propagate ovm_PROPAGATE_UNITARITY
      | Prop_Col_Unitarity ->
          propagate ovm_PROPAGATE_COL_UNITARITY
      | Prop_Feynman ->
          propagate ovm_PROPAGATE_FEYNMAN
      | Prop_Col_Feynman ->
          propagate ovm_PROPAGATE_COL_FEYNMAN
      | Prop_Vectorspinor ->
          propagate ovm_PROPAGATE_VECTORSPINOR
      | Prop_Tensor_2 ->
          propagate ovm_PROPAGATE_TENSOR2
      | Aux_Col_Scalar | Aux_Col_Vector | Aux_Col_Tensor_1 ->
          failwith "print_fusion: Aux_Col_* not implemented!"
      | Aux_Vector | Aux_Tensor_1 | Aux_Scalar | Aux_Spinor | Aux_ConjSpinor
      | Aux_Majorana | Only_Insertion ->
          propagate ovm_PROPAGATE_NONE
      | Prop_Gauge _ ->
          failwith "print_fusion: Prop_Gauge not implemented!"
      | Prop_Tensor_pure ->
          failwith "print_fusion: Prop_Tensor_pure not implemented!"
      | Prop_Vector_pure ->
          failwith "print_fusion: Prop_Vector_pure not implemented!"
      | Prop_Rxi _ ->
          failwith "print_fusion: Prop_Rxi not implemented!"
      | Prop_UFO _ ->
          failwith "print_fusion: Prop_UFO not implemented!"
      end;

(* Since the OVM knows that we want to propagate a wf, we can send the
   necessary fusions now. *)

      List.iter (print_current lookups lhs_wfID amplitude) (F.rhs fusion)

    let print_all_fusions lookups =
      let fusions = CF.fusions lookups.amplitudes in
      let fset = List.fold_left (fun s x -> FSet.add x s) FSet.empty fusions in
      ignore (List.fold_left (fun level (f, amplitude) ->
        let wf = F.lhs f in
        let lhs_momID = mom_ID lookups.pmap wf in
        let level' = List.length (F.momentum_list wf) in
        if (level' > level && level' > 2) then break ();
        print_fusion lookups lhs_momID f amplitude;
        level')
      1 (FSet.elements fset) )

(* \thocwmodulesubsection{Brakets} *)

    let print_braket lookups amplitude braket =
      let bra = F.bra braket
      and ket = F.ket braket in
      let braID = wf_index lookups.wfmap lookups.n_wfs
        (mult_wf lookups.dict amplitude bra) in
      List.iter (print_current lookups braID amplitude) ket

(* \begin{equation}
   \ii T = \ii^{\#\text{vertices}}\ii^{\#\text{propagators}} \cdots
         = \ii^{n-2}\ii^{n-3} \cdots
         = -\ii(-1)^n \cdots
   \end{equation} *)

(* All brakets for one cflow amplitude should be calculated by one
   thread to avoid multiple access on the same memory (amplitude).*)

    let print_brakets lookups (amplitude, i) =
      let n = List.length (F.externals amplitude) in
      let sign = if n mod 2 = 0 then -1 else 1
      and sym = F.symmetry amplitude in
      printi ovm_CALC_BRAKET ~lhs:i ~rhs1:sym ~coupl:sign;
      amplitude |> F.brakets |> List.iter (print_braket lookups amplitude)

(* Fortran arrays/OCaml lists start on 1/0. The amplitude list is sorted by
   [amp_compare] according to their color flows. In this way the amp array
   is sorted in the same way as [table_color_factors]. *)

    let print_all_brakets lookups =
      let g i elt = print_brakets lookups (elt, i+1) in
      lookups.amplitudes |> CF.processes |> List.sort amp_compare
                         |> ThoList.iteri g 0

(* \thocwmodulesubsection{Couplings} *)

(* For now we only care to catch the arrays [gncneu], [gnclep], [gncup] and
   [gncdown] of the SM. This will need an overhaul when it is clear how we store
   the type information of coupling constants. *)

    let strip_array_tag = function
      | Real_Array x -> x
      | Complex_Array x -> x

    let array_constants_list =
      let params = M.parameters()
      and strip_to_constant (lhs, _) = strip_array_tag lhs in
        List.map strip_to_constant params.derived_arrays

    let is_array x = List.mem x array_constants_list

    let constants_map =
      let first = fun (x, _, _) -> x in
      let second = fun (_, y, _) -> y in
      let third = fun (_, _, z) -> z in
      let v3 = List.map third (first (M.vertices () ))
      and v4 = List.map third (second (M.vertices () )) in
      let set = List.fold_left (fun s x -> CSet.add x s) CSet.empty (v3 @ v4) in
      let (arrays, singles) = CSet.partition is_array set in
        (singles |> CSet.elements |> map_of_list,
         arrays  |> CSet.elements |> map_of_list)

(* \thocwmodulesubsection{Output calls} *)

    let amplitudes_to_channel (cmdline : string) (oc : out_channel)
      (diagnostics : (diagnostic * bool) list ) (amplitudes : CF.amplitudes) =

      set_formatter_out_channel oc;
      if (num_particles amplitudes = 0) then begin
        print_description cmdline;
        print_zero_header (); nl ()
      end else begin
        let (wfset, amap) = wfset_amps amplitudes in
        let pset = expand_pset (momenta_set wfset)
        and n_wfs = num_wfs wfset in
        let wfmap = wf_map_of_list (WFSet.elements wfset)
        and pmap = map_of_list (ISet.elements pset)
        and cmap = constants_map in

        let lookups = {pmap = pmap; wfmap = wfmap; cmap = cmap; amap = amap;
          n_wfs = n_wfs; amplitudes = amplitudes;
          dict = CF.dictionary amplitudes} in

        print_description cmdline;
        print_header lookups wfset;
        print_spin_table amplitudes;
        print_flavor_tables amplitudes;
        print_color_tables amplitudes;
        printf "@\n%s" ("OVM instructions for momenta addition," ^
                        " fusions and brakets start here: ");
        break ();
        add_all_mom lookups pset;
        print_ext_amps lookups;
        break ();
        print_all_fusions lookups;
        break ();
        print_all_brakets lookups;
        break (); nl ();
        print_flush ()
      end

    let parameters_to_fortran oc _ =
     (*i The -params options is used as wrapper between OVM and Whizard. Most
       * trouble for the OVM comes from the array dimensionalities of couplings
       * but O'Mega should also know whether a constant is real or complex.
       * Hopefully all will be clearer with the fully general Lorentz structures
       * and UFO support. For now, we stick with this brute-force solution. i*)
      set_formatter_out_channel oc;
      let arrays_to_set = not (IMap.is_empty (snd constants_map)) in
      let set_coupl ty dim cmap = IMap.iter (fun key elt ->
        printf "    %s(%s%d) = %s" ty dim key (M.constant_symbol elt);
        nl () ) cmap in
      let declarations () =
        printf "  complex(%s), dimension(%d) :: ovm_coupl_cmplx"
          !kind (constants_map |> fst |> largest_key); nl ();
        if arrays_to_set then
          printf "  complex(%s), dimension(2, %d) :: ovm_coupl_cmplx2"
            !kind (constants_map |> snd |> largest_key); nl () in
      let print_line str = printf "%s" str; nl() in
      let print_md5sum = function
          | Some s ->
            print_line "  function md5sum ()";
            print_line "    character(len=32) :: md5sum";
            print_line ("    bytecode_file = '" ^ !bytecode_file ^ "'");
            print_line "    call initialize_vm (vm, bytecode_file)";
            print_line "    ! DON'T EVEN THINK of modifying the following line!";
            print_line ("    md5sum = '" ^ s ^ "'");
            print_line "  end function md5sum";
          | None -> ()
      in
      let print_inquiry_function_openmp () = begin
        print_line "  pure function openmp_supported () result (status)";
        print_line "    logical :: status";
        print_line ("    status = " ^ (if !openmp then ".true." else ".false."));
        print_line "  end function openmp_supported";
        nl ()
      end in
      let print_interface whizard =
      if whizard then begin
        print_line "  subroutine init (par, scheme)";
        print_line "    real(kind=default), dimension(*), intent(in) :: par";
        print_line "    integer, intent(in) :: scheme";
        print_line ("     bytecode_file = '" ^ !bytecode_file ^ "'");
        print_line "    call import_from_whizard (par, scheme)";
        print_line "    call initialize_vm (vm, bytecode_file)";
        print_line "  end subroutine init";
        nl ();
        print_line "  subroutine final ()";
        print_line "    call vm%final ()";
        print_line "  end subroutine final";
        nl ();
        print_line "  subroutine update_alpha_s (alpha_s)";
        print_line ("    real(kind=" ^ !kind ^ "), intent(in) :: alpha_s");
        print_line "    call model_update_alpha_s (alpha_s)";
        print_line "  end subroutine update_alpha_s";
        nl ()
      end
      else begin
        print_line "  subroutine init ()";
        print_line ("     bytecode_file = '" ^ !bytecode_file ^ "'");
        print_line "     call init_parameters ()";
        print_line "     call initialize_vm (vm, bytecode_file)";
        print_line "  end subroutine"
      end in
      let print_lookup_functions () = begin
        print_line "  pure function number_particles_in () result (n)";
        print_line "    integer :: n";
        print_line "    n = vm%number_particles_in ()";
        print_line "  end function number_particles_in";
        nl();
        print_line "  pure function number_particles_out () result (n)";
        print_line "    integer :: n";
        print_line "    n = vm%number_particles_out ()";
        print_line "  end function number_particles_out";
        nl();
        print_line "  pure function number_spin_states () result (n)";
        print_line "    integer :: n";
        print_line "    n = vm%number_spin_states ()";
        print_line "  end function number_spin_states";
        nl();
        print_line "  pure subroutine spin_states (a)";
        print_line "    integer, dimension(:,:), intent(out) :: a";
        print_line "    call vm%spin_states (a)";
        print_line "  end subroutine spin_states";
        nl();
        print_line "  pure function number_flavor_states () result (n)";
        print_line "    integer :: n";
        print_line "    n = vm%number_flavor_states ()";
        print_line "  end function number_flavor_states";
        nl();
        print_line "  pure subroutine flavor_states (a)";
        print_line "    integer, dimension(:,:), intent(out) :: a";
        print_line "    call vm%flavor_states (a)";
        print_line "  end subroutine flavor_states";
        nl();
        print_line "  pure function number_color_indices () result (n)";
        print_line "    integer :: n";
        print_line "    n = vm%number_color_indices ()";
        print_line "  end function number_color_indices";
        nl();
        print_line "  pure function number_color_flows () result (n)";
        print_line "    integer :: n";
        print_line "    n = vm%number_color_flows ()";
        print_line "  end function number_color_flows";
        nl();
        print_line "  pure subroutine color_flows (a, g)";
        print_line "    integer, dimension(:,:,:), intent(out) :: a";
        print_line "    logical, dimension(:,:), intent(out) :: g";
        print_line "    call vm%color_flows (a, g)";
        print_line "  end subroutine color_flows";
        nl();
        print_line "  pure function number_color_factors () result (n)";
        print_line "    integer :: n";
        print_line "    n = vm%number_color_factors ()";
        print_line "  end function number_color_factors";
        nl();
        print_line "  pure subroutine color_factors (cf)";
        print_line "    use omega_color";
        print_line "    type(omega_color_factor), dimension(:), intent(out) :: cf";
        print_line "    call vm%color_factors (cf)";
        print_line "  end subroutine color_factors";
        nl();
        print_line "  !pure unless OpenMP";
        print_line "  !pure function color_sum (flv, hel) result (amp2)";
        print_line "  function color_sum (flv, hel) result (amp2)";
        print_line "    use kinds";
        print_line "    integer, intent(in) :: flv, hel";
        print_line "    real(kind=default) :: amp2";
        print_line "    amp2 = vm%color_sum (flv, hel)";
        print_line "  end function color_sum";
        nl();
        print_line "  subroutine new_event (p)";
        print_line "    use kinds";
        print_line "    real(kind=default), dimension(0:3,*), intent(in) :: p";
        print_line "    call vm%new_event (p)";
        print_line "  end subroutine new_event";
        nl();
        print_line "  subroutine reset_helicity_selection (threshold, cutoff)";
        print_line "    use kinds";
        print_line "    real(kind=default), intent(in) :: threshold";
        print_line "    integer, intent(in) :: cutoff";
        print_line "    call vm%reset_helicity_selection (threshold, cutoff)";
        print_line "  end subroutine reset_helicity_selection";
        nl();
        print_line "  pure function is_allowed (flv, hel, col) result (yorn)";
        print_line "    logical :: yorn";
        print_line "    integer, intent(in) :: flv, hel, col";
        print_line "    yorn = vm%is_allowed (flv, hel, col)";
        print_line "  end function is_allowed";
        nl();
        print_line "  pure function get_amplitude (flv, hel, col) result (amp_result)";
        print_line "    use kinds";
        print_line "    complex(kind=default) :: amp_result";
        print_line "    integer, intent(in) :: flv, hel, col";
        print_line "    amp_result = vm%get_amplitude(flv, hel, col)";
        print_line "  end function get_amplitude";
        nl();
      end in
      print_line ("module " ^ !wrapper_module);
      print_line ("  use " ^ !parameter_module_external);
      print_line "  use iso_varying_string, string_t => varying_string";
      print_line "  use kinds";
      print_line "  use omegavm95";
      print_line "  implicit none";
      print_line "  private";
      print_line "  type(vm_t) :: vm";
      print_line "  type(string_t) :: bytecode_file";
      print_line ("  public :: number_particles_in, number_particles_out," ^
          " number_spin_states, &");
      print_line ("    spin_states, number_flavor_states, flavor_states," ^
          " number_color_indices, &");
      print_line ("    number_color_flows, color_flows," ^
          " number_color_factors, color_factors, &");
      print_line ("    color_sum, new_event, reset_helicity_selection," ^
          " is_allowed, get_amplitude, &");
      print_line ("    init, " ^ 
          (match !md5sum with Some _ -> "md5sum, "
                            | None -> "") ^ "openmp_supported");
      if !whizard then
        print_line ("  public :: final, update_alpha_s")
      else
        print_line ("  public :: initialize_vm");
      declarations ();
      print_line "contains";

      print_line "  subroutine setup_couplings ()";
      set_coupl "ovm_coupl_cmplx" "" (fst constants_map);
      if arrays_to_set then
        set_coupl "ovm_coupl_cmplx2" ":," (snd constants_map);
      print_line "  end subroutine setup_couplings";
      print_line "  subroutine initialize_vm (vm, bytecode_file)";
      print_line "    class(vm_t), intent(out) :: vm";
      print_line "    type(string_t), intent(in) :: bytecode_file";
      print_line "    type(string_t) :: version";
      print_line "    type(string_t) :: model";
      print_line ("    version = 'OVM " ^ version ^ "'");
      print_line ("    model = 'Model " ^ model_name ^ "'");
      print_line "    call setup_couplings ()";
      print_line "    call vm%init (bytecode_file, version, model, verbose=.False., &";
      print_line "      coupl_cmplx=ovm_coupl_cmplx, &";
      if arrays_to_set then
        print_line "      coupl_cmplx2=ovm_coupl_cmplx2, &";
      print_line ("      mass=mass, width=width, openmp=" ^ (if !openmp then
        ".true." else ".false.") ^ ")");
      print_line "  end subroutine initialize_vm";
      nl();
      print_md5sum !md5sum;
      print_inquiry_function_openmp ();
      print_interface !whizard;
      print_lookup_functions ();

      print_line ("end module " ^ !wrapper_module)

    let parameters_to_channel oc =
      parameters_to_fortran oc (CM.parameters ())

  end

(* \thocwmodulesection{\texttt{Fortran\,90/95}} *)

(* \thocwmodulesubsection{Dirac Fermions}
   We factor out the code for fermions so that we can use the simpler
   implementation for Dirac fermions if the model contains no Majorana
   fermions. *)

module type Fermions =
  sig
    open Coupling
    val psi_type : string
    val psibar_type : string
    val chi_type : string
    val grav_type : string
    val psi_incoming : string
    val brs_psi_incoming : string
    val psibar_incoming : string
    val brs_psibar_incoming : string
    val chi_incoming : string
    val brs_chi_incoming : string
    val grav_incoming : string
    val psi_outgoing : string
    val brs_psi_outgoing : string
    val psibar_outgoing : string
    val brs_psibar_outgoing : string
    val chi_outgoing : string
    val brs_chi_outgoing : string
    val grav_outgoing : string
    val psi_propagator : string
    val psibar_propagator : string
    val chi_propagator : string
    val grav_propagator : string
    val psi_projector : string
    val psibar_projector : string
    val chi_projector : string
    val grav_projector : string
    val psi_gauss : string
    val psibar_gauss : string
    val chi_gauss : string
    val grav_gauss : string
    val print_current : int * fermionbar * boson * fermion ->
      string -> string -> string -> fuse2 -> unit
    val print_current_mom : int * fermionbar * boson * fermion ->
      string -> string -> string -> string -> string -> string
      -> fuse2 -> unit
    val print_current_p : int * fermion * boson * fermion ->
      string -> string -> string -> fuse2 -> unit
    val print_current_b : int * fermionbar * boson * fermionbar ->
      string -> string -> string -> fuse2 -> unit
    val print_current_g : int * fermionbar * boson * fermion ->
      string -> string -> string -> string -> string -> string
      -> fuse2 -> unit
    val print_current_g4 : int * fermionbar * boson2 * fermion ->
      string -> string -> string -> string -> fuse3 -> unit
    val reverse_braket : bool -> lorentz -> lorentz list -> bool
    val use_module : string
    val require_library : string list
   end

module Fortran_Fermions : Fermions =
  struct
    open Coupling
    open Format

    let psi_type = "spinor"
    let psibar_type = "conjspinor"
    let chi_type = "???"
    let grav_type = "???"

    let psi_incoming = "u"
    let brs_psi_incoming = "brs_u"
    let psibar_incoming = "vbar"
    let brs_psibar_incoming = "brs_vbar"
    let chi_incoming = "???"
    let brs_chi_incoming = "???"
    let grav_incoming = "???"
    let psi_outgoing = "v"
    let brs_psi_outgoing = "brs_v"
    let psibar_outgoing = "ubar"
    let brs_psibar_outgoing = "brs_ubar"
    let chi_outgoing = "???"
    let brs_chi_outgoing = "???"
    let grav_outgoing = "???"

    let psi_propagator = "pr_psi"
    let psibar_propagator = "pr_psibar"
    let chi_propagator = "???"
    let grav_propagator = "???"

    let psi_projector = "pj_psi"
    let psibar_projector = "pj_psibar"
    let chi_projector = "???"
    let grav_projector = "???"

    let psi_gauss = "pg_psi"
    let psibar_gauss = "pg_psibar"
    let chi_gauss = "???"
    let grav_gauss = "???"

    let format_coupling coeff c =
      match coeff with
      | 1 -> c
      | -1 -> "(-" ^ c ^")"
      | coeff -> string_of_int coeff ^ "*" ^ c

    let format_coupling_2 coeff c =
      match coeff with
      | 1 -> c
      | -1 -> "-" ^ c
      | coeff -> string_of_int coeff ^ "*" ^ c

(* \begin{dubious}
     JR's coupling constant HACK, necessitated by tho's bad design descition.
   \end{dubious} *)

    let fastener s i ?p ?q () =
      try
        let offset = (String.index s '(') in
        if ((String.get s (String.length s - 1)) != ')') then
          failwith "fastener: wrong usage of parentheses"
        else
          let func_name = (String.sub s 0 offset) and
              tail =
            (String.sub s (succ offset) (String.length s - offset - 2)) in
          if (String.contains func_name ')') ||
             (String.contains tail '(') ||
             (String.contains tail ')') then
            failwith "fastener: wrong usage of parentheses"
          else
            func_name ^ "(" ^ string_of_int i ^ "," ^ tail ^ ")"
      with
      | Not_found ->
          if (String.contains s ')') then
            failwith "fastener: wrong usage of parentheses"
          else
            match p with
            | None   -> s ^ "(" ^ string_of_int i ^ ")"
            | Some p ->
              match q with
              | None   -> s ^ "(" ^ p ^ "*" ^ p ^ "," ^ string_of_int i ^ ")"
              | Some q -> s ^ "(" ^ p ^ "," ^ q ^ "," ^ string_of_int i ^ ")"

    let print_fermion_current coeff f c wf1 wf2 fusion =
      let c = format_coupling coeff c in
      match fusion with
      | F13 -> printf "%s_ff(%s,%s,%s)" f c wf1 wf2
      | F31 -> printf "%s_ff(%s,%s,%s)" f c wf2 wf1
      | F23 -> printf "f_%sf(%s,%s,%s)" f c wf1 wf2
      | F32 -> printf "f_%sf(%s,%s,%s)" f c wf2 wf1
      | F12 -> printf "f_f%s(%s,%s,%s)" f c wf1 wf2
      | F21 -> printf "f_f%s(%s,%s,%s)" f c wf2 wf1

(* \begin{dubious}
     Using a two element array for the combined vector-axial and scalar-pseudo
     couplings helps to support HELAS as well.  Since we will probably never
     support general boson couplings with HELAS, it might be retired in favor
     of two separate variables.  For this [Model.constant_symbol] has to be
     generalized.
   \end{dubious} *)

(* \begin{dubious}
     NB: passing the array instead of two separate constants would be a
     \emph{bad} idea, because the support for Majorana spinors below will
     have to flip signs!
   \end{dubious} *)

    let print_fermion_current2 coeff f c wf1 wf2 fusion =
      let c = format_coupling_2 coeff c in
      let c1 = fastener c 1 ()
      and c2 = fastener c 2 () in
      match fusion with
      | F13 -> printf "%s_ff(%s,%s,%s,%s)" f c1 c2 wf1 wf2
      | F31 -> printf "%s_ff(%s,%s,%s,%s)" f c1 c2 wf2 wf1
      | F23 -> printf "f_%sf(%s,%s,%s,%s)" f c1 c2 wf1 wf2
      | F32 -> printf "f_%sf(%s,%s,%s,%s)" f c1 c2 wf2 wf1
      | F12 -> printf "f_f%s(%s,%s,%s,%s)" f c1 c2 wf1 wf2
      | F21 -> printf "f_f%s(%s,%s,%s,%s)" f c1 c2 wf2 wf1

    let print_fermion_current_mom_v1 coeff f c wf1 wf2 p1 p2 p12 fusion =
      let c = format_coupling coeff c in
      let c1 = fastener c 1 and
          c2 = fastener c 2 in
      match fusion with
      | F13 -> printf "%s_ff(%s,%s,%s,%s)" f (c1 ~p:p12 ()) (c2 ~p:p12 ()) wf1 wf2
      | F31 -> printf "%s_ff(%s,%s,%s,%s)" f (c1 ~p:p12 ()) (c2 ~p:p12 ()) wf2 wf1
      | F23 -> printf "f_%sf(%s,%s,%s,%s)" f (c1 ~p:p1 ()) (c2 ~p:p1 ()) wf1 wf2
      | F32 -> printf "f_%sf(%s,%s,%s,%s)" f (c1 ~p:p2 ()) (c2 ~p:p2 ()) wf2 wf1
      | F12 -> printf "f_f%s(%s,%s,%s,%s)" f (c1 ~p:p2 ()) (c2 ~p:p2 ()) wf1 wf2
      | F21 -> printf "f_f%s(%s,%s,%s,%s)" f (c1 ~p:p1 ()) (c2 ~p:p1 ()) wf2 wf1

    let print_fermion_current_mom_v2 coeff f c wf1 wf2 p1 p2 p12 fusion =
      let c = format_coupling coeff c in
      let c1 = fastener c 1 and
          c2 = fastener c 2 in
      match fusion with
      | F13 -> printf "%s_ff(%s,%s,@,%s,%s,%s)" f (c1 ~p:p12 ()) (c2 ~p:p12 ()) wf1 wf2 p12
      | F31 -> printf "%s_ff(%s,%s,@,%s,%s,%s)" f (c1 ~p:p12 ()) (c2 ~p:p12 ()) wf2 wf1 p12
      | F23 -> printf "f_%sf(%s,%s,@,%s,%s,%s)" f (c1 ~p:p1 ()) (c2 ~p:p1 ()) wf1 wf2 p1
      | F32 -> printf "f_%sf(%s,%s,@,%s,%s,%s)" f (c1 ~p:p2 ()) (c2 ~p:p2 ()) wf2 wf1 p2
      | F12 -> printf "f_f%s(%s,%s,@,%s,%s,%s)" f (c1 ~p:p2 ()) (c2 ~p:p2 ()) wf1 wf2 p2
      | F21 -> printf "f_f%s(%s,%s,@,%s,%s,%s)" f (c1 ~p:p1 ()) (c2 ~p:p1 ()) wf2 wf1 p1

    let print_fermion_current_mom_ff coeff f c wf1 wf2 p1 p2 p12 fusion =
      let c = format_coupling coeff c in
      let c1 = fastener c 1 and
          c2 = fastener c 2 in
      match fusion with
      | F13 -> printf "%s_ff(%s,%s,%s,%s)" f (c1 ~p:p1 ~q:p2 ()) (c2 ~p:p1 ~q:p2 ()) wf1 wf2
      | F31 -> printf "%s_ff(%s,%s,%s,%s)" f (c1 ~p:p1 ~q:p2 ()) (c2 ~p:p1 ~q:p2 ()) wf2 wf1
      | F23 -> printf "f_%sf(%s,%s,%s,%s)" f (c1 ~p:p12 ~q:p2 ()) (c2 ~p:p12 ~q:p2 ()) wf1 wf2
      | F32 -> printf "f_%sf(%s,%s,%s,%s)" f (c1 ~p:p12 ~q:p1 ()) (c2 ~p:p12 ~q:p1 ()) wf2 wf1
      | F12 -> printf "f_f%s(%s,%s,%s,%s)" f (c1 ~p:p12 ~q:p1 ()) (c2 ~p:p12 ~q:p1 ()) wf1 wf2
      | F21 -> printf "f_f%s(%s,%s,%s,%s)" f (c1 ~p:p12 ~q:p2 ()) (c2 ~p:p12 ~q:p2 ()) wf2 wf1

    let print_current = function
      | coeff, Psibar, VA, Psi -> print_fermion_current2 coeff "va"
      | coeff, Psibar, VA2, Psi -> print_fermion_current coeff "va2"
      | coeff, Psibar, VA3, Psi -> print_fermion_current coeff "va3"
      | coeff, Psibar, V, Psi -> print_fermion_current coeff "v"
      | coeff, Psibar, A, Psi -> print_fermion_current coeff "a"
      | coeff, Psibar, VL, Psi -> print_fermion_current coeff "vl"
      | coeff, Psibar, VR, Psi -> print_fermion_current coeff "vr"
      | coeff, Psibar, VLR, Psi -> print_fermion_current2 coeff "vlr"
      | coeff, Psibar, SP, Psi -> print_fermion_current2 coeff "sp"
      | coeff, Psibar, S, Psi -> print_fermion_current coeff "s"
      | coeff, Psibar, P, Psi -> print_fermion_current coeff "p"
      | coeff, Psibar, SL, Psi -> print_fermion_current coeff "sl"
      | coeff, Psibar, SR, Psi -> print_fermion_current coeff "sr"
      | coeff, Psibar, SLR, Psi -> print_fermion_current2 coeff "slr"
      | _, Psibar, _, Psi -> invalid_arg
            "Targets.Fortran_Fermions: no superpotential here"
      | _, Chibar, _, _ | _, _, _, Chi -> invalid_arg
            "Targets.Fortran_Fermions: Majorana spinors not handled"
      | _, Gravbar, _, _ | _, _, _, Grav -> invalid_arg
            "Targets.Fortran_Fermions: Gravitinos not handled"

    let print_current_mom = function
      | coeff, Psibar, VLRM, Psi -> print_fermion_current_mom_v1 coeff "vlr"
      | coeff, Psibar, VAM, Psi -> print_fermion_current_mom_ff coeff "va"
      | coeff, Psibar, VA3M, Psi -> print_fermion_current_mom_ff coeff "va3"
      | coeff, Psibar, SPM, Psi -> print_fermion_current_mom_v1 coeff "sp"
      | coeff, Psibar, TVA, Psi -> print_fermion_current_mom_v1 coeff "tva"
      | coeff, Psibar, TVAM, Psi -> print_fermion_current_mom_v2 coeff "tvam"
      | coeff, Psibar, TLR, Psi -> print_fermion_current_mom_v1 coeff "tlr"
      | coeff, Psibar, TLRM, Psi -> print_fermion_current_mom_v2 coeff "tlrm"
      | coeff, Psibar, TRL, Psi -> print_fermion_current_mom_v1 coeff "trl"
      | coeff, Psibar, TRLM, Psi -> print_fermion_current_mom_v2 coeff "trlm"
      | _, Psibar, _, Psi -> invalid_arg
            "Targets.Fortran_Fermions: only sigma tensor coupling here"
      | _, Chibar, _, _ | _, _, _, Chi -> invalid_arg
            "Targets.Fortran_Fermions: Majorana spinors not handled"
      | _, Gravbar, _, _ | _, _, _, Grav -> invalid_arg
            "Targets.Fortran_Fermions: Gravitinos not handled"

    let print_current_p = function
      | _, _, _, _ -> invalid_arg
            "Targets.Fortran_Fermions: No clashing arrows here"

    let print_current_b = function
      | _, _, _, _ -> invalid_arg
            "Targets.Fortran_Fermions: No clashing arrows here"

    let print_current_g = function
      | _, _, _, _ -> invalid_arg
            "Targets.Fortran_Fermions: No gravitinos here"

    let print_current_g4 = function
      | _, _, _, _ -> invalid_arg
            "Targets.Fortran_Fermions: No gravitinos here"

    let reverse_braket vintage bra ket =
      match bra with
      | Spinor -> true
      | _ -> false

    let use_module = "omega95"
    let require_library =
      ["omega_spinors_2010_01_A"; "omega_spinor_cpls_2010_01_A"]
  end

(* \thocwmodulesubsection{Main Functor} *)

module Make_Fortran (Fermions : Fermions)
    (Fusion_Maker : Fusion.Maker) (P : Momentum.T) (M : Model.T) =
  struct

    let require_library =
      Fermions.require_library @
      [ "omega_vectors_2010_01_A"; "omega_polarizations_2010_01_A";
        "omega_couplings_2010_01_A"; "omega_color_2010_01_A";
        "omega_utils_2010_01_A" ]

    module CM = Colorize.It(M)
    module F = Fusion_Maker(P)(M)

    module CF = Fusion.Multi(Fusion_Maker)(P)(M)
    type amplitudes = CF.amplitudes

    open Coupling
    open Format

    type output_mode =
      | Single_Function
      | Single_Module of int
      | Single_File of int
      | Multi_File of int

    let line_length = ref 80
    let continuation_lines = ref (-1) (* 255 *)
    let kind = ref "default"
    let fortran95 = ref true
    let module_name = ref "omega_amplitude"
    let output_mode = ref (Single_Module 10)
    let use_modules = ref []
    let whizard = ref false
    let amp_triv = ref false
    let parameter_module = ref ""
    let md5sum = ref None
    let no_write = ref false
    let km_write = ref false
    let km_pure = ref false
    let km_2_write = ref false
    let km_2_pure = ref false
    let openmp = ref false
    let pure_unless_openmp = false

    let options = Options.create
      [ "90", Arg.Clear fortran95,
        "don't use Fortran95 features that are not in Fortran90";
        "kind", Arg.String (fun s -> kind := s),
        "real and complex kind (default: " ^ !kind ^ ")";
        "width", Arg.Int (fun w -> line_length := w), "maximum line length";
        "continuation", Arg.Int (fun l -> continuation_lines := l),
        "maximum # of continuation lines";
        "module", Arg.String (fun s -> module_name := s), "module name";
        "single_function", Arg.Unit (fun () -> output_mode := Single_Function),
        "compute the matrix element(s) in a monolithic function";
        "split_function", Arg.Int (fun n -> output_mode := Single_Module n),
        "split the matrix element(s) into small functions [default, size = 10]";
        "split_module", Arg.Int (fun n -> output_mode := Single_File n),
        "split the matrix element(s) into small modules";
        "split_file", Arg.Int (fun n -> output_mode := Multi_File n),
        "split the matrix element(s) into small files";
        "use", Arg.String (fun s -> use_modules := s :: !use_modules),
        "use module";
        "parameter_module", Arg.String (fun s -> parameter_module := s),
        "parameter_module";
        "md5sum", Arg.String (fun s -> md5sum := Some s),
        "transfer MD5 checksum";
        "whizard", Arg.Set whizard, "include WHIZARD interface";
	"amp_triv", Arg.Set amp_triv, "only print trivial amplitude";
        "no_write", Arg.Set no_write, "no 'write' statements";
        "kmatrix_write", Arg.Set km_2_write, "write K matrix functions";
        "kmatrix_2_write", Arg.Set km_write, "write K matrix 2 functions";
        "kmatrix_write_pure", Arg.Set km_pure, "write K matrix pure functions";
        "kmatrix_2_write_pure", Arg.Set km_2_pure, "write Kmatrix2pure functions";
        "openmp", Arg.Set openmp, "activate OpenMP support in generated code"]

(* Fortran style line continuation: *)
    let nl = Format_Fortran.newline

    let print_list = function
      | [] -> ()
      | a :: rest ->
          print_string a;
          List.iter (fun s -> printf ",@ %s" s) rest

(* \thocwmodulesubsection{Variables and Declarations} *)

    (* ["NC"] is already used up in the module ["constants"]: *)
    let nc_parameter = "N_"
    let omega_color_factor_abbrev = "OCF"
    let openmp_tld_type = "thread_local_data"
    let openmp_tld = "tld"

    let flavors_symbol ?(decl = false) flavors =
      (if !openmp && not decl then openmp_tld ^ "%" else "" ) ^
      "oks_" ^ String.concat "" (List.map CM.flavor_symbol flavors)

    let p2s p =
      if p >= 0 && p <= 9 then
        string_of_int p
      else if p <= 36 then
        String.make 1 (Char.chr (Char.code 'A' + p - 10))
      else
        "_"

    let format_momentum p =
      "p" ^ String.concat "" (List.map p2s p)

    let format_p wf =
      String.concat "" (List.map p2s (F.momentum_list wf))

    let ext_momentum wf =
      match F.momentum_list wf with
      | [n] -> n
      | _ -> invalid_arg "Targets.Fortran.ext_momentum"

    module PSet = Set.Make (struct type t = int list let compare = compare end)
    module WFSet = Set.Make (struct type t = F.wf let compare = compare end)

    let add_tag wf name =
      match F.wf_tag wf with
      | None -> name
      | Some tag -> name ^ "_" ^ tag

    let variable ?(decl = false) wf =
      (if !openmp && not decl then openmp_tld ^ "%" else "")
      ^ add_tag wf ("owf_" ^ CM.flavor_symbol (F.flavor wf) ^ "_" ^ format_p wf)

    let momentum wf = "p" ^ format_p wf
    let spin wf = "s(" ^ string_of_int (ext_momentum wf) ^ ")"

    let format_multiple_variable ?(decl = false) wf i =
      variable ~decl wf ^ "_X" ^ string_of_int i

    let multiple_variable ?(decl = false) amplitude dictionary wf =
      try
        format_multiple_variable ~decl wf (dictionary amplitude wf)
      with
      | Not_found -> variable wf

    let multiple_variables ?(decl = false) multiplicity wf =
      try
        List.map
          (format_multiple_variable ~decl wf)
          (ThoList.range 1 (multiplicity wf))
      with
      | Not_found -> [variable ~decl wf]

    let declaration_chunk_size = 64

    let declare_list_chunk multiplicity t = function
      | [] -> ()
      | wfs ->
          printf "    @[<2>%s :: " t;
          print_list (ThoList.flatmap (multiple_variables ~decl:true multiplicity) wfs); nl ()

    let declare_list multiplicity t = function
      | [] -> ()
      | wfs ->
          List.iter
            (declare_list_chunk multiplicity t)
            (ThoList.chopn declaration_chunk_size wfs)

    type declarations =
        { scalars : F.wf list;
          spinors : F.wf list;
          conjspinors : F.wf list;
          realspinors : F.wf list;
          ghostspinors : F.wf list;
          vectorspinors : F.wf list;
          vectors : F.wf list;
          ward_vectors : F.wf list;
          massive_vectors : F.wf list;
          tensors_1 : F.wf list;
          tensors_2 : F.wf list;
          brs_scalars : F.wf list;
          brs_spinors : F.wf list;
          brs_conjspinors : F.wf list;
          brs_realspinors : F.wf list;
          brs_vectorspinors : F.wf list;
          brs_vectors : F.wf list;
          brs_massive_vectors : F.wf list }

    let rec classify_wfs' acc = function
      | [] -> acc
      | wf :: rest ->
          classify_wfs'
            (match CM.lorentz (F.flavor wf) with
            | Scalar -> {acc with scalars = wf :: acc.scalars}
            | Spinor -> {acc with spinors = wf :: acc.spinors}
            | ConjSpinor -> {acc with conjspinors = wf :: acc.conjspinors}
            | Majorana -> {acc with realspinors = wf :: acc.realspinors}
            | Maj_Ghost -> {acc with ghostspinors = wf :: acc.ghostspinors}
            | Vectorspinor ->
                {acc with vectorspinors = wf :: acc.vectorspinors}
            | Vector -> {acc with vectors = wf :: acc.vectors}
(*i            | Ward_Vector -> {acc with ward_vectors = wf :: acc.ward_vectors}
i*)
            | Massive_Vector ->
                {acc with massive_vectors = wf :: acc.massive_vectors}
            | Tensor_1 -> {acc with tensors_1 = wf :: acc.tensors_1}
            | Tensor_2 -> {acc with tensors_2 = wf :: acc.tensors_2}
            | BRS Scalar -> {acc with brs_scalars = wf :: acc.brs_scalars}
            | BRS Spinor -> {acc with brs_spinors = wf :: acc.brs_spinors}
            | BRS ConjSpinor -> {acc with brs_conjspinors =
                                 wf :: acc.brs_conjspinors}
            | BRS Majorana -> {acc with brs_realspinors =
                               wf :: acc.brs_realspinors}
            | BRS Vectorspinor -> {acc with brs_vectorspinors =
                                   wf :: acc.brs_vectorspinors}
            | BRS Vector -> {acc with brs_vectors = wf :: acc.brs_vectors}
            | BRS Massive_Vector -> {acc with brs_massive_vectors =
                                     wf :: acc.brs_massive_vectors}
            | BRS _ -> invalid_arg "Targets.wfs_classify': not needed here")
            rest

    let classify_wfs wfs = classify_wfs'
        { scalars = []; spinors = []; conjspinors = []; realspinors = [];
          ghostspinors = []; vectorspinors = []; vectors = [];
          ward_vectors = [];
          massive_vectors = []; tensors_1 = []; tensors_2 = [];
          brs_scalars = [] ; brs_spinors = []; brs_conjspinors = [];
          brs_realspinors = []; brs_vectorspinors = [];
          brs_vectors = []; brs_massive_vectors = []}
        wfs

(* \thocwmodulesubsection{Parameters} *)

    type 'a parameters =
        { real_singles : 'a list;
          real_arrays : ('a * int) list;
          complex_singles : 'a list;
          complex_arrays : ('a * int) list }

    let rec classify_singles acc = function
      | [] -> acc
      | Real p :: rest -> classify_singles
            { acc with real_singles = p :: acc.real_singles } rest
      | Complex p :: rest -> classify_singles
            { acc with complex_singles = p :: acc.complex_singles } rest

    let rec classify_arrays acc = function
      | [] -> acc
      | (Real_Array p, rhs) :: rest -> classify_arrays
            { acc with real_arrays =
              (p, List.length rhs) :: acc.real_arrays } rest
      | (Complex_Array p, rhs) :: rest -> classify_arrays
            { acc with complex_arrays =
              (p, List.length rhs) :: acc.complex_arrays } rest

    let classify_parameters params =
      classify_arrays
        (classify_singles
           { real_singles = [];
             real_arrays = [];
             complex_singles = [];
             complex_arrays = [] }
           (List.map fst params.derived)) params.derived_arrays

    let schisma = ThoList.chopn

    let schisma_num i n l =
      ThoList.enumerate i (schisma n l)

    let declare_parameters' t = function
      | [] -> ()
      | plist ->
          printf "  @[<2>%s(kind=%s), public, save :: " t !kind;
          print_list (List.map CM.constant_symbol plist); nl ()

    let declare_parameters t plist =
      List.iter (declare_parameters' t) plist

    let declare_parameter_array t (p, n) =
      printf "  @[<2>%s(kind=%s), dimension(%d), public, save :: %s"
        t !kind n (CM.constant_symbol p); nl ()

    (* NB: we use [string_of_float] to make sure that a decimal
       point is included to make Fortran compilers happy. *)
    let default_parameter (x, v) =
      printf "@ %s = %s_%s" (CM.constant_symbol x) (string_of_float v) !kind

    let declare_default_parameters t = function
      | [] -> ()
      | p :: plist ->
          printf "  @[<2>%s(kind=%s), public, save ::" t !kind;
          default_parameter p;
          List.iter (fun p' -> printf ","; default_parameter p') plist;
          nl ()

    let format_constant = function
      | I -> "(0,1)"
      | Integer c ->
         if c < 0 then
           sprintf "(%d.0_%s)" c !kind
         else
           sprintf "%d.0_%s" c !kind
      | Float x ->
         if x < 0. then
           "(" ^ string_of_float x ^ "_" ^ !kind ^ ")"
         else
           string_of_float x ^ "_" ^ !kind
      | _ -> invalid_arg "format_constant"

    let rec eval_parameter' = function
      | (I | Integer _ | Float _) as c ->
         printf "%s" (format_constant c)
      | Atom x -> printf "%s" (CM.constant_symbol x)
      | Sum [] -> printf "0.0_%s" !kind
      | Sum [x] -> eval_parameter' x
      | Sum (x :: xs) ->
          printf "@,("; eval_parameter' x;
          List.iter (fun x -> printf "@, + "; eval_parameter' x) xs;
          printf ")"
      | Diff (x, y) ->
          printf "@,("; eval_parameter' x;
          printf " - "; eval_parameter' y; printf ")"
      | Neg x -> printf "@,( - "; eval_parameter' x; printf ")"
      | Prod [] -> printf "1.0_%s" !kind
      | Prod [x] -> eval_parameter' x
      | Prod (x :: xs) ->
          printf "@,("; eval_parameter' x;
          List.iter (fun x -> printf " * "; eval_parameter' x) xs;
          printf ")"
      | Quot (x, y) ->
          printf "@,("; eval_parameter' x;
          printf " / "; eval_parameter' y; printf ")"
      | Rec x ->
          printf "@, (1.0_%s / " !kind; eval_parameter' x; printf ")"
      | Pow (x, n) ->
         printf "@,("; eval_parameter' x; 
         if n < 0 then
           printf "**(%d)" n
         else
           printf "**%d" n;
         printf ")"
      | PowX (x, y) ->
          printf "@,("; eval_parameter' x;
           printf "**"; eval_parameter' y; printf ")"
      | Sqrt x -> printf "@,sqrt ("; eval_parameter' x; printf ")"
      | Sin x -> printf "@,sin ("; eval_parameter' x; printf ")"
      | Cos x -> printf "@,cos ("; eval_parameter' x; printf ")"
      | Tan x -> printf "@,tan ("; eval_parameter' x; printf ")"
      | Cot x -> printf "@,cot ("; eval_parameter' x; printf ")"
      | Asin x -> printf "@,asin ("; eval_parameter' x; printf ")"
      | Acos x -> printf "@,acos ("; eval_parameter' x; printf ")"
      | Atan x -> printf "@,atan ("; eval_parameter' x; printf ")"
      | Atan2 (y, x) -> printf "@,atan2 ("; eval_parameter' y;
          printf ",@ "; eval_parameter' x; printf ")"
      | Sinh x -> printf "@,sinh ("; eval_parameter' x; printf ")"
      | Cosh x -> printf "@,cosh ("; eval_parameter' x; printf ")"
      | Tanh x -> printf "@,tanh ("; eval_parameter' x; printf ")"
      | Exp x -> printf "@,exp ("; eval_parameter' x; printf ")"
      | Log x -> printf "@,log ("; eval_parameter' x; printf ")"
      | Log10 x -> printf "@,log10 ("; eval_parameter' x; printf ")"
      | Conj (Integer _ | Float _ as x) -> eval_parameter' x
      | Conj x -> printf "@,cconjg ("; eval_parameter' x; printf ")"
      | Abs x -> printf "@,abs ("; eval_parameter' x; printf ")"

    let strip_single_tag = function
      | Real x -> x
      | Complex x -> x

    let strip_array_tag = function
      | Real_Array x -> x
      | Complex_Array x -> x

    let eval_parameter (lhs, rhs) =
      let x = CM.constant_symbol (strip_single_tag lhs) in
      printf "    @[<2>%s = " x; eval_parameter' rhs; nl ()

    let eval_para_list n l =
      printf "  subroutine setup_parameters_%03d ()" n; nl ();
      List.iter eval_parameter l;
      printf "  end subroutine setup_parameters_%03d" n; nl ()

    let eval_parameter_pair (lhs, rhs) =
      let x = CM.constant_symbol (strip_array_tag lhs) in
      let _ = List.fold_left (fun i rhs' ->
        printf "    @[<2>%s(%d) = " x i; eval_parameter' rhs'; nl ();
        succ i) 1 rhs in
      ()

    let eval_para_pair_list n l =
      printf "  subroutine setup_parameters_%03d ()" n; nl ();
      List.iter eval_parameter_pair l;
      printf "  end subroutine setup_parameters_%03d" n; nl ()

    let print_echo fmt p =
      let s = CM.constant_symbol p in
      printf "    write (unit = *, fmt = fmt_%s) \"%s\", %s"
        fmt s s; nl ()

    let print_echo_array fmt (p, n) =
      let s = CM.constant_symbol p in
      for i = 1 to n do
        printf "    write (unit = *, fmt = fmt_%s_array) " fmt ;
        printf "\"%s\", %d, %s(%d)" s i s i; nl ()
      done

    let contains params couplings =
      List.exists
        (fun (name, _) -> List.mem (CM.constant_symbol name) params)
        couplings.input

    let rec depends_on params = function
      | I | Integer _ | Float _ -> false
      | Atom name -> List.mem (CM.constant_symbol name) params
      | Sum es | Prod es ->
         List.exists (depends_on params) es
      | Diff (e1, e2) | Quot (e1, e2) | PowX (e1, e2) ->
         depends_on params e1 || depends_on params e2
      | Neg e | Rec e | Pow (e, _) ->
         depends_on params e
      | Sqrt e | Exp e | Log e | Log10 e
      | Sin e | Cos e | Tan e | Cot e
      | Asin e | Acos e | Atan e
      | Sinh e | Cosh e | Tanh e
      | Conj e | Abs e ->
         depends_on params e
      | Atan2 (e1, e2) ->
         depends_on params e1 || depends_on params e2

    let dependencies params couplings =
      if contains params couplings then
        List.rev
          (fst (List.fold_left
                  (fun (deps, plist) (param, v) ->
                    match param with
                    | Real name | Complex name ->
                       if depends_on plist v then
                         ((param, v) :: deps, CM.constant_symbol name :: plist)
                       else
                         (deps, plist))
                  ([], params) couplings.derived))
      else
        []

    let dependencies_arrays params couplings =
      if contains params couplings then
        List.rev
          (fst (List.fold_left
                  (fun (deps, plist) (param, vlist) ->
                    match param with
                    | Real_Array name | Complex_Array name ->
                       if List.exists (depends_on plist) vlist then
                         ((param, vlist) :: deps,
                          CM.constant_symbol name :: plist)
                       else
                         (deps, plist))
                  ([], params) couplings.derived_arrays))
      else
        []

    let parameters_to_fortran oc params =
      Format_Fortran.set_formatter_out_channel ~width:!line_length oc;
      let declarations = classify_parameters params in
      printf "module %s" !parameter_module; nl ();
      printf "  use kinds"; nl ();
      printf "  use constants"; nl ();
      printf "  implicit none"; nl ();
      printf "  private"; nl ();
      printf "  @[<2>public :: setup_parameters";
      printf ",@ import_from_whizard";
      printf ",@ model_update_alpha_s";
      if !no_write then begin
        printf "! No print_parameters";
      end else begin
        printf ",@ print_parameters";
      end; nl ();
      declare_default_parameters "real" params.input;
      declare_parameters "real" (schisma 69 declarations.real_singles);
      List.iter (declare_parameter_array "real") declarations.real_arrays;
      declare_parameters "complex" (schisma 69 declarations.complex_singles);
      List.iter (declare_parameter_array "complex") declarations.complex_arrays;
      printf "  interface cconjg"; nl ();
      printf "    module procedure cconjg_real, cconjg_complex"; nl ();
      printf "  end interface"; nl ();
      printf "  private :: cconjg_real, cconjg_complex"; nl ();
      printf "contains"; nl ();
      printf "  function cconjg_real (x) result (xc)"; nl ();
      printf "    real(kind=default), intent(in) :: x"; nl ();
      printf "    real(kind=default) :: xc"; nl ();
      printf "    xc = x"; nl ();
      printf "  end function cconjg_real"; nl ();
      printf "  function cconjg_complex (z) result (zc)"; nl ();
      printf "    complex(kind=default), intent(in) :: z"; nl ();
      printf "    complex(kind=default) :: zc"; nl ();
      printf "    zc = conjg (z)"; nl ();
      printf "  end function cconjg_complex"; nl ();
      printf "  ! derived parameters:"; nl ();
      let shredded = schisma_num 1 120 params.derived in
      let shredded_arrays = schisma_num 1 120 params.derived_arrays in
      let num_sub = List.length shredded in
      let num_sub_arrays = List.length shredded_arrays in
      List.iter (fun (i,l) -> eval_para_list i l) shredded;
      List.iter (fun (i,l) -> eval_para_pair_list (num_sub + i) l)
        shredded_arrays;
      printf "  subroutine setup_parameters ()"; nl ();
      for i = 1 to num_sub + num_sub_arrays do
        printf "    call setup_parameters_%03d ()" i; nl ();
      done;
      printf "  end subroutine setup_parameters"; nl ();
      printf "  subroutine import_from_whizard (par_array, scheme)"; nl ();
      printf
        "    real(%s), dimension(%d), intent(in) :: par_array"
        !kind (List.length params.input); nl ();
      printf "    integer, intent(in) :: scheme"; nl ();
      let i = ref 1 in
      List.iter
        (fun (p, _) ->
          printf "    %s = par_array(%d)" (CM.constant_symbol p) !i; nl ();
          incr i)
        params.input;
      printf "    call setup_parameters ()"; nl ();
      printf "  end subroutine import_from_whizard"; nl ();
      printf "  subroutine model_update_alpha_s (alpha_s)"; nl ();
      printf "    real(%s), intent(in) :: alpha_s" !kind; nl ();
      begin match (dependencies ["aS"] params,
                   dependencies_arrays ["aS"] params) with
      | [], [] ->
         printf "    ! 'aS' not among the input parameters"; nl ();
      | deps, deps_arrays ->
         printf "    aS = alpha_s"; nl ();
         List.iter eval_parameter deps;
         List.iter eval_parameter_pair deps_arrays
      end;
      printf "  end subroutine model_update_alpha_s"; nl ();
      if !no_write then begin
        printf "! No print_parameters"; nl ();
      end else begin
        printf "  subroutine print_parameters ()"; nl ();
        printf "    @[<2>character(len=*), parameter ::";
        printf "@ fmt_real = \"(A12,4X,' = ',E25.18)\",";
        printf "@ fmt_complex = \"(A12,4X,' = ',E25.18,' + i*',E25.18)\",";
        printf "@ fmt_real_array = \"(A12,'(',I2.2,')',' = ',E25.18)\",";
        printf "@ fmt_complex_array = ";
        printf "\"(A12,'(',I2.2,')',' = ',E25.18,' + i*',E25.18)\""; nl ();
        printf "    @[<2>write (unit = *, fmt = \"(A)\") @,";
        printf "\"default values for the input parameters:\""; nl ();
        List.iter (fun (p, _) -> print_echo "real" p) params.input;
        printf "    @[<2>write (unit = *, fmt = \"(A)\") @,";
        printf "\"derived parameters:\""; nl ();
        List.iter (print_echo "real") declarations.real_singles;
        List.iter (print_echo "complex") declarations.complex_singles;
        List.iter (print_echo_array "real") declarations.real_arrays;
        List.iter (print_echo_array "complex") declarations.complex_arrays;
        printf "  end subroutine print_parameters"; nl ();
      end;
      printf "end module %s" !parameter_module; nl ()

(* \thocwmodulesubsection{Run-Time Diagnostics} *)

    type diagnostic = All | Arguments | Momenta | Gauge

    type diagnostic_mode = Off | Warn | Panic

    let warn mode =
      match !mode with
      | Off -> false
      | Warn -> true
      | Panic -> true

    let panic mode =
      match !mode with
      | Off -> false
      | Warn -> false
      | Panic -> true

    let suffix mode =
      if panic mode then
        "panic"
      else
        "warn"

    let diagnose_arguments = ref Off
    let diagnose_momenta = ref Off
    let diagnose_gauge = ref Off

    let rec parse_diagnostic = function
      | All, panic ->
          parse_diagnostic (Arguments, panic);
          parse_diagnostic (Momenta, panic);
          parse_diagnostic (Gauge, panic)
      | Arguments, panic ->
          diagnose_arguments := if panic then Panic else Warn
      | Momenta, panic ->
          diagnose_momenta := if panic then Panic else Warn
      | Gauge, panic ->
          diagnose_gauge := if panic then Panic else Warn

(* If diagnostics are required, we have to switch off
   Fortran95 features like pure functions. *)

    let parse_diagnostics = function
      | [] -> ()
      | diagnostics ->
          fortran95 := false;
          List.iter parse_diagnostic diagnostics

(* \thocwmodulesubsection{Amplitude} *)

    let declare_momenta_chunk = function
      | [] -> ()
      | momenta ->
          printf "    @[<2>type(momentum) :: ";
          print_list (List.map format_momentum momenta); nl ()

    let declare_momenta = function
      | [] -> ()
      | momenta ->
          List.iter
            declare_momenta_chunk
            (ThoList.chopn declaration_chunk_size momenta)

    let declare_wavefunctions multiplicity wfs =
      let wfs' = classify_wfs wfs in
      declare_list multiplicity ("complex(kind=" ^ !kind ^ ")")
        (wfs'.scalars @ wfs'.brs_scalars);
      declare_list multiplicity ("type(" ^ Fermions.psi_type ^ ")")
        (wfs'.spinors @ wfs'.brs_spinors);
      declare_list multiplicity ("type(" ^ Fermions.psibar_type ^ ")")
        (wfs'.conjspinors @ wfs'.brs_conjspinors);
      declare_list multiplicity ("type(" ^ Fermions.chi_type ^ ")")
        (wfs'.realspinors @ wfs'.brs_realspinors @ wfs'.ghostspinors);
      declare_list multiplicity ("type(" ^ Fermions.grav_type ^ ")") wfs'.vectorspinors;
      declare_list multiplicity "type(vector)" (wfs'.vectors @ wfs'.massive_vectors @
         wfs'.brs_vectors @ wfs'.brs_massive_vectors @ wfs'.ward_vectors);
      declare_list multiplicity "type(tensor2odd)" wfs'.tensors_1;
      declare_list multiplicity "type(tensor)" wfs'.tensors_2

    let flavors a = F.incoming a @ F.outgoing a

    let declare_brakets_chunk = function
      | [] -> ()
      | amplitudes ->
          printf "    @[<2>complex(kind=%s) :: " !kind;
          print_list (List.map (fun a -> flavors_symbol ~decl:true (flavors a)) amplitudes); nl ()

    let declare_brakets = function
      | [] -> ()
      | amplitudes ->
          List.iter
            declare_brakets_chunk
            (ThoList.chopn declaration_chunk_size amplitudes)

    let print_variable_declarations amplitudes =
      let multiplicity = CF.multiplicity amplitudes
      and processes = CF.processes amplitudes in
      if not !amp_triv then begin
	declare_momenta
          (PSet.elements
             (List.fold_left
		(fun set a ->
                  PSet.union set (List.fold_right
                                    (fun wf -> PSet.add (F.momentum_list wf))
                                    (F.externals a) PSet.empty))
		PSet.empty processes));
	declare_momenta
          (PSet.elements
             (List.fold_left
		(fun set a ->
                  PSet.union set (List.fold_right
                                    (fun wf -> PSet.add (F.momentum_list wf))
                                    (F.variables a) PSet.empty))
		PSet.empty processes));
	if !openmp then begin
          printf "  type %s@[<2>" openmp_tld_type;
          nl ();
	end ;
	declare_wavefunctions multiplicity
          (WFSet.elements
             (List.fold_left
		(fun set a ->
                  WFSet.union set (List.fold_right WFSet.add (F.externals a) WFSet.empty))
		WFSet.empty processes));
	declare_wavefunctions multiplicity
          (WFSet.elements
             (List.fold_left
		(fun set a ->
                  WFSet.union set (List.fold_right WFSet.add (F.variables a) WFSet.empty))
		WFSet.empty processes));
	declare_brakets processes;
	if !openmp then begin
          printf "@]  end type %s\n" openmp_tld_type;
          printf "  type(%s) :: %s" openmp_tld_type openmp_tld;
          nl ();
	end;
      end

(* [print_current] is the most important function that has to match the functions
   in \verb+omega95+ (see appendix~\ref{sec:fortran}).  It offers plentiful
   opportunities for making mistakes, in particular those related to signs.
   We start with a few auxiliary functions:  *)

    let children2 rhs =
      match F.children rhs with
      | [wf1; wf2] -> (wf1, wf2)
      | _ -> failwith "Targets.children2: can't happen"

    let children3 rhs =
      match F.children rhs with
      | [wf1; wf2; wf3] -> (wf1, wf2, wf3)
      | _ -> invalid_arg "Targets.children3: can't happen"

(* Note that it is (marginally) faster to multiply the two scalar products
   with the coupling constant than the four vector components.
   \begin{dubious}
     This could be part of \verb+omegalib+ as well \ldots
   \end{dubious} *)

    let format_coeff = function
      | 1 -> ""
      | -1 -> "-"
      | coeff -> "(" ^ string_of_int coeff ^ ")*"

    let format_coupling coeff c =
      match coeff with
      | 1 -> c
      | -1 -> "(-" ^ c ^")"
      | coeff -> string_of_int coeff ^ "*" ^ c

(* \begin{dubious}
     The following is error prone and should be generated automagically.
   \end{dubious} *)

    let print_vector4 c wf1 wf2 wf3 fusion (coeff, contraction) =
      match contraction, fusion with
      | C_12_34, (F341|F431|F342|F432|F123|F213|F124|F214)
      | C_13_42, (F241|F421|F243|F423|F132|F312|F134|F314)
      | C_14_23, (F231|F321|F234|F324|F142|F412|F143|F413) ->
          printf "((%s%s)*(%s*%s))*%s" (format_coeff coeff) c wf1 wf2 wf3
      | C_12_34, (F134|F143|F234|F243|F312|F321|F412|F421)
      | C_13_42, (F124|F142|F324|F342|F213|F231|F413|F431)
      | C_14_23, (F123|F132|F423|F432|F214|F241|F314|F341) ->
          printf "((%s%s)*(%s*%s))*%s" (format_coeff coeff) c wf2 wf3 wf1
      | C_12_34, (F314|F413|F324|F423|F132|F231|F142|F241)
      | C_13_42, (F214|F412|F234|F432|F123|F321|F143|F341)
      | C_14_23, (F213|F312|F243|F342|F124|F421|F134|F431) ->
          printf "((%s%s)*(%s*%s))*%s" (format_coeff coeff) c wf1 wf3 wf2
          
    let print_vector4_t_0 c wf1 p1 wf2 p2 wf3 p3 fusion (coeff, contraction) =
      match contraction, fusion with
      | C_12_34, (F234|F243|F134|F143|F421|F321|F412|F312)
      | C_13_42, (F324|F342|F124|F142|F431|F231|F413|F213)
      | C_14_23, (F423|F432|F123|F132|F341|F241|F314|F214) ->
          printf "g_dim8g3_t_0(%s,%s,%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2 wf3 p3
      | C_12_34, (F324|F314|F423|F413|F142|F132|F241|F231)
      | C_13_42, (F234|F214|F432|F412|F143|F123|F341|F321)
      | C_14_23, (F243|F213|F342|F312|F134|F124|F431|F421) ->
          printf "g_dim8g3_t_0(%s,%s,%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1 wf3 p3
      | C_12_34, (F342|F341|F432|F431|F124|F123|F214|F213)
      | C_13_42, (F243|F241|F423|F421|F134|F132|F314|F312)
      | C_14_23, (F234|F231|F324|F321|F143|F142|F413|F412) ->
          printf "g_dim8g3_t_0(%s,%s,%s,%s,%s,%s,%s)" c wf3 p3 wf1 p1 wf2 p2

    let print_vector4_t_1 c wf1 p1 wf2 p2 wf3 p3 fusion (coeff, contraction) =
      match contraction, fusion with
      | C_12_34, (F234|F243|F134|F143|F421|F321|F412|F312)
      | C_13_42, (F324|F342|F124|F142|F431|F231|F413|F213)
      | C_14_23, (F423|F432|F123|F132|F341|F241|F314|F214) ->
          printf "g_dim8g3_t_1(%s,%s,%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2 wf3 p3
      | C_12_34, (F324|F314|F423|F413|F142|F132|F241|F231)
      | C_13_42, (F234|F214|F432|F412|F143|F123|F341|F321)
      | C_14_23, (F243|F213|F342|F312|F134|F124|F431|F421) ->
          printf "g_dim8g3_t_1(%s,%s,%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1 wf3 p3
      | C_12_34, (F342|F341|F432|F431|F124|F123|F214|F213)
      | C_13_42, (F243|F241|F423|F421|F134|F132|F314|F312)
      | C_14_23, (F234|F231|F324|F321|F143|F142|F413|F412) ->
          printf "g_dim8g3_t_1(%s,%s,%s,%s,%s,%s,%s)" c wf3 p3 wf1 p1 wf2 p2          

    let print_vector4_t_2 c wf1 p1 wf2 p2 wf3 p3 fusion (coeff, contraction) =
      match contraction, fusion with
      | C_12_34, (F234|F243|F134|F143|F421|F321|F412|F312)
      | C_13_42, (F324|F342|F124|F142|F431|F231|F413|F213)
      | C_14_23, (F423|F432|F123|F132|F341|F241|F314|F214) ->
          printf "g_dim8g3_t_2(%s,%s,%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2 wf3 p3
      | C_12_34, (F324|F314|F423|F413|F142|F132|F241|F231)
      | C_13_42, (F234|F214|F432|F412|F143|F123|F341|F321)
      | C_14_23, (F243|F213|F342|F312|F134|F124|F431|F421) ->
          printf "g_dim8g3_t_2(%s,%s,%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1 wf3 p3
      | C_12_34, (F342|F341|F432|F431|F124|F123|F214|F213)
      | C_13_42, (F243|F241|F423|F421|F134|F132|F314|F312)
      | C_14_23, (F234|F231|F324|F321|F143|F142|F413|F412) ->
          printf "g_dim8g3_t_2(%s,%s,%s,%s,%s,%s,%s)" c wf3 p3 wf1 p1 wf2 p2

    let print_vector4_m_0 c wf1 p1 wf2 p2 wf3 p3 fusion (coeff, contraction) =
      match contraction, fusion with
      | C_12_34, (F234|F243|F134|F143|F421|F321|F412|F312)
      | C_13_42, (F324|F342|F124|F142|F431|F231|F413|F213)
      | C_14_23, (F423|F432|F123|F132|F341|F241|F314|F214) ->
          printf "g_dim8g3_m_0(%s,%s,%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2 wf3 p3
      | C_12_34, (F324|F314|F423|F413|F142|F132|F241|F231)
      | C_13_42, (F234|F214|F432|F412|F143|F123|F341|F321)
      | C_14_23, (F243|F213|F342|F312|F134|F124|F431|F421) ->
          printf "g_dim8g3_m_0(%s,%s,%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1 wf3 p3
      | C_12_34, (F342|F341|F432|F431|F124|F123|F214|F213)
      | C_13_42, (F243|F241|F423|F421|F134|F132|F314|F312)
      | C_14_23, (F234|F231|F324|F321|F143|F142|F413|F412) ->
          printf "g_dim8g3_m_0(%s,%s,%s,%s,%s,%s,%s)" c wf3 p3 wf1 p1 wf2 p2

    let print_vector4_m_1 c wf1 p1 wf2 p2 wf3 p3 fusion (coeff, contraction) =
      match contraction, fusion with
      | C_12_34, (F234|F243|F134|F143|F421|F321|F412|F312)
      | C_13_42, (F324|F342|F124|F142|F431|F231|F413|F213)
      | C_14_23, (F423|F432|F123|F132|F341|F241|F314|F214) ->
          printf "g_dim8g3_m_1(%s,%s,%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2 wf3 p3
      | C_12_34, (F324|F314|F423|F413|F142|F132|F241|F231)
      | C_13_42, (F234|F214|F432|F412|F143|F123|F341|F321)
      | C_14_23, (F243|F213|F342|F312|F134|F124|F431|F421) ->
          printf "g_dim8g3_m_1(%s,%s,%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1 wf3 p3
      | C_12_34, (F342|F341|F432|F431|F124|F123|F214|F213)
      | C_13_42, (F243|F241|F423|F421|F134|F132|F314|F312)
      | C_14_23, (F234|F231|F324|F321|F143|F142|F413|F412) ->
          printf "g_dim8g3_m_1(%s,%s,%s,%s,%s,%s,%s)" c wf3 p3 wf1 p1 wf2 p2
          
     let print_vector4_m_7 c wf1 p1 wf2 p2 wf3 p3 fusion (coeff, contraction) =
      match contraction, fusion with
      | C_12_34, (F234|F243|F134|F143|F421|F321|F412|F312)
      | C_13_42, (F324|F342|F124|F142|F431|F231|F413|F213)
      | C_14_23, (F423|F432|F123|F132|F341|F241|F314|F214) ->
          printf "g_dim8g3_m_7(%s,%s,%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2 wf3 p3
      | C_12_34, (F324|F314|F423|F413|F142|F132|F241|F231)
      | C_13_42, (F234|F214|F432|F412|F143|F123|F341|F321)
      | C_14_23, (F243|F213|F342|F312|F134|F124|F431|F421) ->
          printf "g_dim8g3_m_7(%s,%s,%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1 wf3 p3
      | C_12_34, (F342|F341|F432|F431|F124|F123|F214|F213)
      | C_13_42, (F243|F241|F423|F421|F134|F132|F314|F312)
      | C_14_23, (F234|F231|F324|F321|F143|F142|F413|F412) ->
          printf "g_dim8g3_m_7(%s,%s,%s,%s,%s,%s,%s)" c wf3 p3 wf1 p1 wf2 p2      

    let print_add_vector4 c wf1 wf2 wf3 fusion (coeff, contraction) =
      printf "@ + ";
      print_vector4 c wf1 wf2 wf3 fusion (coeff, contraction)

    let print_vector4_km c pa pb wf1 wf2 wf3 fusion (coeff, contraction) =
      match contraction, fusion with
      | C_12_34, (F341|F431|F342|F432|F123|F213|F124|F214)
      | C_13_42, (F241|F421|F243|F423|F132|F312|F134|F314)
      | C_14_23, (F231|F321|F234|F324|F142|F412|F143|F413) ->
          printf "((%s%s%s+%s))*(%s*%s))*%s"
            (format_coeff coeff) c pa pb wf1 wf2 wf3
      | C_12_34, (F134|F143|F234|F243|F312|F321|F412|F421)
      | C_13_42, (F124|F142|F324|F342|F213|F231|F413|F431)
      | C_14_23, (F123|F132|F423|F432|F214|F241|F314|F341) ->
          printf "((%s%s%s+%s))*(%s*%s))*%s"
            (format_coeff coeff) c pa pb wf2 wf3 wf1
      | C_12_34, (F314|F413|F324|F423|F132|F231|F142|F241)
      | C_13_42, (F214|F412|F234|F432|F123|F321|F143|F341)
      | C_14_23, (F213|F312|F243|F342|F124|F421|F134|F431) ->
          printf "((%s%s%s+%s))*(%s*%s))*%s"
            (format_coeff coeff) c pa pb wf1 wf3 wf2
            
    let print_vector4_km_t_0 c pa pb wf1 p1 wf2 p2 wf3 p3 fusion (coeff, contraction) =
      match contraction, fusion with
      | C_12_34, (F234|F243|F134|F143|F421|F321|F412|F312)
      | C_13_42, (F324|F342|F124|F142|F431|F231|F413|F213)
      | C_14_23, (F423|F432|F123|F132|F341|F241|F314|F214) ->
          printf "@[(%s%s%s+%s)*g_dim8g3_t_0(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf1 p1 wf2 p2 wf3 p3
      | C_12_34, (F324|F314|F423|F413|F142|F132|F241|F231)
      | C_13_42, (F234|F214|F432|F412|F143|F123|F341|F321)
      | C_14_23, (F243|F213|F342|F312|F134|F124|F431|F421) ->
          printf "@[(%s%s%s+%s)*g_dim8g3_t_0(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf2 p2 wf1 p1 wf3 p3
      | C_12_34, (F342|F341|F432|F431|F124|F123|F214|F213)
      | C_13_42, (F243|F241|F423|F421|F134|F132|F314|F312)
      | C_14_23, (F234|F231|F324|F321|F143|F142|F413|F412) ->
          printf "@[(%s%s%s+%s)*g_dim8g3_t_0(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf3 p3 wf1 p1 wf2 p2 

    let print_vector4_km_t_1 c pa pb wf1 p1 wf2 p2 wf3 p3 fusion (coeff, contraction) =
      match contraction, fusion with
      | C_12_34, (F234|F243|F134|F143|F421|F321|F412|F312)
      | C_13_42, (F324|F342|F124|F142|F431|F231|F413|F213)
      | C_14_23, (F423|F432|F123|F132|F341|F241|F314|F214) ->
          printf "@[(%s%s%s+%s)*g_dim8g3_t_1(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf1 p1 wf2 p2 wf3 p3
      | C_12_34, (F324|F314|F423|F413|F142|F132|F241|F231)
      | C_13_42, (F234|F214|F432|F412|F143|F123|F341|F321)
      | C_14_23, (F243|F213|F342|F312|F134|F124|F431|F421) ->
          printf "@[(%s%s%s+%s)*g_dim8g3_t_1(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf2 p2 wf1 p1 wf3 p3
      | C_12_34, (F342|F341|F432|F431|F124|F123|F214|F213)
      | C_13_42, (F243|F241|F423|F421|F134|F132|F314|F312)
      | C_14_23, (F234|F231|F324|F321|F143|F142|F413|F412) ->
          printf "@[(%s%s%s+%s)*g_dim8g3_t_1(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf3 p3 wf1 p1 wf2 p2
            
    let print_vector4_km_t_2 c pa pb wf1 p1 wf2 p2 wf3 p3 fusion (coeff, contraction) =
      match contraction, fusion with
      | C_12_34, (F234|F243|F134|F143|F421|F321|F412|F312)
      | C_13_42, (F324|F342|F124|F142|F431|F231|F413|F213)
      | C_14_23, (F423|F432|F123|F132|F341|F241|F314|F214) ->
          printf "@[(%s%s%s+%s)*g_dim8g3_t_2(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf1 p1 wf2 p2 wf3 p3
      | C_12_34, (F324|F314|F423|F413|F142|F132|F241|F231)
      | C_13_42, (F234|F214|F432|F412|F143|F123|F341|F321)
      | C_14_23, (F243|F213|F342|F312|F134|F124|F431|F421) ->
          printf "@[(%s%s%s+%s)*g_dim8g3_t_2(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf2 p2 wf1 p1 wf3 p3
      | C_12_34, (F342|F341|F432|F431|F124|F123|F214|F213)
      | C_13_42, (F243|F241|F423|F421|F134|F132|F314|F312)
      | C_14_23, (F234|F231|F324|F321|F143|F142|F413|F412) ->
          printf "@[(%s%s%s+%s)*g_dim8g3_t_2(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf3 p3 wf1 p1 wf2 p2
            
    let print_vector4_km_t_rsi c pa pb pc wf1 p1 wf2 p2 wf3 p3 fusion (coeff, contraction) =
      match contraction, fusion with
      | C_12_34, (F234|F243|F134|F143|F421|F321|F412|F312)
      | C_13_42, (F324|F342|F124|F142|F431|F231|F413|F213)
      | C_14_23, (F423|F432|F123|F132|F341|F241|F314|F214) ->
          printf "@[(%s%s%s+%s)*g_dim8g3_t_0(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf1 p1 wf2 p2 wf3 p3
      | C_12_34, (F324|F314|F423|F413|F142|F132|F241|F231)
      | C_13_42, (F234|F214|F432|F412|F143|F123|F341|F321)
      | C_14_23, (F243|F213|F342|F312|F134|F124|F431|F421) ->
          printf "@[(%s%s%s+%s)*g_dim8g3_t_0(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))*((%s+%s)*(%s+%s)/((%s+%s)*(%s+%s)))@]"
            (format_coeff coeff) c pa pb wf2 p2 wf1 p1 wf3 p3 pa pb pa pb pb pc pb pc
      | C_12_34, (F342|F341|F432|F431|F124|F123|F214|F213)
      | C_13_42, (F243|F241|F423|F421|F134|F132|F314|F312)
      | C_14_23, (F234|F231|F324|F321|F143|F142|F413|F412) ->
          printf "@[(%s%s%s+%s)*g_dim8g3_t_0(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))*((%s+%s)*(%s+%s)/((%s+%s)*(%s+%s)))@]"
            (format_coeff coeff) c pa pb wf3 p3 wf1 p1 wf2 p2 pa pb pa pb pa pc pa pc            

    let print_vector4_km_m_0 c pa pb wf1 p1 wf2 p2 wf3 p3 fusion (coeff, contraction) =
      match contraction, fusion with
      | C_12_34, (F234|F243|F134|F143|F421|F321|F412|F312)
      | C_13_42, (F324|F342|F124|F142|F431|F231|F413|F213)
      | C_14_23, (F423|F432|F123|F132|F341|F241|F314|F214) ->
          if (String.contains c 'w' || String.contains c '4') then
             printf "@[(%s%s%s+%s)*g_dim8g3_m_0(cmplx(1,kind=default),cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
                (format_coeff coeff) c pa pb wf1 p1 wf2 p2 wf3 p3
          else
             printf "@[((%s%s%s+%s))*g_dim8g3_m_0(cmplx(costhw**(-2),kind=default),cmplx(costhw**2,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
                (format_coeff coeff) c pa pb wf1 p1 wf2 p2 wf3 p3
      | C_12_34, (F324|F314|F423|F413|F142|F132|F241|F231)
      | C_13_42, (F234|F214|F432|F412|F143|F123|F341|F321)
      | C_14_23, (F243|F213|F342|F312|F134|F124|F431|F421) ->
          if (String.contains c 'w' || String.contains c '4') then
             printf "@[(%s%s%s+%s)*g_dim8g3_m_0(cmplx(1,kind=default),cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
                (format_coeff coeff) c pa pb  wf2 p2 wf1 p1 wf3 p3
          else
             printf "@[(%s%s%s+%s)*g_dim8g3_m_0(cmplx(costhw**(-2),kind=default),cmplx(costhw**2,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
                (format_coeff coeff) c pa pb wf2 p2 wf1 p1 wf3 p3
      | C_12_34, (F342|F341|F432|F431|F124|F123|F214|F213)
      | C_13_42, (F243|F241|F423|F421|F134|F132|F314|F312)
      | C_14_23, (F234|F231|F324|F321|F143|F142|F413|F412) ->
          if (String.contains c 'w' || String.contains c '4') then
             printf "@[(%s%s%s+%s)*g_dim8g3_m_0(cmplx(1,kind=default),cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
                (format_coeff coeff) c pa pb wf3 p3 wf1 p1 wf2 p2
          else
             printf "@[(%s%s%s+%s)*g_dim8g3_m_0(cmplx(costhw**(-2),kind=default),cmplx(costhw**2,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
                (format_coeff coeff) c pa pb wf3 p3 wf1 p1 wf2 p2

    let print_vector4_km_m_1 c pa pb wf1 p1 wf2 p2 wf3 p3 fusion (coeff, contraction) =
      match contraction, fusion with
      | C_12_34, (F234|F243|F134|F143|F421|F321|F412|F312)
      | C_13_42, (F324|F342|F124|F142|F431|F231|F413|F213)
      | C_14_23, (F423|F432|F123|F132|F341|F241|F314|F214) ->
          if (String.contains c 'w' || String.contains c '4') then
             printf "@[(%s%s%s+%s)*g_dim8g3_m_1(cmplx(1,kind=default),cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
                (format_coeff coeff) c pa pb wf1 p1 wf2 p2 wf3 p3
          else
             printf "@[(%s%s%s+%s)*g_dim8g3_m_1(cmplx(costhw**(-2),kind=default),cmplx(costhw**2,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
                (format_coeff coeff) c pa pb wf1 p1 wf2 p2 wf3 p3
      | C_12_34, (F324|F314|F423|F413|F142|F132|F241|F231)
      | C_13_42, (F234|F214|F432|F412|F143|F123|F341|F321)
      | C_14_23, (F243|F213|F342|F312|F134|F124|F431|F421) ->
          if (String.contains c 'w' || String.contains c '4') then
             printf "@[(%s%s%s+%s)*g_dim8g3_m_1(cmplx(1,kind=default),cmplx(1,kind=default),@  %s,%s,%s,%s,%s,%s))@]"
                (format_coeff coeff) c pa pb wf2 p2 wf1 p1 wf3 p3
          else
             printf "@[(%s%s%s+%s)*g_dim8g3_m_1(cmplx(costhw**(-2),kind=default),cmplx(costhw**2,kind=default),@  %s,%s,%s,%s,%s,%s))@]"
                (format_coeff coeff) c pa pb wf2 p2 wf1 p1 wf3 p3
      | C_12_34, (F342|F341|F432|F431|F124|F123|F214|F213)
      | C_13_42, (F243|F241|F423|F421|F134|F132|F314|F312)
      | C_14_23, (F234|F231|F324|F321|F143|F142|F413|F412) ->
          if (String.contains c 'w' || String.contains c '4') then
             printf "@[(%s%s%s+%s)*g_dim8g3_m_1(cmplx(1,kind=default),cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
                (format_coeff coeff) c pa pb wf3 p3 wf1 p1 wf2 p2
          else
             printf "@[(%s%s%s+%s)*g_dim8g3_m_1(cmplx(costhw**(-2),kind=default),cmplx(costhw**2,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
                (format_coeff coeff) c pa pb wf3 p3 wf1 p1 wf2 p2
                
    let print_vector4_km_m_7 c pa pb wf1 p1 wf2 p2 wf3 p3 fusion (coeff, contraction) =
      match contraction, fusion with
      | C_12_34, (F234|F243|F134|F143|F421|F321|F412|F312)
      | C_13_42, (F324|F342|F124|F142|F431|F231|F413|F213)
      | C_14_23, (F423|F432|F123|F132|F341|F241|F314|F214) ->
          if (String.contains c 'w' || String.contains c '4') then
             printf "@[(%s%s%s+%s)*@ g_dim8g3_m_7(cmplx(1,kind=default),cmplx(1,kind=default),cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
                (format_coeff coeff) c pa pb wf1 p1 wf2 p2 wf3 p3
          else
             printf "@[(%s%s%s+%s)*@ g_dim8g3_m_7(cmplx(costhw**(-2),kind=default),cmplx(1,kind=default),cmplx(costhw**2,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
                (format_coeff coeff) c pa pb wf1 p1 wf2 p2 wf3 p3
      | C_12_34, (F324|F314|F423|F413|F142|F132|F241|F231)
      | C_13_42, (F234|F214|F432|F412|F143|F123|F341|F321)
      | C_14_23, (F243|F213|F342|F312|F134|F124|F431|F421) ->
          if (String.contains c 'w' || String.contains c '4') then
             printf "@[(%s%s%s+%s)*@ g_dim8g3_m_7(cmplx(1,kind=default),cmplx(1,kind=default),cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
                (format_coeff coeff) c pa pb wf2 p2 wf1 p1 wf3 p3
          else
             printf "@[(%s%s%s+%s)*@ g_dim8g3_m_7(cmplx(costhw**(-2),kind=default),cmplx(1,kind=default),cmplx(costhw**2,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
                (format_coeff coeff) c pa pb wf2 p2 wf1 p1 wf3 p3
      | C_12_34, (F342|F341|F432|F431|F124|F123|F214|F213)
      | C_13_42, (F243|F241|F423|F421|F134|F132|F314|F312)
      | C_14_23, (F234|F231|F324|F321|F143|F142|F413|F412) ->
          if (String.contains c 'w' || String.contains c '4') then
             printf "@[(%s%s%s+%s)*@ g_dim8g3_m_7(cmplx(1,kind=default),cmplx(1,kind=default),cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
                (format_coeff coeff) c pa pb wf3 p3 wf1 p1 wf2 p2
          else
             printf "@[(%s%s%s+%s)*@ g_dim8g3_m_7(cmplx(costhw**(-2),kind=default),cmplx(1,kind=default),cmplx(costhw**2,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
                (format_coeff coeff) c pa pb wf3 p3 wf1 p1 wf2 p2        

    let print_add_vector4_km c pa pb wf1 wf2 wf3 fusion (coeff, contraction) =
      printf "@ + ";
      print_vector4_km c pa pb wf1 wf2 wf3 fusion (coeff, contraction)

    let print_dscalar4 c wf1 wf2 wf3 p1 p2 p3 p123
        fusion (coeff, contraction) =
      match contraction, fusion with
      | C_12_34, (F341|F431|F342|F432|F123|F213|F124|F214)
      | C_13_42, (F241|F421|F243|F423|F132|F312|F134|F314)
      | C_14_23, (F231|F321|F234|F324|F142|F412|F143|F413) ->
          printf "((%s%s)*(%s*%s)*(%s*%s)*%s*%s*%s)"
            (format_coeff coeff) c p1 p2 p3 p123 wf1 wf2 wf3
      | C_12_34, (F134|F143|F234|F243|F312|F321|F412|F421)
      | C_13_42, (F124|F142|F324|F342|F213|F231|F413|F431)
      | C_14_23, (F123|F132|F423|F432|F214|F241|F314|F341) ->
          printf "((%s%s)*(%s*%s)*(%s*%s)*%s*%s*%s)"
            (format_coeff coeff) c p2 p3 p1 p123 wf1 wf2 wf3
      | C_12_34, (F314|F413|F324|F423|F132|F231|F142|F241)
      | C_13_42, (F214|F412|F234|F432|F123|F321|F143|F341)
      | C_14_23, (F213|F312|F243|F342|F124|F421|F134|F431) ->
          printf "((%s%s)*(%s*%s)*(%s*%s)*%s*%s*%s)"
            (format_coeff coeff) c p1 p3 p2 p123 wf1 wf2 wf3

    let print_add_dscalar4 c wf1 wf2 wf3 p1 p2 p3 p123
        fusion (coeff, contraction) =
      printf "@ + ";
      print_dscalar4 c wf1 wf2 wf3 p1 p2 p3 p123 fusion (coeff, contraction)

    let print_dscalar2_vector2 c wf1 wf2 wf3 p1 p2 p3 p123 fusion (coeff, contraction) =
      match contraction, fusion with
      | C_12_34, (F123|F213|F124|F214) ->
          printf "(%s%s)*(%s*%s)*(%s*%s)*%s"
            (format_coeff coeff) c p1 p2 wf1 wf2 wf3
      | C_12_34, (F134|F143|F234|F243) ->
          printf "(%s%s)*(%s*%s)*(%s*%s)*%s"
            (format_coeff coeff) c p1 p123 wf2 wf3 wf1
      | C_12_34, (F132|F231|F142|F241) ->
          printf "(%s%s)*(%s*%s)*(%s*%s)*%s"
            (format_coeff coeff) c p1 p3 wf1 wf3 wf2  
      | C_12_34, (F312|F321|F412|F421) ->
          printf "(%s%s)*(%s*%s)*(%s*%s)*%s"
            (format_coeff coeff) c p2 p3 wf2 wf3 wf1 
      | C_12_34, (F314|F413|F324|F423) ->
          printf "(%s%s)*(%s*%s)*(%s*%s)*%s"
            (format_coeff coeff) c p2 p123 wf1 wf3 wf2  
      | C_12_34, (F341|F431|F342|F432) ->
          printf "(%s%s)*(%s*%s)*(%s*%s)*%s"
            (format_coeff coeff) c p3 p123 wf1 wf2 wf3
      | C_13_42, (F123|F214) 
      | C_14_23, (F124|F213) ->
          printf "((%s%s)*(%s*%s*%s)*%s*%s)"
            (format_coeff coeff) c wf1 p1 wf3 wf2 p2
      | C_13_42, (F124|F213) 
      | C_14_23, (F123|F214) ->
          printf "((%s%s)*(%s*%s*%s)*%s*%s)"
            (format_coeff coeff) c wf2 p2 wf3 wf1 p1
      | C_13_42, (F132|F241) 
      | C_14_23, (F142|F231) ->
          printf "((%s%s)*(%s*%s*%s)*%s*%s)"
            (format_coeff coeff) c wf1 p1 wf2 wf3 p3
      | C_13_42, (F142|F231) 
      | C_14_23, (F132|F241) ->
          printf "((%s%s)*(%s*%s*%s)*%s*%s)"
            (format_coeff coeff) c wf3 p3 wf2 wf1 p1
      | C_13_42, (F312|F421) 
      | C_14_23, (F412|F321) ->
          printf "((%s%s)*(%s*%s*%s)*%s*%s)"
            (format_coeff coeff) c wf2 p2 wf1 wf3 p3
      | C_13_42, (F321|F412) 
      | C_14_23, (F421|F312) ->
          printf "((%s%s)*(%s*%s*%s)*%s*%s)"
            (format_coeff coeff) c wf3 p3 wf1 wf2 p2
      | C_13_42, (F134|F243) 
      | C_14_23, (F143|F234) ->
          printf "((%s%s)*(%s*%s)*(%s*%s*%s))"
            (format_coeff coeff) c wf3 p123 wf1 p1 wf2
      | C_13_42, (F143|F234) 
      | C_14_23, (F134|F243) ->
          printf "((%s%s)*(%s*%s)*(%s*%s*%s))"
            (format_coeff coeff) c wf2 p123 wf1 p1 wf3
      | C_13_42, (F314|F423) 
      | C_14_23, (F413|F324) ->
          printf "((%s%s)*(%s*%s)*(%s*%s*%s))"
            (format_coeff coeff) c wf3 p123 wf2 p2 wf1
      | C_13_42, (F324|F413) 
      | C_14_23, (F423|F314) ->
          printf "((%s%s)*(%s*%s)*(%s*%s*%s))"
            (format_coeff coeff) c wf1 p123 wf2 p2 wf3
      | C_13_42, (F341|F432) 
      | C_14_23, (F431|F342) ->
          printf "((%s%s)*(%s*%s)*(%s*%s*%s))"
            (format_coeff coeff) c wf2 p123 wf3 p3 wf1
      | C_13_42, (F342|F431) 
      | C_14_23, (F432|F341) ->
          printf "((%s%s)*(%s*%s)*(%s*%s*%s))"
            (format_coeff coeff) c wf1 p123 wf3 p3 wf2

    let print_add_dscalar2_vector2 c wf1 wf2 wf3 p1 p2 p3 p123
        fusion (coeff, contraction) =
      printf "@ + ";
      print_dscalar2_vector2 c wf1 wf2 wf3 p1 p2 p3 p123
        fusion (coeff, contraction)

    let print_dscalar2_vector2_km c pa pb wf1 wf2 wf3 p1 p2 p3 p123 fusion (coeff, contraction) =
      match contraction, fusion with
      | C_12_34, (F123|F213|F124|F214) ->
          printf "(%s%s%s+%s))*(%s*%s)*(%s*%s)*%s"
            (format_coeff coeff) c pa pb p1 p2 wf1 wf2 wf3
      | C_12_34, (F134|F143|F234|F243) ->
          printf "(%s%s%s+%s))*(%s*%s)*(%s*%s)*%s"
            (format_coeff coeff) c pa pb p1 p123 wf2 wf3 wf1
      | C_12_34, (F132|F231|F142|F241) ->
          printf "(%s%s%s+%s))*(%s*%s)*(%s*%s)*%s"
            (format_coeff coeff) c pa pb p1 p3 wf1 wf3 wf2  
      | C_12_34, (F312|F321|F412|F421) ->
          printf "(%s%s%s+%s))*(%s*%s)*(%s*%s)*%s"
            (format_coeff coeff) c pa pb p2 p3 wf2 wf3 wf1 
      | C_12_34, (F314|F413|F324|F423) ->
          printf "(%s%s%s+%s))*(%s*%s)*(%s*%s)*%s"
            (format_coeff coeff) c pa pb p2 p123 wf1 wf3 wf2  
      | C_12_34, (F341|F431|F342|F432) ->
          printf "(%s%s%s+%s))*(%s*%s)*(%s*%s)*%s"
            (format_coeff coeff) c pa pb p3 p123 wf1 wf2 wf3
      | C_13_42, (F123|F214) 
      | C_14_23, (F124|F213) ->
          printf "((%s%s%s+%s))*(%s*%s*%s)*%s*%s)"
            (format_coeff coeff) c pa pb wf1 p1 wf3 wf2 p2
      | C_13_42, (F124|F213) 
      | C_14_23, (F123|F214) ->
          printf "((%s%s%s+%s))*(%s*%s*%s)*%s*%s)"
            (format_coeff coeff) c pa pb wf2 p2 wf3 wf1 p1
      | C_13_42, (F132|F241) 
      | C_14_23, (F142|F231) ->
          printf "((%s%s%s+%s))*(%s*%s*%s)*%s*%s)"
            (format_coeff coeff) c pa pb wf1 p1 wf2 wf3 p3
      | C_13_42, (F142|F231) 
      | C_14_23, (F132|F241) ->
          printf "((%s%s%s+%s))*(%s*%s*%s)*%s*%s)"
            (format_coeff coeff) c pa pb wf3 p3 wf2 wf1 p1
      | C_13_42, (F312|F421) 
      | C_14_23, (F412|F321) ->
          printf "((%s%s%s+%s))*(%s*%s*%s)*%s*%s)"
            (format_coeff coeff) c pa pb wf2 p2 wf1 wf3 p3
      | C_13_42, (F321|F412) 
      | C_14_23, (F421|F312) ->
          printf "((%s%s%s+%s))*(%s*%s*%s)*%s*%s)"
            (format_coeff coeff) c pa pb wf3 p3 wf1 wf2 p2
      | C_13_42, (F134|F243) 
      | C_14_23, (F143|F234) ->
          printf "((%s%s%s+%s))*(%s*%s)*(%s*%s*%s))"
            (format_coeff coeff) c pa pb wf3 p123 wf1 p1 wf2
      | C_13_42, (F143|F234) 
      | C_14_23, (F134|F243) ->
          printf "((%s%s%s+%s))*(%s*%s)*(%s*%s*%s))"
            (format_coeff coeff) c pa pb wf2 p123 wf1 p1 wf3
      | C_13_42, (F314|F423) 
      | C_14_23, (F413|F324) ->
          printf "((%s%s%s+%s))*(%s*%s)*(%s*%s*%s))"
            (format_coeff coeff) c pa pb wf3 p123 wf2 p2 wf1
      | C_13_42, (F324|F413) 
      | C_14_23, (F423|F314) ->
          printf "((%s%s%s+%s))*(%s*%s)*(%s*%s*%s))"
            (format_coeff coeff) c pa pb wf1 p123 wf2 p2 wf3
      | C_13_42, (F341|F432) 
      | C_14_23, (F431|F342) ->
          printf "((%s%s%s+%s))*(%s*%s)*(%s*%s*%s))"
            (format_coeff coeff) c pa pb wf2 p123 wf3 p3 wf1
      | C_13_42, (F342|F431) 
      | C_14_23, (F432|F341) ->
          printf "((%s%s%s+%s))*(%s*%s)*(%s*%s*%s))"
            (format_coeff coeff) c pa pb wf1 p123 wf3 p3 wf2

    let print_add_dscalar2_vector2_km c pa pb wf1 wf2 wf3 p1 p2 p3 p123 fusion (coeff, contraction) =
      printf "@ + ";
      print_dscalar2_vector2_km c pa pb wf1 wf2 wf3 p1 p2 p3 p123 fusion (coeff, contraction)
      
    let print_dscalar2_vector2_m_0_km c pa pb wf1 wf2 wf3 p1 p2 p3 fusion (coeff, contraction) =
      match contraction, fusion with
      | C_12_34, (F123|F213|F124|F214) ->
          printf "@[((%s%s%s+%s))*v_phi2v_m_0(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf1 p1 wf2 p2 wf3 p3
      | C_12_34, (F134|F143|F234|F243) ->
          printf "@[((%s%s%s+%s))*phi_phi2v_m_0(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf1 p1 wf2 p2 wf3 p3
      | C_12_34, (F132|F231|F142|F241) ->
          printf "@[((%s%s%s+%s))*v_phi2v_m_0(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf1 p1 wf3 p3 wf2 p2
      | C_12_34, (F312|F321|F412|F421) ->
          printf "@[((%s%s%s+%s))*v_phi2v_m_0(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf3 p3 wf2 p2 wf1 p1
      | C_12_34, (F314|F413|F324|F423) ->
          printf "@[((%s%s%s+%s))*phi_phi2v_m_0(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf2 p2 wf1 p1 wf3 p3
      | C_12_34, (F341|F431|F342|F432) ->
          printf "@[((%s%s%s+%s))*phi_phi2v_m_0(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf3 p3 wf2 p2 wf1 p1
      | C_13_42, (F123|F214)
      | C_14_23, (F124|F213) ->
          printf "@[((%s%s%s+%s))*v_phi2v_m_0(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf1 p1 wf2 p3 wf3 p2
      | C_13_42, (F124|F213)
      | C_14_23, (F123|F214) ->
          printf "@[((%s%s%s+%s))*v_phi2v_m_0(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf2 p2 wf1 p3 wf3 p1
      | C_13_42, (F132|F241)
      | C_14_23, (F142|F231) ->
          printf "@[((%s%s%s+%s))*v_phi2v_m_0(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf1 p1 wf3 p2 wf2 p3
      | C_13_42, (F142|F231)
      | C_14_23, (F132|F241) ->
          printf "@[((%s%s%s+%s))*v_phi2v_m_0(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf3 p3 wf1 p2 wf2 p1
      | C_13_42, (F312|F421)
      | C_14_23, (F412|F321) ->
          printf "@[((%s%s%s+%s))*v_phi2v_m_0(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf2 p2 wf3 p1 wf1 p3
      | C_13_42, (F321|F412)
      | C_14_23, (F421|F312) ->
          printf "@[((%s%s%s+%s))*v_phi2v_m_0(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf3 p3 wf2 p1 wf1 p2
      | C_13_42, (F134|F243)
      | C_14_23, (F143|F234) ->
          printf "@[((%s%s%s+%s))*phi_phi2v_m_0(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf1 p3 wf3 p1 wf2 p2
      | C_13_42, (F143|F234)
      | C_14_23, (F134|F243) ->
          printf "@[((%s%s%s+%s))*phi_phi2v_m_0(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf1 p2 wf2 p1 wf3 p3
      | C_13_42, (F314|F423)
      | C_14_23, (F413|F324) ->
          printf "@[((%s%s%s+%s))*phi_phi2v_m_0(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf2 p3 wf3 p2 wf1 p1
      | C_13_42, (F324|F413)
      | C_14_23, (F423|F314) ->
          printf "@[((%s%s%s+%s))*phi_phi2v_m_0(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf2 p1 wf1 p2 wf3 p3
      | C_13_42, (F341|F432)
      | C_14_23, (F431|F342) ->
          printf "@[((%s%s%s+%s))*phi_phi2v_m_0(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf3 p2 wf2 p3 wf1 p1
      | C_13_42, (F342|F431)
      | C_14_23, (F432|F341) ->
          printf "@[((%s%s%s+%s))*phi_phi2v_m_0(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf3 p1 wf1 p3 wf2 p2

    let print_add_dscalar2_vector2_m_0_km c pa pb wf1 wf2 wf3 p1 p2 p3 fusion (coeff, contraction) =
      printf "@ + ";
      print_dscalar2_vector2_m_0_km c pa pb wf1 wf2 wf3 p1 p2 p3 fusion (coeff, contraction)

   let print_dscalar2_vector2_m_1_km c pa pb wf1 wf2 wf3 p1 p2 p3 fusion (coeff, contraction) =
      match contraction, fusion with
      | C_12_34, (F123|F213|F124|F214) ->
          printf "@[((%s%s%s+%s))*v_phi2v_m_1(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf1 p1 wf2 p2 wf3 p3
      | C_12_34, (F134|F143|F234|F243) ->
          printf "@[((%s%s%s+%s))*phi_phi2v_m_1(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf1 p1 wf2 p2 wf3 p3
      | C_12_34, (F132|F231|F142|F241) ->
          printf "@[((%s%s%s+%s))*v_phi2v_m_1(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf1 p1 wf3 p3 wf2 p2
      | C_12_34, (F312|F321|F412|F421) ->
          printf "@[((%s%s%s+%s))*v_phi2v_m_1(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf3 p3 wf2 p2 wf1 p1
      | C_12_34, (F314|F413|F324|F423) ->
          printf "@[((%s%s%s+%s))*phi_phi2v_m_1(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf2 p2 wf1 p1 wf3 p3
      | C_12_34, (F341|F431|F342|F432) ->
          printf "@[((%s%s%s+%s))*phi_phi2v_m_1(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf3 p3 wf2 p2 wf1 p1
      | C_13_42, (F123|F214)
      | C_14_23, (F124|F213) ->
          printf "@[((%s%s%s+%s))*v_phi2v_m_1(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf1 p1 wf2 p3 wf3 p2
      | C_13_42, (F124|F213)
      | C_14_23, (F123|F214) ->
          printf "@[((%s%s%s+%s))*v_phi2v_m_1(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf2 p2 wf1 p3 wf3 p1
      | C_13_42, (F132|F241)
      | C_14_23, (F142|F231) ->
          printf "@[((%s%s%s+%s))*v_phi2v_m_1(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf1 p1 wf3 p2 wf2 p3
      | C_13_42, (F142|F231)
      | C_14_23, (F132|F241) ->
          printf "@[((%s%s%s+%s))*v_phi2v_m_1(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf3 p3 wf1 p2 wf2 p1
      | C_13_42, (F312|F421)
      | C_14_23, (F412|F321) ->
          printf "@[((%s%s%s+%s))*v_phi2v_m_1(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf2 p2 wf3 p1 wf1 p3
      | C_13_42, (F321|F412)
      | C_14_23, (F421|F312) ->
          printf "@[((%s%s%s+%s))*v_phi2v_m_1(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf3 p3 wf2 p1 wf1 p2
      | C_13_42, (F134|F243)
      | C_14_23, (F143|F234) ->
          printf "@[((%s%s%s+%s))*phi_phi2v_m_1(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf1 p3 wf3 p1 wf2 p2
      | C_13_42, (F143|F234)
      | C_14_23, (F134|F243) ->
          printf "@[((%s%s%s+%s))*phi_phi2v_m_1(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf1 p2 wf2 p1 wf3 p3
      | C_13_42, (F314|F423)
      | C_14_23, (F413|F324) ->
          printf "@[((%s%s%s+%s))*phi_phi2v_m_1(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf2 p3 wf3 p2 wf1 p1
      | C_13_42, (F324|F413)
      | C_14_23, (F423|F314) ->
          printf "@[((%s%s%s+%s))*phi_phi2v_m_1(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf2 p1 wf1 p2 wf3 p3
      | C_13_42, (F341|F432)
      | C_14_23, (F431|F342) ->
          printf "@[((%s%s%s+%s))*phi_phi2v_m_1(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf3 p2 wf2 p3 wf1 p1
      | C_13_42, (F342|F431)
      | C_14_23, (F432|F341) ->
          printf "@[((%s%s%s+%s))*phi_phi2v_m_1(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf3 p1 wf1 p3 wf2 p2

    let print_add_dscalar2_vector2_m_1_km c pa pb wf1 wf2 wf3 p1 p2 p3 fusion (coeff, contraction) =
      printf "@ + ";
      print_dscalar2_vector2_m_1_km c pa pb wf1 wf2 wf3 p1 p2 p3 fusion (coeff, contraction)
      
   let print_dscalar2_vector2_m_7_km c pa pb wf1 wf2 wf3 p1 p2 p3 fusion (coeff, contraction) =
      match contraction, fusion with
      | C_12_34, (F123|F213|F124|F214) ->
          printf "@[((%s%s%s+%s))*v_phi2v_m_7(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf1 p1 wf2 p2 wf3 p3
      | C_12_34, (F134|F143|F234|F243) ->
          printf "@[((%s%s%s+%s))*phi_phi2v_m_7(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf1 p1 wf2 p2 wf3 p3
      | C_12_34, (F132|F231|F142|F241) ->
          printf "@[((%s%s%s+%s))*v_phi2v_m_7(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf1 p1 wf3 p3 wf2 p2
      | C_12_34, (F312|F321|F412|F421) ->
          printf "@[((%s%s%s+%s))*v_phi2v_m_7(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf3 p3 wf2 p2 wf1 p1
      | C_12_34, (F314|F413|F324|F423) ->
          printf "@[((%s%s%s+%s))*phi_phi2v_m_7(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf2 p2 wf1 p1 wf3 p3
      | C_12_34, (F341|F431|F342|F432) ->
          printf "@[((%s%s%s+%s))*phi_phi2v_m_7(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf3 p3 wf2 p2 wf1 p1
      | C_13_42, (F123|F214)
      | C_14_23, (F124|F213) ->
          printf "@[((%s%s%s+%s))*v_phi2v_m_7(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf1 p1 wf2 p3 wf3 p2
      | C_13_42, (F124|F213)
      | C_14_23, (F123|F214) ->
          printf "@[((%s%s%s+%s))*v_phi2v_m_7(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf2 p2 wf1 p3 wf3 p1
      | C_13_42, (F132|F241)
      | C_14_23, (F142|F231) ->
          printf "@[((%s%s%s+%s))*v_phi2v_m_7(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf1 p1 wf3 p2 wf2 p3
      | C_13_42, (F142|F231)
      | C_14_23, (F132|F241) ->
          printf "@[((%s%s%s+%s))*v_phi2v_m_7(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf3 p3 wf1 p2 wf2 p1
      | C_13_42, (F312|F421)
      | C_14_23, (F412|F321) ->
          printf "@[((%s%s%s+%s))*v_phi2v_m_7(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf2 p2 wf3 p1 wf1 p3
      | C_13_42, (F321|F412)
      | C_14_23, (F421|F312) ->
          printf "@[((%s%s%s+%s))*v_phi2v_m_7(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf3 p3 wf2 p1 wf1 p2
      | C_13_42, (F134|F243)
      | C_14_23, (F143|F234) ->
          printf "@[((%s%s%s+%s))*phi_phi2v_m_7(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf1 p3 wf3 p1 wf2 p2
      | C_13_42, (F143|F234)
      | C_14_23, (F134|F243) ->
          printf "@[((%s%s%s+%s))*phi_phi2v_m_7(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf1 p2 wf2 p1 wf3 p3
      | C_13_42, (F314|F423)
      | C_14_23, (F413|F324) ->
          printf "@[((%s%s%s+%s))*phi_phi2v_m_7(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf2 p3 wf3 p2 wf1 p1
      | C_13_42, (F324|F413)
      | C_14_23, (F423|F314) ->
          printf "@[((%s%s%s+%s))*phi_phi2v_m_7(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf2 p1 wf1 p2 wf3 p3
      | C_13_42, (F341|F432)
      | C_14_23, (F431|F342) ->
          printf "@[((%s%s%s+%s))*phi_phi2v_m_7(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf3 p2 wf2 p3 wf1 p1
      | C_13_42, (F342|F431)
      | C_14_23, (F432|F341) ->
          printf "@[((%s%s%s+%s))*phi_phi2v_m_7(cmplx(1,kind=default),@ %s,%s,%s,%s,%s,%s))@]"
            (format_coeff coeff) c pa pb wf3 p1 wf1 p3 wf2 p2

    let print_add_dscalar2_vector2_m_7_km c pa pb wf1 wf2 wf3 p1 p2 p3 fusion (coeff, contraction) =
      printf "@ + ";
      print_dscalar2_vector2_m_7_km c pa pb wf1 wf2 wf3 p1 p2 p3 fusion (coeff, contraction)  

    let print_dscalar4_km c pa pb wf1 wf2 wf3 p1 p2 p3 p123 fusion (coeff, contraction) =
      match contraction, fusion with
      | C_12_34, (F341|F431|F342|F432|F123|F213|F124|F214)
      | C_13_42, (F241|F421|F243|F423|F132|F312|F134|F314)
      | C_14_23, (F231|F321|F234|F324|F142|F412|F143|F413) ->
          printf "((%s%s%s+%s))*(%s*%s)*(%s*%s)*%s*%s*%s)"
            (format_coeff coeff) c pa pb p1 p2 p3 p123 wf1 wf2 wf3
      | C_12_34, (F134|F143|F234|F243|F312|F321|F412|F421)
      | C_13_42, (F124|F142|F324|F342|F213|F231|F413|F431)
      | C_14_23, (F123|F132|F423|F432|F214|F241|F314|F341) ->
          printf "((%s%s%s+%s))*(%s*%s)*(%s*%s)*%s*%s*%s)"
            (format_coeff coeff) c pa pb p2 p3 p1 p123 wf1 wf2 wf3
      | C_12_34, (F314|F413|F324|F423|F132|F231|F142|F241)
      | C_13_42, (F214|F412|F234|F432|F123|F321|F143|F341)
      | C_14_23, (F213|F312|F243|F342|F124|F421|F134|F431) ->
          printf "((%s%s%s+%s))*(%s*%s)*(%s*%s)*%s*%s*%s)"
            (format_coeff coeff) c pa pb p1 p3 p2 p123 wf1 wf2 wf3

    let print_add_dscalar4_km c pa pb wf1 wf2 wf3 p1 p2 p3 p123 fusion (coeff, contraction) =
      printf "@ + ";
      print_dscalar4_km c pa pb wf1 wf2 wf3 p1 p2 p3 p123 fusion (coeff, contraction)

    let print_current amplitude dictionary rhs =
      match F.coupling rhs with
      | V3 (vertex, fusion, constant) ->
          let ch1, ch2 = children2 rhs in
          let wf1 = multiple_variable amplitude dictionary ch1
          and wf2 = multiple_variable amplitude dictionary ch2
          and p1 = momentum ch1
          and p2 = momentum ch2
          and m1 = CM.mass_symbol (F.flavor ch1)
          and m2 = CM.mass_symbol (F.flavor ch2) in
          let c = CM.constant_symbol constant in
          printf "@, %s " (if (F.sign rhs) < 0 then "-" else "+");
          begin match vertex with

(* Fermionic currents $\bar\psi\fmslash{A}\psi$ and $\bar\psi\phi\psi$
   are handled by the [Fermions] module, since they depend on the
   choice of Feynman rules: Dirac or Majorana. *)

          | FBF (coeff, fb, b, f) ->
              begin match coeff, fb, b, f with
              | _, _, (VLRM|SPM|VAM|VA3M|TVA|TVAM|TLR|TLRM|TRL|TRLM), _ ->
                  let p12 = Printf.sprintf "(-%s-%s)" p1 p2 in
                  Fermions.print_current_mom (coeff, fb, b, f) c wf1 wf2 p1 p2
                      p12 fusion
              | _, _, _, _ ->
                  Fermions.print_current (coeff, fb, b, f) c wf1 wf2 fusion
              end
          | PBP (coeff, f1, b, f2) ->
              Fermions.print_current_p (coeff, f1, b, f2) c wf1 wf2 fusion
          | BBB (coeff, fb1, b, fb2) ->
              Fermions.print_current_b (coeff, fb1, b, fb2) c wf1 wf2 fusion
          | GBG (coeff, fb, b, f) ->  let p12 =
              Printf.sprintf "(-%s-%s)" p1 p2 in
              Fermions.print_current_g (coeff, fb, b, f) c wf1 wf2 p1 p2
                   p12 fusion

(* Table~\ref{tab:dim4-bosons} is a bit misleading, since if includes
   totally antisymmetric structure constants.  The space-time part alone
   is also totally antisymmetric: *)

          | Gauge_Gauge_Gauge coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | (F23|F31|F12) ->
                  printf "g_gg(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | (F32|F13|F21) ->
                  printf "g_gg(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              end

          | I_Gauge_Gauge_Gauge coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | (F23|F31|F12) ->
                  printf "g_gg((0,1)*(%s),%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | (F32|F13|F21) ->
                  printf "g_gg((0,1)*(%s),%s,%s,%s,%s)" c wf2 p2 wf1 p1
              end

(* In [Aux_Gauge_Gauge], we can not rely on antisymmetry alone, because of the
   different Lorentz representations of the auxialiary and the gauge field.
   Instead we have to provide the sign in
   \begin{equation}
     (V_2 \wedge V_3) \cdot T_1 =
       \begin{cases}
          V_2 \cdot (T_1 \cdot V_3) = - V_2 \cdot (V_3 \cdot T_1) & \\
          V_3 \cdot (V_2 \cdot T_1) = - V_3 \cdot (T_1 \cdot V_2) &
       \end{cases}
   \end{equation}
   ourselves. Alternatively, one could provide \verb+g_xg+ mirroring
   \verb+g_gx+. *)

          | Aux_Gauge_Gauge coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | F23 -> printf "x_gg(%s,%s,%s)" c wf1 wf2
              | F32 -> printf "x_gg(%s,%s,%s)" c wf2 wf1
              | F12 -> printf "g_gx(%s,%s,%s)" c wf2 wf1
              | F21 -> printf "g_gx(%s,%s,%s)" c wf1 wf2
              | F13 -> printf "(-1)*g_gx(%s,%s,%s)" c wf2 wf1
              | F31 -> printf "(-1)*g_gx(%s,%s,%s)" c wf1 wf2
              end

(* These cases are symmetric and we just have to juxtapose the correct fields
   and provide parentheses to minimize the number of multiplications. *)

          | Scalar_Vector_Vector coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | (F23|F32) -> printf "%s*(%s*%s)" c wf1 wf2
              | (F12|F13) -> printf "(%s*%s)*%s" c wf1 wf2
              | (F21|F31) -> printf "(%s*%s)*%s" c wf2 wf1
              end

          | Aux_Vector_Vector coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | (F23|F32) -> printf "%s*(%s*%s)" c wf1 wf2
              | (F12|F13) -> printf "(%s*%s)*%s" c wf1 wf2
              | (F21|F31) -> printf "(%s*%s)*%s" c wf2 wf1
              end

(* Even simpler: *)

          | Scalar_Scalar_Scalar coeff ->
              printf "(%s*%s*%s)" (format_coupling coeff c) wf1 wf2

          | Aux_Scalar_Scalar coeff ->
              printf "(%s*%s*%s)" (format_coupling coeff c) wf1 wf2

          | Aux_Scalar_Vector coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | (F13|F31) -> printf "%s*(%s*%s)" c wf1 wf2
              | (F23|F21) -> printf "(%s*%s)*%s" c wf1 wf2
              | (F32|F12) -> printf "(%s*%s)*%s" c wf2 wf1
              end

          | Vector_Scalar_Scalar coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | F23 -> printf "v_ss(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F32 -> printf "v_ss(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              | F12 -> printf "s_vs(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F21 -> printf "s_vs(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              | F13 -> printf "(-1)*s_vs(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F31 -> printf "(-1)*s_vs(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              end

          | Graviton_Scalar_Scalar coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | F12 -> printf "s_gravs(%s,%s,-(%s+%s),%s,%s,%s)" c m2 p1 p2 p2 wf1 wf2
              | F21 -> printf "s_gravs(%s,%s,-(%s+%s),%s,%s,%s)" c m1 p1 p2 p1 wf2 wf1
              | F13 -> printf "s_gravs(%s,%s,%s,-(%s+%s),%s,%s)" c m2 p2 p1 p2 wf1 wf2
              | F31 -> printf "s_gravs(%s,%s,%s,-(%s+%s),%s,%s)" c m1 p1 p1 p2 wf2 wf1
              | F23 -> printf "grav_ss(%s,%s,%s,%s,%s,%s)" c m1 p1 p2 wf1 wf2
              | F32 -> printf "grav_ss(%s,%s,%s,%s,%s,%s)" c m1 p2 p1 wf2 wf1
              end

(* In producing a vector in the fusion we always contract the rightmost index with the
   vector wavefunction from [rhs]. So the first momentum is always the one of the
   vector boson produced in the fusion, while the second one is that from the [rhs].
   This makes the cases [F12] and [F13] as well as [F21] and [F31] equal. In principle,
   we could have already done this for the [Graviton_Scalar_Scalar] case. *)


          | Graviton_Vector_Vector coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | (F12|F13) -> printf "v_gravv(%s,%s,-(%s+%s),%s,%s,%s)" c m2 p1 p2 p2 wf1 wf2
              | (F21|F31) -> printf "v_gravv(%s,%s,-(%s+%s),%s,%s,%s)" c m1 p1 p2 p1 wf2 wf1
              | F23 -> printf "grav_vv(%s,%s,%s,%s,%s,%s)" c m1 p1 p2 wf1 wf2
              | F32 -> printf "grav_vv(%s,%s,%s,%s,%s,%s)" c m1 p2 p1 wf2 wf1
              end

          | Graviton_Spinor_Spinor coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | F23 -> printf "f_gravf(%s,%s,-(%s+%s),(-%s),%s,%s)" c m2 p1 p2 p2 wf1 wf2
              | F32 -> printf "f_gravf(%s,%s,-(%s+%s),(-%s),%s,%s)" c m1 p1 p2 p1 wf2 wf1
              | F12 -> printf "f_fgrav(%s,%s,%s,%s+%s,%s,%s)" c m1 p1 p1 p2 wf1 wf2
              | F21 -> printf "f_fgrav(%s,%s,%s,%s+%s,%s,%s)" c m2 p2 p1 p2 wf2 wf1
              | F13 -> printf "grav_ff(%s,%s,%s,(-%s),%s,%s)" c m1 p1 p2 wf1 wf2
              | F31 -> printf "grav_ff(%s,%s,%s,(-%s),%s,%s)" c m1 p2 p1 wf2 wf1
              end

          | Dim4_Vector_Vector_Vector_T coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | F23 -> printf "tkv_vv(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F32 -> printf "tkv_vv(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              | F12 -> printf "tv_kvv(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F21 -> printf "tv_kvv(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              | F13 -> printf "(-1)*tv_kvv(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F31 -> printf "(-1)*tv_kvv(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              end

          | Dim4_Vector_Vector_Vector_L coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | F23 -> printf "lkv_vv(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F32 -> printf "lkv_vv(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              | F12 | F13 -> printf "lv_kvv(%s,%s,%s,%s)" c wf1 p1 wf2
              | F21 | F31 -> printf "lv_kvv(%s,%s,%s,%s)" c wf2 p2 wf1
              end

          | Dim6_Gauge_Gauge_Gauge coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | F23 | F31 | F12 ->
                  printf "kg_kgkg(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F32 | F13 | F21 ->
                  printf "kg_kgkg(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              end

          | Dim4_Vector_Vector_Vector_T5 coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | F23 -> printf "t5kv_vv(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F32 -> printf "t5kv_vv(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              | F12 | F13 -> printf "t5v_kvv(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F21 | F31 -> printf "t5v_kvv(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              end

          | Dim4_Vector_Vector_Vector_L5 coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | F23 -> printf "l5kv_vv(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F32 -> printf "l5kv_vv(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              | F12 -> printf "l5v_kvv(%s,%s,%s,%s)" c wf1 p1 wf2
              | F21 -> printf "l5v_kvv(%s,%s,%s,%s)" c wf2 p2 wf1
              | F13 -> printf "(-1)*l5v_kvv(%s,%s,%s,%s)" c wf1 p1 wf2
              | F31 -> printf "(-1)*l5v_kvv(%s,%s,%s,%s)" c wf2 p2 wf1
              end

          | Dim6_Gauge_Gauge_Gauge_5 coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | F23 -> printf "kg5_kgkg(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F32 -> printf "kg5_kgkg(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              | F12 -> printf "kg_kg5kg(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F21 -> printf "kg_kg5kg(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              | F13 -> printf "(-1)*kg_kg5kg(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F31 -> printf "(-1)*kg_kg5kg(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              end

          | Aux_DScalar_DScalar coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | (F23|F32) ->
                  printf "%s*(%s*%s)*(%s*%s)" c p1 p2 wf1 wf2
              | (F12|F13) ->
                  printf "%s*(-((%s+%s)*%s))*(%s*%s)" c p1 p2 p2 wf1 wf2
              | (F21|F31) ->
                  printf "%s*(-((%s+%s)*%s))*(%s*%s)" c p1 p2 p1 wf1 wf2
              end

          | Aux_Vector_DScalar coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | F23 -> printf "%s*(%s*%s)*%s" c wf1 p2 wf2
              | F32 -> printf "%s*(%s*%s)*%s" c wf2 p1 wf1
              | F12 -> printf "%s*(-((%s+%s)*%s))*%s" c p1 p2 wf2 wf1
              | F21 -> printf "%s*(-((%s+%s)*%s))*%s" c p1 p2 wf1 wf2
              | (F13|F31) -> printf "(-(%s+%s))*(%s*%s*%s)" p1 p2 c wf1 wf2
              end

          | Dim5_Scalar_Gauge2 coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | (F23|F32) -> printf "(%s)*((%s*%s)*(%s*%s) - (%s*%s)*(%s*%s))"
                    c p1 wf2 p2 wf1 p1 p2 wf2 wf1
              | (F12|F13) -> printf "(%s)*%s*((-((%s+%s)*%s))*%s - ((-(%s+%s)*%s))*%s)"
                    c wf1 p1 p2 wf2 p2 p1 p2 p2 wf2
              | (F21|F31) -> printf "(%s)*%s*((-((%s+%s)*%s))*%s - ((-(%s+%s)*%s))*%s)"
                    c wf2 p2 p1 wf1 p1 p1 p2 p1 wf1
              end

          | Dim5_Scalar_Gauge2_Skew coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | (F23|F32) -> printf "(- phi_vv (%s, %s, %s, %s, %s))" c p1 p2 wf1 wf2
              | (F12|F13) -> printf "(- v_phiv (%s, %s, %s, %s, %s))" c wf1 p1 p2 wf2
              | (F21|F31) -> printf "v_phiv (%s, %s, %s, %s, %s)" c wf2 p1 p2 wf1
              end

          | Dim5_Scalar_Vector_Vector_T coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | (F23|F32) -> printf "(%s)*(%s*%s)*(%s*%s)" c p1 wf2 p2 wf1
              | (F12|F13) -> printf "(%s)*%s*(-((%s+%s)*%s))*%s" c wf1 p1 p2 wf2 p2
              | (F21|F31) -> printf "(%s)*%s*(-((%s+%s)*%s))*%s" c wf2 p2 p1 wf1 p1
              end

          | Dim5_Scalar_Vector_Vector_U coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | (F23|F32) -> printf "phi_u_vv (%s, %s, %s, %s, %s)" c p1 p2 wf1 wf2
              | (F12|F13) -> printf "v_u_phiv (%s, %s, %s, %s, %s)" c wf1 p1 p2 wf2
              | (F21|F31) -> printf "v_u_phiv (%s, %s, %s, %s, %s)" c wf2 p2 p1 wf1
              end

          | Dim5_Scalar_Vector_Vector_TU coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | F23 -> printf "(%s)*((%s*%s)*(-(%s+%s)*%s) - (-(%s+%s)*%s)*(%s*%s))"
                    c p1 wf2 p1 p2 wf1 p1 p2 p1 wf1 wf2
              | F32 -> printf "(%s)*((%s*%s)*(-(%s+%s)*%s) - (-(%s+%s)*%s)*(%s*%s))"
                    c p2 wf1 p1 p2 wf2 p1 p2 p2 wf1 wf2
              | F12 -> printf "(%s)*%s*((%s*%s)*%s - (%s*%s)*%s)"
                    c wf1 p1 wf2 p2 p1 p2 wf2
              | F21 -> printf "(%s)*%s*((%s*%s)*%s - (%s*%s)*%s)"
                    c wf2 p2 wf1 p1 p1 p2 wf1
              | F13 -> printf "(%s)*%s*((-(%s+%s)*%s)*%s - (-(%s+%s)*%s)*%s)"
                    c wf1 p1 p2 wf2 p1 p1 p2 p1 wf2
              | F31 -> printf "(%s)*%s*((-(%s+%s)*%s)*%s - (-(%s+%s)*%s)*%s)"
                    c wf2 p1 p2 wf1 p2 p1 p2 p2 wf1
              end

          | Dim5_Scalar_Scalar2 coeff->
              let c = format_coupling coeff c in
	      begin match fusion with
	      | (F23|F32) -> printf "phi_dim5s2(%s, %s ,%s, %s, %s)" 
	          c wf1 p1 wf2 p2 
	      | (F12|F13) -> let p12 = Printf.sprintf "(-%s-%s)" p1 p2 in
	          printf "phi_dim5s2(%s,%s,%s,%s,%s)" c wf1 p12 wf2 p2
	      | (F21|F31) -> let p12 = Printf.sprintf "(-%s-%s)" p1 p2 in
	          printf "phi_dim5s2(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p12
	      end

          | Scalar_Vector_Vector_t coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | (F23|F32) -> printf "s_vv_t(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | (F12|F13) -> printf "v_sv_t(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | (F21|F31) -> printf "v_sv_t(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              end

          | Dim6_Vector_Vector_Vector_T coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | F23 -> printf "(%s)*(%s*%s)*(%s*%s)*(%s-%s)" c p2 wf1 p1 wf2 p1 p2
              | F32 -> printf "(%s)*(%s*%s)*(%s*%s)*(%s-%s)" c p1 wf2 p2 wf1 p2 p1
              | (F12|F13) -> printf "(%s)*((%s+2*%s)*%s)*(-((%s+%s)*%s))*%s"
                    c p1 p2 wf1 p1 p2 wf2 p2
              | (F21|F31) -> printf "(%s)*((-((%s+%s)*%s))*(%s+2*%s)*%s)*%s"
                    c p2 p1 wf1 p2 p1 wf2 p1
              end

          | Tensor_2_Vector_Vector coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | (F23|F32) -> printf "t2_vv(%s,%s,%s)" c wf1 wf2
              | (F12|F13) -> printf "v_t2v(%s,%s,%s)" c wf1 wf2
              | (F21|F31) -> printf "v_t2v(%s,%s,%s)" c wf2 wf1
              end

          | Tensor_2_Scalar_Scalar coeff->
              let c = format_coupling coeff c in
              begin match fusion with
              | (F23|F32) -> printf "t2_phi2(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | (F12|F13) -> printf "phi_t2phi(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | (F21|F31) -> printf "phi_t2phi(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              end

          | Tensor_2_Vector_Vector_1 coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | (F23|F32) -> printf "t2_vv_1(%s,%s,%s)" c wf1 wf2
              | (F12|F13) -> printf "v_t2v_1(%s,%s,%s)" c wf1 wf2
              | (F21|F31) -> printf "v_t2v_1(%s,%s,%s)" c wf2 wf1
              end

          | Tensor_2_Vector_Vector_cf coeff->
              let c = format_coupling coeff c in
              begin match fusion with
              | (F23|F32) -> printf "t2_vv_cf(%s,%s,%s)" c wf1 wf2 
              | (F12|F13) -> printf "v_t2v_cf(%s,%s,%s)" c wf1 wf2
              | (F21|F31) -> printf "v_t2v_cf(%s,%s,%s)" c wf2 wf1
              end

	  | Tensor_2_Scalar_Scalar_cf coeff->
              let c = format_coupling coeff c in
              begin match fusion with
              | (F23|F32) -> printf "t2_phi2_cf(%s,%s,%s,%s, %s)" c wf1 p1 wf2 p2 
              | (F12|F13) -> printf "phi_t2phi_cf(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | (F21|F31) -> printf "phi_t2phi_cf(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              end

          | Dim5_Tensor_2_Vector_Vector_1 coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | (F23|F32) -> printf "t2_vv_d5_1(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | (F12|F13) -> printf "v_t2v_d5_1(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | (F21|F31) -> printf "v_t2v_d5_1(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              end

         | Tensor_2_Vector_Vector_t coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | (F23|F32) -> printf "t2_vv_t(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | (F12|F13) -> printf "v_t2v_t(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | (F21|F31) -> printf "v_t2v_t(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              end

          | Dim5_Tensor_2_Vector_Vector_2 coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | F23 -> printf "t2_vv_d5_2(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F32 -> printf "t2_vv_d5_2(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              | (F12|F13) -> printf "v_t2v_d5_2(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | (F21|F31) -> printf "v_t2v_d5_2(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              end

          | TensorVector_Vector_Vector coeff->
              let c = format_coupling coeff c in
              begin match fusion with
              | (F23|F32) -> printf "dv_vv(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2 
              | (F12|F13) -> printf "v_dvv(%s,%s,%s,%s)" c wf1 p1 wf2 
              | (F21|F31) -> printf "v_dvv(%s,%s,%s,%s)" c wf2 p2 wf1 
              end

          | TensorVector_Vector_Vector_cf coeff->
              let c = format_coupling coeff c in
              begin match fusion with
              | (F23|F32) -> printf "dv_vv_cf(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2 
              | (F12|F13) -> printf "v_dvv_cf(%s,%s,%s,%s)" c wf1 p1 wf2
              | (F21|F31) -> printf "v_dvv_cf(%s,%s,%s,%s)" c wf2 p2 wf1
              end

          | TensorVector_Scalar_Scalar coeff->
              let c = format_coupling coeff c in
              begin match fusion with
              | (F23|F32) -> printf "dv_phi2(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2 
              | (F12|F13) -> printf "phi_dvphi(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | (F21|F31) -> printf "phi_dvphi(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              end

          | TensorVector_Scalar_Scalar_cf coeff->
              let c = format_coupling coeff c in
              begin match fusion with
              | (F23|F32) -> printf "dv_phi2_cf(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2 
              | (F12|F13) -> printf "phi_dvphi_cf(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | (F21|F31) -> printf "phi_dvphi_cf(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              end

          | TensorScalar_Vector_Vector coeff->
              let c = format_coupling coeff c in
              begin match fusion with
              | (F23|F32) -> printf "tphi_vv(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2 
              | (F12|F13) -> printf "v_tphiv(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | (F21|F31) -> printf "v_tphiv(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              end

          | TensorScalar_Vector_Vector_cf coeff->
              let c = format_coupling coeff c in
              begin match fusion with
              | (F23|F32) -> printf "tphi_vv_cf(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2 
              | (F12|F13) -> printf "v_tphiv_cf(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | (F21|F31) -> printf "v_tphiv_cf(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              end

          | TensorScalar_Scalar_Scalar coeff->
              let c = format_coupling coeff c in
              begin match fusion with
              | (F23|F32) -> printf "tphi_ss(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2 
              | (F12|F13) -> printf "s_tphis(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | (F21|F31) -> printf "s_tphis(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              end

          | TensorScalar_Scalar_Scalar_cf coeff->
              let c = format_coupling coeff c in
              begin match fusion with
              | (F23|F32) -> printf "tphi_ss_cf(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2 
              | (F12|F13) -> printf "s_tphis_cf(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | (F21|F31) -> printf "s_tphis_cf(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              end

          | Dim7_Tensor_2_Vector_Vector_T coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | F23 -> printf "t2_vv_d7(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F32 -> printf "t2_vv_d7(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              | (F12|F13) -> printf "v_t2v_d7(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | (F21|F31) -> printf "v_t2v_d7(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              end

          | Dim6_Scalar_Vector_Vector_D coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | (F23|F32) -> printf "s_vv_6D(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | (F12|F13) -> printf "v_sv_6D(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | (F21|F31) -> printf "v_sv_6D(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              end

          | Dim6_Scalar_Vector_Vector_DP coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | (F23|F32) -> printf "s_vv_6DP(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | (F12|F13) -> printf "v_sv_6DP(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | (F21|F31) -> printf "v_sv_6DP(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              end

          | Dim6_HAZ_D coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | F23 -> printf "h_az_D(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F32 -> printf "h_az_D(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              | F13 -> printf "a_hz_D(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F31 -> printf "a_hz_D(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              | F12 -> printf "z_ah_D(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              | F21 -> printf "z_ah_D(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              end
          | Dim6_HAZ_DP coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | F23 -> printf "h_az_DP(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F32 -> printf "h_az_DP(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              | F13 -> printf "a_hz_DP(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F31 -> printf "a_hz_DP(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              | F12 -> printf "z_ah_DP(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              | F21 -> printf "z_ah_DP(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              end
          | Gauge_Gauge_Gauge_i coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | F23 -> printf "g_gg_23(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F32 -> printf "g_gg_23(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              | F13 -> printf "g_gg_13(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F31 -> printf "g_gg_13(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              | F12 -> printf "(-1) * g_gg_13(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F21 -> printf "(-1) * g_gg_13(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              end

          | Dim6_GGG coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | F23 -> printf "g_gg_6(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F32 -> printf "g_gg_6(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              | F12 -> printf "g_gg_6(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F21 -> printf "g_gg_6(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              | F13 -> printf "(-1) * g_gg_6(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F31 -> printf "(-1) * g_gg_6(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              end 

          | Dim6_AWW_DP coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | F23 -> printf "a_ww_DP(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F32 -> printf "a_ww_DP(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              | F13 -> printf "w_aw_DP(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F31 -> printf "w_aw_DP(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              | F12 -> printf "(-1) * w_aw_DP(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F21 -> printf "(-1) * w_aw_DP(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              end

          | Dim6_AWW_DW coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | F23 -> printf "a_ww_DW(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F32 -> printf "a_ww_DW(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              | F13 -> printf "(-1) * a_ww_DW(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F31 -> printf "(-1) * a_ww_DW(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              | F12 -> printf "a_ww_DW(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F21 -> printf "a_ww_DW(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              end

          | Dim6_Gauge_Gauge_Gauge_i coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | F23 | F31 | F12 ->
                  printf "kg_kgkg_i(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F32 | F13 | F21 ->
                  printf "kg_kgkg_i(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              end

          | Dim6_HHH coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | (F23|F32|F12|F21|F13|F31) -> 
		printf "h_hh_6(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              end

          | Dim6_WWZ_DPWDW coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | F23 -> printf "w_wz_DPW(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F32 -> printf "w_wz_DPW(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              | F13 -> printf "(-1) * w_wz_DPW(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F31 -> printf "(-1) * w_wz_DPW(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              | F12 -> printf "z_ww_DPW(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F21 -> printf "z_ww_DPW(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              end
          | Dim6_WWZ_DW coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | F23 -> printf "w_wz_DW(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F32 -> printf "w_wz_DW(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              | F13 -> printf "(-1) * w_wz_DW(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F31 -> printf "(-1) * w_wz_DW(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              | F12 -> printf "z_ww_DW(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F21 -> printf "z_ww_DW(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              end
          | Dim6_WWZ_D coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | F23 -> printf "w_wz_D(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F32 -> printf "w_wz_D(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              | F13 -> printf "(-1) * w_wz_D(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F31 -> printf "(-1) * w_wz_D(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              | F12 -> printf "z_ww_D(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | F21 -> printf "z_ww_D(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              end

(*i
          | Dim6_Glu_Glu_Glu coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | (F23|F31|F12) -> 
                   printf "g_gg_glu(%s,%s,%s,%s,%s)" c wf1 p1 wf2 p2
              | (F32|F13|F21) -> 
                   printf "g_gg_glu(%s,%s,%s,%s,%s)" c wf2 p2 wf1 p1
              end   
i*)

          end

(* Flip the sign to account for the~$\mathrm{i}^2$ relative to diagrams
   with only cubic couplings.
   \label{hack:sign(V4)} *)
(* \begin{dubious}
     That's an \emph{slightly dangerous} hack!!!  How do we accnount
     for such signs when treating $n$-ary vertices uniformly?
   \end{dubious} *)

      | V4 (vertex, fusion, constant) ->
          let c = CM.constant_symbol constant
          and ch1, ch2, ch3 = children3 rhs in
          let wf1 = multiple_variable amplitude dictionary ch1
          and wf2 = multiple_variable amplitude dictionary ch2
          and wf3 = multiple_variable amplitude dictionary ch3
          and p1 = momentum ch1
          and p2 = momentum ch2
          and p3 = momentum ch3 in
          printf "@, %s " (if (F.sign rhs) < 0 then "+" else "-");
          begin match vertex with
          | Scalar4 coeff ->
              printf "(%s*%s*%s*%s)" (format_coupling coeff c) wf1 wf2 wf3
          | Scalar2_Vector2 coeff ->
              let c = format_coupling coeff c in
              begin match fusion with
              | F134 | F143 | F234 | F243 ->
                  printf "%s*%s*(%s*%s)" c wf1 wf2 wf3
              | F314 | F413 | F324 | F423 ->
                  printf "%s*%s*(%s*%s)" c wf2 wf1 wf3
              | F341 | F431 | F342 | F432 ->
                  printf "%s*%s*(%s*%s)" c wf3 wf1 wf2
              | F312 | F321 | F412 | F421 ->
                  printf "(%s*%s*%s)*%s" c wf2 wf3 wf1
              | F231 | F132 | F241 | F142 ->
                  printf "(%s*%s*%s)*%s" c wf1 wf3 wf2
              | F123 | F213 | F124 | F214 ->
                  printf "(%s*%s*%s)*%s" c wf1 wf2 wf3
              end
          | Vector4 contractions ->
              begin match contractions with
              | [] -> invalid_arg "Targets.print_current: Vector4 []"
              | head :: tail ->
                  printf "(";
                  print_vector4 c wf1 wf2 wf3 fusion head;
                  List.iter (print_add_vector4 c wf1 wf2 wf3 fusion) tail;
                  printf ")"
              end
          | Dim8_Vector4_t_0 contractions ->
              begin match contractions with
              | [] -> invalid_arg "Targets.print_current: Vector4 []"
              | head :: tail ->
                  print_vector4_t_0 c wf1 p1 wf2 p2 wf3 p3 fusion head;
                  List.iter (print_add_vector4 c wf1 wf2 wf3 fusion) tail;
              end
          | Dim8_Vector4_t_1 contractions ->
              begin match contractions with
              | [] -> invalid_arg "Targets.print_current: Vector4 []"
              | head :: tail ->
                  print_vector4_t_1 c wf1 p1 wf2 p2 wf3 p3 fusion head;
                  List.iter (print_add_vector4 c wf1 wf2 wf3 fusion) tail;
              end              
          | Dim8_Vector4_t_2 contractions ->
              begin match contractions with
              | [] -> invalid_arg "Targets.print_current: Vector4 []"
              | head :: tail ->
                  print_vector4_t_2 c wf1 p1 wf2 p2 wf3 p3 fusion head;
                  List.iter (print_add_vector4 c wf1 wf2 wf3 fusion) tail;
              end
          | Dim8_Vector4_m_0 contractions ->
              begin match contractions with
              | [] -> invalid_arg "Targets.print_current: Vector4 []"
              | head :: tail ->
                  print_vector4_m_0 c wf1 p1 wf2 p2 wf3 p3 fusion head;
                  List.iter (print_add_vector4 c wf1 wf2 wf3 fusion) tail;
              end
          | Dim8_Vector4_m_1 contractions ->
              begin match contractions with
              | [] -> invalid_arg "Targets.print_current: Vector4 []"
              | head :: tail ->
                  print_vector4_m_1 c wf1 p1 wf2 p2 wf3 p3 fusion head;
                  List.iter (print_add_vector4 c wf1 wf2 wf3 fusion) tail;
              end
          | Dim8_Vector4_m_7 contractions ->
              begin match contractions with
              | [] -> invalid_arg "Targets.print_current: Vector4 []"
              | head :: tail ->
                  print_vector4_m_7 c wf1 p1 wf2 p2 wf3 p3 fusion head;
                  List.iter (print_add_vector4 c wf1 wf2 wf3 fusion) tail;
              end    
          | Vector4_K_Matrix_tho (_, poles) ->
              let pa, pb =
                begin match fusion with
                | (F341|F431|F342|F432|F123|F213|F124|F214) -> (p1, p2)
                | (F134|F143|F234|F243|F312|F321|F412|F421) -> (p2, p3)
                | (F314|F413|F324|F423|F132|F231|F142|F241) -> (p1, p3)
                end in
              printf "(%s*(%s*%s)*(%s*%s)*(%s*%s)@,*("
                c p1 wf1 p2 wf2 p3 wf3;
              List.iter (fun (coeff, pole) ->
                printf "+%s/((%s+%s)*(%s+%s)-%s)"
                  (CM.constant_symbol coeff) pa pb pa pb
                  (CM.constant_symbol pole))
                poles;
              printf ")*(-%s-%s-%s))" p1 p2 p3
          | Vector4_K_Matrix_jr (disc, contractions) ->
              let pa, pb =
                begin match disc, fusion with
                | 3, (F143|F413|F142|F412|F321|F231|F324|F234) -> (p1, p2)
                | 3, (F314|F341|F214|F241|F132|F123|F432|F423) -> (p2, p3)
                | 3, (F134|F431|F124|F421|F312|F213|F342|F243) -> (p1, p3)
                | _, (F341|F431|F342|F432|F123|F213|F124|F214) -> (p1, p2)
                | _, (F134|F143|F234|F243|F312|F321|F412|F421) -> (p2, p3)
                | _, (F314|F413|F324|F423|F132|F231|F142|F241) -> (p1, p3)
                end in
              begin match contractions with
              | [] -> invalid_arg "Targets.print_current: Vector4_K_Matrix_jr []"
              | head :: tail ->
                  printf "(";
                  print_vector4_km c pa pb wf1 wf2 wf3 fusion head;
                  List.iter (print_add_vector4_km c pa pb wf1 wf2 wf3 fusion)
                    tail;
                  printf ")"
              end
          | Vector4_K_Matrix_cf_t0 (disc, contractions) ->
              let pa, pb, pc =
                begin match disc, fusion with
                | 3, (F143|F413|F142|F412|F321|F231|F324|F234) -> (p1, p2, p3)
                | 3, (F314|F341|F214|F241|F132|F123|F432|F423) -> (p2, p3, p1)
                | 3, (F134|F431|F124|F421|F312|F213|F342|F243) -> (p1, p3, p2)
                | _, (F341|F431|F342|F432|F123|F213|F124|F214) -> (p1, p2, p3)
                | _, (F134|F143|F234|F243|F312|F321|F412|F421) -> (p2, p3, p1)
                | _, (F314|F413|F324|F423|F132|F231|F142|F241) -> (p1, p3, p2)
                end in
              begin match contractions with
              | [] -> invalid_arg "Targets.print_current: Vector4_K_Matrix_cf_t0 []"
              | head :: tail ->
                  printf "(";
                  print_vector4_km_t_0 c pa pb wf1 p1 wf2 p2 wf3 p3 fusion head;
                  List.iter (print_add_vector4_km c pa pb wf1 wf2 wf3 fusion)
                    tail;
                  printf ")"
              end              
          | Vector4_K_Matrix_cf_t1 (disc, contractions) ->
              let pa, pb =
                begin match disc, fusion with
                | 3, (F143|F413|F142|F412|F321|F231|F324|F234) -> (p1, p2)
                | 3, (F314|F341|F214|F241|F132|F123|F432|F423) -> (p2, p3)
                | 3, (F134|F431|F124|F421|F312|F213|F342|F243) -> (p1, p3)
                | _, (F341|F431|F342|F432|F123|F213|F124|F214) -> (p1, p2)
                | _, (F134|F143|F234|F243|F312|F321|F412|F421) -> (p2, p3)
                | _, (F314|F413|F324|F423|F132|F231|F142|F241) -> (p1, p3)
                end in
              begin match contractions with
              | [] -> invalid_arg "Targets.print_current: Vector4_K_Matrix_cf_t1 []"
              | head :: tail ->
                  printf "(";
                  print_vector4_km_t_1 c pa pb wf1 p1 wf2 p2 wf3 p3 fusion head;
                  List.iter (print_add_vector4_km c pa pb wf1 wf2 wf3 fusion)
                    tail;
                  printf ")"
              end              
          | Vector4_K_Matrix_cf_t2 (disc, contractions) ->
              let pa, pb =
                begin match disc, fusion with
                | 3, (F143|F413|F142|F412|F321|F231|F324|F234) -> (p1, p2)
                | 3, (F314|F341|F214|F241|F132|F123|F432|F423) -> (p2, p3)
                | 3, (F134|F431|F124|F421|F312|F213|F342|F243) -> (p1, p3)
                | _, (F341|F431|F342|F432|F123|F213|F124|F214) -> (p1, p2)
                | _, (F134|F143|F234|F243|F312|F321|F412|F421) -> (p2, p3)
                | _, (F314|F413|F324|F423|F132|F231|F142|F241) -> (p1, p3)
                end in
              begin match contractions with
              | [] -> invalid_arg "Targets.print_current: Vector4_K_Matrix_cf_t2 []"
              | head :: tail ->
                  printf "(";
                  print_vector4_km_t_2 c pa pb wf1 p1 wf2 p2 wf3 p3 fusion head;
                  List.iter (print_add_vector4_km c pa pb wf1 wf2 wf3 fusion)
                    tail;
                  printf ")"
              end
          | Vector4_K_Matrix_cf_t_rsi (disc, contractions) ->
              let pa, pb, pc =
                begin match disc, fusion with
                | 3, (F143|F413|F142|F412|F321|F231|F324|F234) -> (p1, p2, p3)
                | 3, (F314|F341|F214|F241|F132|F123|F432|F423) -> (p2, p3, p1)
                | 3, (F134|F431|F124|F421|F312|F213|F342|F243) -> (p1, p3, p2)
                | _, (F341|F431|F342|F432|F123|F213|F124|F214) -> (p1, p2, p3)
                | _, (F134|F143|F234|F243|F312|F321|F412|F421) -> (p2, p3, p1)
                | _, (F314|F413|F324|F423|F132|F231|F142|F241) -> (p1, p3, p2)
                end in
              begin match contractions with
              | [] -> invalid_arg "Targets.print_current: Vector4_K_Matrix_cf_t_rsi []"
              | head :: tail ->
                  printf "(";
                  print_vector4_km_t_rsi c pa pb pc wf1 p1 wf2 p2 wf3 p3 fusion head;
                  List.iter (print_add_vector4_km c pa pb wf1 wf2 wf3 fusion)
                    tail;
                  printf ")"
              end              
          | Vector4_K_Matrix_cf_m0 (disc, contractions) ->
              let pa, pb =
                begin match disc, fusion with
                | 3, (F143|F413|F142|F412|F321|F231|F324|F234) -> (p1, p2)
                | 3, (F314|F341|F214|F241|F132|F123|F432|F423) -> (p2, p3)
                | 3, (F134|F431|F124|F421|F312|F213|F342|F243) -> (p1, p3)
                | _, (F341|F431|F342|F432|F123|F213|F124|F214) -> (p1, p2)
                | _, (F134|F143|F234|F243|F312|F321|F412|F421) -> (p2, p3)
                | _, (F314|F413|F324|F423|F132|F231|F142|F241) -> (p1, p3)
                end in
              begin match contractions with
              | [] -> invalid_arg "Targets.print_current: Vector4_K_Matrix_cf_m0 []"
              | head :: tail ->
                  printf "(";
                  print_vector4_km_m_0 c pa pb wf1 p1 wf2 p2 wf3 p3 fusion head;
                  List.iter (print_add_vector4_km c pa pb wf1 wf2 wf3 fusion)
                    tail;
                  printf ")"
              end
          | Vector4_K_Matrix_cf_m1 (disc, contractions) ->
              let pa, pb =
                begin match disc, fusion with
                | 3, (F143|F413|F142|F412|F321|F231|F324|F234) -> (p1, p2)
                | 3, (F314|F341|F214|F241|F132|F123|F432|F423) -> (p2, p3)
                | 3, (F134|F431|F124|F421|F312|F213|F342|F243) -> (p1, p3)
                | _, (F341|F431|F342|F432|F123|F213|F124|F214) -> (p1, p2)
                | _, (F134|F143|F234|F243|F312|F321|F412|F421) -> (p2, p3)
                | _, (F314|F413|F324|F423|F132|F231|F142|F241) -> (p1, p3)
                end in
              begin match contractions with
              | [] -> invalid_arg "Targets.print_current: Vector4_K_Matrix_cf_m1 []"
              | head :: tail ->
                  printf "(";
                  print_vector4_km_m_1 c pa pb wf1 p1 wf2 p2 wf3 p3 fusion head;
                  List.iter (print_add_vector4_km c pa pb wf1 wf2 wf3 fusion)
                    tail;
                  printf ")"
              end
          | Vector4_K_Matrix_cf_m7 (disc, contractions) ->
              let pa, pb =
                begin match disc, fusion with
                | 3, (F143|F413|F142|F412|F321|F231|F324|F234) -> (p1, p2)
                | 3, (F314|F341|F214|F241|F132|F123|F432|F423) -> (p2, p3)
                | 3, (F134|F431|F124|F421|F312|F213|F342|F243) -> (p1, p3)
                | _, (F341|F431|F342|F432|F123|F213|F124|F214) -> (p1, p2)
                | _, (F134|F143|F234|F243|F312|F321|F412|F421) -> (p2, p3)
                | _, (F314|F413|F324|F423|F132|F231|F142|F241) -> (p1, p3)
                end in
              begin match contractions with
              | [] -> invalid_arg "Targets.print_current: Vector4_K_Matrix_cf_m7 []"
              | head :: tail ->
                  printf "(";
                  print_vector4_km_m_7 c pa pb wf1 p1 wf2 p2 wf3 p3 fusion head;
                  List.iter (print_add_vector4_km c pa pb wf1 wf2 wf3 fusion)
                    tail;
                  printf ")"
              end    
          | DScalar2_Vector2_K_Matrix_ms (disc, contractions) ->
              let p123 = Printf.sprintf "(-%s-%s-%s)" p1 p2 p3 in
              let pa, pb =
                begin match disc, fusion with
                | 3, (F143|F413|F142|F412|F321|F231|F324|F234) -> (p1, p2)
                | 3, (F314|F341|F214|F241|F132|F123|F432|F423) -> (p2, p3)
                | 3, (F134|F431|F124|F421|F312|F213|F342|F243) -> (p1, p3)
                | 4, (F143|F413|F142|F412|F321|F231|F324|F234) -> (p1, p2)
                | 4, (F314|F341|F214|F241|F132|F123|F432|F423) -> (p2, p3)
                | 4, (F134|F431|F124|F421|F312|F213|F342|F243) -> (p1, p3)
                | 5, (F143|F413|F142|F412|F321|F231|F324|F234) -> (p1, p2)
                | 5, (F314|F341|F214|F241|F132|F123|F432|F423) -> (p2, p3)
                | 5, (F134|F431|F124|F421|F312|F213|F342|F243) -> (p1, p3)
                | 6, (F134|F132|F314|F312|F241|F243|F421|F423) -> (p1, p2)
                | 6, (F213|F413|F231|F431|F124|F324|F142|F342) -> (p2, p3)
                | 6, (F143|F123|F341|F321|F412|F214|F432|F234) -> (p1, p3)
                | 7, (F134|F132|F314|F312|F241|F243|F421|F423) -> (p1, p2)
                | 7, (F213|F413|F231|F431|F124|F324|F142|F342) -> (p2, p3)
                | 7, (F143|F123|F341|F321|F412|F214|F432|F234) -> (p1, p3)
                | 8, (F134|F132|F314|F312|F241|F243|F421|F423) -> (p1, p2)
                | 8, (F213|F413|F231|F431|F124|F324|F142|F342) -> (p2, p3)
                | 8, (F143|F123|F341|F321|F412|F214|F432|F234) -> (p1, p3)
                | _, (F341|F431|F342|F432|F123|F213|F124|F214) -> (p1, p2)
                | _, (F134|F143|F234|F243|F312|F321|F412|F421) -> (p2, p3)
                | _, (F314|F413|F324|F423|F132|F231|F142|F241) -> (p1, p3)
                end in
              begin match contractions with
              | [] -> invalid_arg "Targets.print_current: DScalar2_Vector4_K_Matrix_ms []"
              | head :: tail ->
                  printf "(";
                  print_dscalar2_vector2_km
                    c pa pb wf1 wf2 wf3 p1 p2 p3 p123 fusion head; 
                  List.iter (print_add_dscalar2_vector2_km
                                  c pa pb wf1 wf2 wf3 p1 p2 p3 p123 fusion) 
                    tail;
                  printf ")"
              end
          | DScalar2_Vector2_m_0_K_Matrix_cf (disc, contractions) ->
              let pa, pb =
                begin match disc, fusion with
                | 3, (F143|F413|F142|F412|F321|F231|F324|F234) -> (p1, p2)
                | 3, (F314|F341|F214|F241|F132|F123|F432|F423) -> (p2, p3)
                | 3, (F134|F431|F124|F421|F312|F213|F342|F243) -> (p1, p3)
                | 4, (F143|F413|F142|F412|F321|F231|F324|F234) -> (p1, p2)
                | 4, (F314|F341|F214|F241|F132|F123|F432|F423) -> (p2, p3)
                | 4, (F134|F431|F124|F421|F312|F213|F342|F243) -> (p1, p3)
                | 5, (F143|F413|F142|F412|F321|F231|F324|F234) -> (p1, p2)
                | 5, (F314|F341|F214|F241|F132|F123|F432|F423) -> (p2, p3)
                | 5, (F134|F431|F124|F421|F312|F213|F342|F243) -> (p1, p3)
                | 6, (F134|F132|F314|F312|F241|F243|F421|F423) -> (p1, p2)
                | 6, (F213|F413|F231|F431|F124|F324|F142|F342) -> (p2, p3)
                | 6, (F143|F123|F341|F321|F412|F214|F432|F234) -> (p1, p3)
                | 7, (F134|F132|F314|F312|F241|F243|F421|F423) -> (p1, p2)
                | 7, (F213|F413|F231|F431|F124|F324|F142|F342) -> (p2, p3)
                | 7, (F143|F123|F341|F321|F412|F214|F432|F234) -> (p1, p3)
                | 8, (F134|F132|F314|F312|F241|F243|F421|F423) -> (p1, p2)
                | 8, (F213|F413|F231|F431|F124|F324|F142|F342) -> (p2, p3)
                | 8, (F143|F123|F341|F321|F412|F214|F432|F234) -> (p1, p3)
                | _, (F341|F431|F342|F432|F123|F213|F124|F214) -> (p1, p2)
                | _, (F134|F143|F234|F243|F312|F321|F412|F421) -> (p2, p3)
                | _, (F314|F413|F324|F423|F132|F231|F142|F241) -> (p1, p3)
                end in
              begin match contractions with
              | [] -> invalid_arg "Targets.print_current: DScalar2_Vector4_K_Matrix_cf_m0 []"
              | head :: tail ->
                  printf "(";
                  print_dscalar2_vector2_m_0_km
                    c pa pb wf1 wf2 wf3 p1 p2 p3 fusion head;
                  List.iter (print_add_dscalar2_vector2_m_0_km
                                  c pa pb wf1 wf2 wf3 p1 p2 p3 fusion)
                    tail;
                  printf ")"
              end
          | DScalar2_Vector2_m_1_K_Matrix_cf (disc, contractions) ->
              let pa, pb =
                begin match disc, fusion with
                | 3, (F143|F413|F142|F412|F321|F231|F324|F234) -> (p1, p2)
                | 3, (F314|F341|F214|F241|F132|F123|F432|F423) -> (p2, p3)
                | 3, (F134|F431|F124|F421|F312|F213|F342|F243) -> (p1, p3)
                | 4, (F143|F413|F142|F412|F321|F231|F324|F234) -> (p1, p2)
                | 4, (F314|F341|F214|F241|F132|F123|F432|F423) -> (p2, p3)
                | 4, (F134|F431|F124|F421|F312|F213|F342|F243) -> (p1, p3)
                | 5, (F143|F413|F142|F412|F321|F231|F324|F234) -> (p1, p2)
                | 5, (F314|F341|F214|F241|F132|F123|F432|F423) -> (p2, p3)
                | 5, (F134|F431|F124|F421|F312|F213|F342|F243) -> (p1, p3)
                | 6, (F134|F132|F314|F312|F241|F243|F421|F423) -> (p1, p2)
                | 6, (F213|F413|F231|F431|F124|F324|F142|F342) -> (p2, p3)
                | 6, (F143|F123|F341|F321|F412|F214|F432|F234) -> (p1, p3)
                | 7, (F134|F132|F314|F312|F241|F243|F421|F423) -> (p1, p2)
                | 7, (F213|F413|F231|F431|F124|F324|F142|F342) -> (p2, p3)
                | 7, (F143|F123|F341|F321|F412|F214|F432|F234) -> (p1, p3)
                | 8, (F134|F132|F314|F312|F241|F243|F421|F423) -> (p1, p2)
                | 8, (F213|F413|F231|F431|F124|F324|F142|F342) -> (p2, p3)
                | 8, (F143|F123|F341|F321|F412|F214|F432|F234) -> (p1, p3)
                | _, (F341|F431|F342|F432|F123|F213|F124|F214) -> (p1, p2)
                | _, (F134|F143|F234|F243|F312|F321|F412|F421) -> (p2, p3)
                | _, (F314|F413|F324|F423|F132|F231|F142|F241) -> (p1, p3)
                end in
              begin match contractions with
              | [] -> invalid_arg "Targets.print_current: DScalar2_Vector4_K_Matrix_cf_m1 []"
              | head :: tail ->
                  printf "(";
                  print_dscalar2_vector2_m_1_km
                    c pa pb wf1 wf2 wf3 p1 p2 p3 fusion head;
                  List.iter (print_add_dscalar2_vector2_m_1_km
                                  c pa pb wf1 wf2 wf3 p1 p2 p3 fusion)
                    tail;
                  printf ")"
              end
          | DScalar2_Vector2_m_7_K_Matrix_cf (disc, contractions) ->
              let pa, pb =
                begin match disc, fusion with
                | 3, (F143|F413|F142|F412|F321|F231|F324|F234) -> (p1, p2)
                | 3, (F314|F341|F214|F241|F132|F123|F432|F423) -> (p2, p3)
                | 3, (F134|F431|F124|F421|F312|F213|F342|F243) -> (p1, p3)
                | 4, (F143|F413|F142|F412|F321|F231|F324|F234) -> (p1, p2)
                | 4, (F314|F341|F214|F241|F132|F123|F432|F423) -> (p2, p3)
                | 4, (F134|F431|F124|F421|F312|F213|F342|F243) -> (p1, p3)
                | 5, (F143|F413|F142|F412|F321|F231|F324|F234) -> (p1, p2)
                | 5, (F314|F341|F214|F241|F132|F123|F432|F423) -> (p2, p3)
                | 5, (F134|F431|F124|F421|F312|F213|F342|F243) -> (p1, p3)
                | 6, (F134|F132|F314|F312|F241|F243|F421|F423) -> (p1, p2)
                | 6, (F213|F413|F231|F431|F124|F324|F142|F342) -> (p2, p3)
                | 6, (F143|F123|F341|F321|F412|F214|F432|F234) -> (p1, p3)
                | 7, (F134|F132|F314|F312|F241|F243|F421|F423) -> (p1, p2)
                | 7, (F213|F413|F231|F431|F124|F324|F142|F342) -> (p2, p3)
                | 7, (F143|F123|F341|F321|F412|F214|F432|F234) -> (p1, p3)
                | 8, (F134|F132|F314|F312|F241|F243|F421|F423) -> (p1, p2)
                | 8, (F213|F413|F231|F431|F124|F324|F142|F342) -> (p2, p3)
                | 8, (F143|F123|F341|F321|F412|F214|F432|F234) -> (p1, p3)
                | _, (F341|F431|F342|F432|F123|F213|F124|F214) -> (p1, p2)
                | _, (F134|F143|F234|F243|F312|F321|F412|F421) -> (p2, p3)
                | _, (F314|F413|F324|F423|F132|F231|F142|F241) -> (p1, p3)
                end in
              begin match contractions with
              | [] -> invalid_arg "Targets.print_current: DScalar2_Vector4_K_Matrix_cf_m7 []"
              | head :: tail ->
                  printf "(";
                  print_dscalar2_vector2_m_7_km
                    c pa pb wf1 wf2 wf3 p1 p2 p3 fusion head;
                  List.iter (print_add_dscalar2_vector2_m_7_km
                                  c pa pb wf1 wf2 wf3 p1 p2 p3 fusion)
                    tail;
                  printf ")"
              end    
          | DScalar4_K_Matrix_ms (disc, contractions) ->
              let p123 = Printf.sprintf "(-%s-%s-%s)" p1 p2 p3 in
              let pa, pb =
                begin match disc, fusion with
                | 3, (F143|F413|F142|F412|F321|F231|F324|F234) -> (p1, p2)
                | 3, (F314|F341|F214|F241|F132|F123|F432|F423) -> (p2, p3)
                | 3, (F134|F431|F124|F421|F312|F213|F342|F243) -> (p1, p3)
                | _, (F341|F431|F342|F432|F123|F213|F124|F214) -> (p1, p2)
                | _, (F134|F143|F234|F243|F312|F321|F412|F421) -> (p2, p3)
                | _, (F314|F413|F324|F423|F132|F231|F142|F241) -> (p1, p3)
                end in
              begin match contractions with
              | [] -> invalid_arg "Targets.print_current: DScalar4_K_Matrix_ms []"
              | head :: tail ->
                  printf "(";
                  print_dscalar4_km
                    c pa pb wf1 wf2 wf3 p1 p2 p3 p123 fusion head; 
                  List.iter (print_add_dscalar4_km
                                  c pa pb wf1 wf2 wf3 p1 p2 p3 p123 fusion) 
                    tail;
                  printf ")"
              end
          | Dim8_Scalar2_Vector2_1 coeff ->
              let c = format_coupling coeff c in
                  begin match fusion with
                  | F134 | F143 | F234 | F243 ->
                      printf "phi_phi2v_1(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F314 | F413 | F324 | F423 ->
                      printf "phi_phi2v_1(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F341 | F431 | F342 | F432 ->
                      printf "phi_phi2v_1(%s,%s,%s,%s,%s,%s,%s)" 
                          c wf3 p3 wf2 p2 wf1 p1
                  | F312 | F321 | F412 | F421 ->
	              printf "v_phi2v_1(%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1
                  | F231 | F132 | F241 | F142 ->
	              printf "v_phi2v_1(%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2
                  | F123 | F213 | F124 | F214 ->
	              printf "v_phi2v_1(%s,%s,%s,%s,%s,%s)" 
                          c wf1 p1 wf2 p2 wf3
                  end
          | Dim8_Scalar2_Vector2_2 coeff ->
              let c = format_coupling coeff c in
                  begin match fusion with
                  | F134 | F143 | F234 | F243 ->
                      printf "phi_phi2v_2(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F314 | F413 | F324 | F423 ->
                      printf "phi_phi2v_2(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F341 | F431 | F342 | F432 ->
                      printf "phi_phi2v_2(%s,%s,%s,%s,%s,%s,%s)" 
                          c wf3 p3 wf2 p2 wf1 p1
                  | F312 | F321 | F412 | F421 ->
	              printf "v_phi2v_2(%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1
                  | F231 | F132 | F241 | F142 ->
	              printf "v_phi2v_2(%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2
                  | F123 | F213 | F124 | F214 ->
	              printf "v_phi2v_2(%s,%s,%s,%s,%s,%s)" 
                          c wf1 p1 wf2 p2 wf3
                  end
          | Dim8_Scalar2_Vector2_m_0 coeff ->
              let c = format_coupling coeff c in
                  begin match fusion with
                  | F134 | F143 | F234 | F243 ->
                      printf "phi_phi2v_m_0(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F314 | F413 | F324 | F423 ->
                      printf "phi_phi2v_m_0(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F341 | F431 | F342 | F432 ->
                      printf "phi_phi2v_m_0(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F312 | F321 | F412 | F421 ->
                      printf "v_phi2v_m_0(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F231 | F132 | F241 | F142 ->
                      printf "v_phi2v_m_0(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F123 | F213 | F124 | F214 ->
                      printf "v_phi2v_m_0(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  end
          | Dim8_Scalar2_Vector2_m_1 coeff ->
              let c = format_coupling coeff c in
                  begin match fusion with
                  | F134 | F143 | F234 | F243 ->
                      printf "phi_phi2v_m_1(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F314 | F413 | F324 | F423 ->
                      printf "phi_phi2v_m_1(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F341 | F431 | F342 | F432 ->
                      printf "phi_phi2v_m_1(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F312 | F321 | F412 | F421 ->
                      printf "v_phi2v_m_1(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F231 | F132 | F241 | F142 ->
                      printf "v_phi2v_m_1(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F123 | F213 | F124 | F214 ->
                      printf "v_phi2v_m_1(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  end
          | Dim8_Scalar2_Vector2_m_7 coeff ->
              let c = format_coupling coeff c in
                  begin match fusion with
                  | F134 | F143 | F234 | F243 ->
                      printf "phi_phi2v_m_7(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F314 | F413 | F324 | F423 ->
                      printf "phi_phi2v_m_7(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F341 | F431 | F342 | F432 ->
                      printf "phi_phi2v_m_7(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F312 | F321 | F412 | F421 ->
                      printf "v_phi2v_m_7(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F231 | F132 | F241 | F142 ->
                      printf "v_phi2v_m_7(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F123 | F213 | F124 | F214 ->
                      printf "v_phi2v_m_7(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  end        
          | Dim8_Scalar4 coeff ->
              let c = format_coupling coeff c in
                  begin match fusion with
                      | F134 | F143 | F234 | F243 | F314 | F413 | F324 | F423
                      | F341 | F431 | F342 | F432 | F312 | F321 | F412 | F421
                      | F231 | F132 | F241 | F142 | F123 | F213 | F124 | F214 ->
	                  printf "s_dim8s3 (%s,%s,%s,%s,%s,%s,%s)" 
                              c wf1 p1 wf2 p2 wf3 p3
                  end
          | GBBG (coeff, fb, b, f) ->
              Fermions.print_current_g4 (coeff, fb, b, f) c wf1 wf2 wf3
                   fusion

          | Dim6_H4_P2 coeff ->
              let c = format_coupling coeff c in
                  begin match fusion with
                      | F134 | F143 | F234 | F243 | F314 | F413 | F324 | F423
                      | F341 | F431 | F342 | F432 | F312 | F321 | F412 | F421
                      | F231 | F132 | F241 | F142 | F123 | F213 | F124 | F214 ->
	                  printf "hhhh_p2 (%s,%s,%s,%s,%s,%s,%s)" 
                              c wf1 p1 wf2 p2 wf3 p3
                  end
          | Dim6_AHWW_DPB coeff ->
              let c = format_coupling coeff c in
                  begin match fusion with
                  | F234 -> 
                      printf "a_hww_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F243 -> 
                      printf "a_hww_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F342 -> 
                      printf "a_hww_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F324 -> 
                      printf "a_hww_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F423 -> 
                      printf "a_hww_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F432 -> 
                      printf "a_hww_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F134 ->
                      printf "h_aww_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F143 ->
                      printf "h_aww_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F341 ->
                      printf "h_aww_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F314 ->
                      printf "h_aww_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F413 ->
                      printf "h_aww_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F431 ->
                      printf "h_aww_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F124 ->
                      printf "w_ahw_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F142 ->
                      printf "w_ahw_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F241 ->
                      printf "w_ahw_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F214 ->
                      printf "w_ahw_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F412 ->
                      printf "w_ahw_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F421 ->
                      printf "w_ahw_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F123  ->
                      printf "(-1)*w_ahw_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F132  ->
                      printf "(-1)*w_ahw_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F231  ->
                      printf "(-1)*w_ahw_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F213  ->
                      printf "(-1)*w_ahw_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F312  ->
                      printf "(-1)*w_ahw_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F321  ->
                      printf "(-1)*w_ahw_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
		  end
          | Dim6_AHWW_DPW coeff ->
              let c = format_coupling coeff c in
                  begin match fusion with
                  | F234 -> 
                      printf "a_hww_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F243 -> 
                      printf "a_hww_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F342 -> 
                      printf "a_hww_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F324 -> 
                      printf "a_hww_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F423 -> 
                      printf "a_hww_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F432 -> 
                      printf "a_hww_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F134 ->
                      printf "h_aww_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F143 ->
                      printf "h_aww_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F341 ->
                      printf "h_aww_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F314 ->
                      printf "h_aww_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F413 ->
                      printf "h_aww_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F431 ->
                      printf "h_aww_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F124 ->
                      printf "w_ahw_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F142 ->
                      printf "w_ahw_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F241 ->
                      printf "w_ahw_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F214 ->
                      printf "w_ahw_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F412 ->
                      printf "w_ahw_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F421 ->
                      printf "w_ahw_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F123  ->
                      printf "(-1)*w_ahw_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F132  ->
                      printf "(-1)*w_ahw_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F231  ->
                      printf "(-1)*w_ahw_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F213  ->
                      printf "(-1)*w_ahw_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F312  ->
                      printf "(-1)*w_ahw_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F321  ->
                      printf "(-1)*w_ahw_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
		  end
          | Dim6_AHWW_DW coeff ->
              let c = format_coupling coeff c in
                 begin match fusion with
                  | F234 -> 
                      printf "a_hww_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F243 -> 
                      printf "a_hww_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F342 -> 
                      printf "a_hww_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F324 -> 
                      printf "a_hww_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F423 -> 
                      printf "a_hww_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F432 -> 
                      printf "a_hww_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F134 ->
                      printf "h_aww_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F143 ->
                      printf "h_aww_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F341 ->
                      printf "h_aww_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F314 ->
                      printf "h_aww_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F413 ->
                      printf "h_aww_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F431 ->
                      printf "h_aww_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F124 ->
                      printf "w3_ahw_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F142 ->
                      printf "w3_ahw_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F241 ->
                      printf "w3_ahw_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F214 ->
                      printf "w3_ahw_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F412 ->
                      printf "w3_ahw_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F421 ->
                      printf "w3_ahw_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F123  ->
                      printf "w4_ahw_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F132  ->
                      printf "w4_ahw_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F231  ->
                      printf "w4_ahw_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F213  ->
                      printf "w4_ahw_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F312  ->
                      printf "w4_ahw_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F321  ->
                      printf "w4_ahw_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
(*i               | F234 | F134 | F124 | F123 -> 
                      printf "a_hww_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F243 | F143 | F142 | F132 -> 
                      printf "a_hww_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F342 | F341 | F241 | F231 -> 
                      printf "a_hww_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F324 | F314 | F214 | F213 -> 
                      printf "a_hww_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F423 | F413 | F412 | F312 -> 
                      printf "a_hww_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F432 | F431 | F421 | F321 -> 
                      printf "a_hww_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1 i*)
		  end
            | Dim6_Scalar2_Vector2_D coeff ->
              let c = format_coupling coeff c in
                  begin match fusion with 
                  | F234 | F134 ->
                      printf "h_hww_D(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F243 | F143 ->
                      printf "h_hww_D(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F342 | F341 -> 
                      printf "h_hww_D(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F324 | F314 ->
                      printf "h_hww_D(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F423 | F413 ->
                      printf "h_hww_D(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F432 | F431 ->
                      printf "h_hww_D(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F124 | F123  ->
                      printf "w_hhw_D(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F142 | F132  ->
                      printf "w_hhw_D(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F241 | F231  ->
                      printf "w_hhw_D(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F214 | F213  ->
                      printf "w_hhw_D(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F412 | F312 -> 
                      printf "w_hhw_D(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F421 | F321  ->
                      printf "w_hhw_D(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  end
            | Dim6_Scalar2_Vector2_DP coeff ->
              let c = format_coupling coeff c in
                  begin match fusion with
                  | F234 | F134 -> 
                      printf "h_hww_DP(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F342 | F341  ->
                      printf "h_hww_DP(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F423 | F413  ->
                      printf "h_hww_DP(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F243 | F143  ->
                      printf "h_hww_DP(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F324 | F314  ->
                      printf "h_hww_DP(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F432 | F431  ->
                      printf "h_hww_DP(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F123 | F124 -> 
                      printf "w_hhw_DP(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F231 | F241->
                      printf "w_hhw_DP(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F312 | F412 ->
                      printf "w_hhw_DP(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F132 | F142->
                      printf "w_hhw_DP(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F213 | F214 ->
                      printf "w_hhw_DP(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F321 | F421 ->
                      printf "w_hhw_DP(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
(*i               | F234 ->
                      printf "h_hww_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F243 ->
                      printf "h_hww_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F342  -> 
                      printf "h_hww_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F324 ->
                      printf "h_hww_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F423 ->
                      printf "h_hww_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F432 ->
                      printf "h_hww_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F124 ->
                      printf "w_hhw_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F142 ->
                      printf "w_hhw_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F241 ->
                      printf "w_hhw_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F214 ->
                      printf "w_hhw_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F412 -> 
                      printf "w_hhw_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F421 ->
                      printf "w_hhw_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F134 ->
                      printf "h_hww_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F143 ->
                      printf "h_hww_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F341 -> 
                      printf "h_hww_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F314 ->
                      printf "h_hww_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F413 ->
                      printf "h_hww_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F431 ->
                      printf "h_hww_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F123  ->
                      printf "w_hhw_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F132  ->
                      printf "w_hhw_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F231  ->
                      printf "w_hhw_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F213  ->
                      printf "w_hhw_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F312 -> 
                      printf "w_hhw_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F321  ->
                      printf "w_hhw_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1   i*)
                  end
          | Dim6_Scalar2_Vector2_PB coeff ->
              let c = format_coupling coeff c in
                  begin match fusion with
                  | F234 | F134 -> 
                      printf "h_hvv_PB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F342 | F341  ->
                      printf "h_hvv_PB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F423 | F413  ->
                      printf "h_hvv_PB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F243 | F143  ->
                      printf "h_hvv_PB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F324 | F314  ->
                      printf "h_hvv_PB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F432 | F431  ->
                      printf "h_hvv_PB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F123 | F124 -> 
                      printf "v_hhv_PB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F231 | F241->
                      printf "v_hhv_PB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F312 | F412 ->
                      printf "v_hhv_PB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F132 | F142->
                      printf "v_hhv_PB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F213 | F214 ->
                      printf "v_hhv_PB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F321 | F421 ->
                      printf "v_hhv_PB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  end  

          | Dim6_HHZZ_T coeff ->
              let c = format_coupling coeff c in
                  begin match fusion with
                  | F234 | F134 -> 
                      printf "(%s)*(%s)*(%s)*(%s)"  c wf1 wf2 wf3
                  | F342 | F341  ->
                      printf "(%s)*(%s)*(%s)*(%s)"  c wf3 wf1 wf2
                  | F423 | F413  ->
                      printf "(%s)*(%s)*(%s)*(%s)"  c wf2 wf3 wf1
                  | F243 | F143  ->
                      printf "(%s)*(%s)*(%s)*(%s)"  c wf1 wf3 wf2
                  | F324 | F314  ->
                      printf "(%s)*(%s)*(%s)*(%s)"  c wf2 wf1 wf3 
                  | F432 | F431  ->
                      printf "(%s)*(%s)*(%s)*(%s)"  c wf3 wf2 wf1
                  | F123 | F124 | F231 | F241 | F312 | F412 ->
                      printf "(%s)*(%s)*(%s)*(%s)"  c wf1 wf2 wf3
                  | F132 | F142 | F213 | F214 | F321 | F421 ->
                      printf "(%s)*(%s)*(%s)*(%s)"  c wf1 wf2 wf3
                  end  
          | Dim6_Vector4_DW coeff ->
              let c = format_coupling coeff c in
                  begin match fusion with
                  | F234 | F134 -> 
                      printf "a_aww_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F342 | F341  ->
                      printf "a_aww_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F423 | F413  ->
                      printf "a_aww_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F243 | F143  ->
                      printf "a_aww_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F324 | F314  ->
                      printf "a_aww_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3 
                  | F432 | F431  ->
                      printf "a_aww_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F124 | F123 -> 
                      printf "w_aaw_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F241 | F231 ->
                      printf "w_aaw_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F412 | F312 ->
                      printf "w_aaw_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F142 | F132 ->
                      printf "w_aaw_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F214 | F213 ->
                      printf "w_aaw_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F421 | F321 ->
                      printf "w_aaw_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  end
          | Dim6_Vector4_W coeff ->
              let c = format_coupling coeff c in
                  begin match fusion with
                  | F234 | F134 -> 
                      printf "a_aww_W(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F342 | F341  ->
                      printf "a_aww_W(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F423 | F413  ->
                      printf "a_aww_W(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F243 | F143  ->
                      printf "a_aww_W(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F324 | F314  ->
                      printf "a_aww_W(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F432 | F431  ->
                      printf "a_aww_W(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F123 | F124 -> 
                      printf "w_aaw_W(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F231 | F241->
                      printf "w_aaw_W(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F312 | F412 ->
                      printf "w_aaw_W(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F132 | F142->
                      printf "w_aaw_W(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F213 | F214 ->
                      printf "w_aaw_W(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F321 | F421 ->
                      printf "w_aaw_W(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  end  
            | Dim6_HWWZ_DW coeff ->
              let c = format_coupling coeff c in
                  begin match fusion with 
                  | F234 ->
                      printf "h_wwz_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F243 ->
                      printf "h_wwz_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F342  -> 
                      printf "h_wwz_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F324 ->
                      printf "h_wwz_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F423 ->
                      printf "h_wwz_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F432 ->
                      printf "h_wwz_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F124 ->
                      printf "(-1)*w_hwz_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F142 ->
                      printf "(-1)*w_hwz_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F241 ->
                      printf "(-1)*w_hwz_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F214 ->
                      printf "(-1)*w_hwz_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F412 -> 
                      printf "(-1)*w_hwz_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F421 ->
                      printf "(-1)*w_hwz_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F134 ->
                      printf "w_hwz_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F143 ->
                      printf "w_hwz_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F341 -> 
                      printf "w_hwz_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F314 ->
                      printf "w_hwz_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F413 ->
                      printf "w_hwz_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F431 ->
                      printf "w_hwz_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F123  ->
                      printf "z_hww_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F132  ->
                      printf "z_hww_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F231  ->
                      printf "z_hww_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F213  ->
                      printf "z_hww_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F312 -> 
                      printf "z_hww_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F321  ->
                      printf "z_hww_DW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  end
            | Dim6_HWWZ_DPB coeff ->
              let c = format_coupling coeff c in
                  begin match fusion with 
                  | F234 ->
                      printf "h_wwz_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F243 ->
                      printf "h_wwz_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F342  -> 
                      printf "h_wwz_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F324 ->
                      printf "h_wwz_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F423 ->
                      printf "h_wwz_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F432 ->
                      printf "h_wwz_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F124 ->
                      printf "(-1)*w_hwz_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F142 ->
                      printf "(-1)*w_hwz_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F241 ->
                      printf "(-1)*w_hwz_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F214 ->
                      printf "(-1)*w_hwz_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F412 -> 
                      printf "(-1)*w_hwz_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F421 ->
                      printf "(-1)*w_hwz_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F134 ->
                      printf "w_hwz_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F143 ->
                      printf "w_hwz_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F341 -> 
                      printf "w_hwz_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F314 ->
                      printf "w_hwz_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F413 ->
                      printf "w_hwz_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F431 ->
                      printf "w_hwz_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F123  ->
                      printf "z_hww_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F132  ->
                      printf "z_hww_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F231  ->
                      printf "z_hww_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F213  ->
                      printf "z_hww_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F312 -> 
                      printf "z_hww_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F321  ->
                      printf "z_hww_DPB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  end
            | Dim6_HWWZ_DDPW coeff ->
              let c = format_coupling coeff c in
                  begin match fusion with 
                  | F234 ->
                      printf "h_wwz_DDPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F243 ->
                      printf "h_wwz_DDPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F342  -> 
                      printf "h_wwz_DDPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F324 ->
                      printf "h_wwz_DDPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F423 ->
                      printf "h_wwz_DDPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F432 ->
                      printf "h_wwz_DDPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F124 ->
                      printf "(-1)*w_hwz_DDPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F142 ->
                      printf "(-1)*w_hwz_DDPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F241 ->
                      printf "(-1)*w_hwz_DDPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F214 ->
                      printf "(-1)*w_hwz_DDPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F412 -> 
                      printf "(-1)*w_hwz_DDPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F421 ->
                      printf "(-1)*w_hwz_DDPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F134 ->
                      printf "w_hwz_DDPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F143 ->
                      printf "w_hwz_DDPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F341 -> 
                      printf "w_hwz_DDPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F314 ->
                      printf "w_hwz_DDPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F413 ->
                      printf "w_hwz_DDPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F431 ->
                      printf "w_hwz_DDPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F123  ->
                      printf "z_hww_DDPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F132  ->
                      printf "z_hww_DDPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F231  ->
                      printf "z_hww_DDPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F213  ->
                      printf "z_hww_DDPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F312 -> 
                      printf "z_hww_DDPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F321  ->
                      printf "z_hww_DDPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  end
            | Dim6_HWWZ_DPW coeff ->
              let c = format_coupling coeff c in
                  begin match fusion with 
                  | F234 ->
                      printf "h_wwz_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F243 ->
                      printf "h_wwz_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F342  -> 
                      printf "h_wwz_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F324 ->
                      printf "h_wwz_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F423 ->
                      printf "h_wwz_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F432 ->
                      printf "h_wwz_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F124 ->
                      printf "(-1)*w_hwz_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F142 ->
                      printf "(-1)*w_hwz_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F241 ->
                      printf "(-1)*w_hwz_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F214 ->
                      printf "(-1)*w_hwz_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F412 -> 
                      printf "(-1)*w_hwz_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F421 ->
                      printf "(-1)*w_hwz_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F134 ->
                      printf "w_hwz_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F143 ->
                      printf "w_hwz_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F341 -> 
                      printf "w_hwz_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F314 ->
                      printf "w_hwz_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F413 ->
                      printf "w_hwz_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F431 ->
                      printf "w_hwz_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F123  ->
                      printf "z_hww_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F132  ->
                      printf "z_hww_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F231  ->
                      printf "z_hww_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F213  ->
                      printf "z_hww_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F312 -> 
                      printf "z_hww_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F321  ->
                      printf "z_hww_DPW(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  end
            | Dim6_AHHZ_D coeff ->
              let c = format_coupling coeff c in
                  begin match fusion with 
                  | F234 ->
                      printf "a_hhz_D(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F243 ->
                      printf "a_hhz_D(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F342  -> 
                      printf "a_hhz_D(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F324 ->
                      printf "a_hhz_D(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F423 ->
                      printf "a_hhz_D(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F432 ->
                      printf "a_hhz_D(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F124 ->
                      printf "h_ahz_D(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F142 ->
                      printf "h_ahz_D(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F241 ->
                      printf "h_ahz_D(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F214 ->
                      printf "h_ahz_D(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F412 -> 
                      printf "h_ahz_D(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F421 ->
                      printf "h_ahz_D(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F134 ->
                      printf "h_ahz_D(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F143 ->
                      printf "h_ahz_D(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F341 -> 
                      printf "h_ahz_D(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F314 ->
                      printf "h_ahz_D(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F413 ->
                      printf "h_ahz_D(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F431 ->
                      printf "h_ahz_D(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F123  ->
                      printf "z_ahh_D(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F132  ->
                      printf "z_ahh_D(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F231  ->
                      printf "z_ahh_D(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F213  ->
                      printf "z_ahh_D(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F312 -> 
                      printf "z_ahh_D(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F321  ->
                      printf "z_ahh_D(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  end
            | Dim6_AHHZ_DP coeff ->
              let c = format_coupling coeff c in
                  begin match fusion with 
                  | F234 ->
                      printf "a_hhz_DP(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F243 ->
                      printf "a_hhz_DP(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F342  -> 
                      printf "a_hhz_DP(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F324 ->
                      printf "a_hhz_DP(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F423 ->
                      printf "a_hhz_DP(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F432 ->
                      printf "a_hhz_DP(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F124 ->
                      printf "h_ahz_DP(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F142 ->
                      printf "h_ahz_DP(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F241 ->
                      printf "h_ahz_DP(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F214 ->
                      printf "h_ahz_DP(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F412 -> 
                      printf "h_ahz_DP(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F421 ->
                      printf "h_ahz_DP(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F134 ->
                      printf "h_ahz_DP(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F143 ->
                      printf "h_ahz_DP(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F341 -> 
                      printf "h_ahz_DP(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F314 ->
                      printf "h_ahz_DP(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F413 ->
                      printf "h_ahz_DP(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F431 ->
                      printf "h_ahz_DP(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F123  ->
                      printf "z_ahh_DP(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F132  ->
                      printf "z_ahh_DP(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F231  ->
                      printf "z_ahh_DP(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F213  ->
                      printf "z_ahh_DP(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F312 -> 
                      printf "z_ahh_DP(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F321  ->
                      printf "z_ahh_DP(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  end
            | Dim6_AHHZ_PB coeff ->
              let c = format_coupling coeff c in
                  begin match fusion with 
                  | F234 ->
                      printf "a_hhz_PB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F243 ->
                      printf "a_hhz_PB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F342  -> 
                      printf "a_hhz_PB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F324 ->
                      printf "a_hhz_PB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F423 ->
                      printf "a_hhz_PB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F432 ->
                      printf "a_hhz_PB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F124 ->
                      printf "h_ahz_PB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F142 ->
                      printf "h_ahz_PB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F241 ->
                      printf "h_ahz_PB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F214 ->
                      printf "h_ahz_PB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F412 -> 
                      printf "h_ahz_PB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F421 ->
                      printf "h_ahz_PB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F134 ->
                      printf "h_ahz_PB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F143 ->
                      printf "h_ahz_PB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F341 -> 
                      printf "h_ahz_PB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F314 ->
                      printf "h_ahz_PB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F413 ->
                      printf "h_ahz_PB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F431 ->
                      printf "h_ahz_PB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  | F123  ->
                      printf "z_ahh_PB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf2 p2 wf3 p3
                  | F132  ->
                      printf "z_ahh_PB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf1 p1 wf3 p3 wf2 p2
                  | F231  ->
                      printf "z_ahh_PB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf1 p1 wf2 p2
                  | F213  ->
                      printf "z_ahh_PB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf1 p1 wf3 p3
                  | F312 -> 
                      printf "z_ahh_PB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf2 p2 wf3 p3 wf1 p1
                  | F321  ->
                      printf "z_ahh_PB(%s,%s,%s,%s,%s,%s,%s)"
                          c wf3 p3 wf2 p2 wf1 p1
                  end  

(* \begin{dubious}
     In principle, [p4] could be obtained from the left hand side \ldots
   \end{dubious} *)
          | DScalar4 contractions ->
              let p123 = Printf.sprintf "(-%s-%s-%s)" p1 p2 p3 in
              begin match contractions with
              | [] -> invalid_arg "Targets.print_current: DScalar4 []"
              | head :: tail ->
                  printf "(";
                  print_dscalar4 c wf1 wf2 wf3 p1 p2 p3 p123 fusion head;
                  List.iter (print_add_dscalar4
                               c wf1 wf2 wf3 p1 p2 p3 p123 fusion) tail;
                  printf ")"
              end

          | DScalar2_Vector2 contractions ->
              let p123 = Printf.sprintf "(-%s-%s-%s)" p1 p2 p3 in
              begin match contractions with
              | [] -> invalid_arg "Targets.print_current: DScalar4 []"
              | head :: tail ->
                  printf "(";
                  print_dscalar2_vector2
                    c wf1 wf2 wf3 p1 p2 p3 p123 fusion head;
                  List.iter (print_add_dscalar2_vector2
                               c wf1 wf2 wf3 p1 p2 p3 p123 fusion) tail;
                  printf ")"
              end

          end

      (* \begin{dubious}
           This reproduces the hack on page~\pageref{hack:sign(V4)}
           and gives the correct results up to quartic vertices.
           Make sure that it is also correct in light
           of~\eqref{eq:factors-of-i}, i.\,e.
           \begin{equation*}
             \ii T = \ii^{\#\text{vertices}}\ii^{\#\text{propagators}} \cdots
                   = \ii^{n-2}\ii^{n-3} \cdots
                   = -\ii(-1)^n \cdots
           \end{equation*}
         \end{dubious} *)
      | Vn (UFO (c, v, s, fl, color), fusion, constant) ->
         if Color.Vertex.trivial color then
           let g = CM.constant_symbol constant
           and chn = F.children rhs in
           let wfs = List.map (multiple_variable amplitude dictionary) chn
           and ps = List.map momentum chn in
           let n = List.length fusion in
           let eps = if n mod 2 = 0 then -1 else 1 in
           printf "@, %s " (if (eps * F.sign rhs) < 0 then "-" else "+");
           UFO.Targets.Fortran.fuse c v s fl g wfs ps fusion
         else
           failwith "print_current: nontrivial color structure"

    let print_propagator f p m gamma =
      let minus_third = "(-1.0_" ^ !kind ^ "/3.0_" ^ !kind ^ ")" in
      let w =
        begin match CM.width f with
          | Vanishing | Fudged -> "0.0_" ^ !kind
          | Constant | Complex_Mass -> gamma
          | Timelike -> "wd_tl(" ^ p ^ "," ^ gamma ^ ")"
          | Running -> "wd_run(" ^ p ^ "," ^ m ^ "," ^ gamma ^ ")"
          | Custom f -> f ^ "(" ^ p ^ "," ^ gamma ^ ")"
        end in
      let cms =
	begin match CM.width f with
	  | Complex_Mass -> ".true."
	  | _ -> ".false."
	end in
      match CM.propagator f with
	| Prop_Scalar ->
          printf "pr_phi(%s,%s,%s," p m w
	| Prop_Col_Scalar ->
          printf "%s * pr_phi(%s,%s,%s," minus_third p m w
	| Prop_Ghost -> printf "(0,1) * pr_phi(%s, %s, %s," p m w
	| Prop_Spinor ->
          printf "%s(%s,%s,%s,%s," Fermions.psi_propagator p m w cms
	| Prop_ConjSpinor ->
          printf "%s(%s,%s,%s,%s," Fermions.psibar_propagator p m w cms
	| Prop_Majorana ->
          printf "%s(%s,%s,%s,%s," Fermions.chi_propagator p m w cms
	| Prop_Col_Majorana ->
          printf "%s * %s(%s,%s,%s,%s," minus_third Fermions.chi_propagator p m w cms
	| Prop_Unitarity ->
          printf "pr_unitarity(%s,%s,%s,%s," p m w cms
	| Prop_Col_Unitarity ->
          printf "%s * pr_unitarity(%s,%s,%s,%s," minus_third p m w cms
	| Prop_Feynman ->
          printf "pr_feynman(%s," p
	| Prop_Col_Feynman ->
          printf "%s * pr_feynman(%s," minus_third p
	| Prop_Gauge xi ->
          printf "pr_gauge(%s,%s," p (CM.gauge_symbol xi)
	| Prop_Rxi xi ->
          printf "pr_rxi(%s,%s,%s,%s," p m w (CM.gauge_symbol xi)
	| Prop_Tensor_2 ->
          printf "pr_tensor(%s,%s,%s," p m w
	| Prop_Tensor_pure ->
          printf "pr_tensor_pure(%s,%s,%s," p m w
	| Prop_Vector_pure ->
          printf "pr_vector_pure(%s,%s,%s," p m w
	| Prop_Vectorspinor ->
          printf "pr_grav(%s,%s,%s," p m w
	| Aux_Scalar | Aux_Spinor | Aux_ConjSpinor | Aux_Majorana
	| Aux_Vector | Aux_Tensor_1 -> printf "("
	| Aux_Col_Scalar | Aux_Col_Vector | Aux_Col_Tensor_1 -> printf "%s * (" minus_third
	| Only_Insertion -> printf "("
	| Prop_UFO name ->
          printf "pr_U_%s(%s,%s,%s," name p m w

    let print_projector f p m gamma =
      let minus_third = "(-1.0_" ^ !kind ^ "/3.0_" ^ !kind ^ ")" in
      match CM.propagator f with
      | Prop_Scalar ->
          printf "pj_phi(%s,%s," m gamma
      | Prop_Col_Scalar ->
          printf "%s * pj_phi(%s,%s," minus_third m gamma
      | Prop_Ghost ->
          printf "(0,1) * pj_phi(%s,%s," m gamma
      | Prop_Spinor ->
          printf "%s(%s,%s,%s," Fermions.psi_projector p m gamma
      | Prop_ConjSpinor ->
          printf "%s(%s,%s,%s," Fermions.psibar_projector p m gamma
      | Prop_Majorana ->
          printf "%s(%s,%s,%s," Fermions.chi_projector p m gamma
      | Prop_Col_Majorana ->
          printf "%s * %s(%s,%s,%s," minus_third Fermions.chi_projector p m gamma
      | Prop_Unitarity ->
          printf "pj_unitarity(%s,%s,%s," p m gamma
      | Prop_Col_Unitarity ->
          printf "%s * pj_unitarity(%s,%s,%s," minus_third p m gamma
      | Prop_Feynman | Prop_Col_Feynman ->
          invalid_arg "no on-shell Feynman propagator!"
      | Prop_Gauge _ ->
          invalid_arg "no on-shell massless gauge propagator!"
      | Prop_Rxi _ ->
          invalid_arg "no on-shell Rxi propagator!"
      | Prop_Vectorspinor ->
          printf "pj_grav(%s,%s,%s," p m gamma
      | Prop_Tensor_2 ->
          printf "pj_tensor(%s,%s,%s," p m gamma
      | Prop_Tensor_pure ->
          invalid_arg "no on-shell pure Tensor propagator!"
      | Prop_Vector_pure ->
          invalid_arg "no on-shell pure Vector propagator!"
      | Aux_Scalar | Aux_Spinor | Aux_ConjSpinor | Aux_Majorana
      | Aux_Vector | Aux_Tensor_1 -> printf "("
      | Aux_Col_Scalar | Aux_Col_Vector | Aux_Col_Tensor_1 -> printf "%s * (" minus_third
      | Only_Insertion -> printf "("
      | Prop_UFO name ->
         invalid_arg "no on shell UFO propagator"

    let print_gauss f p m gamma =
      let minus_third = "(-1.0_" ^ !kind ^ "/3.0_" ^ !kind ^ ")" in
      match CM.propagator f with
      | Prop_Scalar ->
          printf "pg_phi(%s,%s,%s," p m gamma
      | Prop_Ghost ->
          printf "(0,1) * pg_phi(%s,%s,%s," p m gamma
      | Prop_Spinor ->
          printf "%s(%s,%s,%s," Fermions.psi_projector p m gamma
      | Prop_ConjSpinor ->
          printf "%s(%s,%s,%s," Fermions.psibar_projector p m gamma
      | Prop_Majorana ->
          printf "%s(%s,%s,%s," Fermions.chi_projector p m gamma
      | Prop_Col_Majorana ->
          printf "%s * %s(%s,%s,%s," minus_third Fermions.chi_projector p m gamma
      | Prop_Unitarity ->
          printf "pg_unitarity(%s,%s,%s," p m gamma
      | Prop_Feynman | Prop_Col_Feynman ->
          invalid_arg "no on-shell Feynman propagator!"
      | Prop_Gauge _ ->
          invalid_arg "no on-shell massless gauge propagator!"
      | Prop_Rxi _ ->
          invalid_arg "no on-shell Rxi propagator!"
      | Prop_Tensor_2 ->
          printf "pg_tensor(%s,%s,%s," p m gamma
      | Prop_Tensor_pure ->
          invalid_arg "no pure tensor propagator!"
      | Prop_Vector_pure ->
          invalid_arg "no pure vector propagator!"
      | Aux_Scalar | Aux_Spinor | Aux_ConjSpinor | Aux_Majorana
      | Aux_Vector | Aux_Tensor_1 -> printf "("
      | Only_Insertion -> printf "("
      | Prop_UFO name ->
         invalid_arg "no UFO gauss insertion"
      | _ -> invalid_arg "targets:print_gauss: not available"

    let print_fusion_diagnostics amplitude dictionary fusion =
      if warn diagnose_gauge then begin
        let lhs = F.lhs fusion in
        let f = F.flavor lhs
        and v = variable lhs
        and p = momentum lhs in
        let mass = CM.mass_symbol f in
        match CM.propagator f with
        | Prop_Gauge _ | Prop_Feynman
        | Prop_Rxi _ | Prop_Unitarity ->
            printf "      @[<2>%s =" v;
            List.iter (print_current amplitude dictionary) (F.rhs fusion); nl ();
            begin match CM.goldstone f with
            | None ->
                printf "      call omega_ward_%s(\"%s\",%s,%s,%s)"
                  (suffix diagnose_gauge) v mass p v; nl ()
            | Some (g, phase) ->
                let gv = add_tag lhs (CM.flavor_symbol g ^ "_" ^ format_p lhs) in
                printf "      call omega_slavnov_%s"
                  (suffix diagnose_gauge);
                printf "(@[\"%s\",%s,%s,%s,@,%s*%s)"
                  v mass p v (format_constant phase) gv; nl ()
            end
        | _ -> ()
      end

    let print_fusion amplitude dictionary fusion =
      let lhs = F.lhs fusion in
      let f = F.flavor lhs in
      printf "      @[<2>%s =@, " (multiple_variable amplitude dictionary lhs);
      if F.on_shell amplitude lhs then
        print_projector f (momentum lhs)
          (CM.mass_symbol f) (CM.width_symbol f)
      else
        if F.is_gauss amplitude lhs then
          print_gauss f (momentum lhs)
            (CM.mass_symbol f) (CM.width_symbol f)
        else
          print_propagator f (momentum lhs)
            (CM.mass_symbol f) (CM.width_symbol f);
      List.iter (print_current amplitude dictionary) (F.rhs fusion);
      printf ")"; nl ()

    let print_momenta seen_momenta amplitude =
      List.fold_left (fun seen f ->
        let wf = F.lhs f in
        let p = F.momentum_list wf in
        if not (PSet.mem p seen) then begin
          let rhs1 = List.hd (F.rhs f) in
          printf "    %s = %s" (momentum wf)
            (String.concat " + "
               (List.map momentum (F.children rhs1))); nl ()
        end;
        PSet.add p seen)
        seen_momenta (F.fusions amplitude)

    let print_fusions dictionary fusions =
      List.iter
        (fun (f, amplitude) ->
          print_fusion_diagnostics amplitude dictionary f;
          print_fusion amplitude dictionary f)
        fusions

(* \begin{dubious}
     The following will need a bit more work, because
     the decision when to [reverse_braket] for UFO models
     with Majorana fermions needs collaboration
     from [UFO.Targets.Fortran.fuse] which is called by
     [print_current].  See the function
     [UFO_targets.Fortran.jrr_print_majorana_current_transposing]
     for illustration (the function is never used and only for
     documentation).
   \end{dubious} *)

    let spins_of_rhs rhs =
      List.map (fun wf -> CM.lorentz (F.flavor wf)) (F.children rhs)

    let spins_of_ket ket =
      match ThoList.uniq (List.map spins_of_rhs ket) with
      | [spins] -> spins
      | [] -> failwith "Targets.Fortran.spins_of_ket: empty"
      | _ -> [] (* HACK! *)

    let print_braket amplitude dictionary name braket =
      let bra = F.bra braket
      and ket = F.ket braket in
      let spin_bra = CM.lorentz (F.flavor bra)
      and spins_ket = spins_of_ket ket in
      let vintage = true (* [F.vintage] *) in
      printf "      @[<2>%s = %s@, + " name name;
      if Fermions.reverse_braket vintage spin_bra spins_ket then
        begin
          printf "@,(";
          List.iter (print_current amplitude dictionary) ket;
          printf ")*%s" (multiple_variable amplitude dictionary bra)
        end
      else
        begin
          printf "%s*@,(" (multiple_variable amplitude dictionary bra);
          List.iter (print_current amplitude dictionary) ket;
          printf ")"
        end;
      nl ()

(* \begin{equation}
   \label{eq:factors-of-i}
     \ii T = \ii^{\#\text{vertices}}\ii^{\#\text{propagators}} \cdots
           = \ii^{n-2}\ii^{n-3} \cdots
           = -\ii(-1)^n \cdots
   \end{equation} *)

(* \begin{dubious}
     [tho:] we write some brakets twice using different names.  Is it useful
     to cache them?
   \end{dubious} *)

    let print_brakets dictionary amplitude =
      let name = flavors_symbol (flavors amplitude) in
      printf "      %s = 0" name; nl ();
      List.iter (print_braket amplitude dictionary name) (F.brakets amplitude);
      let n = List.length (F.externals amplitude) in
      if n mod 2 = 0 then begin
        printf "      @[<2>%s =@, - %s ! %d vertices, %d propagators"
          name name (n - 2) (n - 3); nl ()
      end else begin
        printf "      ! %s = %s ! %d vertices, %d propagators"
          name name (n - 2) (n - 3); nl ()
      end;
      let s = F.symmetry amplitude in
      if s > 1 then
        printf "      @[<2>%s =@, %s@, / sqrt(%d.0_%s) ! symmetry factor" name name s !kind
      else
        printf "      ! unit symmetry factor";
      nl ()

    let print_incoming wf =
      let p = momentum wf
      and s = spin wf
      and f = F.flavor wf in
      let m = CM.mass_symbol f in
      match CM.lorentz f with
      | Scalar -> printf "1"
      | BRS Scalar -> printf "(0,-1) * (%s * %s - %s**2)" p p m
      | Spinor ->
          printf "%s (%s, - %s, %s)" Fermions.psi_incoming m p s
      | BRS Spinor ->
          printf "%s (%s, - %s, %s)" Fermions.brs_psi_incoming m p s
      | ConjSpinor ->
          printf "%s (%s, - %s, %s)" Fermions.psibar_incoming m p s
      | BRS ConjSpinor ->
          printf "%s (%s, - %s, %s)" Fermions.brs_psibar_incoming m p s
      | Majorana ->
          printf "%s (%s, - %s, %s)" Fermions.chi_incoming m p s
      | Maj_Ghost -> printf "ghost (%s, - %s, %s)" m p s
      | BRS Majorana ->
          printf "%s (%s, - %s, %s)" Fermions.brs_chi_incoming m p s
      | Vector | Massive_Vector ->
          printf "eps (%s, - %s, %s)" m p s
(*i   | Ward_Vector -> printf "%s" p   i*)
      | BRS Vector | BRS Massive_Vector -> printf
            "(0,1) * (%s * %s - %s**2) * eps (%s, -%s, %s)" p p m m p s
      | Vectorspinor | BRS Vectorspinor ->
          printf "%s (%s, - %s, %s)" Fermions.grav_incoming m p s
      | Tensor_1 -> invalid_arg "Tensor_1 only internal"
      | Tensor_2 -> printf "eps2 (%s, - %s, %s)" m p s
      | _ -> invalid_arg "no such BRST transformations"

    let print_outgoing wf =
      let p = momentum wf
      and s = spin wf
      and f = F.flavor wf in
      let m = CM.mass_symbol f in
      match CM.lorentz f with
      | Scalar -> printf "1"
      | BRS Scalar -> printf "(0,-1) * (%s * %s - %s**2)" p p m
      | Spinor ->
          printf "%s (%s, %s, %s)" Fermions.psi_outgoing m p s
      | BRS Spinor ->
          printf "%s (%s, %s, %s)" Fermions.brs_psi_outgoing m p s
      | ConjSpinor ->
          printf "%s (%s, %s, %s)" Fermions.psibar_outgoing m p s
      | BRS ConjSpinor ->
          printf "%s (%s, %s, %s)" Fermions.brs_psibar_outgoing m p s
      | Majorana ->
          printf "%s (%s, %s, %s)" Fermions.chi_outgoing m p s
      | BRS Majorana ->
          printf "%s (%s, %s, %s)" Fermions.brs_chi_outgoing m p s
      | Maj_Ghost -> printf "ghost (%s, %s, %s)" m p s
      | Vector | Massive_Vector ->
          printf "conjg (eps (%s, %s, %s))" m p s
(*i   | Ward_Vector -> printf "%s" p   i*)
      | BRS Vector | BRS Massive_Vector -> printf
            "(0,1) * (%s*%s-%s**2) * (conjg (eps (%s, %s, %s)))" p p m m p s
      | Vectorspinor | BRS Vectorspinor ->
          printf "%s (%s, %s, %s)" Fermions.grav_incoming m p s
      | Tensor_1 -> invalid_arg "Tensor_1 only internal"
      | Tensor_2 -> printf "conjg (eps2 (%s, %s, %s))" m p s
      | BRS _ -> invalid_arg "no such BRST transformations"

    (*i unused value
    let twice_spin wf =
      match CM.lorentz (F.flavor wf) with
      | Scalar | BRS Scalar -> "0"
      | Spinor | ConjSpinor | Majorana | Maj_Ghost | Vectorspinor
      | BRS Spinor | BRS ConjSpinor | BRS Majorana | BRS Vectorspinor -> "1"
      | Vector | BRS Vector | Massive_Vector | BRS Massive_Vector -> "2"
      | Tensor_1 -> "2"
      | Tensor_2 -> "4"
      | BRS _ -> invalid_arg "Targets.twice_spin: no such BRST transformation"
     i*)

    (*i unused value
    let print_argument_diagnostics amplitude =
      let externals = (F.externals amplitude) in
      let n = List.length externals
      and masses = List.map (fun wf -> CM.mass_symbol (F.flavor wf)) externals in
      if warn diagnose_arguments then begin
        printf "    call omega_check_arguments_%s (%d, k)"
          (suffix diagnose_arguments) n; nl ()
      end;
      if warn diagnose_momenta then begin
        printf "    @[<2>call omega_check_momenta_%s ((/ "
          (suffix diagnose_momenta);
        print_list masses;
        printf " /), k)"; nl ()
      end
     i*)

    let print_external_momenta amplitude =
      let externals =
        List.combine
          (F.externals amplitude)
          (List.map (fun _ -> true) (F.incoming amplitude) @
           List.map (fun _ -> false) (F.outgoing amplitude)) in
      List.iter (fun (wf, incoming) ->
        if incoming then
          printf "    %s = - k(:,%d) ! incoming"
            (momentum wf) (ext_momentum wf)
        else
          printf "    %s =   k(:,%d) ! outgoing"
            (momentum wf) (ext_momentum wf); nl ()) externals

    let print_externals seen_wfs amplitude =
      let externals =
        List.combine
          (F.externals amplitude)
          (List.map (fun _ -> true) (F.incoming amplitude) @
           List.map (fun _ -> false) (F.outgoing amplitude)) in
      List.fold_left (fun seen (wf, incoming) ->
        if not (WFSet.mem wf seen) then begin
          printf "      @[<2>%s =@, " (variable wf);
          (if incoming then print_incoming else print_outgoing) wf; nl ()
        end;
        WFSet.add wf seen) seen_wfs externals

    (*i unused value
    let flavors_to_string flavors =
      String.concat " " (List.map CM.flavor_to_string flavors)
    i*)

    (*i unused value
    let process_to_string amplitude =
      flavors_to_string (F.incoming amplitude) ^ " -> " ^
      flavors_to_string (F.outgoing amplitude)
    i*)

    let flavors_sans_color_to_string flavors =
      String.concat " " (List.map M.flavor_to_string flavors)

    let process_sans_color_to_string (fin, fout) =
      flavors_sans_color_to_string fin ^ " -> " ^
      flavors_sans_color_to_string fout

    let print_fudge_factor amplitude =
      let name = flavors_symbol (flavors amplitude) in
      List.iter (fun wf ->
        let p = momentum wf
        and f = F.flavor wf in
        match CM.width f with
        | Fudged ->
            let m = CM.mass_symbol f
            and w = CM.width_symbol f in
            printf "      if (%s > 0.0_%s) then" w !kind; nl ();
            printf "        @[<2>%s = %s@ * (%s*%s - %s**2)"
              name name p p m;
            printf "@ / cmplx (%s*%s - %s**2, %s*%s, kind=%s)"
              p p m m w !kind; nl ();
            printf "      end if"; nl ()
        | _ -> ()) (F.s_channel amplitude)

    let num_helicities amplitudes =
      List.length (CF.helicities amplitudes)

(* \thocwmodulesubsection{Spin, Flavor \&\ Color Tables} *)

(* The following abomination is required to keep the number of continuation
   lines as low as possible.  FORTRAN77-style \texttt{DATA} statements
   are actually a bit nicer here, but they are nor available for
   \emph{constant} arrays. *)

(* \begin{dubious}
     We used to have a more elegant design with a sentinel~0 added to each
     initializer, but some revisions of the Compaq/Digital Compiler have a
     bug that causes it to reject this variant.
   \end{dubious} *)

(* \begin{dubious}
     The actual table writing code using \texttt{reshape} should be factored,
     since it's the same algorithm every time.
   \end{dubious} *)

    let print_integer_parameter name value =
      printf "  @[<2>integer, parameter :: %s = %d" name value; nl ()

    let print_real_parameter name value =
      printf "  @[<2>real(kind=%s), parameter :: %s = %d"
        !kind name value; nl ()

    let print_logical_parameter name value =
      printf "  @[<2>logical, parameter :: %s = .%s."
        name (if value then "true" else "false"); nl ()

    let num_particles_in amplitudes =
      match CF.flavors amplitudes with
      | [] -> 0
      | (fin, _) :: _ -> List.length fin

    let num_particles_out amplitudes =
      match CF.flavors amplitudes with
      | [] -> 0
      | (_, fout) :: _ -> List.length fout

    let num_particles amplitudes =
      match CF.flavors amplitudes with
      | [] -> 0
      | (fin, fout) :: _ -> List.length fin + List.length fout

    module CFlow = Color.Flow

    let num_color_flows amplitudes =
      if !amp_triv then
        1
      else
        List.length (CF.color_flows amplitudes)

    let num_color_indices_default = 2 (* Standard model *)

    let num_color_indices amplitudes =
      try CFlow.rank (List.hd (CF.color_flows amplitudes)) with _ -> num_color_indices_default

    let color_to_string c =
      "(" ^ (String.concat "," (List.map (Printf.sprintf "%3d") c)) ^ ")"

    let cflow_to_string cflow =
      String.concat " " (List.map color_to_string (CFlow.in_to_lists cflow)) ^ " -> " ^
      String.concat " " (List.map color_to_string (CFlow.out_to_lists cflow))

    let protected = ", protected" (* Fortran 2003! *)

    (*i unused value
    let print_spin_table_old abbrev name = function
      | [] ->
          printf "  @[<2>integer, dimension(n_prt,0) ::";
          printf "@ table_spin_%s" name; nl ()
      | _ :: tuples' as tuples ->
          ignore (List.fold_left (fun i (tuple1, tuple2) ->
            printf "  @[<2>integer, dimension(n_prt), parameter, private ::";
            printf "@ %s%04d = (/ %s /)" abbrev i
              (String.concat ", " (List.map (Printf.sprintf "%2d") (tuple1 @ tuple2)));
            nl (); succ i) 1 tuples);
          printf
            "  @[<2>integer, dimension(n_prt,n_hel), parameter ::";
          printf "@ table_spin_%s =@ reshape ( (/" name;
          printf "@ %s%04d" abbrev 1;
          ignore (List.fold_left (fun i tuple ->
            printf ",@ %s%04d" abbrev i; succ i) 2 tuples');
          printf "@ /), (/ n_prt, n_hel /) )"; nl ()
    i*)

    let print_spin_table name tuples =
      printf "  @[<2>integer, dimension(n_prt,n_hel), save%s :: table_spin_%s"
        protected name; nl ();
      match tuples with
      | [] -> ()
      | _ ->
          ignore (List.fold_left (fun i (tuple1, tuple2) ->
            printf "  @[<2>data table_spin_%s(:,%4d) / %s /" name i
              (String.concat ", " (List.map (Printf.sprintf "%2d") (tuple1 @ tuple2)));
            nl (); succ i) 1 tuples)

    let print_spin_tables amplitudes =
      (* [print_spin_table_old "s" "states_old" (CF.helicities amplitudes);] *)
      print_spin_table "states" (CF.helicities amplitudes);
      nl ()

    (*i unused value
    let print_flavor_table_old n abbrev name = function
      | [] ->
          printf "  @[<2>integer, dimension(n_prt,0) ::";
          printf "@ table_flavor_%s" name; nl ()
      | _ :: tuples' as tuples ->
          ignore (List.fold_left (fun i tuple ->
            printf
              "  @[<2>integer, dimension(n_prt), parameter, private ::";
            printf "@ %s%04d = (/ %s /) ! %s" abbrev i
              (String.concat ", "
                 (List.map (fun f -> Printf.sprintf "%3d" (M.pdg f)) tuple))
              (String.concat " " (List.map M.flavor_to_string tuple));
            nl (); succ i) 1 tuples);
          printf
            "  @[<2>integer, dimension(n_prt,n_flv), parameter ::";
          printf "@ table_flavor_%s =@ reshape ( (/" name;
          printf "@ %s%04d" abbrev 1;
          ignore (List.fold_left (fun i tuple ->
            printf ",@ %s%04d" abbrev i; succ i) 2 tuples');
          printf "@ /), (/ n_prt, n_flv /) )"; nl ()
    i*)

    let print_flavor_table name tuples =
      printf "  @[<2>integer, dimension(n_prt,n_flv), save%s :: table_flavor_%s"
        protected name; nl ();
      match tuples with
      | [] -> ()
      | _ ->
          ignore (List.fold_left (fun i tuple ->
            printf "  @[<2>data table_flavor_%s(:,%4d) / %s / ! %s" name i
              (String.concat ", "
                 (List.map (fun f -> Printf.sprintf "%3d" (M.pdg f)) tuple))
              (String.concat " " (List.map M.flavor_to_string tuple));
            nl (); succ i) 1 tuples)

    let print_flavor_tables amplitudes =
      (* [let n = num_particles amplitudes in] *)
      (* [print_flavor_table_old n "f" "states_old"
        (List.map (fun (fin, fout) -> fin @ fout) (CF.flavors amplitudes));] *)
      print_flavor_table "states"
        (List.map (fun (fin, fout) -> fin @ fout) (CF.flavors amplitudes));
      nl ()

    let num_flavors amplitudes =
      List.length (CF.flavors amplitudes)

    (*i unused value
    let print_color_flows_table_old abbrev = function
      | [] ->
          printf "  @[<2>integer, dimension(n_cindex, n_prt, n_cflow) ::";
          printf "@ table_color_flows"; nl ()
      | _ :: tuples' as tuples ->
          ignore (List.fold_left (fun i tuple ->
            printf
              "  @[<2>integer, dimension(n_cindex, n_prt), parameter, private ::";
            printf "@ %s%04d = reshape ( (/ " abbrev i;
            begin match CFlow.to_lists tuple with
            | [] -> ()
            | cf1 :: cfn ->
                printf "@ %s" (String.concat "," (List.map string_of_int cf1));
                List.iter (function cf ->
                  printf ",@  %s" (String.concat "," (List.map string_of_int cf))) cfn
            end;
            printf "@ /),@ (/ n_cindex, n_prt /) )";
            nl (); succ i) 1 tuples);
          printf
            "  @[<2>integer, dimension(n_cindex, n_prt, n_cflow), parameter ::";
          printf "@ table_color_flows_old =@ reshape ( (/";
          printf "@ %s%04d" abbrev 1;
          ignore (List.fold_left (fun i tuple ->
            printf ",@ %s%04d" abbrev i; succ i) 2 tuples');
          printf "@ /),@ (/ n_cindex, n_prt, n_cflow /) )"; nl ()
    i*)

    (*i unused value
    let print_ghost_flags_table_old abbrev = function
      | [] ->
          printf "  @[<2>logical, dimension(n_prt, n_cflow) ::";
          printf "@ table_ghost_flags"; nl ()
      | _ :: tuples' as tuples ->
          ignore (List.fold_left (fun i tuple ->
            printf
              "  @[<2>logical, dimension(n_prt), parameter, private ::";
            printf "@ %s%04d = (/ " abbrev i;
            begin match CFlow.ghost_flags tuple with
            | [] -> ()
            | gf1 :: gfn ->
                printf "@ %s" (if gf1 then "T" else "F");
                List.iter (function gf -> printf ",@  %s" (if gf then "T" else "F")) gfn
            end;
            printf "@ /)";
            nl (); succ i) 1 tuples);
          printf
            "  @[<2>logical, dimension(n_prt, n_cflow), parameter ::";
          printf "@ table_ghost_flags_old =@ reshape ( (/";
          printf "@ %s%04d" abbrev 1;
          ignore (List.fold_left (fun i tuple ->
            printf ",@ %s%04d" abbrev i; succ i) 2 tuples');
          printf "@ /),@ (/ n_prt, n_cflow /) )"; nl ()
    i*)

    let print_color_flows_table tuples =
      if !amp_triv then begin
        printf
          "  @[<2>integer, dimension(n_cindex,n_prt,n_cflow), save%s :: table_color_flows = 0"
          protected; nl ();
	end
      else begin
        printf
          "  @[<2>integer, dimension(n_cindex,n_prt,n_cflow), save%s :: table_color_flows"
          protected; nl ();
      end;
      if not !amp_triv then begin
        match tuples with
        | [] -> ()
        | _ :: _ as tuples ->
            ignore (List.fold_left (fun i tuple ->
              begin match CFlow.to_lists tuple with
              | [] -> ()
              | cf1 :: cfn ->
                  printf "  @[<2>data table_color_flows(:,:,%4d) /" i;
                  printf "@ %s" (String.concat "," (List.map string_of_int cf1));
                  List.iter (function cf ->
                    printf ",@  %s" (String.concat "," (List.map string_of_int cf))) cfn;
                  printf "@ /"; nl ()
              end;
              succ i) 1 tuples)
      end

    let print_ghost_flags_table tuples =
      if !amp_triv then begin
        printf
          "  @[<2>logical, dimension(n_prt,n_cflow), save%s :: table_ghost_flags = F"
          protected; nl ();
	end
      else begin
        printf
          "  @[<2>logical, dimension(n_prt,n_cflow), save%s :: table_ghost_flags"
          protected; nl ();
        match tuples with
        | [] -> ()
        | _ ->
            ignore (List.fold_left (fun i tuple ->
              begin match CFlow.ghost_flags tuple with
              | [] -> ()
              | gf1 :: gfn ->
                  printf "  @[<2>data table_ghost_flags(:,%4d) /" i;
                  printf "@ %s" (if gf1 then "T" else "F");
                  List.iter (function gf -> printf ",@  %s" (if gf then "T" else "F")) gfn;
                  printf " /";
                  nl ()
              end;
              succ i) 1 tuples)
      end

    let format_power_of x
        { Color.Flow.num = num; Color.Flow.den = den; Color.Flow.power = pwr } =
      match num, den, pwr with
      | _, 0, _ -> invalid_arg "format_power_of: zero denominator"
      | 0, _, _ -> "+zero"
      | 1, 1, 0 | -1, -1, 0 -> "+one"
      | -1, 1, 0 | 1, -1, 0 -> "-one"
      | 1, 1, 1 | -1, -1, 1 -> "+" ^ x
      | -1, 1, 1 | 1, -1, 1 -> "-" ^ x
      | 1, 1, -1 | -1, -1, -1 -> "+1/" ^ x
      | -1, 1, -1 | 1, -1, -1 -> "-1/" ^ x
      | 1, 1, p | -1, -1, p ->
          "+" ^ (if p > 0 then "" else "1/") ^ x ^ "**" ^ string_of_int (abs p)
      | -1, 1, p | 1, -1, p ->
          "-" ^ (if p > 0 then "" else "1/") ^ x ^ "**" ^ string_of_int (abs p)
      | n, 1, 0 ->
          (if n < 0 then "-" else "+") ^ string_of_int (abs n) ^ ".0_" ^ !kind
      | n, d, 0 ->
          (if n * d < 0 then "-" else "+") ^
          string_of_int (abs n) ^ ".0_" ^ !kind ^ "/" ^
          string_of_int (abs d)
      | n, 1, 1 ->
          (if n < 0 then "-" else "+") ^ string_of_int (abs n) ^ "*" ^ x
      | n, 1, -1 ->
          (if n < 0 then "-" else "+") ^ string_of_int (abs n) ^ "/" ^ x
      | n, d, 1 ->
          (if n * d < 0 then "-" else "+") ^
          string_of_int (abs n) ^ ".0_" ^ !kind ^ "/" ^
          string_of_int (abs d) ^ "*" ^ x
      | n, d, -1 ->
          (if n * d < 0 then "-" else "+") ^
          string_of_int (abs n) ^ ".0_" ^ !kind ^ "/" ^
          string_of_int (abs d) ^ "/" ^ x
      | n, 1, p ->
          (if n < 0 then "-" else "+") ^ string_of_int (abs n) ^
          (if p > 0 then "*" else "/") ^ x ^ "**" ^ string_of_int (abs p)
      | n, d, p ->
          (if n * d < 0 then "-" else "+") ^
          string_of_int (abs n) ^ ".0_" ^ !kind ^ "/" ^
          string_of_int (abs d) ^
          (if p > 0 then "*" else "/") ^ x ^ "**" ^ string_of_int (abs p)

    let format_powers_of x = function
      | [] -> "zero"
      | powers -> String.concat "" (List.map (format_power_of x) powers)

    (*i unused value
    let print_color_factor_table_old table =
      let n_cflow = Array.length table in
      let n_cfactors = ref 0 in
      for c1 = 0 to pred n_cflow do
        for c2 = 0 to pred n_cflow do
          match table.(c1).(c2) with
          | [] -> ()
          | _ -> incr n_cfactors
        done
      done;
      print_integer_parameter "n_cfactors"  !n_cfactors;
      if n_cflow <= 0 then begin
        printf "  @[<2>type(%s), dimension(n_cfactors) ::"
          omega_color_factor_abbrev;
        printf "@ table_color_factors"; nl ()
      end else begin
        printf
          "  @[<2>type(%s), dimension(n_cfactors), parameter ::"
          omega_color_factor_abbrev;
        printf "@ table_color_factors = (/@ ";
        let comma = ref "" in
        for c1 = 0 to pred n_cflow do
          for c2 = 0 to pred n_cflow do
            match table.(c1).(c2) with
            | [] -> ()
            | cf ->
                printf "%s@ %s(%d,%d,%s)" !comma omega_color_factor_abbrev
                  (succ c1) (succ c2) (format_powers_of nc_parameter cf);
                comma := ","
          done
        done;
        printf "@ /)"; nl ()
      end
    i*)

(* \begin{dubious}
     We can optimize the following slightly by reusing common color factor [parameter]s.
   \end{dubious} *)

    let print_color_factor_table table =
      let n_cflow = Array.length table in
      let n_cfactors = ref 0 in
      for c1 = 0 to pred n_cflow do
        for c2 = 0 to pred n_cflow do
          match table.(c1).(c2) with
          | [] -> ()
          | _ -> incr n_cfactors
        done
      done;
      print_integer_parameter "n_cfactors"  !n_cfactors;
      printf "  @[<2>type(%s), dimension(n_cfactors), save%s ::"
        omega_color_factor_abbrev protected;
      printf "@ table_color_factors"; nl ();
      if not !amp_triv then begin
        let i = ref 1 in
        if n_cflow > 0 then begin
          for c1 = 0 to pred n_cflow do
            for c2 = 0 to pred n_cflow do
              match table.(c1).(c2) with
              | [] -> ()
              | cf ->
                  printf "  @[<2>real(kind=%s), parameter, private :: color_factor_%06d = %s"
                    !kind !i (format_powers_of nc_parameter cf);
                  nl ();
                  printf "  @[<2>data table_color_factors(%6d) / %s(%d,%d,color_factor_%06d) /"
                    !i omega_color_factor_abbrev (succ c1) (succ c2) !i;
                  incr i;
                  nl ();
            done
          done
        end;
      end

    let print_color_tables amplitudes =
      let cflows =  CF.color_flows amplitudes
      and cfactors = CF.color_factors amplitudes in
      (* [print_color_flows_table_old "c" cflows; nl ();] *)
      print_color_flows_table cflows; nl ();
      (* [print_ghost_flags_table_old "g" cflows; nl ();] *)
      print_ghost_flags_table cflows; nl ();
      (* [print_color_factor_table_old cfactors; nl ();] *)
      print_color_factor_table cfactors; nl ()

    let option_to_logical = function
      | Some _ -> "T"
      | None -> "F"

    (*i unused value
    let print_flavor_color_table_old abbrev n_flv n_cflow table =
      if n_flv <= 0 || n_cflow <= 0 then begin
        printf "  @[<2>logical, dimension(n_flv, n_cflow) ::";
        printf "@ flv_col_is_allowed"; nl ()
      end else begin
        for c = 0 to pred n_cflow do
          printf
            "  @[<2>logical, dimension(n_flv), parameter, private ::";
          printf "@ %s%04d = (/@ %s" abbrev (succ c) (option_to_logical table.(0).(c));
          for f = 1 to pred n_flv do
            printf ",@ %s" (option_to_logical table.(f).(c))
          done;
          printf "@ /)"; nl ()
        done;
        printf
          "  @[<2>logical, dimension(n_flv, n_cflow), parameter ::";
        printf "@ flv_col_is_allowed_old =@ reshape ( (/@ %s%04d" abbrev 1;
        for c = 1 to pred n_cflow do
          printf ",@ %s%04d" abbrev (succ c)
        done;
        printf "@ /),@ (/ n_flv, n_cflow /) )"; nl ()
      end
    i*)

    let print_flavor_color_table n_flv n_cflow table =
      if !amp_triv then begin
        printf
          "  @[<2>logical, dimension(n_flv, n_cflow), save%s :: @ flv_col_is_allowed = T"
        protected; nl ();
	end
      else begin
        printf
          "  @[<2>logical, dimension(n_flv, n_cflow), save%s :: @ flv_col_is_allowed"
        protected; nl ();
        if n_flv > 0 then begin
          for c = 0 to pred n_cflow do
            printf
              "  @[<2>data flv_col_is_allowed(:,%4d) /" (succ c);
            printf "@ %s" (option_to_logical table.(0).(c));
            for f = 1 to pred n_flv do
              printf ",@ %s" (option_to_logical table.(f).(c))
            done;
            printf "@ /"; nl ()
          done;
	end;
      end

    let print_amplitude_table a =
      (* [print_flavor_color_table_old "a"
        (num_flavors a) (List.length (CF.color_flows a)) (CF.process_table a);
      nl ();] *)
      print_flavor_color_table
        (num_flavors a) (List.length (CF.color_flows a)) (CF.process_table a);
      nl ();
      printf
        "  @[<2>complex(kind=%s), dimension(n_flv, n_cflow, n_hel), save :: amp" !kind;
      nl ();
      nl ()

    let print_helicity_selection_table () =
      printf "  @[<2>logical, dimension(n_hel), save :: ";
      printf "hel_is_allowed = T"; nl ();
      printf "  @[<2>real(kind=%s), dimension(n_hel), save :: " !kind;
      printf "hel_max_abs = 0"; nl ();
      printf "  @[<2>real(kind=%s), save :: " !kind;
      printf "hel_sum_abs = 0, ";
      printf "hel_threshold = 1E10_%s" !kind; nl ();
      printf "  @[<2>integer, save :: ";
      printf "hel_count = 0, ";
      printf "hel_cutoff = 100"; nl ();
      printf "  @[<2>integer :: ";
      printf "i"; nl ();
      printf "  @[<2>integer, save, dimension(n_hel) :: ";
      printf "hel_map = (/(i, i = 1, n_hel)/)"; nl ();
      printf "  @[<2>integer, save :: hel_finite = n_hel"; nl ();
      nl ()

(* \thocwmodulesubsection{Optional MD5 sum function} *)

    let print_md5sum_functions = function
      | Some s ->
          printf "  @[<5>"; if !fortran95 then printf "pure ";
          printf "function md5sum ()"; nl ();
          printf "    character(len=32) :: md5sum"; nl ();
          printf "    ! DON'T EVEN THINK of modifying the following line!"; nl ();
          printf "    md5sum = \"%s\"" s; nl ();
          printf "  end function md5sum"; nl ();
          nl ()
      | None -> ()

(* \thocwmodulesubsection{Maintenance \&\ Inquiry Functions} *)

    let print_maintenance_functions () =
      if !whizard then begin
        printf "  subroutine init (par, scheme)"; nl ();
        printf "    real(kind=%s), dimension(*), intent(in) :: par" !kind; nl ();
        printf "    integer, intent(in) :: scheme"; nl ();
        printf "    call import_from_whizard (par, scheme)"; nl ();
        printf "  end subroutine init"; nl ();
        nl ();
        printf "  subroutine final ()"; nl ();
        printf "  end subroutine final"; nl ();
        nl ();
        printf "  subroutine update_alpha_s (alpha_s)"; nl ();
        printf "    real(kind=%s), intent(in) :: alpha_s" !kind; nl ();
        printf "    call model_update_alpha_s (alpha_s)"; nl ();
        printf "  end subroutine update_alpha_s"; nl ();
        nl ()
      end

    let print_inquiry_function_openmp () = begin
      printf "  pure function openmp_supported () result (status)"; nl ();
      printf "    logical :: status"; nl ();
      printf "    status = %s" (if !openmp then ".true." else ".false."); nl ();
      printf "  end function openmp_supported"; nl ();
      nl ()
    end

    (*i unused value
    let print_inquiry_function_declarations name =
      printf "  @[<2>public :: number_%s,@ %s" name name;
      nl ()
    i*)

    (*i unused value
    let print_numeric_inquiry_functions () =
      printf "  @[<5>"; if !fortran95 then printf "pure ";
      printf "function number_particles_in () result (n)"; nl ();
      printf "    integer :: n"; nl ();
      printf "    n = n_in"; nl ();
      printf "  end function number_particles_in"; nl ();
      nl ();
      printf "  @[<5>"; if !fortran95 then printf "pure ";
      printf "function number_particles_out () result (n)"; nl ();
      printf "    integer :: n"; nl ();
      printf "    n = n_out"; nl ();
      printf "  end function number_particles_out"; nl ();
      nl ()
    i*)

    let print_external_mass_case flv (fin, fout) =
      printf "    case (%3d)" (succ flv); nl ();
      List.iteri
        (fun i f ->
          printf "      m(%2d) = %s" (succ i) (M.mass_symbol f); nl ())
        (fin @ fout)

    let print_external_masses amplitudes =
      printf "  @[<5>"; if !fortran95 then printf "pure ";
      printf "subroutine external_masses (m, flv)"; nl ();
      printf "    real(kind=%s), dimension(:), intent(out) :: m" !kind; nl ();
      printf "    integer, intent(in) :: flv"; nl ();
      printf "    select case (flv)"; nl ();
      List.iteri print_external_mass_case (CF.flavors amplitudes);
      printf "    end select"; nl ();
      printf "  end subroutine external_masses"; nl ();
      nl ()

    let print_numeric_inquiry_functions (f, v) =
      printf "  @[<5>"; if !fortran95 then printf "pure ";
      printf "function %s () result (n)" f; nl ();
      printf "    integer :: n"; nl ();
      printf "    n = %s" v; nl ();
      printf "  end function %s" f; nl ();
      nl ()

    let print_inquiry_functions name =
      printf "  @[<5>"; if !fortran95 then printf "pure ";
      printf "function number_%s () result (n)" name; nl ();
      printf "    integer :: n"; nl ();
      printf "    n = size (table_%s, dim=2)" name; nl ();
      printf "  end function number_%s" name; nl ();
      nl ();
      printf "  @[<5>"; if !fortran95 then printf "pure ";
      printf "subroutine %s (a)" name; nl ();
      printf "    integer, dimension(:,:), intent(out) :: a"; nl ();
      printf "    a = table_%s" name; nl ();
      printf "  end subroutine %s" name; nl ();
      nl ()

    let print_color_flows () =
      printf "  @[<5>"; if !fortran95 then printf "pure ";
      printf "function number_color_indices () result (n)"; nl ();
      printf "    integer :: n"; nl ();
      if !amp_triv then begin
        printf "    n = n_cindex"; nl ();
	end
      else begin
        printf "    n = size (table_color_flows, dim=1)"; nl ();
      end;
      printf "  end function number_color_indices"; nl ();
      nl ();
      printf "  @[<5>"; if !fortran95 then printf "pure ";
      printf "function number_color_flows () result (n)"; nl ();
      printf "    integer :: n"; nl ();
      if !amp_triv then begin
        printf "    n = n_cflow"; nl ();
	end
      else begin
        printf "    n = size (table_color_flows, dim=3)"; nl ();
      end;
      printf "  end function number_color_flows"; nl ();
      nl ();
      printf "  @[<5>"; if !fortran95 then printf "pure ";
      printf "subroutine color_flows (a, g)"; nl ();
      printf "    integer, dimension(:,:,:), intent(out) :: a"; nl ();
      printf "    logical, dimension(:,:), intent(out) :: g"; nl ();
      printf "    a = table_color_flows"; nl ();
      printf "    g = table_ghost_flags"; nl ();
      printf "  end subroutine color_flows"; nl ();
      nl ()

    let print_color_factors () =
      printf "  @[<5>"; if !fortran95 then printf "pure ";
      printf "function number_color_factors () result (n)"; nl ();
      printf "    integer :: n"; nl ();
      printf "    n = size (table_color_factors)"; nl ();
      printf "  end function number_color_factors"; nl ();
      nl ();
      printf "  @[<5>"; if !fortran95 then printf "pure ";
      printf "subroutine color_factors (cf)"; nl ();
      printf "    type(%s), dimension(:), intent(out) :: cf"
        omega_color_factor_abbrev; nl ();
      printf "    cf = table_color_factors"; nl ();
      printf "  end subroutine color_factors"; nl ();
      nl ();
      printf "  @[<5>"; if !fortran95 && pure_unless_openmp then printf "pure ";
      printf "function color_sum (flv, hel) result (amp2)"; nl ();
      printf "    integer, intent(in) :: flv, hel"; nl ();
      printf "    real(kind=%s) :: amp2" !kind; nl ();
      printf "    amp2 = real (omega_color_sum (flv, hel, amp, table_color_factors))"; nl ();
      printf "  end function color_sum"; nl ();
      nl ()

    let print_dispatch_functions () =
      printf "  @[<5>";
      printf "subroutine new_event (p)"; nl ();
      printf "    real(kind=%s), dimension(0:3,*), intent(in) :: p" !kind; nl ();
      printf "    logical :: mask_dirty"; nl ();
      printf "    integer :: hel"; nl ();
      printf "    call calculate_amplitudes (amp, p, hel_is_allowed)"; nl ();
      printf "    if ((hel_threshold .gt. 0) .and. (hel_count .le. hel_cutoff)) then"; nl ();
      printf "      call @[<3>omega_update_helicity_selection@ (hel_count,@ amp,@ ";
      printf "hel_max_abs,@ hel_sum_abs,@ hel_is_allowed,@ hel_threshold,@ hel_cutoff,@ mask_dirty)"; nl ();
      printf "      if (mask_dirty) then"; nl ();
      printf "        hel_finite = 0"; nl ();
      printf "        do hel = 1, n_hel"; nl ();
      printf "          if (hel_is_allowed(hel)) then"; nl ();
      printf "            hel_finite = hel_finite + 1"; nl ();
      printf "            hel_map(hel_finite) = hel"; nl ();
      printf "          end if"; nl ();
      printf "        end do"; nl ();
      printf "      end if"; nl ();
      printf "    end if"; nl ();
      printf "  end subroutine new_event"; nl ();
      nl ();
      printf "  @[<5>";
      printf "subroutine reset_helicity_selection (threshold, cutoff)"; nl ();
      printf "    real(kind=%s), intent(in) :: threshold" !kind; nl ();
      printf "    integer, intent(in) :: cutoff"; nl ();
      printf "    integer :: i"; nl ();
      printf "    hel_is_allowed = T"; nl ();
      printf "    hel_max_abs = 0"; nl ();
      printf "    hel_sum_abs = 0"; nl ();
      printf "    hel_count = 0"; nl ();
      printf "    hel_threshold = threshold"; nl ();
      printf "    hel_cutoff = cutoff"; nl ();
      printf "    hel_map = (/(i, i = 1, n_hel)/)"; nl ();
      printf "    hel_finite = n_hel"; nl ();
      printf "  end subroutine reset_helicity_selection"; nl ();
      nl ();
      printf "  @[<5>"; if !fortran95 then printf "pure ";
      printf "function is_allowed (flv, hel, col) result (yorn)"; nl ();
      printf "    logical :: yorn"; nl ();
      printf "    integer, intent(in) :: flv, hel, col"; nl ();
      if !amp_triv then begin
         printf "    ! print *, 'inside is_allowed'"; nl ();
      end;
      if not !amp_triv then begin
         printf "    yorn = hel_is_allowed(hel) .and. ";
         printf "flv_col_is_allowed(flv,col)"; nl ();
         end
      else begin
         printf "    yorn = .false."; nl ();
      end;
      printf "  end function is_allowed"; nl ();
      nl ();
      printf "  @[<5>"; if !fortran95 then printf "pure ";
      printf "function get_amplitude (flv, hel, col) result (amp_result)"; nl ();
      printf "    complex(kind=%s) :: amp_result" !kind; nl ();
      printf "    integer, intent(in) :: flv, hel, col"; nl ();
      printf "    amp_result = amp(flv, col, hel)"; nl ();
      printf "  end function get_amplitude"; nl ();
      nl ()

(* \thocwmodulesubsection{Main Function} *)

    let format_power_of_nc
        { Color.Flow.num = num; Color.Flow.den = den; Color.Flow.power = pwr } =
      match num, den, pwr with
      | _, 0, _ -> invalid_arg "format_power_of_nc: zero denominator"
      | 0, _, _ -> ""
      | 1, 1, 0 | -1, -1, 0 -> "+ 1"
      | -1, 1, 0 | 1, -1, 0 -> "- 1"
      | 1, 1, 1 | -1, -1, 1 -> "+ N"
      | -1, 1, 1 | 1, -1, 1 -> "- N"
      | 1, 1, -1 | -1, -1, -1 -> "+ 1/N"
      | -1, 1, -1 | 1, -1, -1 -> "- 1/N"
      | 1, 1, p | -1, -1, p ->
          "+ " ^ (if p > 0 then "" else "1/") ^ "N^" ^ string_of_int (abs p)
      | -1, 1, p | 1, -1, p ->
          "- " ^ (if p > 0 then "" else "1/") ^ "N^" ^ string_of_int (abs p)
      | n, 1, 0 ->
          (if n < 0 then "- " else "+ ") ^ string_of_int (abs n)
      | n, d, 0 ->
          (if n * d < 0 then "- " else "+ ") ^
          string_of_int (abs n) ^ "/" ^ string_of_int (abs d)
      | n, 1, 1 ->
          (if n < 0 then "- " else "+ ") ^ string_of_int (abs n) ^ "N"
      | n, 1, -1 ->
          (if n < 0 then "- " else "+ ") ^ string_of_int (abs n) ^ "/N"
      | n, d, 1 ->
          (if n * d < 0 then "- " else "+ ") ^
          string_of_int (abs n) ^ "/" ^ string_of_int (abs d) ^ "N"
      | n, d, -1 ->
          (if n * d < 0 then "- " else "+ ") ^
          string_of_int (abs n) ^ "/" ^ string_of_int (abs d) ^ "/N"
      | n, 1, p ->
          (if n < 0 then "- " else "+ ") ^ string_of_int (abs n) ^
          (if p > 0 then "*" else "/") ^ "N^" ^ string_of_int (abs p)
      | n, d, p ->
          (if n * d < 0 then "- " else "+ ") ^ string_of_int (abs n) ^ "/" ^
          string_of_int (abs d) ^ (if p > 0 then "*" else "/") ^ "N^" ^ string_of_int (abs p)

    let format_powers_of_nc = function
      | [] -> "0"
      | powers -> String.concat " " (List.map format_power_of_nc powers)

    let print_description cmdline amplitudes () =
      printf
        "! File generated automatically by O'Mega %s %s %s"
        Config.version Config.status Config.date; nl ();
      List.iter (fun s -> printf "! %s" s; nl ()) (M.caveats ());
      printf "!"; nl ();
      printf "!   %s" cmdline; nl ();
      printf "!"; nl ();
      printf "! with all scattering amplitudes for the process(es)"; nl ();
      printf "!"; nl ();
      printf "!   flavor combinations:"; nl ();
      printf "!"; nl ();
      ThoList.iteri
        (fun i process ->
          printf "!     %3d: %s" i (process_sans_color_to_string process); nl ())
        1 (CF.flavors amplitudes);
      printf "!"; nl ();
      printf "!   color flows:"; nl ();
      if not !amp_triv then begin
	printf "!"; nl ();
	ThoList.iteri
          (fun i cflow ->
            printf "!     %3d: %s" i (cflow_to_string cflow); nl ())
          1 (CF.color_flows amplitudes);
	printf "!"; nl ();
	printf "!     NB: i.g. not all color flows contribute to all flavor"; nl ();
	printf "!     combinations.  Consult the array FLV_COL_IS_ALLOWED"; nl ();
	printf "!     below for the allowed combinations."; nl ();
      end;
      printf "!"; nl ();
      printf "!   Color Factors:"; nl ();
      printf "!"; nl ();
      if not !amp_triv then begin
	let cfactors = CF.color_factors amplitudes in
	for c1 = 0 to pred (Array.length cfactors) do
          for c2 = 0 to c1 do
            match cfactors.(c1).(c2) with
            | [] -> ()
            | cfactor ->
               printf "!     (%3d,%3d): %s"
                 (succ c1) (succ c2) (format_powers_of_nc cfactor); nl ()
          done
	done;
      end;
      if not !amp_triv then begin
         printf "!"; nl ();
         printf "!   vanishing or redundant flavor combinations:"; nl ();
         printf "!"; nl ();
         List.iter (fun process ->
           printf "!          %s" (process_sans_color_to_string process); nl ())
           (CF.vanishing_flavors amplitudes);
         printf "!"; nl ();
      end;
      begin
        match CF.constraints amplitudes with
        | None -> ()
        | Some s ->
            printf
              "!   diagram selection (MIGHT BREAK GAUGE INVARIANCE!!!):"; nl ();
            printf "!"; nl ();
            printf "!     %s" s; nl ();
            printf "!"; nl ()
      end;
      printf "!"; nl ()

(* \thocwmodulesubsection{Printing Modules} *)

    type accessibility =
      | Public
      | Private
      | Protected (* Fortran 2003 *)

    let accessibility_to_string = function
      | Public -> "public"
      | Private -> "private"
      | Protected -> "protected"

    type used_symbol =
      | As_Is of string
      | Aliased of string * string

    let print_used_symbol = function
      | As_Is name -> printf "%s" name
      | Aliased (orig, alias) -> printf "%s => %s" alias orig

    type used_module =
      | Full of string
      | Full_Aliased of string * (string * string) list
      | Subset of string * used_symbol list

    let print_used_module = function
      | Full name
      | Full_Aliased (name, [])
      | Subset (name, []) ->
          printf "  use %s" name;
          nl ()
      | Full_Aliased (name, aliases) ->
          printf "  @[<5>use %s" name;
          List.iter
            (fun (orig, alias) -> printf ", %s => %s" alias orig)
            aliases;
          nl ()
      | Subset (name, used_symbol :: used_symbols) ->
          printf "  @[<5>use %s, only: " name;
          print_used_symbol used_symbol;
          List.iter (fun s -> printf ", "; print_used_symbol s) used_symbols;
          nl ()

    type fortran_module =
        { module_name : string;
          default_accessibility : accessibility;
          used_modules : used_module list;
          public_symbols : string list;
          print_declarations : (unit -> unit) list;
          print_implementations : (unit -> unit) list }

    let print_public = function
      | name1 :: names ->
          printf "  @[<2>public :: %s" name1;
          List.iter (fun n -> printf ",@ %s" n) names; nl ()
      | [] -> ()

    (*i unused value
    let print_public_interface generic procedures =
      printf "  public :: %s" generic; nl ();
      begin match procedures with
      | name1 :: names ->
          printf "  interface %s" generic; nl ();
          printf "     @[<2>module procedure %s" name1;
          List.iter (fun n -> printf ",@ %s" n) names; nl ();
          printf "  end interface"; nl ();
          print_public procedures
      | [] -> ()
      end
    i*)

    let print_module m =
      printf "module %s" m.module_name; nl ();
      List.iter print_used_module m.used_modules;
      printf "  implicit none"; nl ();
      printf "  %s" (accessibility_to_string m.default_accessibility); nl ();
      print_public m.public_symbols; nl ();
      begin match m.print_declarations with
      | [] -> ()
      | print_declarations ->
          List.iter (fun f -> f ()) print_declarations; nl ()
      end;
      begin match m.print_implementations with
      | [] -> ()
      | print_implementations ->
          printf "contains"; nl (); nl ();
          List.iter (fun f -> f ()) print_implementations; nl ();
      end;
      printf "end module %s" m.module_name; nl ()

    let print_modules modules =
      List.iter print_module modules;
      print_flush ()

    let module_to_file line_length oc prelude m =
      output_string oc (m.module_name ^ "\n");
      let filename = m.module_name ^ ".f90" in
      let channel = open_out filename in
      Format_Fortran.set_formatter_out_channel ~width:line_length channel;
      prelude ();
      print_modules [m];
      close_out channel

    let modules_to_file line_length oc prelude = function
      | [] -> ()
      | m :: mlist ->
          module_to_file line_length oc prelude m;
          List.iter (module_to_file line_length oc (fun () -> ())) mlist

(* \thocwmodulesubsection{Chopping Up Amplitudes} *)

    let num_fusions_brakets size amplitudes =
      let num_fusions =
        max 1 size in
      let count_brakets =
        List.fold_left
          (fun sum process -> sum + List.length (F.brakets process))
          0 (CF.processes amplitudes)
      and count_processes =
        List.length (CF.processes amplitudes) in
      if count_brakets > 0 then
        let num_brakets =
          max 1 ((num_fusions * count_processes) / count_brakets) in
        (num_fusions, num_brakets)
      else
        (num_fusions, 1)

    let chop_amplitudes size amplitudes =
      let num_fusions, num_brakets = num_fusions_brakets size amplitudes in
      (ThoList.enumerate 1 (ThoList.chopn num_fusions (CF.fusions amplitudes)),
       ThoList.enumerate 1 (ThoList.chopn num_brakets (CF.processes amplitudes)))

    let print_compute_fusions1 dictionary (n, fusions) =
      if not !amp_triv then begin
	if !openmp then begin
          printf "  subroutine compute_fusions_%04d (%s)" n openmp_tld; nl ();
          printf "  @[<5>type(%s), intent(inout) :: %s" openmp_tld_type openmp_tld; nl ();
	end else begin
          printf "  @[<5>subroutine compute_fusions_%04d ()" n; nl ();
	end;
	print_fusions dictionary fusions;
	printf "  end subroutine compute_fusions_%04d" n; nl ();
      end

    and print_compute_brakets1 dictionary (n, processes) =
      if not !amp_triv then begin
	if !openmp then begin
          printf "  subroutine compute_brakets_%04d (%s)" n openmp_tld; nl ();
          printf "  @[<5>type(%s), intent(inout) :: %s" openmp_tld_type openmp_tld; nl ();
	end else begin
          printf "  @[<5>subroutine compute_brakets_%04d ()" n; nl ();
	end;
	List.iter (print_brakets dictionary) processes;
	printf "  end subroutine compute_brakets_%04d" n; nl ();
      end

(* \thocwmodulesubsection{Common Stuff} *)

    let omega_public_symbols =
      ["number_particles_in"; "number_particles_out";
       "number_color_indices";
       "reset_helicity_selection"; "new_event";
       "is_allowed"; "get_amplitude"; "color_sum";
       "external_masses"; "openmp_supported"] @
      ThoList.flatmap
        (fun n -> ["number_" ^ n; n])
        ["spin_states"; "flavor_states"; "color_flows"; "color_factors"]

    let whizard_public_symbols md5sum =
      ["init"; "final"; "update_alpha_s"] @
      (match md5sum with Some _ -> ["md5sum"] | None -> [])

    let used_modules () =
      [Full "kinds";
       Full Fermions.use_module;
       Full_Aliased ("omega_color", ["omega_color_factor", omega_color_factor_abbrev])] @
      List.map
        (fun m -> Full m)
        (match !parameter_module with
         | "" -> !use_modules
         | pm -> pm :: !use_modules)

    let public_symbols () =
      if !whizard then
        omega_public_symbols @ (whizard_public_symbols !md5sum)
      else
        omega_public_symbols

    let print_constants amplitudes =

      printf "  ! DON'T EVEN THINK of removing the following!"; nl ();
      printf "  ! If the compiler complains about undeclared"; nl ();
      printf "  ! or undefined variables, you are compiling"; nl ();
      printf "  ! against an incompatible omega95 module!"; nl ();
      printf "  @[<2>integer, dimension(%d), parameter, private :: "
        (List.length require_library);
      printf "require =@ (/ @[";
      print_list require_library;
      printf " /)"; nl (); nl ();

      (* Using these parameters makes sense for documentation, but in
         practice, there is no need to ever change them. *)
      List.iter
        (function name, value -> print_integer_parameter name (value amplitudes))
        [ ("n_prt", num_particles);
          ("n_in", num_particles_in);
          ("n_out", num_particles_out);
          ("n_cflow", num_color_flows); (* Number of different color amplitudes. *)
          ("n_cindex", num_color_indices);  (* Maximum rank of color tensors. *)
          ("n_flv", num_flavors); (* Number of different flavor amplitudes. *)
          ("n_hel", num_helicities)  (* Number of different helicty amplitudes. *) ];
      nl ();

      (* Abbreviations.  *)
      printf "  ! NB: you MUST NOT change the value of %s here!!!" nc_parameter;
      nl ();
      printf "  !     It is defined here for convenience only and must be"; nl ();
      printf "  !     compatible with hardcoded values in the amplitude!"; nl ();
      print_real_parameter nc_parameter (CM.nc ()); (* $N_C$ *)
      List.iter
        (function name, value -> print_logical_parameter name value)
        [ ("F", false); ("T", true) ]; nl ();

      print_spin_tables amplitudes;
      print_flavor_tables amplitudes;
      print_color_tables amplitudes;
      print_amplitude_table amplitudes;
      print_helicity_selection_table ()

    let print_interface amplitudes =
      print_md5sum_functions !md5sum;
      print_maintenance_functions ();
      List.iter print_numeric_inquiry_functions
        [("number_particles_in", "n_in");
         ("number_particles_out", "n_out")];
      List.iter print_inquiry_functions
        ["spin_states"; "flavor_states"];
      print_external_masses amplitudes;
      print_inquiry_function_openmp ();
      print_color_flows ();
      print_color_factors ();
      print_dispatch_functions ();
      nl ();
      (* Is this really necessary? *)
      Format_Fortran.switch_line_continuation false;
      if !km_write || !km_pure then (Targets_Kmatrix.Fortran.print !km_pure);
      if !km_2_write || !km_2_pure then (Targets_Kmatrix_2.Fortran.print !km_2_pure);
      Format_Fortran.switch_line_continuation true;
      nl ()

    let print_calculate_amplitudes declarations computations amplitudes =
      printf "  @[<5>subroutine calculate_amplitudes (amp, k, mask)"; nl ();
      printf "    complex(kind=%s), dimension(:,:,:), intent(out) :: amp" !kind; nl ();
      printf "    real(kind=%s), dimension(0:3,*), intent(in) :: k" !kind; nl ();
      printf "    logical, dimension(:), intent(in) :: mask"; nl ();
      printf "    integer, dimension(n_prt) :: s"; nl ();
      printf "    integer :: h, hi"; nl ();
      declarations ();
      if not !amp_triv then begin
	begin match CF.processes amplitudes with
	| p :: _ -> print_external_momenta p
	|  _ -> ()
	end;
	ignore (List.fold_left print_momenta PSet.empty (CF.processes amplitudes));
      end;
      printf "    amp = 0"; nl ();
      if not !amp_triv then begin
	if num_helicities amplitudes > 0 then begin
          printf "    if (hel_finite == 0) return"; nl ();
          if !openmp then begin
            printf "!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(s, h, %s) SCHEDULE(STATIC)" openmp_tld; nl ();
          end;
          printf "    do hi = 1, hel_finite"; nl ();
          printf "      h = hel_map(hi)"; nl ();
          printf "      s = table_spin_states(:,h)"; nl ();
          ignore (List.fold_left print_externals WFSet.empty (CF.processes amplitudes));
          computations ();
          List.iter print_fudge_factor (CF.processes amplitudes);
        (* This sorting should slightly improve cache locality. *)
          let triple_snd = fun (_,  x, _) -> x
          in let triple_fst = fun (x, _, _) -> x
             in let rec builder1 flvi flowi flows = match flows with
             | (Some a) :: tl -> (flvi, flowi, flavors_symbol (flavors a)) :: (builder1 flvi (flowi + 1) tl)
             | None :: tl -> builder1 flvi (flowi + 1) tl
             | [] -> []
		in let rec builder2 flvi flvs = match flvs with
		| flv :: tl -> (builder1 flvi 1 flv) @ (builder2 (flvi + 1) tl)
		| [] -> []
		   in let unsorted = builder2 1 (List.map Array.to_list (Array.to_list (CF.process_table amplitudes)))
		      in let sorted = List.sort (fun a b ->
			if (triple_snd a != triple_snd b) then triple_snd a - triple_snd b else (triple_fst a - triple_fst b))
			   unsorted
			 in List.iter (fun (flvi, flowi, flv) ->
			   (printf "      amp(%d,%d,h) = %s" flvi flowi flv; nl ();)) sorted;

	(*i     printf "     else"; nl ();
          printf "      amp(:,h,:) = 0"; nl (); i*)
			 printf "    end do"; nl ();
			 if !openmp then begin
			   printf "!$OMP END PARALLEL DO"; nl ();
			 end;
	end;
      end;
      printf "  end subroutine calculate_amplitudes"; nl ()

    let print_compute_chops chopped_fusions chopped_brakets () =
      List.iter
        (fun (i, _) -> printf "      call compute_fusions_%04d (%s)" i
           (if !openmp then openmp_tld else ""); nl ())
        chopped_fusions;
      List.iter
        (fun (i, _) -> printf "      call compute_brakets_%04d (%s)" i
           (if !openmp then openmp_tld else ""); nl ())
        chopped_brakets

    (* \thocwmodulesubsection{UFO Fusions} *)

    module VSet =
      Set.Make (struct type t = F.constant Coupling.t let compare = compare end)

    let ufo_fusions_used amplitudes =
      let couplings =
        List.fold_left
          (fun acc p ->
            let fusions = ThoList.flatmap F.rhs (F.fusions p)
            and brakets = ThoList.flatmap F.ket (F.brakets p) in
            let couplings =
              VSet.of_list (List.map F.coupling (fusions @ brakets)) in
            VSet.union acc couplings)
          VSet.empty (CF.processes amplitudes) in
      VSet.fold
        (fun v acc ->
          match v with
          | Coupling.Vn (Coupling.UFO (_, v, _, _, _), _, _) ->
             Sets.String.add v acc
          | _ -> acc)
        couplings Sets.String.empty

(* \thocwmodulesubsection{Single Function} *)

    let amplitudes_to_channel_single_function cmdline oc amplitudes =

      let print_declarations () =
        print_constants amplitudes

      and print_implementations () =
        print_interface amplitudes;
        print_calculate_amplitudes
          (fun () -> print_variable_declarations amplitudes)
          (fun () ->
            print_fusions (CF.dictionary amplitudes) (CF.fusions amplitudes);
            List.iter
              (print_brakets (CF.dictionary amplitudes))
              (CF.processes amplitudes))
          amplitudes in

      let fortran_module =
        { module_name = !module_name;
          used_modules = used_modules ();
          default_accessibility = Private;
          public_symbols = public_symbols ();
          print_declarations = [print_declarations];
          print_implementations = [print_implementations] } in

      Format_Fortran.set_formatter_out_channel ~width:!line_length oc;
      print_description cmdline amplitudes ();
      print_modules [fortran_module]

(* \thocwmodulesubsection{Single Module} *)

    let amplitudes_to_channel_single_module cmdline oc size amplitudes =

      let print_declarations () =
        print_constants amplitudes;
        print_variable_declarations amplitudes

      and print_implementations () =
        print_interface amplitudes in

      let chopped_fusions, chopped_brakets =
        chop_amplitudes size amplitudes in

      let dictionary = CF.dictionary amplitudes in

      let print_compute_amplitudes () =
        print_calculate_amplitudes
          (fun () -> ())
          (print_compute_chops chopped_fusions chopped_brakets)
          amplitudes

      and print_compute_fusions () =
        List.iter (print_compute_fusions1 dictionary) chopped_fusions

      and print_compute_brakets () =
        List.iter (print_compute_brakets1 dictionary) chopped_brakets in

      let fortran_module =
        { module_name = !module_name;
          used_modules = used_modules ();
          default_accessibility = Private;
          public_symbols = public_symbols ();
          print_declarations = [print_declarations];
          print_implementations = [print_implementations;
                                   print_compute_amplitudes;
                                   print_compute_fusions;
                                   print_compute_brakets] } in

      Format_Fortran.set_formatter_out_channel ~width:!line_length oc;
      print_description cmdline amplitudes ();
      print_modules [fortran_module]

(* \thocwmodulesubsection{Multiple Modules} *)

    let modules_of_amplitudes _ _ size amplitudes =

      let name = !module_name in

      let print_declarations () =
        print_constants amplitudes
      and print_variables () =
        print_variable_declarations amplitudes in

      let constants_module =
        { module_name = name ^ "_constants";
          used_modules = used_modules ();
          default_accessibility = Public;
          public_symbols = [];
          print_declarations = [print_declarations];
          print_implementations = [] } in

      let variables_module =
        { module_name = name ^ "_variables";
          used_modules = used_modules ();
          default_accessibility = Public;
          public_symbols = [];
          print_declarations = [print_variables];
          print_implementations = [] } in

      let dictionary = CF.dictionary amplitudes in

      let print_compute_fusions (n, fusions) () =
	if not !amp_triv then begin
          if !openmp then begin
            printf "  subroutine compute_fusions_%04d (%s)" n openmp_tld; nl ();
            printf "  @[<5>type(%s), intent(inout) :: %s" openmp_tld_type openmp_tld; nl ();
          end else begin
            printf "  @[<5>subroutine compute_fusions_%04d ()" n; nl ();
          end;
          print_fusions dictionary fusions;
          printf "  end subroutine compute_fusions_%04d" n; nl ();
	end in

      let print_compute_brakets (n, processes) () =
	if not !amp_triv then begin
          if !openmp then begin
            printf "  subroutine compute_brakets_%04d (%s)" n openmp_tld; nl ();
            printf "  @[<5>type(%s), intent(inout) :: %s" openmp_tld_type openmp_tld; nl ();
          end else begin
            printf "  @[<5>subroutine compute_brakets_%04d ()" n; nl ();
          end;
          List.iter (print_brakets dictionary) processes;
          printf "  end subroutine compute_brakets_%04d" n; nl ();
	end in

      let fusions_module (n, _ as fusions) =
        let tag = Printf.sprintf "_fusions_%04d" n in
        { module_name = name ^ tag;
          used_modules = (used_modules () @
                          [Full constants_module.module_name;
                           Full variables_module.module_name]);
          default_accessibility = Private;
          public_symbols = ["compute" ^ tag];
          print_declarations = [];
          print_implementations = [print_compute_fusions fusions] } in

      let brakets_module (n, _ as processes) =
        let tag = Printf.sprintf "_brakets_%04d" n in
        { module_name = name ^ tag;
          used_modules = (used_modules () @
                          [Full constants_module.module_name;
                           Full variables_module.module_name]);
          default_accessibility = Private;
          public_symbols = ["compute" ^ tag];
          print_declarations = [];
          print_implementations = [print_compute_brakets processes] } in

      let chopped_fusions, chopped_brakets =
        chop_amplitudes size amplitudes in

      let fusions_modules =
        List.map fusions_module chopped_fusions in

      let brakets_modules =
        List.map brakets_module chopped_brakets in

      let print_implementations () =
        print_interface amplitudes;
        print_calculate_amplitudes
          (fun () -> ())
          (print_compute_chops chopped_fusions chopped_brakets)
          amplitudes in

      let public_module =
        { module_name = name;
           used_modules = (used_modules () @
                           [Full constants_module.module_name;
                            Full variables_module.module_name ] @
                           List.map
                             (fun m -> Full m.module_name)
                             (fusions_modules @ brakets_modules));
          default_accessibility = Private;
          public_symbols = public_symbols ();
          print_declarations = [];
          print_implementations = [print_implementations] }
      and private_modules =
        [constants_module; variables_module] @
          fusions_modules @ brakets_modules in
      (public_module, private_modules)

    let amplitudes_to_channel_single_file cmdline oc size amplitudes =
      let public_module, private_modules =
        modules_of_amplitudes cmdline oc size amplitudes in
      Format_Fortran.set_formatter_out_channel ~width:!line_length oc;
      print_description cmdline amplitudes ();
      print_modules (private_modules @ [public_module])

    let amplitudes_to_channel_multi_file cmdline oc size amplitudes =
      let public_module, private_modules =
        modules_of_amplitudes cmdline oc size amplitudes in
      modules_to_file !line_length oc
        (print_description cmdline amplitudes)
        (public_module :: private_modules)

(* \thocwmodulesubsection{Dispatch} *)

    let amplitudes_to_channel cmdline oc diagnostics amplitudes =
      parse_diagnostics diagnostics;
      let ufo_fusions =
        let ufo_fusions_set = ufo_fusions_used amplitudes in
        if Sets.String.is_empty ufo_fusions_set then
          None
        else
          Some ufo_fusions_set in
      begin match ufo_fusions with
      | Some only ->
         let name = !module_name ^ "_ufo"
         and fortran_module = Fermions.use_module in
         use_modules := name :: !use_modules;
         UFO.Targets.Fortran.lorentz_module
           ~only ~name ~fortran_module ~parameter_module:!parameter_module
           (Format_Fortran.formatter_of_out_channel oc) ()
      | None -> ()
      end;
      match !output_mode with
      | Single_Function ->
          amplitudes_to_channel_single_function cmdline oc amplitudes
      | Single_Module size ->
          amplitudes_to_channel_single_module cmdline oc size amplitudes
      | Single_File size ->
          amplitudes_to_channel_single_file cmdline oc size amplitudes
      | Multi_File size ->
          amplitudes_to_channel_multi_file cmdline oc size amplitudes

    let parameters_to_channel oc =
      parameters_to_fortran oc (CM.parameters ())

  end

module Fortran = Make_Fortran(Fortran_Fermions)

(* \thocwmodulesubsection{Majorana Fermions} *)

(* \begin{JR}
   For this function we need a different approach due to our aim of
   implementing the fermion vertices with the right line as ingoing (in a
   calculational sense) and the left line in a fusion as outgoing. In
   defining all external lines and the fermionic wavefunctions built out of
   them as ingoing we have to invert the left lines to make them outgoing.
   This happens by multiplying them with the inverse charge conjugation
   matrix in an appropriate representation and then transposing it. We must
   distinguish whether the direction of calculation and the physical direction
   of the fermion number flow are parallel or antiparallel. In the first case
   we can use the "normal" Feynman rules for Dirac particles, while in the
   second, according to the paper of Denner et al., we have to reverse the
   sign of the vector and antisymmetric bilinears of the Dirac spinors, cf.
   the [Coupling] module.

   Note the subtlety for the left- and righthanded couplings: Only the vector
   part of these couplings changes in the appropriate cases its sign,
   changing the chirality to the negative of the opposite.
   \end{JR} *)

module Fortran_Majorana_Fermions : Fermions =
  struct
    open Coupling
    open Format

    let psi_type = "bispinor"
    let psibar_type = "bispinor"
    let chi_type = "bispinor"
    let grav_type = "vectorspinor"

(* \begin{JR}
   Because of our rules for fermions we are going to give all incoming fermions
   a [u] spinor and all outgoing fermions a [v] spinor, no matter whether they
   are Dirac fermions, antifermions or Majorana fermions.
   \end{JR} *)

    let psi_incoming = "u"
    let brs_psi_incoming = "brs_u"
    let psibar_incoming = "u"
    let brs_psibar_incoming = "brs_u"
    let chi_incoming = "u"
    let brs_chi_incoming = "brs_u"
    let grav_incoming = "ueps"

    let psi_outgoing = "v"
    let brs_psi_outgoing = "brs_v"
    let psibar_outgoing = "v"
    let brs_psibar_outgoing = "brs_v"
    let chi_outgoing = "v"
    let brs_chi_outgoing = "brs_v"
    let grav_outgoing = "veps"

    let psi_propagator = "pr_psi"
    let psibar_propagator = "pr_psi"
    let chi_propagator = "pr_psi"
    let grav_propagator = "pr_grav"

    let psi_projector = "pj_psi"
    let psibar_projector = "pj_psi"
    let chi_projector = "pj_psi"
    let grav_projector = "pj_grav"

    let psi_gauss = "pg_psi"
    let psibar_gauss = "pg_psi"
    let chi_gauss = "pg_psi"
    let grav_gauss = "pg_grav"

    let format_coupling coeff c =
      match coeff with
      | 1 -> c
      | -1 -> "(-" ^ c ^")"
      | coeff -> string_of_int coeff ^ "*" ^ c

    let format_coupling_2 coeff c =
      match coeff with
      | 1 -> c
      | -1 -> "-" ^ c
      | coeff -> string_of_int coeff ^ "*" ^ c

(* \begin{dubious}
     JR's coupling constant HACK, necessitated by tho's bad design descition.
   \end{dubious} *)

    let fastener s i =
      try
        let offset = (String.index s '(') in
        if ((String.get s (String.length s - 1)) != ')') then
          failwith "fastener: wrong usage of parentheses"
        else
          let func_name = (String.sub s 0 offset) and
              tail =
                (String.sub s (succ offset) (String.length s - offset - 2)) in
          if (String.contains func_name ')') ||
             (String.contains tail '(') ||
             (String.contains tail ')') then
            failwith "fastener: wrong usage of parentheses"
          else
            func_name ^ "(" ^ string_of_int i ^ "," ^ tail ^ ")"
      with
      | Not_found ->
          if (String.contains s ')') then
            failwith "fastener: wrong usage of parentheses"
          else
            s ^ "(" ^ string_of_int i ^ ")"

    let print_fermion_current coeff f c wf1 wf2 fusion =
      let c = format_coupling coeff c in
      match fusion with
      | F13 | F31 -> printf "%s_ff(%s,%s,%s)" f c wf1 wf2
      | F23 | F21 -> printf "f_%sf(%s,%s,%s)" f c wf1 wf2
      | F32 | F12 -> printf "f_%sf(%s,%s,%s)" f c wf2 wf1

    let print_fermion_current2 coeff f c wf1 wf2 fusion =
      let c = format_coupling_2 coeff c in
      let c1 = fastener c 1 and
          c2 = fastener c 2 in
      match fusion with
      | F13 | F31 -> printf "%s_ff(%s,%s,%s,%s)" f c1 c2 wf1 wf2
      | F23 | F21 -> printf "f_%sf(%s,%s,%s,%s)" f c1 c2 wf1 wf2
      | F32 | F12 -> printf "f_%sf(%s,%s,%s,%s)" f c1 c2 wf2 wf1

    let print_fermion_current_mom_v1 coeff f c wf1 wf2 p1 p2 p12 fusion =
      let c = format_coupling coeff c in
      let c1 = fastener c 1 and
          c2 = fastener c 2 in
      match fusion with
      | F13 -> printf "%s_ff(%s,%s,%s,%s)" f c1 c2 wf1 wf2
      | F31 -> printf "%s_ff(-(%s),%s,%s,%s)" f c1 c2 wf1 wf2
      | F23 -> printf "f_%sf(%s,%s,%s,%s)" f c1 c2 wf1 wf2
      | F32 -> printf "f_%sf(%s,%s,%s,%s)" f c1 c2 wf2 wf1
      | F12 -> printf "f_f%s(-(%s),%s,%s,%s)" f c1 c2 wf2 wf1
      | F21 -> printf "f_f%s(-(%s),%s,%s,%s)" f c1 c2 wf1 wf2

    let print_fermion_current_mom_v1_chiral coeff f c wf1 wf2 p1 p2 p12 fusion =
      let c = format_coupling coeff c in
      let c1 = fastener c 1 and
          c2 = fastener c 2 in
      match fusion with
      | F13 -> printf "%s_ff(%s,%s,%s,%s)" f c1 c2 wf1 wf2
      | F31 -> printf "%s_ff(-(%s),-(%s),%s,%s)" f c2 c1 wf1 wf2
      | F23 -> printf "f_%sf(%s,%s,%s,%s)" f c1 c2 wf1 wf2
      | F32 -> printf "f_%sf(%s,%s,%s,%s)" f c1 c2 wf2 wf1
      | F12 -> printf "f_f%s(-(%s),-(%s),%s,%s)" f c2 c1 wf2 wf1
      | F21 -> printf "f_f%s(-(%s),-(%s),%s,%s)" f c2 c1 wf2 wf1

    let print_fermion_current_mom_v2 coeff f c wf1 wf2 p1 p2 p12 fusion =
      let c = format_coupling coeff c in
      let c1 = fastener c 1 and
          c2 = fastener c 2 in
      match fusion with
      | F13 -> printf "%s_ff(%s,%s,%s,%s,%s)" f c1 c2 wf1 wf2 p12
      | F31 -> printf "%s_ff(-(%s),%s,%s,%s,%s)" f c1 c2 wf1 wf2 p12
      | F23 -> printf "f_%sf(%s,%s,%s,%s,%s)" f c1 c2 wf1 wf2 p1
      | F32 -> printf "f_%sf(%s,%s,%s,%s,%s)" f c1 c2 wf2 wf1 p2
      | F12 -> printf "f_f%s(-(%s),%s,%s,%s,%s)" f c1 c2 wf2 wf1 p2
      | F21 -> printf "f_f%s(-(%s),%s,%s,%s,%s)" f c1 c2 wf1 wf2 p1

    let print_fermion_current_mom_v2_chiral coeff f c wf1 wf2 p1 p2 p12 fusion =
      let c = format_coupling coeff c in
      let c1 = fastener c 1 and
          c2 = fastener c 2 in
      match fusion with
      | F13 -> printf "%s_ff(%s,%s,%s,%s,%s)" f c1 c2 wf1 wf2 p12
      | F31 -> printf "%s_ff(-(%s),-(%s),%s,%s,%s)" f c2 c1 wf2 wf1 p12
      | F23 -> printf "f_%sf(%s,%s,%s,%s,%s)" f c1 c2 wf1 wf2 p1
      | F32 -> printf "f_%sf(%s,%s,%s,%s,%s)" f c1 c2 wf2 wf1 p2
      | F12 -> printf "f_f%s(-(%s),-(%s),%s,%s,%s)" f c2 c1 wf1 wf2 p2
      | F21 -> printf "f_f%s(-(%s),-(%s),%s,%s,%s)" f c2 c1 wf2 wf1 p1

    let print_fermion_current_vector coeff f c wf1 wf2 fusion =
      let c = format_coupling coeff c in
      match fusion with
      | F13 -> printf "%s_ff(%s,%s,%s)" f c wf1 wf2
      | F31 -> printf "%s_ff(-%s,%s,%s)" f c wf1 wf2
      | F23 -> printf "f_%sf(%s,%s,%s)" f c wf1 wf2
      | F32 -> printf "f_%sf(%s,%s,%s)" f c wf2 wf1
      | F12 -> printf "f_%sf(-%s,%s,%s)" f c wf2 wf1
      | F21 -> printf "f_%sf(-%s,%s,%s)" f c wf1 wf2

    let print_fermion_current2_vector coeff f c wf1 wf2 fusion =
      let c  = format_coupling_2 coeff c in
      let c1 = fastener c 1 and
          c2 = fastener c 2 in
      match fusion with
      | F13 -> printf "%s_ff(%s,%s,%s,%s)" f c1 c2 wf1 wf2
      | F31 -> printf "%s_ff(-(%s),%s,%s,%s)" f c1 c2 wf1 wf2
      | F23 -> printf "f_%sf(%s,%s,%s,%s)" f c1 c2 wf1 wf2
      | F32 -> printf "f_%sf(%s,%s,%s,%s)" f c1 c2 wf2 wf1
      | F12 -> printf "f_%sf(-(%s),%s,%s,%s)" f c1 c2 wf2 wf1
      | F21 -> printf "f_%sf(-(%s),%s,%s,%s)" f c1 c2 wf1 wf2

    let print_fermion_current_chiral coeff f1 f2 c wf1 wf2 fusion =
      let c = format_coupling coeff c in
      match fusion with
      | F13 -> printf "%s_ff(%s,%s,%s)" f1 c wf1 wf2
      | F31 -> printf "%s_ff(-%s,%s,%s)" f2 c wf1 wf2
      | F23 -> printf "f_%sf(%s,%s,%s)" f1 c wf1 wf2
      | F32 -> printf "f_%sf(%s,%s,%s)" f1 c wf2 wf1
      | F12 -> printf "f_%sf(-%s,%s,%s)" f2 c wf2 wf1
      | F21 -> printf "f_%sf(-%s,%s,%s)" f2 c wf1 wf2

    let print_fermion_current2_chiral coeff f c wf1 wf2 fusion =
      let c = format_coupling_2 coeff c in
      let c1 = fastener c 1 and
          c2 = fastener c 2 in
      match fusion with
      | F13 -> printf "%s_ff(%s,%s,%s,%s)" f c1 c2 wf1 wf2
      | F31 -> printf "%s_ff(-(%s),-(%s),%s,%s)" f c2 c1 wf1 wf2
      | F23 -> printf "f_%sf(%s,%s,%s,%s)" f c1 c2 wf1 wf2
      | F32 -> printf "f_%sf(%s,%s,%s,%s)" f c1 c2 wf2 wf1
      | F12 -> printf "f_%sf(-(%s),-(%s),%s,%s)" f c2 c1 wf2 wf1
      | F21 -> printf "f_%sf(-(%s),-(%s),%s,%s)" f c2 c1 wf1 wf2

    let print_current = function
      | coeff, _, VA, _ -> print_fermion_current2_vector coeff "va"
      | coeff, _, V, _ -> print_fermion_current_vector coeff "v"
      | coeff, _, A, _ -> print_fermion_current coeff "a"
      | coeff, _, VL, _ -> print_fermion_current_chiral coeff "vl" "vr"
      | coeff, _, VR, _ -> print_fermion_current_chiral coeff "vr" "vl"
      | coeff, _, VLR, _ -> print_fermion_current2_chiral coeff "vlr"
      | coeff, _, SP, _ -> print_fermion_current2 coeff "sp"
      | coeff, _, S, _ -> print_fermion_current coeff "s"
      | coeff, _, P, _ -> print_fermion_current coeff "p"
      | coeff, _, SL, _ -> print_fermion_current coeff "sl"
      | coeff, _, SR, _ -> print_fermion_current coeff "sr"
      | coeff, _, SLR, _ -> print_fermion_current2 coeff "slr"
      | coeff, _, POT, _ -> print_fermion_current_vector coeff "pot"
      | _, _, _, _ -> invalid_arg
            "Targets.Fortran_Majorana_Fermions: Not needed in the models"

    let print_current_p = function
      | coeff, Psi, SL, Psi -> print_fermion_current coeff "sl"
      | coeff, Psi, SR, Psi -> print_fermion_current coeff "sr"
      | coeff, Psi, SLR, Psi -> print_fermion_current2 coeff "slr"
      | _, _, _, _ -> invalid_arg
            "Targets.Fortran_Majorana_Fermions: Not needed in the used models"

    let print_current_b = function
      | coeff, Psibar, SL, Psibar -> print_fermion_current coeff "sl"
      | coeff, Psibar, SR, Psibar -> print_fermion_current coeff "sr"
      | coeff, Psibar, SLR, Psibar -> print_fermion_current2 coeff "slr"
      | _, _, _, _  -> invalid_arg
            "Targets.Fortran_Majorana_Fermions: Not needed in the used models"

(* This function is for the vertices with three particles including two
   fermions but also a momentum, therefore with a dimensionful coupling
   constant, e.g. the gravitino vertices. One has to dinstinguish between
   the two kinds of canonical orders in the string of gamma matrices. Of
   course, the direction of the string of gamma matrices is reversed if one
   goes from the [Gravbar, _, Psi] to the [Psibar, _, Grav] vertices, and
   the same is true for the couplings of the gravitino to the Majorana
   fermions. For more details see the tables in the [coupling]
   implementation. *)

(* We now have to fix the directions of the momenta. For making the compiler
   happy and because we don't want to make constructions of infinite
   complexity we list the momentum including vertices without gravitinos
   here; the pattern matching says that's better. Perhaps we have to find a
   better name now.

   For the cases of $MOM$, $MOM5$, $MOML$ and $MOMR$ which arise only in
   BRST transformations we take the mass as a coupling constant. For
   $VMOM$ we don't need a mass either. These vertices are like kinetic terms
   and so need not have a coupling constant. By this we avoid a strange and
   awful construction with a new variable. But be careful with a
   generalization if you want to use these vertices for other purposes.
*)

    let format_coupling_mom coeff c =
      match coeff with
      | 1 -> c
      | -1 -> "(-" ^ c ^")"
      | coeff -> string_of_int coeff ^ "*" ^ c

    let commute_proj f =
      match f with
      | "moml" -> "lmom"
      | "momr" -> "rmom"
      | "lmom" -> "moml"
      | "rmom" -> "momr"
      | "svl"  -> "svr"
      | "svr"  -> "svl"
      | "sl" -> "sr"
      | "sr" -> "sl"
      | "s" -> "s"
      | "p" -> "p"
      | _ -> invalid_arg "Targets:Fortran_Majorana_Fermions: wrong case"

    let print_fermion_current_mom coeff f c wf1 wf2 p1 p2 p12 fusion =
      let c = format_coupling_mom coeff c in
      let c1 = fastener c 1 and
          c2 = fastener c 2 in
      match fusion with
      | F13 -> printf "%s_ff(%s,%s,%s,%s,%s)" f c1 c2 wf1 wf2 p12
      | F31 -> printf "%s_ff(%s,%s,%s,%s,%s)" f c1 c2 wf1 wf2 p12
      | F23 -> printf "f_%sf(%s,%s,%s,%s,%s)" f c1 c2 wf1 wf2 p1
      | F32 -> printf "f_%sf(%s,%s,%s,%s,%s)" f c1 c2 wf2 wf1 p2
      | F12 -> printf "f_%sf(%s,%s,%s,%s,%s)" f c1 c2 wf2 wf1 p2
      | F21 -> printf "f_%sf(%s,%s,%s,%s,%s)" f c1 c2 wf1 wf2 p1

    (*i unused value
    let print_fermion_current_mom_vector coeff f c wf1 wf2 p1 p2 p12 fusion =
      let c = format_coupling_mom coeff c in
      let c1 = fastener c 1 and
          c2 = fastener c 2 in
      match fusion with
      | F13 -> printf "%s_ff(%s,%s,%s,%s,%s)" f c1 c2 wf1 wf2 p12
      | F31 -> printf "%s_ff(-%s,%s,%s,%s,%s)" f c1 c2 wf1 wf2 p12
      | F23 -> printf "f_%sf(%s,%s,%s,%s,%s)" f c1 c2 wf1 wf2 p1
      | F32 -> printf "f_%sf(%s,%s,%s,%s,%s)" f c1 c2 wf2 wf1 p2
      | F12 -> printf "f_%sf(-%s,%s,%s,%s,%s)" f c1 c2 wf2 wf1 p2
      | F21 -> printf "f_%sf(-%s,%s,%s,%s,%s)" f c1 c2 wf1 wf2 p1
      i*)

    let print_fermion_current_mom_sign coeff f c wf1 wf2 p1 p2 p12 fusion =
      let c = format_coupling_mom coeff c in
      let c1 = fastener c 1 and
          c2 = fastener c 2 in
      match fusion with
      | F13 -> printf "%s_ff(%s,%s,%s,%s,%s)" f c1 c2 wf1 wf2 p12
      | F31 -> printf "%s_ff(%s,%s,%s,%s,-(%s))" f c1 c2 wf1 wf2 p12
      | F23 -> printf "f_%sf(%s,%s,%s,%s,%s)" f c1 c2 wf1 wf2 p1
      | F32 -> printf "f_%sf(%s,%s,%s,%s,%s)" f c1 c2 wf2 wf1 p2
      | F12 -> printf "f_%sf(%s,%s,%s,%s,-(%s))" f c1 c2 wf2 wf1 p2
      | F21 -> printf "f_%sf(%s,%s,%s,%s,-(%s))" f c1 c2 wf1 wf2 p1

    let print_fermion_current_mom_sign_1 coeff f c wf1 wf2 p1 p2 p12 fusion =
      let c = format_coupling coeff c in
      match fusion with
      | F13 -> printf "%s_ff(%s,%s,%s,%s)" f c wf1 wf2 p12
      | F31 -> printf "%s_ff(%s,%s,%s,-(%s))" f c wf1 wf2 p12
      | F23 -> printf "f_%sf(%s,%s,%s,%s)" f c wf1 wf2 p1
      | F32 -> printf "f_%sf(%s,%s,%s,%s)" f c wf2 wf1 p2
      | F12 -> printf "f_%sf(%s,%s,%s,-(%s))" f c wf2 wf1 p2
      | F21 -> printf "f_%sf(%s,%s,%s,-(%s))" f c wf1 wf2 p1

    let print_fermion_current_mom_chiral coeff f c wf1 wf2 p1 p2 p12 fusion =
      let c  = format_coupling_mom coeff c and
          cf = commute_proj f in
      let c1 = fastener c 1 and
          c2 = fastener c 2 in
      match fusion with
      | F13 -> printf "%s_ff(%s,%s,%s,%s,%s)" f c1 c2 wf1 wf2 p12
      | F31 -> printf "%s_ff(%s,%s,%s, %s,-(%s))" cf c1 c2 wf1 wf2 p12
      | F23 -> printf "f_%sf(%s,%s,%s,%s,%s)" f c1 c2 wf1 wf2 p1
      | F32 -> printf "f_%sf(%s,%s,%s,%s,%s)" f c1 c2 wf2 wf1 p2
      | F12 -> printf "f_%sf(%s,%s,%s,%s,-(%s))" cf c1 c2 wf2 wf1 p2
      | F21 -> printf "f_%sf(%s,%s,%s,%s,-(%s))" cf c1 c2 wf1 wf2 p1

    let print_fermion_g_current coeff f c wf1 wf2 p1 p2 p12 fusion =
      let c = format_coupling coeff c in
      match fusion with
      | F13 -> printf "%s_grf(%s,%s,%s,%s)" f c wf1 wf2 p12
      | F31 -> printf "%s_fgr(%s,%s,%s,%s)" f c wf1 wf2 p12
      | F23 -> printf "gr_%sf(%s,%s,%s,%s)" f c wf1 wf2 p1
      | F32 -> printf "gr_%sf(%s,%s,%s,%s)" f c wf2 wf1 p2
      | F12 -> printf "f_%sgr(%s,%s,%s,%s)" f c wf2 wf1 p2
      | F21 -> printf "f_%sgr(%s,%s,%s,%s)" f c wf1 wf2 p1

    let print_fermion_g_2_current coeff f c wf1 wf2 p1 p2 p12 fusion =
      let c = format_coupling coeff c in
      match fusion with
      | F13 -> printf "%s_grf(%s(1),%s(2),%s,%s,%s)" f c c wf1 wf2 p12
      | F31 -> printf "%s_fgr(%s(1),%s(2),%s,%s,%s)" f c c wf1 wf2 p12
      | F23 -> printf "gr_%sf(%s(1),%s(2),%s,%s,%s)" f c c wf1 wf2 p1
      | F32 -> printf "gr_%sf(%s(1),%s(2),%s,%s,%s)" f c c wf2 wf1 p2
      | F12 -> printf "f_%sgr(%s(1),%s(2),%s,%s,%s)" f c c wf2 wf1 p2
      | F21 -> printf "f_%sgr(%s(1),%s(2),%s,%s,%s)" f c c wf1 wf2 p1

    let print_fermion_g_current_rev coeff f c wf1 wf2 p1 p2 p12 fusion =
      let c = format_coupling coeff c in
      match fusion with
      | F13 -> printf "%s_fgr(%s,%s,%s,%s)" f c wf1 wf2 p12
      | F31 -> printf "%s_grf(%s,%s,%s,%s)" f c wf1 wf2 p12
      | F23 -> printf "f_%sgr(%s,%s,%s,%s)" f c wf1 wf2 p1
      | F32 -> printf "f_%sgr(%s,%s,%s,%s)" f c wf2 wf1 p2
      | F12 -> printf "gr_%sf(%s,%s,%s,%s)" f c wf2 wf1 p2
      | F21 -> printf "gr_%sf(%s,%s,%s,%s)" f c wf1 wf2 p1

    let print_fermion_g_2_current_rev coeff f c wf1 wf2 p1 p2 p12 fusion =
      let c = format_coupling coeff c in
      match fusion with
      | F13 -> printf "%s_fgr(%s(1),%s(2),%s,%s,%s)" f c c wf1 wf2 p12
      | F31 -> printf "%s_grf(%s(1),%s(2),%s,%s,%s)" f c c wf1 wf2 p12
      | F23 -> printf "f_%sgr(%s(1),%s(2),%s,%s,%s)" f c c wf1 wf2 p1
      | F32 -> printf "f_%sgr(%s(1),%s(2),%s,%s,%s)" f c c wf2 wf1 p2
      | F12 -> printf "gr_%sf(%s(1),%s(2),%s,%s,%s)" f c c wf2 wf1 p2
      | F21 -> printf "gr_%sf(%s(1),%s(2),%s,%s,%s)" f c c wf1 wf2 p1

    let print_fermion_g_current_vector coeff f c wf1 wf2 _ _ _ fusion =
      let c = format_coupling coeff c in
      match fusion with
      | F13 -> printf "%s_grf(%s,%s,%s)" f c wf1 wf2
      | F31 -> printf "%s_fgr(-%s,%s,%s)" f c wf1 wf2
      | F23 -> printf "gr_%sf(%s,%s,%s)" f c wf1 wf2
      | F32 -> printf "gr_%sf(%s,%s,%s)" f c wf2 wf1
      | F12 -> printf "f_%sgr(-%s,%s,%s)" f c wf2 wf1
      | F21 -> printf "f_%sgr(-%s,%s,%s)" f c wf1 wf2

    let print_fermion_g_current_vector_rev coeff f c wf1 wf2 _ _ _ fusion =
      let c = format_coupling coeff c in
      match fusion with
      | F13 -> printf "%s_fgr(%s,%s,%s)" f c wf1 wf2
      | F31 -> printf "%s_grf(-%s,%s,%s)" f c wf1 wf2
      | F23 -> printf "f_%sgr(%s,%s,%s)" f c wf1 wf2
      | F32 -> printf "f_%sgr(%s,%s,%s)" f c wf2 wf1
      | F12 -> printf "gr_%sf(-%s,%s,%s)" f c wf2 wf1
      | F21 -> printf "gr_%sf(-%s,%s,%s)" f c wf1 wf2

    let print_current_g = function
      | coeff, _, MOM, _ -> print_fermion_current_mom_sign coeff "mom"
      | coeff, _, MOM5, _ -> print_fermion_current_mom coeff "mom5"
      | coeff, _, MOML, _ -> print_fermion_current_mom_chiral coeff "moml"
      | coeff, _, MOMR, _ -> print_fermion_current_mom_chiral coeff "momr"
      | coeff, _, LMOM, _ -> print_fermion_current_mom_chiral coeff "lmom"
      | coeff, _, RMOM, _ -> print_fermion_current_mom_chiral coeff "rmom"
      | coeff, _, VMOM, _ -> print_fermion_current_mom_sign_1 coeff "vmom"
      | coeff, Gravbar, S, _ -> print_fermion_g_current coeff "s"
      | coeff, Gravbar, SL, _ -> print_fermion_g_current coeff "sl"
      | coeff, Gravbar, SR, _ -> print_fermion_g_current coeff "sr"
      | coeff, Gravbar, SLR, _ -> print_fermion_g_2_current coeff "slr"
      | coeff, Gravbar, P, _ -> print_fermion_g_current coeff "p"
      | coeff, Gravbar, V, _ -> print_fermion_g_current coeff "v"
      | coeff, Gravbar, VLR, _ -> print_fermion_g_2_current coeff "vlr"
      | coeff, Gravbar, POT, _ -> print_fermion_g_current_vector coeff "pot"
      | coeff, _, S, Grav -> print_fermion_g_current_rev coeff "s"
      | coeff, _, SL, Grav -> print_fermion_g_current_rev coeff "sl"
      | coeff, _, SR, Grav -> print_fermion_g_current_rev coeff "sr"
      | coeff, _, SLR, Grav -> print_fermion_g_2_current_rev coeff "slr"
      | coeff, _, P, Grav -> print_fermion_g_current_rev (-coeff) "p"
      | coeff, _, V, Grav -> print_fermion_g_current_rev coeff "v"
      | coeff, _, VLR, Grav -> print_fermion_g_2_current_rev coeff "vlr"
      | coeff, _, POT, Grav -> print_fermion_g_current_vector_rev coeff "pot"
      | _, _, _, _ -> invalid_arg
          "Targets.Fortran_Majorana_Fermions: not used in the models"

    let print_current_mom = function
      | coeff, _, TVA, _ -> print_fermion_current_mom_v1 coeff "tva"
      | coeff, _, TVAM, _ -> print_fermion_current_mom_v2 coeff "tvam"
      | coeff, _, TLR, _ -> print_fermion_current_mom_v1_chiral coeff "tlr"
      | coeff, _, TLRM, _ -> print_fermion_current_mom_v2_chiral coeff "tlrm"
      | _, _, _, _ -> invalid_arg
            "Targets.Fortran_Majorana_Fermions: Not needed in the models"

(* We need support for dimension-5 vertices with two fermions and two
   bosons, appearing in theories of supergravity and also together with in
   insertions of the supersymmetric current. There is a canonical order
   [fermionbar], [boson_1], [boson_2], [fermion], so what one has to do is a
   mapping from the fusions [F123] etc. to the order of the three wave
   functions [wf1], [wf2] and [wf3]. *)

(* The function [d_p] (for distinct the particle) distinguishes which particle
   (scalar or vector) must be fused to in the special functions. *)

    let d_p = function
      | 1, ("sv"|"pv"|"svl"|"svr"|"slrv") -> "1"
      | 1, _ -> ""
      | 2, ("sv"|"pv"|"svl"|"svr"|"slrv") -> "2"
      | 2, _ -> ""
      | _, _ -> invalid_arg "Targets.Fortran_Majorana_Fermions: not used"

    let wf_of_f wf1 wf2 wf3 f =
      match f with
      | (F123|F423) -> [wf2; wf3; wf1]
      | (F213|F243|F143|F142|F413|F412) -> [wf1; wf3; wf2]
      | (F132|F432) -> [wf3; wf2; wf1]
      | (F231|F234|F134|F124|F431|F421) -> [wf1; wf2; wf3]
      | (F312|F342) -> [wf3; wf1; wf2]
      | (F321|F324|F314|F214|F341|F241) -> [wf2; wf1; wf3]

    let print_fermion_g4_brs_vector_current coeff f c wf1 wf2 wf3 fusion =
      let cf = commute_proj f and
          cp = format_coupling coeff c and
          cm = if f = "pv" then
            format_coupling coeff c
          else
            format_coupling (-coeff) c
      and
          d1 = d_p (1,f) and
          d2 = d_p (2,f) and
          f1 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 0) and
          f2 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 1) and
          f3 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 2) in
      match fusion with
      | (F123|F213|F132|F231|F312|F321) ->
          printf "f_%sf(%s,%s,%s,%s)" cf cm f1 f2 f3
      | (F423|F243|F432|F234|F342|F324) ->
          printf "f_%sf(%s,%s,%s,%s)" f cp f1 f2 f3
      | (F134|F143|F314) -> printf "%s%s_ff(%s,%s,%s,%s)" f d1 cp f1 f2 f3
      | (F124|F142|F214) -> printf "%s%s_ff(%s,%s,%s,%s)" f d2 cp f1 f2 f3
      | (F413|F431|F341) -> printf "%s%s_ff(%s,%s,%s,%s)" cf d1 cm f1 f2 f3
      | (F241|F412|F421) -> printf "%s%s_ff(%s,%s,%s,%s)" cf d2 cm f1 f2 f3

    let print_fermion_g4_svlr_current coeff _ c wf1 wf2 wf3 fusion =
      let c = format_coupling_2 coeff c and
          f1 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 0) and
          f2 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 1) and
          f3 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 2) in
      let c1 = fastener c 1 and
          c2 = fastener c 2 in
      match fusion with
      | (F123|F213|F132|F231|F312|F321) ->
          printf "f_svlrf(-(%s),-(%s),%s,%s,%s)" c2 c1 f1 f2 f3
      | (F423|F243|F432|F234|F342|F324) ->
          printf "f_svlrf(%s,%s,%s,%s,%s)" c1 c2 f1 f2 f3
      | (F134|F143|F314) ->
          printf "svlr2_ff(%s,%s,%s,%s,%s)" c1 c2 f1 f2 f3
      | (F124|F142|F214) ->
          printf "svlr1_ff(%s,%s,%s,%s,%s)" c1 c2 f1 f2 f3
      | (F413|F431|F341) ->
          printf "svlr2_ff(-(%s),-(%s),%s,%s,%s)" c2 c1 f1 f2 f3
      | (F241|F412|F421) ->
          printf "svlr1_ff(-(%s),-(%s),%s,%s,%s)" c2 c1 f1 f2 f3

    let print_fermion_s2_current coeff f c wf1 wf2 wf3 fusion =
      let cp = format_coupling coeff c and
          cm = if f = "p" then
            format_coupling (-coeff) c
          else
            format_coupling coeff c
      and
          cf = commute_proj f and
          f1 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 0) and
          f2 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 1) and
          f3 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 2) in
      match fusion with
      | (F123|F213|F132|F231|F312|F321) ->
          printf "%s * f_%sf(%s,%s,%s)" f1 cf cm f2 f3
      | (F423|F243|F432|F234|F342|F324) ->
          printf "%s * f_%sf(%s,%s,%s)" f1 f cp f2 f3
      | (F134|F143|F314) ->
          printf "%s * %s_ff(%s,%s,%s)" f2 f cp f1 f3
      | (F124|F142|F214) ->
          printf "%s * %s_ff(%s,%s,%s)" f2 f cp f1 f3
      | (F413|F431|F341) ->
          printf "%s * %s_ff(%s,%s,%s)" f2 cf cm f1 f3
      | (F241|F412|F421) ->
          printf "%s * %s_ff(%s,%s,%s)" f2 cf cm f1 f3

    let print_fermion_s2p_current coeff f c wf1 wf2 wf3 fusion =
      let c = format_coupling_2 coeff c and
          f1 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 0) and
          f2 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 1) and
          f3 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 2) in
      let c1 = fastener c 1 and
          c2 = fastener c 2 in
      match fusion with
      | (F123|F213|F132|F231|F312|F321) ->
          printf "%s * f_%sf(%s,-(%s),%s,%s)" f1 f c1 c2 f2 f3
      | (F423|F243|F432|F234|F342|F324) ->
          printf "%s * f_%sf(%s,%s,%s,%s)" f1 f c1 c2 f2 f3
      | (F134|F143|F314) ->
          printf "%s * %s_ff(%s,%s,%s,%s)" f2 f c1 c2 f1 f3
      | (F124|F142|F214) ->
          printf "%s * %s_ff(%s,%s,%s,%s)" f2 f c1 c2 f1 f3
      | (F413|F431|F341) ->
          printf "%s * %s_ff(%s,-(%s),%s,%s)" f2 f c1 c2 f1 f3
      | (F241|F412|F421) ->
          printf "%s * %s_ff(%s,-(%s),%s,%s)" f2 f c1 c2 f1 f3

    let print_fermion_s2lr_current coeff f c wf1 wf2 wf3 fusion =
      let c = format_coupling_2 coeff c and
          f1 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 0) and
          f2 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 1) and
          f3 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 2) in
      let c1 = fastener c 1 and
          c2 = fastener c 2 in
      match fusion with
      | (F123|F213|F132|F231|F312|F321) ->
          printf "%s * f_%sf(%s,%s,%s,%s)" f1 f c2 c1 f2 f3
      | (F423|F243|F432|F234|F342|F324) ->
          printf "%s * f_%sf(%s,%s,%s,%s)" f1 f c1 c2 f2 f3
      | (F134|F143|F314) ->
          printf "%s * %s_ff(%s,%s,%s,%s)" f2 f c1 c2 f1 f3
      | (F124|F142|F214) ->
          printf "%s * %s_ff(%s,%s,%s,%s)" f2 f c1 c2 f1 f3
      | (F413|F431|F341) ->
          printf "%s * %s_ff(%s,%s,%s,%s)" f2 f c2 c1 f1 f3
      | (F241|F412|F421) ->
          printf "%s * %s_ff(%s,%s,%s,%s)" f2 f c2 c1 f1 f3

    let print_fermion_g4_current coeff f c wf1 wf2 wf3 fusion =
      let c = format_coupling coeff c and
          f1 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 0) and
          f2 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 1) and
          f3 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 2) in
      match fusion with
      | (F123|F213|F132|F231|F312|F321) ->
          printf "f_%sgr(-%s,%s,%s,%s)" f c f1 f2 f3
      | (F423|F243|F432|F234|F342|F324) ->
          printf "gr_%sf(%s,%s,%s,%s)" f c f1 f2 f3
      | (F134|F143|F314|F124|F142|F214) ->
          printf "%s_grf(%s,%s,%s,%s)" f c f1 f2 f3
      | (F413|F431|F341|F241|F412|F421) ->
          printf "%s_fgr(-%s,%s,%s,%s)" f c f1 f2 f3

    (*i unused value
    let print_fermion_2_g4_current coeff f c wf1 wf2 wf3 fusion =
      let f1 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 0) and
          f2 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 1) and
          f3 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 2) in
      let c = format_coupling_2 coeff c in
      let c1 = fastener c 1 and
          c2 = fastener c 2 in
      match fusion with
      | (F123|F213|F132|F231|F312|F321) ->
          printf "f_%sgr(-(%s),-(%s),%s,%s,%s)" f c2 c1 f1 f2 f3
      | (F423|F243|F432|F234|F342|F324) ->
          printf "gr_%sf(%s,%s,%s,%s,%s)" f c1 c2 f1 f2 f3
      | (F134|F143|F314|F124|F142|F214) ->
          printf "%s_grf(%s,%s,%s,%s,%s)" f c1 c2 f1 f2 f3
      | (F413|F431|F341|F241|F412|F421) ->
          printf "%s_fgr(-(%s),-(%s),%s,%s,%s)" f c2 c1 f1 f2 f3
    i*)

    let print_fermion_2_g4_current coeff f c wf1 wf2 wf3 fusion =
      let f1 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 0) and
          f2 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 1) and
          f3 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 2) in
      let c = format_coupling_2 coeff c in
      let c1 = fastener c 1 and
          c2 = fastener c 2 in
      match fusion with
      | (F123|F213|F132|F231|F312|F321) ->
          printf "f_%sgr(-(%s),-(%s),%s,%s,%s)" f c2 c1 f1 f2 f3
      | (F423|F243|F432|F234|F342|F324) ->
          printf "gr_%sf(%s,%s,%s,%s,%s)" f c1 c2 f1 f2 f3
      | (F134|F143|F314|F124|F142|F214) ->
          printf "%s_grf(%s,%s,%s,%s,%s)" f c1 c2 f1 f2 f3
      | (F413|F431|F341|F241|F412|F421) ->
          printf "%s_fgr(-(%s),-(%s),%s,%s,%s)" f c2 c1 f1 f2 f3


    let print_fermion_g4_current_rev coeff f c wf1 wf2 wf3 fusion =
      let c = format_coupling coeff c and
          f1 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 0) and
          f2 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 1) and
          f3 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 2) in
      match fusion with
      | (F123|F213|F132|F231|F312|F321) ->
          printf "f_%sgr(%s,%s,%s,%s)" f c f1 f2 f3
      | (F423|F243|F432|F234|F342|F324) ->
          printf "gr_%sf(-%s,%s,%s,%s)" f c f1 f2 f3
      | (F134|F143|F314|F124|F142|F214) ->
          printf "%s_grf(-%s,%s,%s,%s)" f c f1 f2 f3
      | (F413|F431|F341|F241|F412|F421) ->
          printf "%s_fgr(%s,%s,%s,%s)" f c f1 f2 f3

(* Here we have to distinguish which of the two bosons is produced in the
   fusion of three particles which include both fermions. *)

    let print_fermion_g4_vector_current coeff f c wf1 wf2 wf3 fusion =
      let c = format_coupling coeff c and
          d1 = d_p (1,f) and
          d2 = d_p (2,f) and
          f1 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 0) and
          f2 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 1) and
          f3 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 2) in
      match fusion with
      | (F123|F213|F132|F231|F312|F321) ->
          printf "f_%sgr(%s,%s,%s,%s)" f c f1 f2 f3
      | (F423|F243|F432|F234|F342|F324) ->
          printf "gr_%sf(%s,%s,%s,%s)" f c f1 f2 f3
      | (F134|F143|F314) -> printf "%s%s_grf(%s,%s,%s,%s)" f d1 c f1 f2 f3
      | (F124|F142|F214) -> printf "%s%s_grf(%s,%s,%s,%s)" f d2 c f1 f2 f3
      | (F413|F431|F341) -> printf "%s%s_fgr(%s,%s,%s,%s)" f d1 c f1 f2 f3
      | (F241|F412|F421) -> printf "%s%s_fgr(%s,%s,%s,%s)" f d2 c f1 f2 f3

    let print_fermion_2_g4_vector_current coeff f c wf1 wf2 wf3 fusion =
      let d1 = d_p (1,f) and
          d2 = d_p (2,f) and
          f1 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 0) and
          f2 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 1) and
          f3 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 2) in
      let c = format_coupling_2 coeff c in
      let c1 = fastener c 1 and
          c2 = fastener c 2 in
      match fusion with
      | (F123|F213|F132|F231|F312|F321) ->
          printf "f_%sgr(%s,%s,%s,%s,%s)" f c1 c2 f1 f2 f3
      | (F423|F243|F432|F234|F342|F324) ->
          printf "gr_%sf(%s,%s,%s,%s,%s)" f c1 c2 f1 f2 f3
      | (F134|F143|F314) -> printf "%s%s_grf(%s,%s,%s,%s,%s)" f d1 c1 c2 f1 f2 f3
      | (F124|F142|F214) -> printf "%s%s_grf(%s,%s,%s,%s,%s)" f d2 c1 c2 f1 f2 f3
      | (F413|F431|F341) -> printf "%s%s_fgr(%s,%s,%s,%s,%s)" f d1 c1 c2 f1 f2 f3
      | (F241|F412|F421) -> printf "%s%s_fgr(%s,%s,%s,%s,%s)" f d2 c1 c2 f1 f2 f3

    let print_fermion_g4_vector_current_rev coeff f c wf1 wf2 wf3 fusion =
      let c = format_coupling coeff c and
          d1 = d_p (1,f) and
          d2 = d_p (2,f) and
          f1 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 0) and
          f2 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 1) and
          f3 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 2) in
      match fusion with
      | (F123|F213|F132|F231|F312|F321) ->
          printf "gr_%sf(%s,%s,%s,%s)" f c f1 f2 f3
      | (F423|F243|F432|F234|F342|F324) ->
          printf "f_%sgr(%s,%s,%s,%s)" f c f1 f2 f3
      | (F134|F143|F314) -> printf "%s%s_fgr(%s,%s,%s,%s)" f d1 c f1 f2 f3
      | (F124|F142|F214) -> printf "%s%s_fgr(%s,%s,%s,%s)" f d2 c f1 f2 f3
      | (F413|F431|F341) -> printf "%s%s_grf(%s,%s,%s,%s)" f d1 c f1 f2 f3
      | (F241|F412|F421) -> printf "%s%s_grf(%s,%s,%s,%s)" f d2 c f1 f2 f3

    let print_fermion_2_g4_current_rev coeff f c wf1 wf2 wf3 fusion =
      let c = format_coupling_2 coeff c in
      let c1 = fastener c 1 and
          c2 = fastener c 2  and
          d1 = d_p (1,f) and
          d2 = d_p (2,f) in
      let f1 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 0) and
          f2 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 1) and
          f3 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 2) in
      match fusion with
      | (F123|F213|F132|F231|F312|F321) ->
          printf "gr_%sf(%s,%s,%s,%s,%s)" f c1 c2 f1 f2 f3
      | (F423|F243|F432|F234|F342|F324) ->
          printf "f_%sgr(-(%s),-(%s),%s,%s,%s)" f c1 c2 f1 f2 f3
      | (F134|F143|F314) ->
          printf "%s%s_fgr(-(%s),-(%s),%s,%s,%s)" f d1 c1 c2 f1 f2 f3
      | (F124|F142|F214) ->
          printf "%s%s_fgr(-(%s),-(%s),%s,%s,%s)" f d2 c1 c2 f1 f2 f3
      | (F413|F431|F341) ->
          printf "%s%s_grf(%s,%s,%s,%s,%s)" f d1 c1 c2 f1 f2 f3
      | (F241|F412|F421) ->
          printf "%s%s_grf(%s,%s,%s,%s,%s)" f d2 c1 c2 f1 f2 f3

    let print_fermion_2_g4_vector_current_rev coeff f c wf1 wf2 wf3 fusion =
      (* Here we put in the extra minus sign from the coeff. *)
      let c = format_coupling coeff c in
      let c1 = fastener c 1 and
          c2 = fastener c 2 in
      let d1 = d_p (1,f) and
          d2 = d_p (2,f) and
          f1 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 0) and
          f2 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 1) and
          f3 = (List.nth (wf_of_f wf1 wf2 wf3 fusion) 2) in
      match fusion with
      | (F123|F213|F132|F231|F312|F321) ->
          printf "gr_%sf(%s,%s,%s,%s,%s)" f c1 c2 f1 f2 f3
      | (F423|F243|F432|F234|F342|F324) ->
          printf "f_%sgr(%s,%s,%s,%s,%s)" f c1 c2 f1 f2 f3
      | (F134|F143|F314) -> printf "%s%s_fgr(%s,%s,%s,%s,%s)" f d1 c1 c2 f1 f2 f3
      | (F124|F142|F214) -> printf "%s%s_fgr(%s,%s,%s,%s,%s)" f d2 c1 c2 f1 f2 f3
      | (F413|F431|F341) -> printf "%s%s_grf(%s,%s,%s,%s,%s)" f d1 c1 c2 f1 f2 f3
      | (F241|F412|F421) -> printf "%s%s_grf(%s,%s,%s,%s,%s)" f d2 c1 c2 f1 f2 f3


    let print_current_g4 = function
      | coeff, Gravbar, S2, _ -> print_fermion_g4_current coeff "s2"
      | coeff, Gravbar, SV, _ -> print_fermion_g4_vector_current coeff "sv"
      | coeff, Gravbar, SLV, _ -> print_fermion_g4_vector_current coeff "slv"
      | coeff, Gravbar, SRV, _ -> print_fermion_g4_vector_current coeff "srv"
      | coeff, Gravbar, SLRV, _ -> print_fermion_2_g4_vector_current coeff "slrv"
      | coeff, Gravbar, PV, _ -> print_fermion_g4_vector_current coeff "pv"
      | coeff, Gravbar, V2, _ -> print_fermion_g4_current coeff "v2"
      | coeff, Gravbar, V2LR, _ -> print_fermion_2_g4_current coeff "v2lr"
      | _, Gravbar, _, _ -> invalid_arg "print_current_g4: not implemented"
      | coeff, _, S2, Grav -> print_fermion_g4_current_rev coeff "s2"
      | coeff, _, SV, Grav -> print_fermion_g4_vector_current_rev (-coeff) "sv"
      | coeff, _, SLV, Grav -> print_fermion_g4_vector_current_rev (-coeff) "slv"
      | coeff, _, SRV, Grav -> print_fermion_g4_vector_current_rev (-coeff) "srv"
      | coeff, _, SLRV, Grav -> print_fermion_2_g4_vector_current_rev coeff "slrv"
      | coeff, _, PV, Grav -> print_fermion_g4_vector_current_rev coeff "pv"
      | coeff, _, V2, Grav -> print_fermion_g4_vector_current_rev coeff "v2"
      | coeff, _, V2LR, Grav -> print_fermion_2_g4_current_rev coeff "v2lr"
      | _, _, _, Grav -> invalid_arg "print_current_g4: not implemented"
      | coeff, _, S2, _ -> print_fermion_s2_current coeff "s"
      | coeff, _, P2, _ -> print_fermion_s2_current coeff "p"
      | coeff, _, S2P, _ -> print_fermion_s2p_current coeff "sp"
      | coeff, _, S2L, _ -> print_fermion_s2_current coeff "sl"
      | coeff, _, S2R, _ -> print_fermion_s2_current coeff "sr"
      | coeff, _, S2LR, _ -> print_fermion_s2lr_current coeff "slr"
      | coeff, _, V2, _ -> print_fermion_g4_brs_vector_current coeff "v2"
      | coeff, _, SV, _ -> print_fermion_g4_brs_vector_current coeff "sv"
      | coeff, _, PV, _ -> print_fermion_g4_brs_vector_current coeff "pv"
      | coeff, _, SLV, _ -> print_fermion_g4_brs_vector_current coeff "svl"
      | coeff, _, SRV, _ -> print_fermion_g4_brs_vector_current coeff "svr"
      | coeff, _, SLRV, _ -> print_fermion_g4_svlr_current coeff "svlr"
      | _, _, V2LR, _ -> invalid_arg "Targets.print_current: not available"

    let reverse_braket vintage bra ket =
      if vintage then
        false
      else
        match bra, ket with
        | Majorana, Majorana :: _ -> true
        | _, _ -> false

    let use_module = "omega95_bispinors"
    let require_library =
      ["omega_bispinors_2010_01_A"; "omega_bispinor_cpls_2010_01_A"]
   end

module Fortran_Majorana = Make_Fortran(Fortran_Majorana_Fermions)

(* \thocwmodulesubsection{\texttt{FORTRAN\,77}} *)

module Fortran77 = Dummy

(* \thocwmodulesection{\texttt{C}} *)

module C = Dummy

(* \thocwmodulesubsection{\texttt{C++}} *)

module Cpp = Dummy

(* \thocwmodulesubsection{Java} *)

module Java = Dummy

(* \thocwmodulesection{O'Caml} *)

module Ocaml = Dummy

(* \thocwmodulesection{\LaTeX} *)

module LaTeX = Dummy
