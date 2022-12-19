(* omega.ml --

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

let (<<) f g x = f (g x)
let (>>) f g x = g (f x)

module P = Momentum.Default
module P_Whizard = Momentum.DefaultW

module type T =
  sig
    val main : unit -> unit
    type flavor
    val diagrams : flavor -> flavor -> flavor list ->
      ((flavor * Momentum.Default.t) *
         (flavor * Momentum.Default.t,
          flavor * Momentum.Default.t) Tree.t) list
  end

module Make (Fusion_Maker : Fusion.Maker) (PHS_Maker : Fusion.Maker)
         (Target_Maker : Target.Maker) (M : Model.T) =
  struct

    module CM = Colorize.It(M)

    type flavor = M.flavor

    module Proc = Process.Make(M)

(* \begin{dubious}
     We must have initialized the vertices \emph{before}
     applying [Fusion_Maker], at least if we want to continue
     using the vertex cache!
   \end{dubious} *)

(* \begin{dubious}
     NB: this causes the constant initializers in [Fusion_Maker] more than once.
     Such side effects must be avoided if the initializers involve expensive
     computations.   \emph{Relying on the fact that the functor will be
     called only once is not a good idea!}
   \end{dubious} *)
    module F = Fusion_Maker(P)(M)
    module CF = Fusion.Multi(Fusion_Maker)(P)(M)
    module T = Target_Maker(Fusion_Maker)(P)(M)
    module W = Whizard.Make(Fusion_Maker)(P)(P_Whizard)(M)
    module C = Cascade.Make(M)(P)

    module VSet =
      Set.Make (struct type t = F.constant Coupling.t let compare = compare end)

(* For the phase space, we need asymmetric DAGs.

   Since we will not use this to compute amplitudes, there's
   no need to supply the proper statistics module and we may
   always use Majorana fermions to be as general as possible.
   In principle, we could expose in [Fusion.T] the [Fusion.Stat_Maker]
   used by [Fusion_Maker] to construct it, but that is just not
   worth the effort.

   \begin{dubious}
     For the phase space, we should be able to work on the
     uncolored model.
   \end{dubious} *)

    module MT = Modeltools.Topology3(M)
    module PHS = PHS_Maker(P)(MT)
    module CT = Cascade.Make(MT)(P)

(* Form a ['a list] from a ['a option array], containing
   the elements that are not [None] in order. *)

    let opt_array_to_list a =
      let rec opt_array_to_list' acc i a =
	if i < 0 then
	  acc
	else
	  begin match a.(i) with
	  | None -> opt_array_to_list' acc (pred i) a
	  | Some x -> opt_array_to_list' (x :: acc) (pred i) a
	  end in
      opt_array_to_list' [] (Array.length a - 1) a
  
(* Return a list of [CF.amplitude list]s, corresponig to the diagrams
   for a specific color flow for each flavor combination. *)

    let amplitudes_by_flavor amplitudes =
      List.map opt_array_to_list (Array.to_list (CF.process_table amplitudes))

(* \begin{dubious}
     If we plan to distiguish different couplings later on,
     we can no long map all instances of [coupling option]
     in the tree to [None].  In this case, we will need
     to normalize different fusion orders [Coupling.fuse2],
     [Coupling.fuse3] or [Coupling.fusen], because they would
     otherwise lead to inequivalent diagrams. Unfortunately, this
     stuff packaged deep in [Fusion.Tagged_Coupling].
   \end{dubious} *)

(*i
    let strip_fuse' = function
      | Coupling.V3 (v, f, c) -> Coupling.V3 (v, Coupling.F12, c)
      | Coupling.V4 (v, f, c) -> Coupling.V4 (v, Coupling.F123, c)
      | Coupling.Vn (v, f, c) -> Coupling.Vn (v, [], c)

    let strip_fuse = function
      | Some c -> Some (strip_fuse' c)
      | None -> None
i*)

(* \begin{dubious}
     The [Tree.canonicalize] below should be necessary to remove
     topologically equivalent duplicates.
   \end{dubious} *)

(* Take a [CF.amplitude list] assumed to correspond to the same
   external states after stripping the color and return a
   pair of the list of external particles and the corresponding
   Feynman diagrams without color. *)

    let wf1 amplitude = 
      match F.externals amplitude with
      | wf :: _ -> wf
      | [] -> failwith "Omega.forest_sans_color: no external particles"

    let uniq l =
      ThoList.uniq (List.sort compare l)

    let forest_sans_color = function
      | amplitude :: _ as amplitudes ->
	let externals = F.externals amplitude in
	let prune_color wf =
	  (F.flavor_sans_color wf, F.momentum_list wf) in
	let prune_color_and_couplings (wf, c) =
	  (prune_color wf, None) in
	(List.map prune_color externals,
	 uniq
	   (List.map
	      (fun t ->
		Tree.canonicalize
		  (Tree.map prune_color_and_couplings prune_color t))
	      (ThoList.flatmap (fun a -> F.forest (wf1 a) a) amplitudes)))
      | [] -> ([], [])

    let dag_sans_color = function
      | amplitude :: _ as amplitudes ->
        let prune a = a in
        List.map prune amplitudes
      | [] -> []

    let p2s p =
      if p >= 0 && p <= 9 then
        string_of_int p
      else if p <= 36 then
        String.make 1 (Char.chr (Char.code 'A' + p - 10))
      else
        "_"

    let format_p wf =
      String.concat "" (List.map p2s (F.momentum_list wf))

    let variable wf =
      M.flavor_to_string (F.flavor_sans_color wf) ^ "[" ^ format_p wf ^ "]"

    let variable' wf =
      CM.flavor_to_TeX (F.flavor wf) ^ "(" ^ format_p wf ^ ")"

    let feynmf_style propagator color =
      { Tree.style =
          begin match propagator with
          | Coupling.Prop_Feynman
          | Coupling.Prop_Gauge _ ->
            begin match color with
            | Color.AdjSUN _ -> Some ("gluon", "")
            | _ -> Some ("boson", "")
            end
          | Coupling.Prop_Col_Feynman -> Some ("gluon", "")
          | Coupling.Prop_Unitarity
          | Coupling.Prop_Rxi _ -> Some ("dbl_wiggly", "")
          | Coupling.Prop_Spinor
          | Coupling.Prop_ConjSpinor -> Some ("fermion", "")
          | _ -> None
          end;
        Tree.rev =
          begin match propagator with
          | Coupling.Prop_Spinor -> true
          | Coupling.Prop_ConjSpinor -> false
          | _ -> false
          end;
        Tree.label = None;
        Tree.tension = None }

    let header incoming outgoing =
      "$ " ^
      String.concat " "
	(List.map (CM.flavor_to_TeX << F.flavor) incoming) ^
      " \\to " ^
      String.concat " "
	(List.map (CM.flavor_to_TeX << CM.conjugate << F.flavor) outgoing) ^
      " $"

    let header_sans_color incoming outgoing =
      "$ " ^
      String.concat " "
	(List.map (M.flavor_to_TeX << fst) incoming) ^
      " \\to " ^
      String.concat " "
	(List.map (M.flavor_to_TeX << M.conjugate << fst) outgoing) ^
      " $"
	
    let diagram incoming tree =
      let fmf wf =
	let f = F.flavor wf in
	feynmf_style (CM.propagator f) (CM.color f) in
      Tree.map
        (fun (n, _) ->
	  let n' = fmf n in
	  if List.mem n incoming then
            { n' with Tree.rev = not n'.Tree.rev }
	  else
            n')
        (fun l ->
	  if List.mem l incoming then
            l
	  else
            F.conjugate l)
	tree

    let diagram_sans_color incoming (tree) =
      let fmf (f, p) =
	feynmf_style (M.propagator f) (M.color f) in
      Tree.map
	(fun (n, c) ->
	  let n' = fmf n in
	  if List.mem n incoming then
	    { n' with Tree.rev = not n'.Tree.rev }
	  else
	    n')
	(fun (f, p) ->
	  if List.mem (f, p) incoming then
	    (f, p)
	  else
	    (M.conjugate f, p))
	tree

    let feynmf_set amplitude =
      match F.externals amplitude with
      | wf1 :: wf2 :: wfs ->
	let incoming = [wf1; wf2] in
    	{ Tree.header = header incoming wfs;
	  Tree.incoming = incoming;
	  Tree.diagrams =
	    List.map (diagram incoming) (F.forest wf1 amplitude) }
      | _ -> failwith "less than two external particles"

    let feynmf_set_sans_color (externals, trees) =
      match externals with
      | wf1 :: wf2 :: wfs ->
	let incoming = [wf1; wf2] in
	{ Tree.header = header_sans_color incoming wfs;
	  Tree.incoming = incoming;
	  Tree.diagrams =
	    List.map (diagram_sans_color incoming) trees }
      | _ -> failwith "less than two external particles"

    let feynmf_set_sans_color_empty (externals, trees) =
      match externals with
      | wf1 :: wf2 :: wfs ->
	let incoming = [wf1; wf2] in
	{ Tree.header = header_sans_color incoming wfs;
	  Tree.incoming = incoming;
	  Tree.diagrams = [] }
      | _ -> failwith "less than two external particles"

    let uncolored_colored amplitudes =
      { Tree.outer = feynmf_set_sans_color (forest_sans_color amplitudes);
	Tree.inner = List.map feynmf_set amplitudes }

    let uncolored_only amplitudes =
      { Tree.outer = feynmf_set_sans_color (forest_sans_color amplitudes);
	Tree.inner = [] }

    let colored_only amplitudes =
      { Tree.outer = feynmf_set_sans_color_empty (forest_sans_color amplitudes);
	Tree.inner = List.map feynmf_set amplitudes }

    let momentum_to_TeX (_, p) =
      String.concat "" (List.map p2s p)

    let wf_to_TeX (f, _ as wf) =
      M.flavor_to_TeX f ^ "(" ^ momentum_to_TeX wf ^ ")"

    let amplitudes_to_feynmf latex name amplitudes =
	Tree.feynmf_sets_wrapped latex name
	  wf_to_TeX momentum_to_TeX variable' format_p
	  (List.map uncolored_colored (amplitudes_by_flavor amplitudes))
	
    let amplitudes_to_feynmf_sans_color latex name amplitudes =
	Tree.feynmf_sets_wrapped latex name
	  wf_to_TeX momentum_to_TeX variable' format_p
	  (List.map uncolored_only (amplitudes_by_flavor amplitudes))

    let amplitudes_to_feynmf_color_only latex name amplitudes =
	Tree.feynmf_sets_wrapped latex name
	  wf_to_TeX momentum_to_TeX variable' format_p
	  (List.map colored_only (amplitudes_by_flavor amplitudes))

    let debug (str, descr, opt, var) =
      [ "-warning:" ^ str, Arg.Unit (fun () -> var := (opt, false):: !var),
        "         check " ^ descr ^ " and print warning on error";
        "-error:" ^ str, Arg.Unit (fun () -> var := (opt, true):: !var),
        "           check " ^ descr ^ " and terminate on error" ]

    let rec include_goldstones = function
      | [] -> false
      | (T.Gauge, _) :: _ -> true
      | _ :: rest -> include_goldstones rest
      
    let read_lines_rev file =
      let ic = open_in file in
      let rev_lines = ref [] in
      let rec slurp () =
        rev_lines := input_line ic :: !rev_lines;
        slurp () in
      try
        slurp ()
      with
      | End_of_file ->
          close_in ic;
          !rev_lines

    let read_lines file = 
      List.rev (read_lines_rev file)

(*i
    type cache_mode =
      | Cache_Default
      | Cache_Initialize of string

    let cache_option =
      ref Cache_Default
i*)

    let unphysical_polarization = ref None

(* \thocwmodulesection{Main Program} *)

    let main () =
      (* Delay evaluation of [M.external_flavors ()]! *)
      let usage () =
        "usage: " ^ Sys.argv.(0) ^
        " [options] [" ^
	  String.concat "|" (List.map M.flavor_to_string 
			       (ThoList.flatmap snd
				  (M.external_flavors ()))) ^ "]"
      and rev_scatterings = ref []
      and rev_decays = ref []
      and cascades = ref []
      and checks = ref []
      and output_file = ref None
      and print_forest = ref false
      and template = ref false
      and diagrams_all = ref None
      and diagrams_sans_color = ref None
      and diagrams_color_only = ref None
      and diagrams_LaTeX = ref false
      and quiet = ref false
      and write = ref true
      and params = ref false
      and poles = ref false
      and dag_out = ref None
      and dag0_out = ref None
      and phase_space_out = ref None in
      Options.parse
        (Options.cmdline "-target:" T.options @
         Options.cmdline "-model:" M.options @
         Options.cmdline "-fusion:" CF.options @
         ThoList.flatmap debug
           ["a", "arguments", T.All, checks;
            "n", "# of input arguments", T.Arguments, checks;
            "m", "input momenta", T.Momenta, checks;
            "g", "internal Ward identities", T.Gauge, checks] @
         [("-o", Arg.String (fun s -> output_file := Some s),
           "file             write to given file instead of /dev/stdout");
          ("-scatter",
           Arg.String (fun s -> rev_scatterings := s :: !rev_scatterings),
           "expr       in1 in2 -> out1 out2 ...");
          ("-scatter_file",
           Arg.String (fun s -> rev_scatterings := read_lines_rev s @ !rev_scatterings),
           "name  each line: in1 in2 -> out1 out2 ...");
          ("-decay", Arg.String (fun s -> rev_decays := s :: !rev_decays),
           "expr         in -> out1 out2 ...");
          ("-decay_file",
           Arg.String (fun s -> rev_decays := read_lines_rev s @ !rev_decays),
           "name    each line: in -> out1 out2 ...");
          ("-cascade", Arg.String (fun s -> cascades := s :: !cascades),
           "expr       select diagrams");
(*i
          ("-initialize",
           Arg.String (fun s -> cache_option := Cache_Initialize s),
           "dir     precompute lookup tables and store them in directory");
i*)
          ("-unphysical", Arg.Int (fun i -> unphysical_polarization := Some i),
           "n       use unphysical polarization for n-th particle / test WIs");
          ("-template", Arg.Set template,
           "          write a template for handwritten amplitudes");
          ("-forest", Arg.Set print_forest,
           "            Diagrammatic expansion");
          ("-diagrams", Arg.String (fun s -> diagrams_sans_color := Some s),
           "file      produce FeynMP output for Feynman diagrams");
          ("-diagrams:c", Arg.String (fun s -> diagrams_color_only := Some s),
           "file    produce FeynMP output for color flow diagrams");
          ("-diagrams:C", Arg.String (fun s -> diagrams_all := Some s),
           "file    produce FeynMP output for Feynman and color flow diagrams");
          ("-diagrams_LaTeX", Arg.Set diagrams_LaTeX,
           "    enclose FeynMP output in LaTeX wrapper");
          ("-quiet", Arg.Set quiet,
           "             don't print a summary");
          ("-summary", Arg.Clear write,
           "           print only a summary");
          ("-params", Arg.Set params,
           "            print the model parameters");
          ("-poles", Arg.Set poles,
           "             print the Monte Carlo poles");
          ("-dag", Arg.String (fun s -> dag_out := Some s),
           "               print minimal DAG");
          ("-full_dag", Arg.String (fun s -> dag0_out := Some s),
           "          print complete DAG");
          ("-phase_space", Arg.String (fun s -> phase_space_out := Some s),
           "       print minimal DAG for phase space")])
(*i       ("-T", Arg.Int Topology.Binary.debug_triplet, "");
          ("-P", Arg.Int Topology.Binary.debug_partition, "")])
i*)
        (fun _ -> prerr_endline (usage ()); exit 1)
        usage;

      let cmdline =
        String.concat " " (List.map ThoString.quote (Array.to_list Sys.argv)) in
        
      let output_channel, close_output_channel =
        match !output_file with
        | None ->
           (stdout, fun () -> ())
        | Some name ->
           let oc = open_out name in
           (oc, fun () -> close_out oc) in

      let processes =
        try
          ThoList.uniq
            (List.sort compare 
               (match List.rev !rev_scatterings, List.rev !rev_decays with
               | [], [] -> []
               | scatterings, [] ->
                   Proc.expand_scatterings (List.map Proc.parse_scattering scatterings)
               | [], decays ->
                   Proc.expand_decays (List.map Proc.parse_decay decays)
               | scatterings, decays -> 
                   invalid_arg "mixed scattering and decay!"))
        with
        | Invalid_argument s ->
            begin 
              Printf.eprintf "O'Mega: invalid process specification: %s!\n" s;
              flush stderr;
              []
            end in

(* \begin{dubious}
     This is still crude.  Eventually, we want to catch \emph{all} exceptions
     and write an empty (but compilable) amplitude unless one of the special options
     is selected.
   \end{dubious} *)

(*i
      begin match processes, !cache_option, !params with
      | [], Cache_Initialize dir, false ->
         (* [F.initialize_cache dir;] *)
          exit 0
      | _, _, true ->
          if !write then
            T.parameters_to_channel output_channel;
          exit 0
      | [], _, false ->
         if !write then
           T.amplitudes_to_channel cmdline output_channel !checks CF.empty;
         exit 0
      | _, _, false ->
i*)

      begin match processes, !params with
      | _, true ->
          if !write then
            T.parameters_to_channel output_channel;
          exit 0
      | [], false ->
         if !write then
           T.amplitudes_to_channel cmdline output_channel !checks CF.empty;
         exit 0
      | _, false ->

        let selectors =
          let fin, fout = List.hd processes in
          C.to_selectors (C.of_string_list (List.length fin + List.length fout) !cascades) in

        let amplitudes =
          try
            begin match F.check_charges () with
            | [] -> ()
            | violators ->
                let violator_strings =
                  String.concat ", "
                    (List.map
                       (fun flist ->
                         "(" ^ String.concat "," (List.map M.flavor_to_string flist) ^ ")")
                       violators) in
                failwith ("charge violating vertices: " ^ violator_strings)
            end;
            CF.amplitudes (include_goldstones !checks) !unphysical_polarization
	      CF.no_exclusions selectors processes
          with
          | Fusion.Majorana ->
             begin
               Printf.eprintf
                 "O'Mega: found Majorana fermions, switching representation!\n";
               flush stderr;
               close_output_channel ();
               Arg.current := 0;
               raise Fusion.Majorana
             end
          | exc ->
              begin 
                Printf.eprintf
                  "O'Mega: exception %s in amplitude construction!\n"
                  (Printexc.to_string exc);
                flush stderr;
                CF.empty;
              end in

        if !write then
          T.amplitudes_to_channel cmdline output_channel !checks amplitudes;

        if not !quiet then begin
          List.iter
            (fun amplitude ->
              Printf.eprintf "SUMMARY: %d fusions, %d propagators"
                (F.count_fusions amplitude) (F.count_propagators amplitude);
              flush stderr;
              Printf.eprintf ", %d diagrams" (F.count_diagrams amplitude);
              Printf.eprintf "\n")
            (CF.processes amplitudes);
          let couplings =
            List.fold_left
              (fun acc p ->
                let fusions = ThoList.flatmap F.rhs (F.fusions p)
                and brakets = ThoList.flatmap F.ket (F.brakets p) in
                let couplings =
                  VSet.of_list (List.map F.coupling (fusions @ brakets)) in
                VSet.union acc couplings)
              VSet.empty (CF.processes amplitudes) in
          Printf.eprintf "SUMMARY: %d vertices\n" (VSet.cardinal couplings);
          let ufo_couplings =
            VSet.fold
              (fun v acc ->
                match v with
                | Coupling.Vn (Coupling.UFO (_, v, _, _, _), _, _) ->
                   Sets.String.add v acc
                | _ -> acc)
              couplings Sets.String.empty in
          if not (Sets.String.is_empty ufo_couplings) then
            Printf.eprintf
              "SUMMARY: %d UFO vertices: %s\n"
              (Sets.String.cardinal ufo_couplings)
              (String.concat ", " (Sets.String.elements ufo_couplings))
        end;

        if !poles then begin
          List.iter
            (fun amplitude ->
              W.write output_channel "omega" (W.merge (W.trees amplitude)))
            (CF.processes amplitudes)
        end;

        begin match !dag0_out with
        | Some name ->
            let ch = open_out name in
            List.iter (F.tower_to_dot ch) (CF.processes amplitudes);
            close_out ch
        | None -> ()
        end;

        begin match !dag_out with
        | Some name ->
            let ch = open_out name in
            List.iter (F.amplitude_to_dot ch) (CF.processes amplitudes);
            close_out ch
        | None -> ()
        end;

        begin match !phase_space_out with
        | Some name ->
           let selectors =
             let fin, fout = List.hd processes in
             CT.to_selectors (CT.of_string_list (List.length fin + List.length fout) !cascades) in
           let ch = open_out name in
           begin try
             List.iter
               (fun (fin, fout) ->
                 Printf.fprintf
                   ch "%s -> %s ::\n"
                   (String.concat " " (List.map M.flavor_to_string fin))
                   (String.concat " " (List.map M.flavor_to_string fout));
                 match fin with
                 | [] ->
                    failwith "Omega(): phase space: no incoming particles"
                 | [f] ->
                    PHS.phase_space_channels
                      ch
                      (PHS.amplitude_sans_color
                         false PHS.no_exclusions selectors fin fout)
                 | [f1; f2] ->
                    PHS.phase_space_channels
                      ch
                      (PHS.amplitude_sans_color
                         false PHS.no_exclusions selectors fin fout);
                    PHS.phase_space_channels_flipped
                      ch
                      (PHS.amplitude_sans_color
                         false PHS.no_exclusions selectors [f2; f1] fout)
                 | _ ->
                    failwith "Omega(): phase space: 3 or more incoming particles")
               processes;
             close_out ch
           with
           | exc ->
              begin
                close_out ch;
                Printf.eprintf
                  "O'Mega: exception %s in phase space construction!\n"
                  (Printexc.to_string exc);
                flush stderr
              end
           end
        | None -> ()
        end;

        if !print_forest then
          List.iter
            (fun amplitude ->
              List.iter (fun t -> Printf.eprintf "%s\n"
                  (Tree.to_string
                     (Tree.map (fun (wf, _) -> variable wf) (fun _ -> "") t)))
                (F.forest (List.hd (F.externals amplitude)) amplitude))
            (CF.processes amplitudes);

        begin match !diagrams_all with
        | Some name ->
	  amplitudes_to_feynmf !diagrams_LaTeX name amplitudes
        | None -> ()
        end;

        begin match !diagrams_sans_color with
        | Some name ->
	  amplitudes_to_feynmf_sans_color !diagrams_LaTeX name amplitudes
        | None -> ()
        end;

        begin match !diagrams_color_only with
        | Some name ->
	  amplitudes_to_feynmf_color_only !diagrams_LaTeX name amplitudes
        | None -> ()
        end;

        close_output_channel ();

        exit 0

      end

(* \begin{dubious}
     This was only intended for debugging O'Giga \ldots
   \end{dubious} *)

    let decode wf =
      (F.flavor wf, (F.momentum wf : Momentum.Default.t))

    let diagrams in1 in2 out =
      match F.amplitudes false F.no_exclusions C.no_cascades [in1; in2] out with
      | a :: _ ->
          let wf1 = List.hd (F.externals a)
          and wf2 = List.hd (List.tl (F.externals a)) in
          let wf2 = decode wf2 in
          List.map (fun t ->
            (wf2,
             Tree.map (fun (wf, _) -> decode wf) decode t))
            (F.forest wf1 a)
      | [] -> []

    let diagrams in1 in2 out =
      failwith "Omega().diagrams: disabled"

  end

module Binary (TM : Target.Maker) (M : Model.T) =
  Make(Fusion.Binary)(Fusion.Helac_Binary)(TM)(M)
module Binary_Majorana (TM : Target.Maker) (M : Model.T) =
  Make(Fusion.Binary_Majorana)(Fusion.Helac_Binary_Majorana)(TM)(M)
module Mixed23 (TM : Target.Maker) (M : Model.T) =
  Make(Fusion.Mixed23)(Fusion.Helac_Mixed23)(TM)(M)
module Mixed23_Majorana (TM : Target.Maker) (M : Model.T) =
  Make(Fusion.Mixed23_Majorana)(Fusion.Helac_Mixed23_Majorana)(TM)(M)
module Mixed23_Majorana_vintage (TM : Target.Maker) (M : Model.T) =
  Make(Fusion_vintage.Mixed23_Majorana)(Fusion.Helac_Mixed23_Majorana)(TM)(M)

module Bound (M : Model.T) : Tuple.Bound =
  struct
    (* \begin{dubious}
         Above [max_degree = 6], the performance drops \emph{dramatically}!
       \end{dubious} *)
    let max_arity () =
      pred (M.max_degree ())
  end

module Nary (TM : Target.Maker) (M : Model.T) =
  Make(Fusion.Nary(Bound(M)))(Fusion.Helac(Bound(M)))(TM)(M)
module Nary_Majorana (TM : Target.Maker) (M : Model.T) =
  Make(Fusion.Nary_Majorana(Bound(M)))(Fusion.Helac_Majorana(Bound(M)))(TM)(M)

