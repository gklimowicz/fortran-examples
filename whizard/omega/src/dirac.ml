(* Dirac.ml --

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

(* \thocwmodulesection{Dirac $\gamma$-matrices} *)

module type T =
  sig
    type qc = Algebra.QC.t
    type t = qc array array
    val zero : qc
    val one : qc
    val minus_one : qc
    val i : qc
    val minus_i : qc
    val unit : t
    val null : t
    val gamma0 : t
    val gamma1 : t
    val gamma2 : t
    val gamma3 : t
    val gamma5 : t
    val gamma : t array
    val cc : t
    val neg : t -> t
    val add : t -> t -> t
    val sub : t -> t -> t
    val mul : t -> t -> t
    val times : qc -> t -> t
    val transpose : t -> t
    val adjoint : t -> t
    val conj : t -> t
    val product : t list -> t
    val pp : Format.formatter -> t -> unit
    val test_suite : OUnit.test
  end

(* \thocwmodulesubsection{Matrices with complex rational entries} *)

module Q = Algebra.Q
module QC = Algebra.QC

type complex_rational = QC.t

let zero = QC.null
let one = QC.unit
let minus_one = QC.neg one
let i = QC.make Q.null Q.unit
let minus_i = QC.conj i

type matrix = complex_rational array array

(* \thocwmodulesubsection{Dirac $\gamma$-matrices} *)

module type R =
  sig
    type qc = complex_rational
    type t = matrix
    val gamma0 : t
    val gamma1 : t
    val gamma2 : t
    val gamma3 : t
    val gamma5 : t
    val cc : t
    val cc_is_i_gamma2_gamma_0 : bool
  end

module Make (R : R) : T =
  struct

    type qc = complex_rational
    type t = matrix

    let zero = zero
    let one = one
    let minus_one = minus_one
    let i = i
    let minus_i = minus_i

    let null =
      [| [| zero; zero; zero; zero |];
         [| zero; zero; zero; zero |];
         [| zero; zero; zero; zero |];
         [| zero; zero; zero; zero |] |]

    let unit =
      [| [| one;  zero; zero; zero |];
         [| zero; one;  zero; zero |];
         [| zero; zero; one;  zero |];
         [| zero; zero; zero; one  |] |]

    let gamma0 = R.gamma0
    let gamma1 = R.gamma1
    let gamma2 = R.gamma2
    let gamma3 = R.gamma3
    let gamma5 = R.gamma5
    let gamma = [| gamma0; gamma1; gamma2; gamma3 |]
    let cc = R.cc

    let neg g =
      let g' = Array.make_matrix 4 4 zero in
      for i = 0 to 3 do
        for j = 0 to 3 do
          g'.(i).(j) <- QC.neg g.(i).(j)
        done
      done;
      g'

    let add g1 g2 =
      let g12 = Array.make_matrix 4 4 zero in
      for i = 0 to 3 do
        for j = 0 to 3 do
          g12.(i).(j) <- QC.add g1.(i).(j) g2.(i).(j)
        done
      done;
      g12

    let sub g1 g2 =
      let g12 = Array.make_matrix 4 4 zero in
      for i = 0 to 3 do
        for j = 0 to 3 do
          g12.(i).(j) <- QC.sub g1.(i).(j) g2.(i).(j)
        done
      done;
      g12

    let mul g1 g2 =
      let g12 = Array.make_matrix 4 4 zero in
      for i = 0 to 3 do
        for k = 0 to 3 do
          for j = 0 to 3 do
            g12.(i).(k) <- QC.add g12.(i).(k) (QC.mul g1.(i).(j) g2.(j).(k))
          done
        done
      done;
      g12

    let times q g =
      let g' = Array.make_matrix 4 4 zero in
      for i = 0 to 3 do
        for j = 0 to 3 do
          g'.(i).(j) <- QC.mul q g.(i).(j)
        done
      done;
      g'

    let transpose g =
      let g' = Array.make_matrix 4 4 zero in
      for i = 0 to 3 do
        for j = 0 to 3 do
          g'.(i).(j) <- g.(j).(i)
        done
      done;
      g'

    let adjoint g =
      let g' = Array.make_matrix 4 4 zero in
      for i = 0 to 3 do
        for j = 0 to 3 do
          g'.(i).(j) <- QC.conj g.(j).(i)
        done
      done;
      g'

    let conj g =
      let g' = Array.make_matrix 4 4 zero in
      for i = 0 to 3 do
        for j = 0 to 3 do
          g'.(i).(j) <- QC.conj g.(i).(j)
        done
      done;
      g'

    let product glist =
      List.fold_right mul glist unit

    let pp fmt g =
      let pp_row i =
        for j = 0 to 3 do
          Format.fprintf fmt " %8s" (QC.to_string g.(i).(j))
        done in
      Format.fprintf fmt "\n /";
      pp_row 0;
      Format.fprintf fmt " \\\n";
      for i = 1 to 2 do
        Format.fprintf fmt " |";
        pp_row i;
        Format.fprintf fmt " |\n"
      done;
      Format.fprintf fmt " \\";
      pp_row 3;
      Format.fprintf fmt " /\n"

    open OUnit

    let two = QC.make (Q.make 2 1) Q.null
    let half = QC.make (Q.make 1 2) Q.null
    let two_unit = times two unit

    let ac_lhs mu nu =
      add (mul gamma.(mu) gamma.(nu)) (mul gamma.(nu) gamma.(mu))

    let ac_rhs mu nu =
      if mu = nu then
        if mu = 0 then
          two_unit
        else
          neg two_unit
      else
        null

    let test_ac mu nu =
      (ac_lhs mu nu) = (ac_rhs mu nu)

    let ac_lhs_all =
      let lhs = Array.make_matrix 4 4 null in
      for mu = 0 to 3 do
        for nu = 0 to 3 do
          lhs.(mu).(nu) <- ac_lhs mu nu
        done
      done;
      lhs
                                                                   
    let ac_rhs_all =
      let rhs = Array.make_matrix 4 4 null in
      for mu = 0 to 3 do
        for nu = 0 to 3 do
          rhs.(mu).(nu) <- ac_rhs mu nu
        done
      done;
      rhs

    let dump2 lhs rhs =
      for i = 0 to 3 do
        for j = 0 to 3 do
          Printf.printf
            "   i = %d, j =%d: %s + %s*I | %s + %s*I\n"
            i j
            (Q.to_string (QC.real lhs.(i).(j)))
            (Q.to_string (QC.imag lhs.(i).(j)))
            (Q.to_string (QC.real rhs.(i).(j)))
            (Q.to_string (QC.imag rhs.(i).(j)))
        done
      done

    let dump2_all lhs rhs =
      for mu = 0 to 3 do
        for nu = 0 to 3 do
          Printf.printf "mu = %d, nu =%d: \n" mu nu;
          dump2 lhs.(mu).(nu) rhs.(mu).(nu)
        done
      done

    let anticommute =
      "anticommutation relations" >::
        (fun () ->
          assert_bool
            ""
            (if ac_lhs_all = ac_rhs_all then
               true
             else
               begin
                 dump2_all ac_lhs_all ac_rhs_all;
                 false
               end))

    let equal_or_dump2 lhs rhs =
      if lhs = rhs then
        true
      else
        begin
          dump2 lhs rhs;
          false
        end

    let gamma5_def =
      "gamma5" >::
        (fun () ->
          assert_bool
            "definition"
            (equal_or_dump2
               gamma5
               (times i (product [gamma0; gamma1; gamma2; gamma3]))))

    let self_adjoint =
      "(anti)selfadjointness" >:::
        [ "gamma0" >::
            (fun () ->
              assert_bool "self" (equal_or_dump2 gamma0 (adjoint gamma0)));
          "gamma1" >::
            (fun () ->
              assert_bool "anti" (equal_or_dump2 gamma1 (neg (adjoint gamma1))));
          "gamma2" >::
            (fun () ->
              assert_bool "anti" (equal_or_dump2 gamma2 (neg (adjoint gamma2))));
          "gamma3" >::
            (fun () ->
              assert_bool "anti" (equal_or_dump2 gamma3 (neg (adjoint gamma3))));
          "gamma5" >::
            (fun () ->
              assert_bool "self" (equal_or_dump2 gamma5 (adjoint gamma5))) ]

    (* $C^2=-\mathbf{1}$ is \emph{not} true in all realizations, but
       we assume it at several points in [UFO_Lorentz].  Therefore we
       must test it here for all realizations that are implemented. *)
    let cc_inv = neg cc

    (* Verify that $\Gamma^T= - C\Gamma C^{-1}$ using the actual
       matrix transpose: *)
    let cc_gamma g =
      equal_or_dump2 (neg (transpose g)) (product [cc; g; cc_inv])

    (* Of course, $C=\ii\gamma^2\gamma^0$ is also not true in \emph{all}
       realizations.  But it is true in the chiral representation
       used here and we can test it. *)
    let charge_conjugation =
      "charge conjugation" >:::
        [ "inverse" >::
            (fun () ->
              assert_bool "" (equal_or_dump2 (mul cc cc_inv) unit));

          "gamma0" >:: (fun () -> assert_bool "" (cc_gamma gamma0));
          "gamma1" >:: (fun () -> assert_bool "" (cc_gamma gamma1));
          "gamma2" >:: (fun () -> assert_bool "" (cc_gamma gamma2));
          "gamma3" >:: (fun () -> assert_bool "" (cc_gamma gamma3));

          "gamma5" >::
            (fun () ->
              assert_bool "" (equal_or_dump2 (transpose gamma5)
                                             (product [cc; gamma5; cc_inv])));
          "=i*g2*g0" >::
            (fun () ->
              skip_if (not R.cc_is_i_gamma2_gamma_0)
                "representation dependence";
              assert_bool "" (equal_or_dump2 cc (times i (mul gamma2 gamma0))))
        ]

    let test_suite =
      "Dirac Matrices" >:::
        [anticommute;
         gamma5_def;
         self_adjoint;
         charge_conjugation]

  end

module Chiral_R : R =
  struct

    type qc = complex_rational
    type t = matrix

    let gamma0 =
      [| [| zero; zero; one;  zero |];
         [| zero; zero; zero; one  |];
         [| one;  zero; zero; zero |];
         [| zero; one;  zero; zero |] |]

    let gamma1 =
      [| [| zero;      zero;      zero; one  |];
         [| zero;      zero;      one;  zero |];
         [| zero;      minus_one; zero; zero |];
         [| minus_one; zero;      zero; zero |] |]

    let gamma2 =
      [| [| zero;    zero; zero; minus_i |];
         [| zero;    zero; i;    zero    |];
         [| zero;    i;    zero; zero    |];
         [| minus_i; zero; zero; zero    |] |]

    let gamma3 =
      [| [| zero;      zero; one;  zero      |];
         [| zero;      zero; zero; minus_one |];
         [| minus_one; zero; zero; zero      |];
         [| zero;      one;  zero; zero      |] |]

    let gamma5 =
      [| [| minus_one; zero;      zero; zero |];
         [| zero;      minus_one; zero; zero |];
         [| zero;      zero;      one;  zero |];
         [| zero;      zero;      zero; one  |] |]

    let cc =
      [| [| zero;      one;  zero; zero      |];
         [| minus_one; zero; zero; zero      |];
         [| zero;      zero; zero; minus_one |];
         [| zero;      zero; one;  zero      |] |]

    let cc_is_i_gamma2_gamma_0 = true

  end

module Dirac_R : R =
  struct

    type qc = complex_rational
    type t = matrix

    let gamma0 =
      [| [| one;  zero; zero;      zero |];
         [| zero; one;  zero;      zero  |];
         [| zero; zero; minus_one; zero |];
         [| zero; zero; zero;      minus_one |] |]

    let gamma1 = Chiral_R.gamma1
    let gamma2 = Chiral_R.gamma2
    let gamma3 = Chiral_R.gamma3

    let gamma5 =
      [| [| zero; zero; one;  zero |];
         [| zero; zero; zero; one  |];
         [| one;  zero; zero; zero |];
         [| zero; one;  zero; zero |] |]

    let cc =
      [| [| zero; zero;      zero; minus_one  |];
         [| zero; zero;      one;  zero       |];
         [| zero; minus_one; zero; zero       |];
         [| one;  zero;      zero; zero       |] |]

    let cc_is_i_gamma2_gamma_0 = true

  end

module Majorana_R : R =
  struct

    type qc = complex_rational
    type t = matrix

    let gamma0 =
      [| [| zero; zero;    zero; minus_i |];
         [| zero; zero;    i;    zero    |];
         [| zero; minus_i; zero; zero    |];
         [| i;    zero;    zero; zero    |] |]

    let gamma1 =
      [| [| i;    zero;    zero; zero    |];
         [| zero; minus_i; zero; zero    |];
         [| zero; zero;    i;    zero    |];
         [| zero; zero;    zero; minus_i |] |]

    let gamma2 =
      [| [| zero; zero;    zero;    i    |];
         [| zero; zero;    minus_i; zero |];
         [| zero; minus_i; zero;    zero |];
         [| i;    zero;    zero;    zero |] |]

    let gamma3 =
      [| [| zero;    minus_i; zero;    zero    |];
         [| minus_i; zero;    zero;    zero    |];
         [| zero;    zero;    zero;    minus_i |];
         [| zero;    zero;    minus_i; zero    |] |]

    let gamma5 =
      [| [| zero; minus_i; zero;    zero |];
         [| i;    zero;    zero;    zero |];
         [| zero; zero;    zero;    i    |];
         [| zero; zero;    minus_i; zero |] |]

    let cc =
      [| [| zero; zero;      zero; minus_one  |];
         [| zero; zero;      one;  zero       |];
         [| zero; minus_one; zero; zero       |];
         [| one;  zero;      zero; zero       |] |]

    let cc_is_i_gamma2_gamma_0 = false

  end

module Chiral = Make (Chiral_R)
module Dirac = Make (Dirac_R)
module Majorana = Make (Majorana_R)
