(* momentum.ml --

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

module type T =
  sig
    type t
    val of_ints : int -> int list -> t
    exception Duplicate of int
    exception Range of int
    exception Mismatch of string * t * t
    exception Negative
    val to_ints : t -> int list
    val dim : t -> int
    val rank : t -> int
    val singleton : int -> int -> t
    val zero : int -> t
    val compare : t -> t -> int
    val neg : t -> t
    val abs : t -> t
    val add : t -> t -> t
    val sub : t -> t -> t
    val try_add : t -> t -> t option
    val try_sub : t -> t -> t option
    val less : t -> t -> bool
    val lesseq : t -> t -> bool
    val try_fusion : t -> t -> t -> (bool * bool) option
    val to_string : t -> string
    val split : int -> int -> t -> t
    module Scattering :
        sig    
          val incoming : t -> bool
          val outgoing : t -> bool
          val timelike : t -> bool
          val spacelike : t -> bool
          val s_channel_in : t -> bool
          val s_channel_out : t -> bool
          val s_channel : t -> bool
          val flip_s_channel_in : t -> t
        end
    module Decay :
        sig
          val incoming : t -> bool
          val outgoing : t -> bool
          val timelike : t -> bool
          val spacelike : t -> bool
        end
  end

(* \thocwmodulesection{Lists of Integers} *)

(* The first implementation (as part of [Fusion]) was based on sorted
   lists, because I did not want to preclude the use of more general
   indices that integers.  However, there's probably not much use for
   this generality (the indices are typically generated automatically
   and integer are the most natural choice) and it is no longer supported.
   by the current signature.   Thus one can also use the
   more efficient implementation based on bitvectors below. *)

module Lists =
  struct

    type t = { d : int; r : int; p : int list }

    exception Range of int 
    exception Duplicate of int 

    let rec check d = function
      | p1 :: p2 :: _ when p2 <= p1 -> raise (Duplicate p1)
      | p1 :: (p2 :: _ as rest) -> check d rest
      | [p] when p < 1 || p > d -> raise (Range p)
      | [p] -> ()
      | [] -> ()

    let of_ints d p =
      let p' = List.sort compare p in
      check d p';
      { d = d; r = List.length p; p = p' }

    let to_ints p = p.p
    let dim p = p.d
    let rank p = p.r
    let zero d = { d = d; r = 0; p = [] }
    let singleton d p = { d = d; r = 1; p = [p] }

    let to_string p =
      "[" ^ String.concat "," (List.map string_of_int p.p) ^
      "/" ^ string_of_int p.r ^ "/" ^ string_of_int p.d ^ "]"

    exception Mismatch of string * t * t
    let mismatch s p1 p2 = raise (Mismatch (s, p1, p2))

    let matching f s p1 p2 =
      if p1.d = p2.d then
        f p1 p2
      else
        mismatch s p1 p2
     
    let compare p1 p2 =
      if p1.d = p2.d then begin
        let c = compare p1.r p2.r in
        if c <> 0 then
          c
        else
          compare p1.p p2.p
      end else
        mismatch "compare" p1 p2

    let rec neg' d i = function
      | [] ->
          if i <= d then
            i :: neg' d (succ i) []
          else
            []
      | i' :: rest as p ->
          if i' > d then
            failwith "Integer_List.neg: internal error"
          else if i' = i then
            neg' d (succ i) rest
          else
            i :: neg' d (succ i) p

    let neg p = { d = p.d; r = p.d - p.r; p = neg' p.d 1 p.p }

    let abs p =
      if 2 * p.r > p.d then
        neg p
      else
        p

    let rec add' p1 p2 =
      match p1, p2 with
      | [], p -> p
      | p, [] -> p
      | x1 :: p1', x2 :: p2' ->
          if x1 < x2 then
            x1 :: add' p1' p2
          else if x2 < x1 then
            x2 :: add' p1 p2'
          else
            raise (Duplicate x1)

    let add p1 p2 =
      if p1.d = p2.d then
        { d = p1.d; r = p1.r + p2.r; p = add' p1.p p2.p }
      else
        mismatch "add" p1 p2

    let rec try_add' d r acc p1 p2 =
      match p1, p2 with
      | [], p -> Some ({ d = d; r = r; p = List.rev_append acc p })
      | p, [] -> Some ({ d = d; r = r; p = List.rev_append acc p })
      | x1 :: p1', x2 :: p2' ->
          if x1 < x2 then
            try_add' d r (x1 :: acc) p1' p2
          else if x2 < x1 then
            try_add' d r (x2 :: acc) p1 p2'
          else
            None

    let try_add p1 p2 =
      if p1.d = p2.d then
        try_add' p1.d (p1.r + p2.r) [] p1.p p2.p
      else 
        mismatch "try_add" p1 p2

    exception Negative

    let rec sub' p1 p2 =
      match p1, p2 with
      | p, [] -> p
      | [], _ -> raise Negative
      | x1 :: p1', x2 :: p2' ->
          if x1 < x2 then
            x1 :: sub' p1' p2
          else if x1 = x2 then
            sub' p1' p2'
          else
            raise Negative

    let rec sub p1 p2 =
      if p1.d = p2.d then begin
        if p1.r >= p2.r then
          { d = p1.d; r = p1.r - p2.r; p = sub' p1.p p2.p }
        else
          neg (sub p2 p1)
       end else
        mismatch "sub" p1 p2

    let rec try_sub' d r acc p1 p2 =
      match p1, p2 with
      | p, [] -> Some ({ d = d; r = r; p = List.rev_append acc p })
      | [], _ -> None
      | x1 :: p1', x2 :: p2' ->
          if x1 < x2 then
            try_sub' d r (x1 :: acc) p1' p2
          else if x1 = x2 then
            try_sub' d r acc p1' p2'
          else
            None

    let try_sub p1 p2 =
      if p1.d = p2.d then begin
        if p1.r >= p2.r then
          try_sub' p1.d (p1.r - p2.r) [] p1.p p2.p
        else
          match try_sub' p1.d (p2.r - p1.r) [] p2.p p1.p with
          | None -> None
          | Some p -> Some (neg p)
      end else
        mismatch "try_sub" p1 p2

    let rec less' equal p1 p2 =
      match p1, p2 with
      | [], [] -> not equal
      | [], _ -> true
      | x1 :: _ , [] -> false
      | x1 :: p1', x2 :: p2' when x1 = x2 -> less' equal p1' p2'
      | x1 :: p1', x2 :: p2' -> less' false p1 p2'

    let less p1 p2 =
      if p1.d = p2.d then
        less' true p1.p p2.p
      else
        mismatch "sub" p1 p2

    let rec lesseq' p1 p2 =
      match p1, p2 with
      | [], _ -> true
      | x1 :: _ , [] -> false
      | x1 :: p1', x2 :: p2' when x1 = x2 -> lesseq' p1' p2'
      | x1 :: p1', x2 :: p2' -> lesseq' p1 p2'
            
    let lesseq p1 p2 =
      if p1.d = p2.d then
        lesseq' p1.p p2.p
      else
        mismatch "lesseq" p1 p2

    module Scattering =
      struct

        let incoming p =
          if p.r = 1 then
            match p.p with
            | [1] | [2] -> true
            | _ -> false
          else
            false

        let outgoing p =
          if p.r = 1 then
            match p.p with
            | [1] | [2] -> false
            | _ -> true
          else
            false

        let s_channel_in p =
          match p.p with
          | [1; 2] -> true
          | _ -> false

        let rec s_channel_out' d i = function
          | [] -> i = succ d
          | i' :: p when i' = i -> s_channel_out' d (succ i) p
          | _ -> false

        let s_channel_out p =
          match p.p with
          | 3 :: p' -> s_channel_out' p.d 4 p'
          | _ -> false

        let s_channel p = s_channel_in p || s_channel_out p

        let timelike p =
          match p.p with
          | p1 :: p2 :: _ -> p1 > 2 || (p1 = 1 && p2 = 2)
          | p1 :: _ -> p1 > 2
          | [] -> false

        let spacelike p = not (timelike p)

        let flip_s_channel_in p =
          if s_channel_in p then
            neg (of_ints p.d [1;2])
          else
            p

      end

    module Decay =
      struct

        let incoming p =
          if p.r = 1 then
            match p.p with
            | [1] -> true
            | _ -> false
          else
            false

        let outgoing p =
          if p.r = 1 then
            match p.p with
            | [1] -> false
            | _ -> true
          else
            false

        let timelike p =
          match p.p with
          | [1] -> true
          | p1 :: _ -> p1 > 1
          | [] -> false

        let spacelike p = not (timelike p)

      end

    let test_sum p inv1 p1 inv2 p2 =
      if p.d = p1.d then begin
        if p.d = p2.d then begin
          match (if inv1 then try_add else try_sub) p p1 with
          | None -> false
          | Some p' ->
              begin match (if inv2 then try_add else try_sub) p' p2 with
              | None -> false
              | Some p'' -> p''.r = 0 || p''.r = p.d
              end
        end else
          mismatch "test_sum" p p2
      end else
        mismatch "test_sum" p p1

    let try_fusion p p1 p2 =
      if test_sum p false p1 false p2 then
        Some (false, false)
      else if test_sum p true p1 false p2 then
        Some (true, false)
      else if test_sum p false p1 true p2 then
        Some (false, true)
      else if test_sum p true p1 true p2 then
        Some (true, true)
      else
        None
          
    let split i n p =
      let n' = n - 1 in
      let rec split' head = function
        | [] -> (p.r, List.rev head)
        | i1 :: ilist ->
            if i1 < i then
              split' (i1 :: head) ilist
            else if i1 > i then
              (p.r, List.rev_append head (List.map ((+) n') (i1 :: ilist)))
            else
              (p.r + n',
               List.rev_append head
                 ((ThoList.range i1 (i1 + n')) @ (List.map ((+) n') ilist))) in
      let r', p' = split' [] p.p in
      { d = p.d + n'; r = r'; p = p' }

  end

(* \thocwmodulesection{Bit Fiddlings} *)

(* Bit vectors are popular in Fortran based
   implementations~\cite{ALPHA:1997,HELAC:2000,Kilian:WHIZARD} and
   can be more efficient. In particular, when all infomation is
   packed into a single integer, much of the memory overhead is
   reduced. *)

module Bits =
  struct

    type t = int

(* Bits $1\ldots21$ are used as a bitvector, indicating whether a
   particular momentum is included. Bits $22\ldots26$ represent the
   numbers of bits set in bits $1\ldots21$ and bits $27\ldots31$
   denote the maximum number of momenta.  *)
    let mask n = (1 lsl n) - 1
    let mask2 = mask 2
    let mask5 = mask 5
    let mask21 = mask 21

    let maskd = mask5 lsl 26
    let maskr = mask5 lsl 21
    let maskb = mask21

    let dim0 p = p land maskd
    let rank0 p = p land maskr
    let bits0 p = p land maskb

    let dim p = (dim0 p) lsr 26
    let rank p = (rank0 p) lsr 21
    let bits p = bits0 p

    let drb0 d r b = d lor r lor b
    let drb d r b = d lsl 26 lor r lsl 21 lor b

(* For a 64-bit architecture, the corresponding sizes could
   be increased to $1\ldots51$, $52\ldots57$, and $58\ldots63$.
   However, the combinatorical complexity will have killed
   us long before we can reach these values. *)

    exception Range of int 
    exception Duplicate of int 

    exception Mismatch of string * t * t
    let mismatch s p1 p2 = raise (Mismatch (s, p1, p2))

    let of_ints d p =
      let r = List.length p in
      if d <= 21 && r <= 21 then begin
        List.fold_left (fun b p' ->
          if p' <= d then
            b lor (1 lsl (pred p'))
          else
            raise (Range p')) (drb d r 0) p
      end else
        raise (Range r)

    let zero d = drb d 0 0

    let singleton d p = drb d 1 (1 lsl (pred p))

    let rec to_ints' acc p b =
      if b = 0 then
        List.rev acc
      else if (b land 1) = 1 then
        to_ints' (p :: acc) (succ p) (b lsr 1)
      else
        to_ints' acc (succ p) (b lsr 1)

    let to_ints p = to_ints' [] 1 (bits p)

    let to_string p =
      "[" ^ String.concat "," (List.map string_of_int (to_ints p)) ^
      "/" ^ string_of_int (rank p) ^ "/" ^ string_of_int (dim p) ^ "]"

    let compare p1 p2 =
      if dim0 p1 = dim0 p2 then begin
        let c = compare (rank0 p1) (rank0 p2) in
        if c <> 0 then
          c
        else
          compare (bits p1) (bits p2)
      end else
        mismatch "compare" p1 p2

    let neg p =
      let d = dim p and r = rank p in
      drb d (d - r) ((mask d) land (lnot p))

    let abs p =
      if 2 * (rank p) > dim p then
        neg p
      else
        p

    let add p1 p2 =
      let d1 = dim0 p1 and d2 = dim0 p2 in
      if d1 = d2 then begin
        let b1 = bits p1 and b2 = bits p2 in
        if b1 land b2 = 0 then
          drb0 d1 (rank0 p1 + rank0 p2) (b1 lor b2)
        else
          raise (Duplicate 0)
      end else
        mismatch "add" p1 p2

    exception Negative

    let rec sub p1 p2 =
      let d1 = dim0 p1 and d2 = dim0 p2 in
      if d1 = d2 then begin
        let r1 = rank0 p1 and r2 = rank0 p2 in
        if r1 >= r2 then begin
          let b1 = bits p1 and b2 = bits p2 in
          if b1 lor b2 = b1 then
            drb0 d1 (r1 - r2) (b1 lxor b2)
          else
            raise Negative
        end else
          neg (sub p2 p1)
      end else
        mismatch "sub" p1 p2

    let try_add p1 p2 =
      let d1 = dim0 p1 and d2 = dim0 p2 in
      if d1 = d2 then begin
        let b1 = bits p1 and b2 = bits p2 in
        if b1 land b2 = 0 then
          Some (drb0 d1 (rank0 p1 + rank0 p2) (b1 lor b2))
        else
          None
      end else
        mismatch "try_add" p1 p2

    let rec try_sub p1 p2 =
      let d1 = dim0 p1 and d2 = dim0 p2 in
      if d1 = d2 then begin
        let r1 = rank0 p1 and r2 = rank0 p2 in
        if r1 >= r2 then begin
          let b1 = bits p1 and b2 = bits p2 in
          if b1 lor b2 = b1 then
            Some (drb0 d1 (r1 - r2) (b1 lxor b2))
          else
            None
        end else
          begin match try_sub p2 p1 with
          | Some p -> Some (neg p)
          | None -> None
          end
      end else
        mismatch "sub" p1 p2

    let lesseq p1 p2 = 
      let d1 = dim0 p1 and d2 = dim0 p2 in
      if d1 = d2 then begin
        let r1 = rank0 p1 and r2 = rank0 p2 in
        if r1 <= r2 then begin
          let b1 = bits p1 and b2 = bits p2 in
          b1 lor b2 = b2
        end else
          false
      end else
        mismatch "less" p1 p2

    let less p1 p2 = p1 <> p2 && lesseq p1 p2

    let mask_in1 = 1
    let mask_in2 = 2
    let mask_in = mask_in1 lor mask_in2

    module Scattering =
      struct

        let incoming p =
          rank p = 1 && (mask_in land p <> 0)

        let outgoing p =
          rank p = 1 && (mask_in land p = 0)

        let timelike p =
          (rank p > 0 && (mask_in land p = 0)) || (bits p = mask_in)

        let spacelike p = 
          (rank p > 0) && not (timelike p)

        let s_channel_in p =
          bits p = mask_in

        let s_channel_out p = 
          rank p > 0 && (mask_in lxor p = 0)

        let s_channel p =
          s_channel_in p || s_channel_out p

        let flip_s_channel_in p =
          if s_channel_in p then
            neg p
          else
            p

      end

    module Decay =
      struct

        let incoming p =
          rank p = 1 && (mask_in1 land p = mask_in1)

        let outgoing p =
          rank p = 1 && (mask_in1 land p = 0)

        let timelike p =
          incoming p || (rank p > 0 && mask_in1 land p = 0)

        let spacelike p =
          not (timelike p)

      end

    let test_sum p inv1 p1 inv2 p2 =
      let d = dim p in
      if d = dim p1 then begin
        if d = dim p2 then begin
          match (if inv1 then try_add else try_sub) p p1 with
          | None -> false
          | Some p' ->
              begin match (if inv2 then try_add else try_sub) p' p2 with
              | None -> false
              | Some p'' ->
                  let r = rank p'' in
                  r = 0 || r = d
              end
        end else
          mismatch "test_sum" p p2
      end else
        mismatch "test_sum" p p1

    let try_fusion p p1 p2 =
      if test_sum p false p1 false p2 then
        Some (false, false)
      else if test_sum p true p1 false p2 then
        Some (true, false)
      else if test_sum p false p1 true p2 then
        Some (false, true)
      else if test_sum p true p1 true p2 then
        Some (true, true)
      else
        None

(* First create a gap of size~$n-1$ and subsequently fill it if and only if
   the bit~$i$ was set.  *)
    let split i n p =
      let delta_d = n - 1
      and b = bits p in
      let mask_low = mask (pred i)
      and mask_i = 1 lsl (pred i)
      and mask_high = lnot (mask i) in
      let b_low = mask_low land b
      and b_med, delta_r =
        if mask_i land b <> 0 then
          ((mask n) lsl (pred i), delta_d)
        else
          (0, 0)
      and b_high =
        if delta_d > 0 then
          (mask_high land b) lsl delta_d
        else if delta_d = 0 then
          mask_high land b
        else
          (mask_high land b) lsr (-delta_d) in
      drb (dim p + delta_d) (rank p + delta_r) (b_low lor b_med lor b_high)

  end

(* \thocwmodulesection{Whizard} *)

module type Whizard =
  sig
    type t
    val of_momentum : t -> int
    val to_momentum : int -> int -> t
  end

module BitsW =
  struct
    type t = Bits.t
    open Bits (* NB: this includes the internal functions not in [T]! *)

    let of_momentum p =
      let d = dim p in
      let bit_in1 = 1 land p
      and bit_in2 = 1 land (p lsr 1)
      and bits_out = ((mask d) land p) lsr 2 in
      bits_out lor (bit_in1 lsl (d - 1)) lor (bit_in2 lsl (d - 2))

    let rec count_non_zero' acc i last b =
      if i > last then
        acc
      else if (1 lsl (pred i)) land b = 0 then
        count_non_zero' acc (succ i) last b
      else
        count_non_zero' (succ acc) (succ i) last b

    let count_non_zero first last b =
      count_non_zero' 0 first last b

    let to_momentum d w =
      let bit_in1 = 1 land (w lsr (d - 1))
      and bit_in2 = 1 land (w lsr (d - 2))
      and bits_out = (mask (d - 2)) land w in
      let b = (bits_out lsl 2) lor bit_in1 lor (bit_in2 lsl 1) in
      drb d (count_non_zero 1 d b) b

  end

(* The following would be a tad more efficient, if coded directly, but
   there's no point in wasting effort on this. *)

module ListsW =
  struct
    type t = Lists.t
    let of_momentum p =
      BitsW.of_momentum (Bits.of_ints p.Lists.d p.Lists.p)
    let to_momentum d w =
      Lists.of_ints d (Bits.to_ints (BitsW.to_momentum d w))
  end

(* \thocwmodulesection{Suggesting a Default Implementation} *)

(* [Lists] is better tested, but the more recent [Bits] appears to
   work as well and is \emph{much} more efficient, resulting in a
   relative factor of better than 2.  This performance ratio
   is larger than I had expected and we are not likely to
   reach its limit of 21 independent vectors anyway.  *)

module Default = Bits
module DefaultW = BitsW

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
