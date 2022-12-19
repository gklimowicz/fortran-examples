(* circe2/thoMatrix.ml --  *)
(* Copyright (C) 2001-2022 by Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
   Circe2 is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by 
   the Free Software Foundation; either version 2, or (at your option)
   any later version.
   Circe2 is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
   GNU General Public License for more details.
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.  *)  

let map ?inf1 ?sup1 ?inf2 ?sup2 f a =
  ThoArray.map ?inf:inf1 ?sup:sup1
    (ThoArray.map ?inf:inf2 ?sup:sup2 f) a

let copy ?inf1 ?sup1 ?inf2 ?sup2 a =
  map ?inf1 ?sup1 ?inf2 ?sup2 (fun x -> x) a

let iter ?inf1 ?sup1 ?inf2 ?sup2 f a =
  ThoArray.iter ?inf:inf1 ?sup:sup1
    (ThoArray.iter ?inf:inf2 ?sup:sup2 f) a

let fold_left ?inf1 ?sup1 ?inf2 ?sup2 f x a =
  ThoArray.fold_left ?inf:inf1 ?sup:sup1
    (ThoArray.fold_left ?inf:inf2 ?sup:sup2 f) x a

let sum_float ?inf1 ?sup1 ?inf2 ?sup2 a =
  fold_left ?inf1 ?sup1 ?inf2 ?sup2 (+.) 0.0 a

let size a =
  Array.fold_left (fun acc v -> Array.length v + acc) 0 a

let transpose a =
  let n1 = Array.length a
  and n2 = Array.length a.(0) in
  let a' = Array.make_matrix n2 n1 a.(0).(0) in
  for i1 = 0 to pred n1 do
    for i2 = 0 to pred n2 do
      a'.(i2).(i1) <- a.(i1).(i2)
    done
  done;
  a'

let suite =
  let open OUnit in
  "Matrix" >:::
  []

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
