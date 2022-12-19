(* circe2/commands.mli --  *)
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

(* An example for a command file:
\begin{verbatim}
{ file = "tesla.circe"
   { design = "TESLA" roots = 500
     { pid/1 = electron pid/2 = positron
       events = "tesla_500.electron_positron" }
     { pid = photon
       events = "tesla_500.gamma_gamma" }
     { pid/1 = photon pid/2 = positron
       events = "tesla_500.gamma_positron" }
     { pid/1 = electron pid/2 = photon
       events = "tesla_500.electron_gamma" } }
   { design = "TESLA" roots = 800
     { pid/1 = electron pid/2 = positron
       events = "tesla_800.electron_positron" } }
   { design = "TESLA" roots = 500
     { pid = photon
       events = "tesla_gg_500.gamma_gamma" } }
   { design = "TESLA" roots = 500
     { pid = electron
       events = "tesla_ee_500.electron_electron" } } }
\end{verbatim}
*)

exception Invalid_interval of float * float

type t
val parse_file : string -> t
val parse_string : string -> t
val execute : t -> unit

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
