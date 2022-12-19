(* cache.ml --

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


let search_path =
  [ Filename.current_dir_name ]

(*i
let search_path =
  [ Filename.current_dir_name;
    ThoFilename.expand_home Config.user_cache_dir;
    Config.system_cache_dir ]
i*)

module type T =
  sig

    type key
    type hash = string
    type value

    type 'a result = 
      | Hit of 'a
      | Miss 
      | Stale of string

    exception Mismatch of string * string * string

    val hash : key -> hash
    val exists : hash -> string -> bool
    val find : hash -> string -> string option
    val write : hash -> string -> value -> unit
    val write_dir : hash -> string -> string -> value -> unit
    val read : hash -> string -> value
    val maybe_read : hash -> string -> value result

  end

module type Key =
  sig
    type t
  end

module type Value =
  sig
    type t
  end

module Make (Key : Key) (Value : Value) =
  struct

    type key = Key.t
    type hash = string
    type value = Value.t
          
    type tagged =
        { tag : hash;
          value : value; }

    let hash value =
      Digest.string (Marshal.to_string value [])

    let find_first path name =
      let rec find_first' = function
        | [] -> raise Not_found
        | dir :: path ->
            let f = Filename.concat dir name in
            if Sys.file_exists f then
              f
            else
              find_first' path
      in
      find_first' path

    let find hash name =
      try Some (find_first search_path name) with Not_found -> None

    let exists hash name =
      match find hash name with
      | None -> false
      | Some _ -> true

    let try_first f path name =
      let rec try_first' = function
        | [] -> raise Not_found
        | dir :: path ->
            try (f (Filename.concat dir name), dir) with _ -> try_first' path
      in
      try_first' path

    let open_in_bin_first = try_first open_in_bin
    let open_out_bin_last path = try_first open_out_bin (List.rev path)

    let write hash name value =
      let oc, _ = open_out_bin_last search_path name in
      Marshal.to_channel oc { tag = hash; value = value } [];
      close_out oc

    let write_dir hash dir name value =
      let oc = open_out_bin (Filename.concat dir name) in
      Marshal.to_channel oc { tag = hash; value = value } [];
      close_out oc

    type 'a result = 
      | Hit of 'a
      | Miss 
      | Stale of string

    exception Mismatch of string * string * string

    let read hash name =
      let ic, dir = open_in_bin_first search_path name in
      let { tag = tag; value = value } = Marshal.from_channel ic in
      close_in ic;
      if tag = hash then
        value
      else
        raise (Mismatch (Filename.concat dir name, hash, tag))

    let maybe_read hash name =
      try
        Hit (read hash name)
      with
      | Not_found -> Miss
      | Mismatch (file, _, _) -> Stale file
      
  end

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)





