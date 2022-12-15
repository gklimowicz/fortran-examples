! Copyright (c) 2011-2016 Harald Klimach <harald@klimachs.de>
! Copyright (c) 2012, 2014 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012 James Spencer <j.spencer@imperial.ac.uk>
! Copyright (c) 2014-2015 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2015 Verena Krupp
!
! Parts of this file were written by Harald Klimach and Simon Zimny
! German Research School of Simulation Sciences.
! Parts of this file were written by Kannan Masilamani, Harald Klimach,
! Peter Vitt and Verena Krupp for University of Siegen.
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in
! all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
! IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
! DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
! OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
! OR OTHER DEALINGS IN THE SOFTWARE.
! **************************************************************************** !

!> The aot_path can be used to track the position of a Lua entity in nested
!! tables.
!!
!! @warning This is obsolete! Use [[aot_references_module]] instead.
!!          Please note that this module might be removed in future versions
!!          of Aotus.
!!
!! The idea is to initialize the path in the very beginning and then append a
!! node whenever a table is opened. Thus you pass down the growing path object
!! and store at in the level, to which you might need to return later.
module aot_path_module
  use flu_binding, only: flu_State
  use aotus_module, only: open_config_file, close_config
  use aot_table_module, only: aot_table_open, aot_table_close
  use aot_fun_module, only: aot_fun_type, aot_fun_open, aot_fun_close

  implicit none

  private

  !> This data structure describes a node in the path through nested tables.
  type aot_path_node_type
    !> What type of node is this?
    !! Currently supported are function and table
    character(len=16) :: NodeType

    !> How to look up this node, by key or position?
    character(len=16) :: ID_kind

    !> Identifying key
    character(len=80) :: key

    !> Identifying position
    integer :: pos

    !> Link to possible child of this node
    type(aot_path_node_type), pointer :: child => NULL()
  end type

  !> This type is the main data structure of the module and describes the path.
  !!
  !! It contains a linked list of all nodes, as well as the name of the Lua
  !! script where this path is recorded in.
  type aot_path_type
    private
    !> Name of the file where this path object is found in.
    character(len=256) :: LuaFilename

    !> Handle to the topmost table opened for the path.
    integer :: rootHandle

    !> Entry level of the path on the global scope of the Lua script.
    type(aot_path_node_type), pointer :: GlobalNode => NULL()

    !> Moving head through the linked list of path nodes.
    type(aot_path_node_type), pointer :: head => NULL()
  end type

  !> Taking care of the linked list in a copying routine for the assignment of
  !! aot_path_type.
  interface assignment(=)
    module procedure aot_path_copy
  end interface

  public :: aot_path_type
  public :: aot_init_path, aot_fin_path
  public :: aot_path_addNode, aot_path_delNode
  public :: assignment(=)
  public :: aot_path_open, aot_path_close
  public :: aot_path_toString
  public :: aot_path_dump

  !> Re-open a previously recorded path through nested Lua tables.
  !!
  !! This opens all the tables recursively down to the last node in the path.
  !! It might be used to open a table, or a function.
  interface aot_path_open
    module procedure aot_path_open_fun
    module procedure aot_path_open_table
  end interface aot_path_open

  ! Close all tables, that were opened for the given path.
  interface aot_path_close
    module procedure aot_path_close_fun
    module procedure aot_path_close_table
  end interface aot_path_close

contains

  !> This subroutine initializes a path object.
  !!
  !! This is done by setting the given file name as reference to the script,
  !! to look the path up in and emptying the path completely.
  subroutine aot_init_path(me, Filename)
    !> Path object to initialize
    type(aot_path_type), intent(out) :: me
    !> Filename of the Lua script, this path is located in
    character(len=*), optional, intent(in) :: Filename

    ! Finalize the path first, just in case it might have had any entries.
    call aot_fin_path(me)
    if (present(Filename)) then
      me%LuaFilename = adjustl(trim(Filename))
    else
      me%LuaFilename = ''
    end if
    me%rootHandle = 0
  end subroutine aot_init_path


  !> This subroutine finalizes a path object and deallocates
  !! all its nodes.
  subroutine aot_fin_path(me)
    !> Path to destroy
    type(aot_path_type), intent(inout) :: me

    logical :: emptied

    emptied = .false.

    do while (.not. emptied)
      call aot_path_delNode(me, emptied)
    end do
    me%LuaFilename = ''
    me%rootHandle = 0
  end subroutine aot_fin_path


  !> With this subroutine a node is appended to the end of
  !! the list of nodes of the given path.
  !!
  !! You need to provide a NodeType (table or function),
  !! and either its position or key to identify it in the
  !! parent object.
  subroutine aot_path_addNode(me, NodeType, pos, key)
    !> Path to append the node to
    type(aot_path_type), intent(inout) :: me
    !> Type of the node (table of function)
    character(len=*), intent(in) :: NodeType
    !> Position in the parenting table
    integer, intent(in), optional :: pos
    !> Key within the parenting table
    character(len=*), intent(in), optional :: key

    if (.not.associated(me%GlobalNode)) then
      ! New list without any nodes so far
      allocate(me%GlobalNode)
      me%head => me%GlobalNode
    else
      ! Existing list, append at the end
      allocate(me%head%child)
      me%head => me%head%child
    end if

    if (present(pos)) then
      me%head%ID_kind = 'position'
      me%head%pos = pos
    end if

    ! Specified keys overwrite positions
    if (present(key)) then
      me%head%ID_kind = 'key'
      me%head%key = key
    end if

    me%head%NodeType = NodeType

  end subroutine aot_path_addNode


  !> The delNode removes the last node from the list of nodes of the given path.
  !!
  !! With the optional isEmpty argument, it can be tested, if the list
  !! is completely empty after this operation.
  subroutine aot_path_delNode(me, isEmpty)
    !> Path to delet the last node from
    type(aot_path_type), intent(inout) :: me
    !> Flag, if resulting path is empty (contains no nodes anymore)
    logical, intent(out), optional :: isEmpty

    type(aot_path_node_type), pointer :: curNode => NULL()
    logical :: emptyList

    emptyList = .true.

    if (associated(me%GlobalNode)) then
      curNode => me%GlobalNode

      do
        if (associated(curNode%child)) then
          if (associated(curNode%child, me%head)) then
            ! Found second Last Node (its child is the head)
            nullify(curNode%child)
            deallocate(me%head)
            me%head => curNode
            ! The list is not empty, there is at least one
            ! node remaining.
            emptyList = .false.
            ! Leave the loop
            exit
          end if
        else
          ! There is just the global node, no childs yet
          nullify(me%globalNode)
          deallocate(me%head)
          ! Leave the loop
          exit
        end if
        curNode => curNode%child
      end do

    end if

    if (present(isEmpty)) then
      isEmpty = emptyList
    end if

  end subroutine aot_path_delNode


  !> Copy a given path object, this is the implementation of the
  !! assignment left = right.
  subroutine aot_path_copy(left, right)
    !> Object to assign a path to
    type(aot_path_type), intent(inout) :: left
    !> Path to be copied
    type(aot_path_type), intent(in) :: right

    type(aot_path_node_type), pointer :: curNode

    call aot_fin_path(left)

    left%LuaFilename = right%LuaFilename
    left%roothandle = right%roothandle

    if (associated(right%globalNode)) then
      allocate(left%globalNode)
      left%globalNode%NodeType = right%globalNode%NodeType
      left%globalNode%ID_kind = right%globalNode%ID_kind
      left%globalNode%key = right%globalNode%key
      left%globalNode%pos = right%globalNode%pos
      left%head => left%globalNode

      curNode => right%globalNode
      do while(associated(curNode%child))
        allocate(left%head%child)
        curNode => curNode%child
        left%head => left%head%child

        left%head%NodeType = curNode%NodeType
        left%head%ID_kind = curNode%ID_kind
        left%head%key = curNode%key
        left%head%pos = curNode%pos
      end do
    end if

  end subroutine aot_path_copy


  !> This subroutine opens all the tables on the way to the final head node,
  !! which ought to be a function.
  !!
  !! The given fun object is then filled by an aot_fun_open
  !! on the head of the given path.
  !! The handle can be either passed in, to be used for the
  !! look up of the path, or, when specifying the optional
  !! openLua argument as true, it will return the handle to
  !! the newly opened Lua script.
  subroutine aot_path_open_fun(me, conf, fun, openLua)
    !> The path object to open as a function
    type(aot_path_type), intent(inout) :: me
    !> The flu_state handle, which is either opened according to
    !! the path, or used to open the path in.
    type(flu_state) :: conf
    !> The opened function
    type(aot_fun_type), intent(out) :: fun
    !> A flag to indicate, wether to open the Lua script, default
    !! is false, in which case the conf argument has to link to
    !! an actual Lua state handle.
    logical, intent(in), optional :: openLua

    integer :: myHandle = 0

    ! open the table until it reaches the final head node
    call aot_path_open_table( me, conf, myHandle, openLua )

    if (me%head%NodeType == 'function') then
      select case(me%head%ID_kind)
      case('key')
        if (associated(me%head, me%GlobalNode)) then
          call aot_fun_open(L=conf, fun=fun, key=me%head%key)
        else
          call aot_fun_open(L=conf, parent=myHandle, fun=fun, key=me%head%key)
        end if
      case('position')
        call aot_fun_open(L=conf, parent=myHandle, fun=fun, pos=me%head%pos)
      end select
    end if

  end subroutine aot_path_open_fun


  !> This subroutine opens all the tables on the way to the final head node of
  !! the given path.
  !!
  !! The handle can be either passed in, to be used for the
  !! look up of the path, or, when specifying the optional
  !! openLua argument as true, it will return the handle to
  !! the newly opened Lua script.
  subroutine aot_path_open_table(me, conf, thandle, openLua)
    !> The path object to open as a function
    type(aot_path_type), intent(inout) :: me
    !> The flu_state handle, which is either opened according to
    !! the path, or used to open the path in.
    type(flu_state) :: conf
    !> return handle of the last opened table
    integer, intent(out) :: thandle
    !> A flag to indicate, wether to open the Lua script, default
    !! is false, in which case the conf argument has to link to
    !! an actual Lua state handle.
    logical, intent(in), optional :: openLua

    logical :: new_conf
    type(aot_path_node_type), pointer :: curNode => NULL()
    integer :: myHandle = 0
    integer :: prevHandle = -1

    if (present(openLua)) then
      new_conf = openLua
    else
      new_conf = .false.
    end if

    if (new_conf) then
      call open_config_file(conf, me%LuaFilename)
    end if

    curNode => me%GlobalNode

    if (curNode%NodeType == 'table') then
      select case(curNode%ID_kind)
      case('key')
        call aot_table_open(L=conf, thandle=me%roothandle, key=curNode%key)
      end select
      if (associated(curNode%child)) then
        curNode => curNode%child
        myHandle = me%rootHandle
      end if
    end if

    do while(associated(curNode%child))
      prevHandle = myHandle

      select case(curNode%ID_kind)
      case('key')
        call aot_table_open(L=conf, thandle=myHandle, parent=prevHandle, &
          &                 key=curNode%key)
      case('position')
        call aot_table_open(L=conf, thandle=myHandle, parent=prevHandle, &
          &                 pos=curNode%pos)
      end select
      curNode => curNode%child

    end do

    thandle = myHandle

  end subroutine aot_path_open_table


  !> This routine closes function and all other tables opened along the path.
  subroutine aot_path_close_fun(me, conf, fun, closeLua)
    !> The path object to open as a function
    type(aot_path_type), intent(inout) :: me
    !> The flu_state handle, which is either opened according to
    !! the path, or used to open the path in.
    type(flu_state) :: conf
    !> The opened function
    type(aot_fun_type), intent(inout) :: fun
    !> A flag to indicate, wether to close the Lua script, default
    !! is false.
    logical, intent(in), optional :: closeLua

    ! close function
    call aot_fun_close(L=conf, fun=fun)
    ! close tables
    call aot_path_close_table( me, conf, closeLua )

  end subroutine aot_path_close_fun


  !> This routine closes all the table opened in aot_path_open_table.
  subroutine aot_path_close_table(me, conf, closeLua)
    !> The path object to open as a function
    type(aot_path_type), intent(inout) :: me
    !> The flu_state handle, which is either opened according to
    !! the path, or used to open the path in.
    type(flu_state) :: conf
    !> A flag to indicate, wether to close the Lua script, default
    !! is false.
    logical, intent(in), optional :: closeLua

    if (me%roothandle /= 0) then
      call aot_table_close(L=conf, thandle=me%roothandle)
    end if

    if (present(closeLua)) then
      if (closeLua) then
        call close_config(conf)
      end if
    end if

  end subroutine aot_path_close_table

  !> Dumps the complete path into a string.
  !!
  !! This routine transforms a given path into a special notation. Each element
  !! is added to the string, separated by a . char.
  !! If the resulting string is to long for the provided buffer /ref
  !! pathAsString, the buffer will stay empty to not have the caller proceed
  !! with incomplete results.
  subroutine aot_path_toString( path, pathAsString )
    !> The path which information should be printed
    type(aot_path_type), intent(in) :: path
    !> The path represented as string
    character(len=*), intent(out) :: pathAsString

    type(aot_path_node_type), pointer :: current
    integer :: pathLength
    integer :: stringLength
    character(len=10) :: posstr

    stringLength = len(pathAsString)
    pathLength = 0

    ! First we measure the size of the result
    if (associated(path%globalNode)) then
      current => path%globalNode
      do while(associated(current))
        if (associated(current,path%globalNode)) then
          ! Add the length of the first node
          pathLength = len_trim(adjustl(current%key))
        else
          if (trim(current%ID_kind) == 'key') then
            ! Add the length of a following node and the delimiter char
            pathLength = pathLength + len_trim(adjustl(current%key)) + 1
          else
            ! Length of the position number and 2 places for brackets.
            write(posstr,'(i0)') current%pos
            pathLength = pathLength + len_trim(posstr) + 2
          end if
        end if
        current => current%child
      end do
    end if

    ! If the result fits into the buffer, we create it
    if (pathLength <= stringLength .and. pathLength > 0) then
      current => path%globalNode
      do while(associated(current))
        if (associated(current,path%globalNode)) then
          pathAsString = trim(adjustl(current%key))
        else
          if (trim(current%ID_kind) == 'key') then
            pathAsString = trim(pathAsString) // '.' &
              &            // trim(adjustl(current%key))
          else
            write(posstr,'(i0)') current%pos
            pathAsString = trim(pathAsString) // '[' &
              &            // trim(posstr) // ']'
          end if
        end if
        current => current%child
      end do
    else
      ! Either the result is empty or too long, thus we clear the buffer
      pathAsString = ''
    end if

  end subroutine aot_path_toString


  !> Dumps the complete path to the given output unit.
  !!
  !! This routine is for debugging purposes. It takes the path and, beginning
  !! with the global node, dumps all following nodes to the output unit provided
  !! by the caller.
  subroutine aot_path_dump( path, outputUnit )
    !> The path which information should be printed
    type(aot_path_type), intent(in) :: path
    !> The unit to use to write the path data
    integer, intent(in) :: outputUnit

    type(aot_path_node_type), pointer :: current

    write(outputUnit,*) 'Path:'
    write(outputUnit,*) '  Filename: ', path%LuaFilename
    write(outputUnit,'(A,I10)') '   root handle: ', path%rootHandle
    if (associated(path%globalNode)) then
      current => path%globalNode
      do while(associated(current))
        if(associated(current,path%globalNode)) then
          write(outputUnit,*) '  Global node: '
        else
          write(outputUnit,*) '  next: '
        end if
        write(outputUnit,*) '    NodeType: ', current%NodeType
        write(outputUnit,*) '    ID_Kind: ', current%ID_Kind
        if (trim(current%ID_Kind) == 'key') then
          write(outputUnit,*) '    key: ', current%key
        else
          write(outputUnit,'(A,I10)') '     pos: ', current%pos
        end if
        current => current%child
      end do
    end if

  end subroutine aot_path_dump

end module aot_path_module
