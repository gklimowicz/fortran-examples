MODULE OpenGL_glut

!  Derived from new_freeglut.h using glut_interfaces.pl
!  new_freeglut.h is /usr/include/GC/freeglut.h with the included files
!  freeglut_std.h and freeglut_ext.h pasted in explicitly.

USE OpenGL_kinds
IMPLICIT NONE
PRIVATE


!  freeglut.h
!  
!  The freeglut library include file
!  
!  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
!  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
!  PAWEL W. OLSZTA BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
!  IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
!  CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

!  ** #include "freeglut_std.h" **

!  freeglut_std.h
!  
!  The GLUT-compatible part of the freeglut library include file
!  
!  Copyright (c) 1999-2000 Pawel W. Olszta. All Rights Reserved.
!  Written by Pawel W. Olszta, <olszta@sourceforge.net>
!  Creation date: Thu Dec 2 1999
!  
!  Permission is hereby granted, free of charge, to any person obtaining a
!  copy of this software and associated documentation files (the "Software"),
!  to deal in the Software without restriction, including without limitation
!  the rights to use, copy, modify, merge, publish, distribute, sublicense,
!  and/or sell copies of the Software, and to permit persons to whom the
!  Software is furnished to do so, subject to the following conditions:
!  
!  The above copyright notice and this permission notice shall be included
!  in all copies or substantial portions of the Software.
!  
!  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
!  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
!  PAWEL W. OLSZTA BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
!  IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
!  CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


!  Under windows, we have to differentiate between static and dynamic libraries
!  Discussion by FreeGLUT developers suggests that
!  Visual C++ specific code involving pragmas may
!  need to move to a separate header.  24th Dec 2003

!  pragmas or to 1 to exclude library pragmas.
!  The default behavior depends on the compiler/platform.


!  Windows static library 



!  Windows shared library (DLL) 





!  Drag in other Windows libraries as required by FreeGLUT 


!  Non-Windows definition of FGAPI and FGAPIENTRY  


!  The freeglut and GLUT API versions
INTEGER(GLenum), PARAMETER, PUBLIC :: FREEGLUT = 1
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_API_VERSION = 4
INTEGER(GLenum), PARAMETER, PUBLIC :: FREEGLUT_VERSION_2_0 = 1
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_XLIB_IMPLEMENTATION = 13

!  Always include OpenGL and GLU headers

!  GLUT API macro definitions -- the special key codes:
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_F1 = int(z'0001', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_F2 = int(z'0002', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_F3 = int(z'0003', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_F4 = int(z'0004', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_F5 = int(z'0005', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_F6 = int(z'0006', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_F7 = int(z'0007', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_F8 = int(z'0008', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_F9 = int(z'0009', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_F10 = int(z'000A', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_F11 = int(z'000B', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_F12 = int(z'000C', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_LEFT = int(z'0064', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_UP = int(z'0065', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_RIGHT = int(z'0066', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_DOWN = int(z'0067', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_PAGE_UP = int(z'0068', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_PAGE_DOWN = int(z'0069', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_HOME = int(z'006A', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_END = int(z'006B', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_INSERT = int(z'006C', kind=GLenum)

!  GLUT API macro definitions -- mouse state definitions
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_LEFT_BUTTON = int(z'0000', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_MIDDLE_BUTTON = int(z'0001', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_RIGHT_BUTTON = int(z'0002', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_DOWN = int(z'0000', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_UP = int(z'0001', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_LEFT = int(z'0000', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_ENTERED = int(z'0001', kind=GLenum)

!  GLUT API macro definitions -- the display mode definitions
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_RGB = int(z'0000', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_RGBA = int(z'0000', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_INDEX = int(z'0001', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_SINGLE = int(z'0000', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_DOUBLE = int(z'0002', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_ACCUM = int(z'0004', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_ALPHA = int(z'0008', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_DEPTH = int(z'0010', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_STENCIL = int(z'0020', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_MULTISAMPLE = int(z'0080', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_STEREO = int(z'0100', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_LUMINANCE = int(z'0200', kind=GLenum)

!  GLUT API macro definitions -- windows and menu related definitions
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_MENU_NOT_IN_USE = int(z'0000', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_MENU_IN_USE = int(z'0001', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_NOT_VISIBLE = int(z'0000', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_VISIBLE = int(z'0001', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_HIDDEN = int(z'0000', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_FULLY_RETAINED = int(z'0001', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_PARTIALLY_RETAINED = int(z'0002', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_FULLY_COVERED = int(z'0003', kind=GLenum)

!  GLUT API macro definitions -- fonts definitions
!  
!  Steve Baker suggested to make it binary compatible with GLUT:


!  GLUT API macro definitions -- the glutGet parameters
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_X = int(z'0064', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_Y = int(z'0065', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_WIDTH = int(z'0066', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_HEIGHT = int(z'0067', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_BUFFER_SIZE = int(z'0068', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_STENCIL_SIZE = int(z'0069', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_DEPTH_SIZE = int(z'006A', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_RED_SIZE = int(z'006B', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_GREEN_SIZE = int(z'006C', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_BLUE_SIZE = int(z'006D', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_ALPHA_SIZE = int(z'006E', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_ACCUM_RED_SIZE = int(z'006F', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_ACCUM_GREEN_SIZE = int(z'0070', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_ACCUM_BLUE_SIZE = int(z'0071', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_ACCUM_ALPHA_SIZE = int(z'0072', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_DOUBLEBUFFER = int(z'0073', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_RGBA = int(z'0074', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_PARENT = int(z'0075', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_NUM_CHILDREN = int(z'0076', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_COLORMAP_SIZE = int(z'0077', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_NUM_SAMPLES = int(z'0078', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_STEREO = int(z'0079', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_CURSOR = int(z'007A', kind=GLenum)

INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_SCREEN_WIDTH = int(z'00C8', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_SCREEN_HEIGHT = int(z'00C9', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_SCREEN_WIDTH_MM = int(z'00CA', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_SCREEN_HEIGHT_MM = int(z'00CB', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_MENU_NUM_ITEMS = int(z'012C', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_DISPLAY_MODE_POSSIBLE = int(z'0190', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_INIT_WINDOW_X = int(z'01F4', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_INIT_WINDOW_Y = int(z'01F5', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_INIT_WINDOW_WIDTH = int(z'01F6', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_INIT_WINDOW_HEIGHT = int(z'01F7', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_INIT_DISPLAY_MODE = int(z'01F8', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_ELAPSED_TIME = int(z'02BC', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_FORMAT_ID = int(z'007B', kind=GLenum)

!  GLUT API macro definitions -- the glutDeviceGet parameters
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_HAS_KEYBOARD = int(z'0258', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_HAS_MOUSE = int(z'0259', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_HAS_SPACEBALL = int(z'025A', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_HAS_DIAL_AND_BUTTON_BOX = int(z'025B', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_HAS_TABLET = int(z'025C', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_NUM_MOUSE_BUTTONS = int(z'025D', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_NUM_SPACEBALL_BUTTONS = int(z'025E', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_NUM_BUTTON_BOX_BUTTONS = int(z'025F', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_NUM_DIALS = int(z'0260', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_NUM_TABLET_BUTTONS = int(z'0261', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_DEVICE_IGNORE_KEY_REPEAT = int(z'0262', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_DEVICE_KEY_REPEAT = int(z'0263', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_HAS_JOYSTICK = int(z'0264', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_OWNS_JOYSTICK = int(z'0265', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_JOYSTICK_BUTTONS = int(z'0266', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_JOYSTICK_AXES = int(z'0267', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_JOYSTICK_POLL_RATE = int(z'0268', kind=GLenum)

!  GLUT API macro definitions -- the glutLayerGet parameters
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_OVERLAY_POSSIBLE = int(z'0320', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_LAYER_IN_USE = int(z'0321', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_HAS_OVERLAY = int(z'0322', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_TRANSPARENT_INDEX = int(z'0323', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_NORMAL_DAMAGED = int(z'0324', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_OVERLAY_DAMAGED = int(z'0325', kind=GLenum)

!  GLUT API macro definitions -- the glutVideoResizeGet parameters
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_VIDEO_RESIZE_POSSIBLE = int(z'0384', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_VIDEO_RESIZE_IN_USE = int(z'0385', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_VIDEO_RESIZE_X_DELTA = int(z'0386', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_VIDEO_RESIZE_Y_DELTA = int(z'0387', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_VIDEO_RESIZE_WIDTH_DELTA = int(z'0388', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_VIDEO_RESIZE_HEIGHT_DELTA = int(z'0389', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_VIDEO_RESIZE_X = int(z'038A', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_VIDEO_RESIZE_Y = int(z'038B', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_VIDEO_RESIZE_WIDTH = int(z'038C', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_VIDEO_RESIZE_HEIGHT = int(z'038D', kind=GLenum)

!  GLUT API macro definitions -- the glutUseLayer parameters
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_NORMAL = int(z'0000', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_OVERLAY = int(z'0001', kind=GLenum)

!  GLUT API macro definitions -- the glutGetModifiers parameters
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_ACTIVE_SHIFT = int(z'0001', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_ACTIVE_CTRL = int(z'0002', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_ACTIVE_ALT = int(z'0004', kind=GLenum)

!  GLUT API macro definitions -- the glutSetCursor parameters
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_RIGHT_ARROW = int(z'0000', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_LEFT_ARROW = int(z'0001', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_INFO = int(z'0002', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_DESTROY = int(z'0003', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_HELP = int(z'0004', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_CYCLE = int(z'0005', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_SPRAY = int(z'0006', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_WAIT = int(z'0007', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_TEXT = int(z'0008', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_CROSSHAIR = int(z'0009', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_UP_DOWN = int(z'000A', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_LEFT_RIGHT = int(z'000B', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_TOP_SIDE = int(z'000C', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_BOTTOM_SIDE = int(z'000D', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_LEFT_SIDE = int(z'000E', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_RIGHT_SIDE = int(z'000F', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_TOP_LEFT_CORNER = int(z'0010', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_TOP_RIGHT_CORNER = int(z'0011', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_BOTTOM_RIGHT_CORNER = int(z'0012', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_BOTTOM_LEFT_CORNER = int(z'0013', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_INHERIT = int(z'0064', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_NONE = int(z'0065', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_FULL_CROSSHAIR = int(z'0066', kind=GLenum)

!  GLUT API macro definitions -- RGB color component specification definitions
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_RED = int(z'0000', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_GREEN = int(z'0001', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_BLUE = int(z'0002', kind=GLenum)

!  GLUT API macro definitions -- additional keyboard and joystick definitions
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_REPEAT_OFF = int(z'0000', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_REPEAT_ON = int(z'0001', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_REPEAT_DEFAULT = int(z'0002', kind=GLenum)

INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_JOYSTICK_BUTTON_A = int(z'0001', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_JOYSTICK_BUTTON_B = int(z'0002', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_JOYSTICK_BUTTON_C = int(z'0004', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_JOYSTICK_BUTTON_D = int(z'0008', kind=GLenum)

!  GLUT API macro definitions -- game mode definitions
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_GAME_MODE_ACTIVE = int(z'0000', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_GAME_MODE_POSSIBLE = int(z'0001', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_GAME_MODE_WIDTH = int(z'0002', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_GAME_MODE_HEIGHT = int(z'0003', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_GAME_MODE_PIXEL_DEPTH = int(z'0004', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_GAME_MODE_REFRESH_RATE = int(z'0005', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_GAME_MODE_DISPLAY_CHANGED = int(z'0006', kind=GLenum)

!  Initialization functions, see fglut_init.c
!  void    glutInit( int* pargc, char** argv )
PUBLIC glutInit
INTERFACE glutInit
MODULE PROCEDURE glutInit_f03
SUBROUTINE glutInit_gl(pargc, argv) BIND(C,NAME="glutInit")
IMPORT
! INTEGER(GLint) :: pargc
INTEGER(GLint), DIMENSION(*) :: pargc
TYPE(C_PTR), INTENT(IN) :: argv
END SUBROUTINE glutInit_gl
END INTERFACE

!  void    glutInitWindowPosition( int x, int y )
PUBLIC glutInitWindowPosition
INTERFACE
SUBROUTINE glutInitWindowPosition(x, y) BIND(C,NAME="glutInitWindowPosition")
IMPORT
INTEGER(GLint), VALUE :: x, y
END SUBROUTINE glutInitWindowPosition
END INTERFACE

!  void    glutInitWindowSize( int width, int height )
PUBLIC glutInitWindowSize
INTERFACE
SUBROUTINE glutInitWindowSize(width, height) BIND(C,NAME="glutInitWindowSize")
IMPORT
INTEGER(GLint), VALUE :: width, height
END SUBROUTINE glutInitWindowSize
END INTERFACE

!  void    glutInitDisplayMode( unsigned int displayMode )
PUBLIC glutInitDisplayMode
INTERFACE
SUBROUTINE glutInitDisplayMode(displayMode) BIND(C,NAME="glutInitDisplayMode")
IMPORT
INTEGER(GLuint), VALUE :: displayMode
END SUBROUTINE glutInitDisplayMode
END INTERFACE

!  void    glutInitDisplayString( const char* displayMode )
PUBLIC glutInitDisplayString
INTERFACE
SUBROUTINE glutInitDisplayString(displayMode) BIND(C,NAME="glutInitDisplayString")
IMPORT
! CHARACTER(C_CHAR), INTENT(IN) :: displayMode
CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: displayMode
END SUBROUTINE glutInitDisplayString
END INTERFACE


!  Process loop function, see freeglut_main.c
!  void    glutMainLoop( void )
PUBLIC glutMainLoop
INTERFACE
SUBROUTINE glutMainLoop() BIND(C,NAME="glutMainLoop")
IMPORT
END SUBROUTINE glutMainLoop
END INTERFACE


!  Window management functions, see freeglut_window.c
!  int     glutCreateWindow( const char* title )
PUBLIC glutCreateWindow
INTERFACE
FUNCTION glutCreateWindow(title) BIND(C,NAME="glutCreateWindow")
IMPORT
INTEGER(GLint) :: glutCreateWindow
! CHARACTER(C_CHAR), INTENT(IN) :: title
CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: title
END FUNCTION glutCreateWindow
END INTERFACE

!  int     glutCreateSubWindow( int window, int x, int y, int width, int height )
PUBLIC glutCreateSubWindow
INTERFACE
FUNCTION glutCreateSubWindow(window, x, y, width, height) BIND(C,NAME="glutCreateSubWindow")
IMPORT
INTEGER(GLint) :: glutCreateSubWindow
INTEGER(GLint), VALUE :: window, x, y, width, height
END FUNCTION glutCreateSubWindow
END INTERFACE

!  void    glutDestroyWindow( int window )
PUBLIC glutDestroyWindow
INTERFACE
SUBROUTINE glutDestroyWindow(window) BIND(C,NAME="glutDestroyWindow")
IMPORT
INTEGER(GLint), VALUE :: window
END SUBROUTINE glutDestroyWindow
END INTERFACE

!  void    glutSetWindow( int window )
PUBLIC glutSetWindow
INTERFACE
SUBROUTINE glutSetWindow(window) BIND(C,NAME="glutSetWindow")
IMPORT
INTEGER(GLint), VALUE :: window
END SUBROUTINE glutSetWindow
END INTERFACE

!  int     glutGetWindow( void )
PUBLIC glutGetWindow
INTERFACE
FUNCTION glutGetWindow() BIND(C,NAME="glutGetWindow")
IMPORT
INTEGER(GLint) :: glutGetWindow
END FUNCTION glutGetWindow
END INTERFACE

!  void    glutSetWindowTitle( const char* title )
PUBLIC glutSetWindowTitle
INTERFACE
SUBROUTINE glutSetWindowTitle(title) BIND(C,NAME="glutSetWindowTitle")
IMPORT
! CHARACTER(C_CHAR), INTENT(IN) :: title
CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: title
END SUBROUTINE glutSetWindowTitle
END INTERFACE

!  void    glutSetIconTitle( const char* title )
PUBLIC glutSetIconTitle
INTERFACE
SUBROUTINE glutSetIconTitle(title) BIND(C,NAME="glutSetIconTitle")
IMPORT
! CHARACTER(C_CHAR), INTENT(IN) :: title
CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: title
END SUBROUTINE glutSetIconTitle
END INTERFACE

!  void    glutReshapeWindow( int width, int height )
PUBLIC glutReshapeWindow
INTERFACE
SUBROUTINE glutReshapeWindow(width, height) BIND(C,NAME="glutReshapeWindow")
IMPORT
INTEGER(GLint), VALUE :: width, height
END SUBROUTINE glutReshapeWindow
END INTERFACE

!  void    glutPositionWindow( int x, int y )
PUBLIC glutPositionWindow
INTERFACE
SUBROUTINE glutPositionWindow(x, y) BIND(C,NAME="glutPositionWindow")
IMPORT
INTEGER(GLint), VALUE :: x, y
END SUBROUTINE glutPositionWindow
END INTERFACE

!  void    glutShowWindow( void )
PUBLIC glutShowWindow
INTERFACE
SUBROUTINE glutShowWindow() BIND(C,NAME="glutShowWindow")
IMPORT
END SUBROUTINE glutShowWindow
END INTERFACE

!  void    glutHideWindow( void )
PUBLIC glutHideWindow
INTERFACE
SUBROUTINE glutHideWindow() BIND(C,NAME="glutHideWindow")
IMPORT
END SUBROUTINE glutHideWindow
END INTERFACE

!  void    glutIconifyWindow( void )
PUBLIC glutIconifyWindow
INTERFACE
SUBROUTINE glutIconifyWindow() BIND(C,NAME="glutIconifyWindow")
IMPORT
END SUBROUTINE glutIconifyWindow
END INTERFACE

!  void    glutPushWindow( void )
PUBLIC glutPushWindow
INTERFACE
SUBROUTINE glutPushWindow() BIND(C,NAME="glutPushWindow")
IMPORT
END SUBROUTINE glutPushWindow
END INTERFACE

!  void    glutPopWindow( void )
PUBLIC glutPopWindow
INTERFACE
SUBROUTINE glutPopWindow() BIND(C,NAME="glutPopWindow")
IMPORT
END SUBROUTINE glutPopWindow
END INTERFACE

!  void    glutFullScreen( void )
PUBLIC glutFullScreen
INTERFACE
SUBROUTINE glutFullScreen() BIND(C,NAME="glutFullScreen")
IMPORT
END SUBROUTINE glutFullScreen
END INTERFACE


!  Display-connected functions, see freeglut_display.c
!  void    glutPostWindowRedisplay( int window )
PUBLIC glutPostWindowRedisplay
INTERFACE
SUBROUTINE glutPostWindowRedisplay(window) BIND(C,NAME="glutPostWindowRedisplay")
IMPORT
INTEGER(GLint), VALUE :: window
END SUBROUTINE glutPostWindowRedisplay
END INTERFACE

!  void    glutPostRedisplay( void )
PUBLIC glutPostRedisplay
INTERFACE
SUBROUTINE glutPostRedisplay() BIND(C,NAME="glutPostRedisplay")
IMPORT
END SUBROUTINE glutPostRedisplay
END INTERFACE

!  void    glutSwapBuffers( void )
PUBLIC glutSwapBuffers
INTERFACE
SUBROUTINE glutSwapBuffers() BIND(C,NAME="glutSwapBuffers")
IMPORT
END SUBROUTINE glutSwapBuffers
END INTERFACE


!  Mouse cursor functions, see freeglut_cursor.c
!  void    glutWarpPointer( int x, int y )
PUBLIC glutWarpPointer
INTERFACE
SUBROUTINE glutWarpPointer(x, y) BIND(C,NAME="glutWarpPointer")
IMPORT
INTEGER(GLint), VALUE :: x, y
END SUBROUTINE glutWarpPointer
END INTERFACE

!  void    glutSetCursor( int cursor )
PUBLIC glutSetCursor
INTERFACE
SUBROUTINE glutSetCursor(cursor) BIND(C,NAME="glutSetCursor")
IMPORT
INTEGER(GLint), VALUE :: cursor
END SUBROUTINE glutSetCursor
END INTERFACE


!  Overlay stuff, see freeglut_overlay.c
!  void    glutEstablishOverlay( void )
PUBLIC glutEstablishOverlay
INTERFACE
SUBROUTINE glutEstablishOverlay() BIND(C,NAME="glutEstablishOverlay")
IMPORT
END SUBROUTINE glutEstablishOverlay
END INTERFACE

!  void    glutRemoveOverlay( void )
PUBLIC glutRemoveOverlay
INTERFACE
SUBROUTINE glutRemoveOverlay() BIND(C,NAME="glutRemoveOverlay")
IMPORT
END SUBROUTINE glutRemoveOverlay
END INTERFACE

!  void    glutUseLayer( GLenum layer )
PUBLIC glutUseLayer
INTERFACE
SUBROUTINE glutUseLayer(layer) BIND(C,NAME="glutUseLayer")
IMPORT
INTEGER(GLenum), VALUE :: layer
END SUBROUTINE glutUseLayer
END INTERFACE

!  void    glutPostOverlayRedisplay( void )
PUBLIC glutPostOverlayRedisplay
INTERFACE
SUBROUTINE glutPostOverlayRedisplay() BIND(C,NAME="glutPostOverlayRedisplay")
IMPORT
END SUBROUTINE glutPostOverlayRedisplay
END INTERFACE

!  void    glutPostWindowOverlayRedisplay( int window )
PUBLIC glutPostWindowOverlayRedisplay
INTERFACE
SUBROUTINE glutPostWindowOverlayRedisplay(window) BIND(C,NAME="glutPostWindowOverlayRedisplay")
IMPORT
INTEGER(GLint), VALUE :: window
END SUBROUTINE glutPostWindowOverlayRedisplay
END INTERFACE

!  void    glutShowOverlay( void )
PUBLIC glutShowOverlay
INTERFACE
SUBROUTINE glutShowOverlay() BIND(C,NAME="glutShowOverlay")
IMPORT
END SUBROUTINE glutShowOverlay
END INTERFACE

!  void    glutHideOverlay( void )
PUBLIC glutHideOverlay
INTERFACE
SUBROUTINE glutHideOverlay() BIND(C,NAME="glutHideOverlay")
IMPORT
END SUBROUTINE glutHideOverlay
END INTERFACE


!  Menu stuff, see freeglut_menu.c
!  int     glutCreateMenu( void (* callback)( int menu ) )
PUBLIC glutCreateMenu
INTERFACE
FUNCTION glutCreateMenu(proc) BIND(C,NAME="glutCreateMenu")
IMPORT
INTEGER(GLint) :: glutCreateMenu
!  void proc( int menu )
INTERFACE
SUBROUTINE proc(menu) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: menu
END SUBROUTINE proc
END INTERFACE
END FUNCTION glutCreateMenu
END INTERFACE

!  void    glutDestroyMenu( int menu )
PUBLIC glutDestroyMenu
INTERFACE
SUBROUTINE glutDestroyMenu(menu) BIND(C,NAME="glutDestroyMenu")
IMPORT
INTEGER(GLint), VALUE :: menu
END SUBROUTINE glutDestroyMenu
END INTERFACE

!  int     glutGetMenu( void )
PUBLIC glutGetMenu
INTERFACE
FUNCTION glutGetMenu() BIND(C,NAME="glutGetMenu")
IMPORT
INTEGER(GLint) :: glutGetMenu
END FUNCTION glutGetMenu
END INTERFACE

!  void    glutSetMenu( int menu )
PUBLIC glutSetMenu
INTERFACE
SUBROUTINE glutSetMenu(menu) BIND(C,NAME="glutSetMenu")
IMPORT
INTEGER(GLint), VALUE :: menu
END SUBROUTINE glutSetMenu
END INTERFACE

!  void    glutAddMenuEntry( const char* label, int value )
PUBLIC glutAddMenuEntry
INTERFACE
SUBROUTINE glutAddMenuEntry(label, value) BIND(C,NAME="glutAddMenuEntry")
IMPORT
! CHARACTER(C_CHAR), INTENT(IN) :: label
CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: label
INTEGER(GLint), VALUE :: value
END SUBROUTINE glutAddMenuEntry
END INTERFACE

!  void    glutAddSubMenu( const char* label, int subMenu )
PUBLIC glutAddSubMenu
INTERFACE
SUBROUTINE glutAddSubMenu(label, subMenu) BIND(C,NAME="glutAddSubMenu")
IMPORT
! CHARACTER(C_CHAR), INTENT(IN) :: label
CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: label
INTEGER(GLint), VALUE :: subMenu
END SUBROUTINE glutAddSubMenu
END INTERFACE

!  void    glutChangeToMenuEntry( int item, const char* label, int value )
PUBLIC glutChangeToMenuEntry
INTERFACE
SUBROUTINE glutChangeToMenuEntry(item, label, value) BIND(C,NAME="glutChangeToMenuEntry")
IMPORT
! CHARACTER(C_CHAR), INTENT(IN) :: label
CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: label
INTEGER(GLint), VALUE :: item, value
END SUBROUTINE glutChangeToMenuEntry
END INTERFACE

!  void    glutChangeToSubMenu( int item, const char* label, int value )
PUBLIC glutChangeToSubMenu
INTERFACE
SUBROUTINE glutChangeToSubMenu(item, label, value) BIND(C,NAME="glutChangeToSubMenu")
IMPORT
! CHARACTER(C_CHAR), INTENT(IN) :: label
CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: label
INTEGER(GLint), VALUE :: item, value
END SUBROUTINE glutChangeToSubMenu
END INTERFACE

!  void    glutRemoveMenuItem( int item )
PUBLIC glutRemoveMenuItem
INTERFACE
SUBROUTINE glutRemoveMenuItem(item) BIND(C,NAME="glutRemoveMenuItem")
IMPORT
INTEGER(GLint), VALUE :: item
END SUBROUTINE glutRemoveMenuItem
END INTERFACE

!  void    glutAttachMenu( int button )
PUBLIC glutAttachMenu
INTERFACE
SUBROUTINE glutAttachMenu(button) BIND(C,NAME="glutAttachMenu")
IMPORT
INTEGER(GLint), VALUE :: button
END SUBROUTINE glutAttachMenu
END INTERFACE

!  void    glutDetachMenu( int button )
PUBLIC glutDetachMenu
INTERFACE
SUBROUTINE glutDetachMenu(button) BIND(C,NAME="glutDetachMenu")
IMPORT
INTEGER(GLint), VALUE :: button
END SUBROUTINE glutDetachMenu
END INTERFACE


!  Global callback functions, see freeglut_callbacks.c
!  void    glutTimerFunc( unsigned int time, void (* callback)( int ), int value )
PUBLIC glutTimerFunc
INTERFACE
SUBROUTINE glutTimerFunc(time, proc, value) BIND(C,NAME="glutTimerFunc")
IMPORT
INTEGER(GLint), VALUE :: value
INTEGER(GLuint), VALUE :: time
!  void proc( int )
INTERFACE
SUBROUTINE proc(arg1) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1
END SUBROUTINE proc
END INTERFACE
END SUBROUTINE glutTimerFunc
END INTERFACE

!  void    glutIdleFunc( void (* callback)( void ) )
PUBLIC glutIdleFunc
INTERFACE
SUBROUTINE glutIdleFunc(proc) BIND(C,NAME="glutIdleFunc")
IMPORT
!  void proc( void )
INTERFACE
SUBROUTINE proc() BIND(C)
IMPORT
END SUBROUTINE proc
END INTERFACE
END SUBROUTINE glutIdleFunc
END INTERFACE


!  Window-specific callback functions, see freeglut_callbacks.c
!  void    glutKeyboardFunc( void (* callback)( unsigned char, int, int ) )
PUBLIC glutKeyboardFunc
INTERFACE
SUBROUTINE glutKeyboardFunc(proc) BIND(C,NAME="glutKeyboardFunc")
IMPORT
!  void proc( unsigned char, int, int )
INTERFACE
SUBROUTINE proc(arg1, arg2, arg3) BIND(C)
IMPORT
INTEGER(GLbyte), VALUE :: arg1
INTEGER(GLint), VALUE :: arg2, arg3
END SUBROUTINE proc
END INTERFACE
END SUBROUTINE glutKeyboardFunc
END INTERFACE

!  void    glutSpecialFunc( void (* callback)( int, int, int ) )
PUBLIC glutSpecialFunc
INTERFACE
SUBROUTINE glutSpecialFunc(proc) BIND(C,NAME="glutSpecialFunc")
IMPORT
!  void proc( int, int, int )
INTERFACE
SUBROUTINE proc(arg1, arg2, arg3) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1, arg2, arg3
END SUBROUTINE proc
END INTERFACE
END SUBROUTINE glutSpecialFunc
END INTERFACE

!  void    glutReshapeFunc( void (* callback)( int, int ) )
PUBLIC glutReshapeFunc
INTERFACE
SUBROUTINE glutReshapeFunc(proc) BIND(C,NAME="glutReshapeFunc")
IMPORT
!  void proc( int, int )
INTERFACE
SUBROUTINE proc(arg1, arg2) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1, arg2
END SUBROUTINE proc
END INTERFACE
END SUBROUTINE glutReshapeFunc
END INTERFACE

!  void    glutVisibilityFunc( void (* callback)( int ) )
PUBLIC glutVisibilityFunc
INTERFACE
SUBROUTINE glutVisibilityFunc(proc) BIND(C,NAME="glutVisibilityFunc")
IMPORT
!  void proc( int )
INTERFACE
SUBROUTINE proc(arg1) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1
END SUBROUTINE proc
END INTERFACE
END SUBROUTINE glutVisibilityFunc
END INTERFACE

!  void    glutDisplayFunc( void (* callback)( void ) )
PUBLIC glutDisplayFunc
INTERFACE
SUBROUTINE glutDisplayFunc(proc) BIND(C,NAME="glutDisplayFunc")
IMPORT
!  void proc( void )
INTERFACE
SUBROUTINE proc() BIND(C)
IMPORT
END SUBROUTINE proc
END INTERFACE
END SUBROUTINE glutDisplayFunc
END INTERFACE

!  void    glutMouseFunc( void (* callback)( int, int, int, int ) )
PUBLIC glutMouseFunc
INTERFACE
SUBROUTINE glutMouseFunc(proc) BIND(C,NAME="glutMouseFunc")
IMPORT
!  void proc( int, int, int, int )
INTERFACE
SUBROUTINE proc(arg1, arg2, arg3, arg4) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1, arg2, arg3, arg4
END SUBROUTINE proc
END INTERFACE
END SUBROUTINE glutMouseFunc
END INTERFACE

!  void    glutMotionFunc( void (* callback)( int, int ) )
PUBLIC glutMotionFunc
INTERFACE
SUBROUTINE glutMotionFunc(proc) BIND(C,NAME="glutMotionFunc")
IMPORT
!  void proc( int, int )
INTERFACE
SUBROUTINE proc(arg1, arg2) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1, arg2
END SUBROUTINE proc
END INTERFACE
END SUBROUTINE glutMotionFunc
END INTERFACE

!  void    glutPassiveMotionFunc( void (* callback)( int, int ) )
PUBLIC glutPassiveMotionFunc
INTERFACE
SUBROUTINE glutPassiveMotionFunc(proc) BIND(C,NAME="glutPassiveMotionFunc")
IMPORT
!  void proc( int, int )
INTERFACE
SUBROUTINE proc(arg1, arg2) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1, arg2
END SUBROUTINE proc
END INTERFACE
END SUBROUTINE glutPassiveMotionFunc
END INTERFACE

!  void    glutEntryFunc( void (* callback)( int ) )
PUBLIC glutEntryFunc
INTERFACE
SUBROUTINE glutEntryFunc(proc) BIND(C,NAME="glutEntryFunc")
IMPORT
!  void proc( int )
INTERFACE
SUBROUTINE proc(arg1) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1
END SUBROUTINE proc
END INTERFACE
END SUBROUTINE glutEntryFunc
END INTERFACE


!  void    glutKeyboardUpFunc( void (* callback)( unsigned char, int, int ) )
PUBLIC glutKeyboardUpFunc
INTERFACE
SUBROUTINE glutKeyboardUpFunc(proc) BIND(C,NAME="glutKeyboardUpFunc")
IMPORT
!  void proc( unsigned char, int, int )
INTERFACE
SUBROUTINE proc(arg1, arg2, arg3) BIND(C)
IMPORT
INTEGER(GLbyte), VALUE :: arg1
INTEGER(GLint), VALUE :: arg2, arg3
END SUBROUTINE proc
END INTERFACE
END SUBROUTINE glutKeyboardUpFunc
END INTERFACE

!  void    glutSpecialUpFunc( void (* callback)( int, int, int ) )
PUBLIC glutSpecialUpFunc
INTERFACE
SUBROUTINE glutSpecialUpFunc(proc) BIND(C,NAME="glutSpecialUpFunc")
IMPORT
!  void proc( int, int, int )
INTERFACE
SUBROUTINE proc(arg1, arg2, arg3) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1, arg2, arg3
END SUBROUTINE proc
END INTERFACE
END SUBROUTINE glutSpecialUpFunc
END INTERFACE

!  void    glutJoystickFunc( void (* callback)( unsigned int, int, int, int ), int pollInterval )
PUBLIC glutJoystickFunc
INTERFACE
SUBROUTINE glutJoystickFunc(proc, pollInterval) BIND(C,NAME="glutJoystickFunc")
IMPORT
INTEGER(GLint), VALUE :: pollInterval
!  void proc( unsigned int, int, int, int )
INTERFACE
SUBROUTINE proc(arg1, arg2, arg3, arg4) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg2, arg3, arg4
INTEGER(GLuint), VALUE :: arg1
END SUBROUTINE proc
END INTERFACE
END SUBROUTINE glutJoystickFunc
END INTERFACE

!  void    glutMenuStateFunc( void (* callback)( int ) )
PUBLIC glutMenuStateFunc
INTERFACE
SUBROUTINE glutMenuStateFunc(proc) BIND(C,NAME="glutMenuStateFunc")
IMPORT
!  void proc( int )
INTERFACE
SUBROUTINE proc(arg1) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1
END SUBROUTINE proc
END INTERFACE
END SUBROUTINE glutMenuStateFunc
END INTERFACE

!  void    glutMenuStatusFunc( void (* callback)( int, int, int ) )
PUBLIC glutMenuStatusFunc
INTERFACE
SUBROUTINE glutMenuStatusFunc(proc) BIND(C,NAME="glutMenuStatusFunc")
IMPORT
!  void proc( int, int, int )
INTERFACE
SUBROUTINE proc(arg1, arg2, arg3) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1, arg2, arg3
END SUBROUTINE proc
END INTERFACE
END SUBROUTINE glutMenuStatusFunc
END INTERFACE

!  void    glutOverlayDisplayFunc( void (* callback)( void ) )
PUBLIC glutOverlayDisplayFunc
INTERFACE
SUBROUTINE glutOverlayDisplayFunc(proc) BIND(C,NAME="glutOverlayDisplayFunc")
IMPORT
!  void proc( void )
INTERFACE
SUBROUTINE proc() BIND(C)
IMPORT
END SUBROUTINE proc
END INTERFACE
END SUBROUTINE glutOverlayDisplayFunc
END INTERFACE

!  void    glutWindowStatusFunc( void (* callback)( int ) )
PUBLIC glutWindowStatusFunc
INTERFACE
SUBROUTINE glutWindowStatusFunc(proc) BIND(C,NAME="glutWindowStatusFunc")
IMPORT
!  void proc( int )
INTERFACE
SUBROUTINE proc(arg1) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1
END SUBROUTINE proc
END INTERFACE
END SUBROUTINE glutWindowStatusFunc
END INTERFACE


!  void    glutSpaceballMotionFunc( void (* callback)( int, int, int ) )
PUBLIC glutSpaceballMotionFunc
INTERFACE
SUBROUTINE glutSpaceballMotionFunc(proc) BIND(C,NAME="glutSpaceballMotionFunc")
IMPORT
!  void proc( int, int, int )
INTERFACE
SUBROUTINE proc(arg1, arg2, arg3) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1, arg2, arg3
END SUBROUTINE proc
END INTERFACE
END SUBROUTINE glutSpaceballMotionFunc
END INTERFACE

!  void    glutSpaceballRotateFunc( void (* callback)( int, int, int ) )
PUBLIC glutSpaceballRotateFunc
INTERFACE
SUBROUTINE glutSpaceballRotateFunc(proc) BIND(C,NAME="glutSpaceballRotateFunc")
IMPORT
!  void proc( int, int, int )
INTERFACE
SUBROUTINE proc(arg1, arg2, arg3) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1, arg2, arg3
END SUBROUTINE proc
END INTERFACE
END SUBROUTINE glutSpaceballRotateFunc
END INTERFACE

!  void    glutSpaceballButtonFunc( void (* callback)( int, int ) )
PUBLIC glutSpaceballButtonFunc
INTERFACE
SUBROUTINE glutSpaceballButtonFunc(proc) BIND(C,NAME="glutSpaceballButtonFunc")
IMPORT
!  void proc( int, int )
INTERFACE
SUBROUTINE proc(arg1, arg2) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1, arg2
END SUBROUTINE proc
END INTERFACE
END SUBROUTINE glutSpaceballButtonFunc
END INTERFACE

!  void    glutButtonBoxFunc( void (* callback)( int, int ) )
PUBLIC glutButtonBoxFunc
INTERFACE
SUBROUTINE glutButtonBoxFunc(proc) BIND(C,NAME="glutButtonBoxFunc")
IMPORT
!  void proc( int, int )
INTERFACE
SUBROUTINE proc(arg1, arg2) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1, arg2
END SUBROUTINE proc
END INTERFACE
END SUBROUTINE glutButtonBoxFunc
END INTERFACE

!  void    glutDialsFunc( void (* callback)( int, int ) )
PUBLIC glutDialsFunc
INTERFACE
SUBROUTINE glutDialsFunc(proc) BIND(C,NAME="glutDialsFunc")
IMPORT
!  void proc( int, int )
INTERFACE
SUBROUTINE proc(arg1, arg2) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1, arg2
END SUBROUTINE proc
END INTERFACE
END SUBROUTINE glutDialsFunc
END INTERFACE

!  void    glutTabletMotionFunc( void (* callback)( int, int ) )
PUBLIC glutTabletMotionFunc
INTERFACE
SUBROUTINE glutTabletMotionFunc(proc) BIND(C,NAME="glutTabletMotionFunc")
IMPORT
!  void proc( int, int )
INTERFACE
SUBROUTINE proc(arg1, arg2) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1, arg2
END SUBROUTINE proc
END INTERFACE
END SUBROUTINE glutTabletMotionFunc
END INTERFACE

!  void    glutTabletButtonFunc( void (* callback)( int, int, int, int ) )
PUBLIC glutTabletButtonFunc
INTERFACE
SUBROUTINE glutTabletButtonFunc(proc) BIND(C,NAME="glutTabletButtonFunc")
IMPORT
!  void proc( int, int, int, int )
INTERFACE
SUBROUTINE proc(arg1, arg2, arg3, arg4) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1, arg2, arg3, arg4
END SUBROUTINE proc
END INTERFACE
END SUBROUTINE glutTabletButtonFunc
END INTERFACE


!  State setting and retrieval functions, see freeglut_state.c
!  int     glutGet( GLenum query )
PUBLIC glutGet
INTERFACE
FUNCTION glutGet(query) BIND(C,NAME="glutGet")
IMPORT
INTEGER(GLint) :: glutGet
INTEGER(GLenum), VALUE :: query
END FUNCTION glutGet
END INTERFACE

!  int     glutDeviceGet( GLenum query )
PUBLIC glutDeviceGet
INTERFACE
FUNCTION glutDeviceGet(query) BIND(C,NAME="glutDeviceGet")
IMPORT
INTEGER(GLint) :: glutDeviceGet
INTEGER(GLenum), VALUE :: query
END FUNCTION glutDeviceGet
END INTERFACE

!  int     glutGetModifiers( void )
PUBLIC glutGetModifiers
INTERFACE
FUNCTION glutGetModifiers() BIND(C,NAME="glutGetModifiers")
IMPORT
INTEGER(GLint) :: glutGetModifiers
END FUNCTION glutGetModifiers
END INTERFACE

!  int     glutLayerGet( GLenum query )
PUBLIC glutLayerGet
INTERFACE
FUNCTION glutLayerGet(query) BIND(C,NAME="glutLayerGet")
IMPORT
INTEGER(GLint) :: glutLayerGet
INTEGER(GLenum), VALUE :: query
END FUNCTION glutLayerGet
END INTERFACE


!  Font stuff, see freeglut_font.c
!  void    glutBitmapCharacter( void* font, int character )
PUBLIC glutBitmapCharacter
INTERFACE
SUBROUTINE glutBitmapCharacter(font, character) BIND(C,NAME="glutBitmapCharacter")
IMPORT
INTEGER(GLint), VALUE :: character
TYPE(C_PTR), VALUE :: font
END SUBROUTINE glutBitmapCharacter
END INTERFACE

!  int     glutBitmapWidth( void* font, int character )
PUBLIC glutBitmapWidth
INTERFACE
FUNCTION glutBitmapWidth(font, character) BIND(C,NAME="glutBitmapWidth")
IMPORT
INTEGER(GLint) :: glutBitmapWidth
INTEGER(GLint), VALUE :: character
TYPE(C_PTR), VALUE :: font
END FUNCTION glutBitmapWidth
END INTERFACE

!  void    glutStrokeCharacter( void* font, int character )
PUBLIC glutStrokeCharacter
INTERFACE
SUBROUTINE glutStrokeCharacter(font, character) BIND(C,NAME="glutStrokeCharacter")
IMPORT
INTEGER(GLint), VALUE :: character
TYPE(C_PTR), VALUE :: font
END SUBROUTINE glutStrokeCharacter
END INTERFACE

!  int     glutStrokeWidth( void* font, int character )
PUBLIC glutStrokeWidth
INTERFACE
FUNCTION glutStrokeWidth(font, character) BIND(C,NAME="glutStrokeWidth")
IMPORT
INTEGER(GLint) :: glutStrokeWidth
INTEGER(GLint), VALUE :: character
TYPE(C_PTR), VALUE :: font
END FUNCTION glutStrokeWidth
END INTERFACE

!  int     glutBitmapLength( void* font, const unsigned char* string )
PUBLIC glutBitmapLength
INTERFACE
FUNCTION glutBitmapLength(font, string) BIND(C,NAME="glutBitmapLength")
IMPORT
INTEGER(GLint) :: glutBitmapLength
! CHARACTER(C_CHAR), INTENT(IN) :: string
CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: string
TYPE(C_PTR), VALUE :: font
END FUNCTION glutBitmapLength
END INTERFACE

!  int     glutStrokeLength( void* font, const unsigned char* string )
PUBLIC glutStrokeLength
INTERFACE
FUNCTION glutStrokeLength(font, string) BIND(C,NAME="glutStrokeLength")
IMPORT
INTEGER(GLint) :: glutStrokeLength
! CHARACTER(C_CHAR), INTENT(IN) :: string
CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: string
TYPE(C_PTR), VALUE :: font
END FUNCTION glutStrokeLength
END INTERFACE


!  Geometry functions, see freeglut_geometry.c
!  void    glutWireCube( GLdouble size )
PUBLIC glutWireCube
INTERFACE
SUBROUTINE glutWireCube(size) BIND(C,NAME="glutWireCube")
IMPORT
REAL(GLdouble), VALUE :: size
END SUBROUTINE glutWireCube
END INTERFACE

!  void    glutSolidCube( GLdouble size )
PUBLIC glutSolidCube
INTERFACE
SUBROUTINE glutSolidCube(size) BIND(C,NAME="glutSolidCube")
IMPORT
REAL(GLdouble), VALUE :: size
END SUBROUTINE glutSolidCube
END INTERFACE

!  void    glutWireSphere( GLdouble radius, GLint slices, GLint stacks )
PUBLIC glutWireSphere
INTERFACE
SUBROUTINE glutWireSphere(radius, slices, stacks) BIND(C,NAME="glutWireSphere")
IMPORT
INTEGER(GLint), VALUE :: slices, stacks
REAL(GLdouble), VALUE :: radius
END SUBROUTINE glutWireSphere
END INTERFACE

!  void    glutSolidSphere( GLdouble radius, GLint slices, GLint stacks )
PUBLIC glutSolidSphere
INTERFACE
SUBROUTINE glutSolidSphere(radius, slices, stacks) BIND(C,NAME="glutSolidSphere")
IMPORT
INTEGER(GLint), VALUE :: slices, stacks
REAL(GLdouble), VALUE :: radius
END SUBROUTINE glutSolidSphere
END INTERFACE

!  void    glutWireCone( GLdouble base, GLdouble height, GLint slices, GLint stacks )
PUBLIC glutWireCone
INTERFACE
SUBROUTINE glutWireCone(base, height, slices, stacks) BIND(C,NAME="glutWireCone")
IMPORT
INTEGER(GLint), VALUE :: slices, stacks
REAL(GLdouble), VALUE :: base, height
END SUBROUTINE glutWireCone
END INTERFACE

!  void    glutSolidCone( GLdouble base, GLdouble height, GLint slices, GLint stacks )
PUBLIC glutSolidCone
INTERFACE
SUBROUTINE glutSolidCone(base, height, slices, stacks) BIND(C,NAME="glutSolidCone")
IMPORT
INTEGER(GLint), VALUE :: slices, stacks
REAL(GLdouble), VALUE :: base, height
END SUBROUTINE glutSolidCone
END INTERFACE


!  void    glutWireTorus( GLdouble innerRadius, GLdouble outerRadius, GLint sides, GLint rings )
PUBLIC glutWireTorus
INTERFACE
SUBROUTINE glutWireTorus(innerRadius, outerRadius, sides, rings) BIND(C,NAME="glutWireTorus")
IMPORT
INTEGER(GLint), VALUE :: sides, rings
REAL(GLdouble), VALUE :: innerRadius, outerRadius
END SUBROUTINE glutWireTorus
END INTERFACE

!  void    glutSolidTorus( GLdouble innerRadius, GLdouble outerRadius, GLint sides, GLint rings )
PUBLIC glutSolidTorus
INTERFACE
SUBROUTINE glutSolidTorus(innerRadius, outerRadius, sides, rings) BIND(C,NAME="glutSolidTorus")
IMPORT
INTEGER(GLint), VALUE :: sides, rings
REAL(GLdouble), VALUE :: innerRadius, outerRadius
END SUBROUTINE glutSolidTorus
END INTERFACE

!  void    glutWireDodecahedron( void )
PUBLIC glutWireDodecahedron
INTERFACE
SUBROUTINE glutWireDodecahedron() BIND(C,NAME="glutWireDodecahedron")
IMPORT
END SUBROUTINE glutWireDodecahedron
END INTERFACE

!  void    glutSolidDodecahedron( void )
PUBLIC glutSolidDodecahedron
INTERFACE
SUBROUTINE glutSolidDodecahedron() BIND(C,NAME="glutSolidDodecahedron")
IMPORT
END SUBROUTINE glutSolidDodecahedron
END INTERFACE

!  void    glutWireOctahedron( void )
PUBLIC glutWireOctahedron
INTERFACE
SUBROUTINE glutWireOctahedron() BIND(C,NAME="glutWireOctahedron")
IMPORT
END SUBROUTINE glutWireOctahedron
END INTERFACE

!  void    glutSolidOctahedron( void )
PUBLIC glutSolidOctahedron
INTERFACE
SUBROUTINE glutSolidOctahedron() BIND(C,NAME="glutSolidOctahedron")
IMPORT
END SUBROUTINE glutSolidOctahedron
END INTERFACE

!  void    glutWireTetrahedron( void )
PUBLIC glutWireTetrahedron
INTERFACE
SUBROUTINE glutWireTetrahedron() BIND(C,NAME="glutWireTetrahedron")
IMPORT
END SUBROUTINE glutWireTetrahedron
END INTERFACE

!  void    glutSolidTetrahedron( void )
PUBLIC glutSolidTetrahedron
INTERFACE
SUBROUTINE glutSolidTetrahedron() BIND(C,NAME="glutSolidTetrahedron")
IMPORT
END SUBROUTINE glutSolidTetrahedron
END INTERFACE

!  void    glutWireIcosahedron( void )
PUBLIC glutWireIcosahedron
INTERFACE
SUBROUTINE glutWireIcosahedron() BIND(C,NAME="glutWireIcosahedron")
IMPORT
END SUBROUTINE glutWireIcosahedron
END INTERFACE

!  void    glutSolidIcosahedron( void )
PUBLIC glutSolidIcosahedron
INTERFACE
SUBROUTINE glutSolidIcosahedron() BIND(C,NAME="glutSolidIcosahedron")
IMPORT
END SUBROUTINE glutSolidIcosahedron
END INTERFACE


!  Teapot rendering functions, found in freeglut_teapot.c
!  void    glutWireTeapot( GLdouble size )
PUBLIC glutWireTeapot
INTERFACE
SUBROUTINE glutWireTeapot(size) BIND(C,NAME="glutWireTeapot")
IMPORT
REAL(GLdouble), VALUE :: size
END SUBROUTINE glutWireTeapot
END INTERFACE

!  void    glutSolidTeapot( GLdouble size )
PUBLIC glutSolidTeapot
INTERFACE
SUBROUTINE glutSolidTeapot(size) BIND(C,NAME="glutSolidTeapot")
IMPORT
REAL(GLdouble), VALUE :: size
END SUBROUTINE glutSolidTeapot
END INTERFACE


!  Game mode functions, see freeglut_gamemode.c
!  void    glutGameModeString( const char* string )
PUBLIC glutGameModeString
INTERFACE
SUBROUTINE glutGameModeString(string) BIND(C,NAME="glutGameModeString")
IMPORT
! CHARACTER(C_CHAR), INTENT(IN) :: string
CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: string
END SUBROUTINE glutGameModeString
END INTERFACE

!  int     glutEnterGameMode( void )
PUBLIC glutEnterGameMode
INTERFACE
FUNCTION glutEnterGameMode() BIND(C,NAME="glutEnterGameMode")
IMPORT
INTEGER(GLint) :: glutEnterGameMode
END FUNCTION glutEnterGameMode
END INTERFACE

!  void    glutLeaveGameMode( void )
PUBLIC glutLeaveGameMode
INTERFACE
SUBROUTINE glutLeaveGameMode() BIND(C,NAME="glutLeaveGameMode")
IMPORT
END SUBROUTINE glutLeaveGameMode
END INTERFACE

!  int     glutGameModeGet( GLenum query )
PUBLIC glutGameModeGet
INTERFACE
FUNCTION glutGameModeGet(query) BIND(C,NAME="glutGameModeGet")
IMPORT
INTEGER(GLint) :: glutGameModeGet
INTEGER(GLenum), VALUE :: query
END FUNCTION glutGameModeGet
END INTERFACE


!  Video resize functions, see freeglut_videoresize.c
!  int     glutVideoResizeGet( GLenum query )
PUBLIC glutVideoResizeGet
INTERFACE
FUNCTION glutVideoResizeGet(query) BIND(C,NAME="glutVideoResizeGet")
IMPORT
INTEGER(GLint) :: glutVideoResizeGet
INTEGER(GLenum), VALUE :: query
END FUNCTION glutVideoResizeGet
END INTERFACE

!  void    glutSetupVideoResizing( void )
PUBLIC glutSetupVideoResizing
INTERFACE
SUBROUTINE glutSetupVideoResizing() BIND(C,NAME="glutSetupVideoResizing")
IMPORT
END SUBROUTINE glutSetupVideoResizing
END INTERFACE

!  void    glutStopVideoResizing( void )
PUBLIC glutStopVideoResizing
INTERFACE
SUBROUTINE glutStopVideoResizing() BIND(C,NAME="glutStopVideoResizing")
IMPORT
END SUBROUTINE glutStopVideoResizing
END INTERFACE

!  void    glutVideoResize( int x, int y, int width, int height )
PUBLIC glutVideoResize
INTERFACE
SUBROUTINE glutVideoResize(x, y, width, height) BIND(C,NAME="glutVideoResize")
IMPORT
INTEGER(GLint), VALUE :: x, y, width, height
END SUBROUTINE glutVideoResize
END INTERFACE

!  void    glutVideoPan( int x, int y, int width, int height )
PUBLIC glutVideoPan
INTERFACE
SUBROUTINE glutVideoPan(x, y, width, height) BIND(C,NAME="glutVideoPan")
IMPORT
INTEGER(GLint), VALUE :: x, y, width, height
END SUBROUTINE glutVideoPan
END INTERFACE


!  Colormap functions, see freeglut_misc.c
!  void    glutSetColor( int color, GLfloat red, GLfloat green, GLfloat blue )
PUBLIC glutSetColor
INTERFACE
SUBROUTINE glutSetColor(color, red, green, blue) BIND(C,NAME="glutSetColor")
IMPORT
INTEGER(GLint), VALUE :: color
REAL(GLfloat), VALUE :: red, green, blue
END SUBROUTINE glutSetColor
END INTERFACE

!  GLfloat glutGetColor( int color, int component )
PUBLIC glutGetColor
INTERFACE
FUNCTION glutGetColor(color, component) BIND(C,NAME="glutGetColor")
IMPORT
REAL(GLfloat) :: glutGetColor
INTEGER(GLint), VALUE :: color, component
END FUNCTION glutGetColor
END INTERFACE

!  void    glutCopyColormap( int window )
PUBLIC glutCopyColormap
INTERFACE
SUBROUTINE glutCopyColormap(window) BIND(C,NAME="glutCopyColormap")
IMPORT
INTEGER(GLint), VALUE :: window
END SUBROUTINE glutCopyColormap
END INTERFACE


!  Misc keyboard and joystick functions, see freeglut_misc.c
!  void    glutIgnoreKeyRepeat( int ignore )
PUBLIC glutIgnoreKeyRepeat
INTERFACE
SUBROUTINE glutIgnoreKeyRepeat(ignore) BIND(C,NAME="glutIgnoreKeyRepeat")
IMPORT
INTEGER(GLint), VALUE :: ignore
END SUBROUTINE glutIgnoreKeyRepeat
END INTERFACE

!  void    glutSetKeyRepeat( int repeatMode )
PUBLIC glutSetKeyRepeat
INTERFACE
SUBROUTINE glutSetKeyRepeat(repeatMode) BIND(C,NAME="glutSetKeyRepeat")
IMPORT
INTEGER(GLint), VALUE :: repeatMode
END SUBROUTINE glutSetKeyRepeat
END INTERFACE

!  void    glutForceJoystickFunc( void )
PUBLIC glutForceJoystickFunc
INTERFACE
SUBROUTINE glutForceJoystickFunc() BIND(C,NAME="glutForceJoystickFunc")
IMPORT
END SUBROUTINE glutForceJoystickFunc
END INTERFACE


!  Misc functions, see freeglut_misc.c
!  int     glutExtensionSupported( const char* extension )
PUBLIC glutExtensionSupported
INTERFACE
FUNCTION glutExtensionSupported(extension) BIND(C,NAME="glutExtensionSupported")
IMPORT
INTEGER(GLint) :: glutExtensionSupported
! CHARACTER(C_CHAR), INTENT(IN) :: extension
CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: extension
END FUNCTION glutExtensionSupported
END INTERFACE

!  void    glutReportErrors( void )
PUBLIC glutReportErrors
INTERFACE
SUBROUTINE glutReportErrors() BIND(C,NAME="glutReportErrors")
IMPORT
END SUBROUTINE glutReportErrors
END INTERFACE


! #if defined(_WIN32) && defined(GLUT_DISABLE_ATEXIT_HACK) && defined(__WATCOMC__)
!  
!     Win32 has an annoying issue where there are multiple C run-time
!     libraries (CRTs).  If the executable is linked with a different CRT
!     from the GLUT DLL, the GLUT DLL will not share the same CRT static
!     data seen by the executable.  In particular, atexit callbacks registered
!     in the executable will not be called if GLUT calls its (different)
!     exit routine).  GLUT is typically built with the
!     "/MD" option (the CRT with multithreading DLL support), but the Visual
!     C++ linker default is "/ML" (the single threaded CRT).
!  
!     One workaround to this issue is requiring users to always link with
!     the same CRT as GLUT is compiled with.  That requires users supply a
!     non-standard option.  GLUT 3.7 has its own built-in workaround where
!     the executable's "exit" function pointer is covertly passed to GLUT.
!     GLUT then calls the executable's exit function pointer to ensure that
!     any "atexit" calls registered by the application are called if GLUT
!     needs to exit.
!  
!     Note that the __glut*WithExit routines should NEVER be called directly.

!  to get the prototype for exit() 

!  Interfaces omitted for the following:
!  void __glutInitWithExit(int *argcp, char **argv, void (__cdecl *exitfunc)(int))
!  int __glutCreateWindowWithExit(const char *title, void (__cdecl *exitfunc)(int))
!  int __glutCreateMenuWithExit(void (* func)(int), void (__cdecl *exitfunc)(int))
! ???, , PUBLIC :: FGUNUSED=__attribute__((unused))


!  ** END OF FILE **


!  ** End of freeglut_std.h **

!  ** #include "freeglut_ext.h" **

!  freeglut_ext.h
!  
!  The non-GLUT-compatible extensions to the freeglut library include file
!  
!  Copyright (c) 1999-2000 Pawel W. Olszta. All Rights Reserved.
!  Written by Pawel W. Olszta, <olszta@sourceforge.net>
!  Creation date: Thu Dec 2 1999
!  
!  Permission is hereby granted, free of charge, to any person obtaining a
!  copy of this software and associated documentation files (the "Software"),
!  to deal in the Software without restriction, including without limitation
!  the rights to use, copy, modify, merge, publish, distribute, sublicense,
!  and/or sell copies of the Software, and to permit persons to whom the
!  Software is furnished to do so, subject to the following conditions:
!  
!  The above copyright notice and this permission notice shall be included
!  in all copies or substantial portions of the Software.
!  
!  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
!  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
!  PAWEL W. OLSZTA BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
!  IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
!  CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


!  Additional GLUT Key definitions for the Special key function
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_NUM_LOCK = int(z'006D', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_BEGIN = int(z'006E', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_DELETE = int(z'006F', kind=GLenum)

!  GLUT API Extension macro definitions -- behaviour when the user clicks on an "x" to close a window
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_ACTION_EXIT = 0
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_ACTION_GLUTMAINLOOP_RETURNS = 1
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_ACTION_CONTINUE_EXECUTION = 2

!  Create a new rendering context when the user opens a new window?
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CREATE_NEW_CONTEXT = 0
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_USE_CURRENT_CONTEXT = 1

!  Direct/Indirect rendering context options (has meaning only in Unix/X11)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_FORCE_INDIRECT_CONTEXT = 0
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_ALLOW_DIRECT_CONTEXT = 1
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_TRY_DIRECT_CONTEXT = 2
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_FORCE_DIRECT_CONTEXT = 3

!  GLUT API Extension macro definitions -- the glutGet parameters
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_INIT_STATE = int(z'007C', kind=GLenum)

INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_ACTION_ON_WINDOW_CLOSE = int(z'01F9', kind=GLenum)

INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_BORDER_WIDTH = int(z'01FA', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_HEADER_HEIGHT = int(z'01FB', kind=GLenum)

INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_VERSION = int(z'01FC', kind=GLenum)

INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_RENDERING_CONTEXT = int(z'01FD', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_DIRECT_RENDERING = int(z'01FE', kind=GLenum)

INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_FULL_SCREEN = int(z'01FF', kind=GLenum)

!  New tokens for glutInitDisplayMode.
!  Only one GLUT_AUXn bit may be used at a time.
!  Value 0x0400 is defined in OpenGLUT.
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_AUX = int(z'1000', kind=GLenum)

INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_AUX1 = int(z'1000', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_AUX2 = int(z'2000', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_AUX3 = int(z'4000', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_AUX4 = int(z'8000', kind=GLenum)

!  Context-related flags, see freeglut_state.c
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_INIT_MAJOR_VERSION = int(z'0200', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_INIT_MINOR_VERSION = int(z'0201', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_INIT_FLAGS = int(z'0202', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_INIT_PROFILE = int(z'0203', kind=GLenum)

!  Flags for glutInitContextFlags, see freeglut_init.c
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_DEBUG = int(z'0001', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_FORWARD_COMPATIBLE = int(z'0002', kind=GLenum)


!  Flags for glutInitContextProfile, see freeglut_init.c
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CORE_PROFILE = int(z'0001', kind=GLenum)

!  Process loop function, see freeglut_main.c
!  void    glutMainLoopEvent( void )
PUBLIC glutMainLoopEvent
INTERFACE
SUBROUTINE glutMainLoopEvent() BIND(C,NAME="glutMainLoopEvent")
IMPORT
END SUBROUTINE glutMainLoopEvent
END INTERFACE

!  void    glutLeaveMainLoop( void )
PUBLIC glutLeaveMainLoop
INTERFACE
SUBROUTINE glutLeaveMainLoop() BIND(C,NAME="glutLeaveMainLoop")
IMPORT
END SUBROUTINE glutLeaveMainLoop
END INTERFACE

!  void    glutExit         ( void )
PUBLIC glutExit
INTERFACE
SUBROUTINE glutExit() BIND(C,NAME="glutExit")
IMPORT
END SUBROUTINE glutExit
END INTERFACE


!  Window management functions, see freeglut_window.c
!  void    glutFullScreenToggle( void )
PUBLIC glutFullScreenToggle
INTERFACE
SUBROUTINE glutFullScreenToggle() BIND(C,NAME="glutFullScreenToggle")
IMPORT
END SUBROUTINE glutFullScreenToggle
END INTERFACE


!  Window-specific callback functions, see freeglut_callbacks.c
!  void    glutMouseWheelFunc( void (* callback)( int, int, int, int ) )
PUBLIC glutMouseWheelFunc
INTERFACE
SUBROUTINE glutMouseWheelFunc(proc) BIND(C,NAME="glutMouseWheelFunc")
IMPORT
!  void proc( int, int, int, int )
INTERFACE
SUBROUTINE proc(arg1, arg2, arg3, arg4) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1, arg2, arg3, arg4
END SUBROUTINE proc
END INTERFACE
END SUBROUTINE glutMouseWheelFunc
END INTERFACE

!  void    glutCloseFunc( void (* callback)( void ) )
PUBLIC glutCloseFunc
INTERFACE
SUBROUTINE glutCloseFunc(proc) BIND(C,NAME="glutCloseFunc")
IMPORT
!  void proc( void )
INTERFACE
SUBROUTINE proc() BIND(C)
IMPORT
END SUBROUTINE proc
END INTERFACE
END SUBROUTINE glutCloseFunc
END INTERFACE

!  void    glutWMCloseFunc( void (* callback)( void ) )
PUBLIC glutWMCloseFunc
INTERFACE
SUBROUTINE glutWMCloseFunc(proc) BIND(C,NAME="glutWMCloseFunc")
IMPORT
!  void proc( void )
INTERFACE
SUBROUTINE proc() BIND(C)
IMPORT
END SUBROUTINE proc
END INTERFACE
END SUBROUTINE glutWMCloseFunc
END INTERFACE

!  A. Donev: Also a destruction callback for menus 
!  void    glutMenuDestroyFunc( void (* callback)( void ) )
PUBLIC glutMenuDestroyFunc
INTERFACE
SUBROUTINE glutMenuDestroyFunc(proc) BIND(C,NAME="glutMenuDestroyFunc")
IMPORT
!  void proc( void )
INTERFACE
SUBROUTINE proc() BIND(C)
IMPORT
END SUBROUTINE proc
END INTERFACE
END SUBROUTINE glutMenuDestroyFunc
END INTERFACE


!  State setting and retrieval functions, see freeglut_state.c
!  void    glutSetOption ( GLenum option_flag, int value )
PUBLIC glutSetOption
INTERFACE
SUBROUTINE glutSetOption(option_flag, value) BIND(C,NAME="glutSetOption")
IMPORT
INTEGER(GLenum), VALUE :: option_flag
INTEGER(GLint), VALUE :: value
END SUBROUTINE glutSetOption
END INTERFACE

!  int *   glutGetModeValues(GLenum mode, int * size)
PUBLIC glutGetModeValues
INTERFACE
FUNCTION glutGetModeValues(mode, size) BIND(C,NAME="glutGetModeValues")
IMPORT
!  int* :: glutGetModeValues
INTEGER(GLenum), VALUE :: mode
! INTEGER(GLint) :: size
INTEGER(GLint), DIMENSION(*) :: size
END FUNCTION glutGetModeValues
END INTERFACE

!  A.Donev: User-data manipulation 
!  void*   glutGetWindowData( void )
PUBLIC glutGetWindowData
INTERFACE
FUNCTION glutGetWindowData() BIND(C,NAME="glutGetWindowData")
IMPORT
TYPE(C_PTR) :: glutGetWindowData
END FUNCTION glutGetWindowData
END INTERFACE

!  void    glutSetWindowData(void* data)
PUBLIC glutSetWindowData
INTERFACE
SUBROUTINE glutSetWindowData(data) BIND(C,NAME="glutSetWindowData")
IMPORT
TYPE(C_PTR), VALUE :: data
END SUBROUTINE glutSetWindowData
END INTERFACE

!  void*   glutGetMenuData( void )
PUBLIC glutGetMenuData
INTERFACE
FUNCTION glutGetMenuData() BIND(C,NAME="glutGetMenuData")
IMPORT
TYPE(C_PTR) :: glutGetMenuData
END FUNCTION glutGetMenuData
END INTERFACE

!  void    glutSetMenuData(void* data)
PUBLIC glutSetMenuData
INTERFACE
SUBROUTINE glutSetMenuData(data) BIND(C,NAME="glutSetMenuData")
IMPORT
TYPE(C_PTR), VALUE :: data
END SUBROUTINE glutSetMenuData
END INTERFACE


!  Font stuff, see freeglut_font.c
!  int     glutBitmapHeight( void* font )
PUBLIC glutBitmapHeight
INTERFACE
FUNCTION glutBitmapHeight(font) BIND(C,NAME="glutBitmapHeight")
IMPORT
INTEGER(GLint) :: glutBitmapHeight
TYPE(C_PTR), VALUE :: font
END FUNCTION glutBitmapHeight
END INTERFACE

!  GLfloat glutStrokeHeight( void* font )
PUBLIC glutStrokeHeight
INTERFACE
FUNCTION glutStrokeHeight(font) BIND(C,NAME="glutStrokeHeight")
IMPORT
REAL(GLfloat) :: glutStrokeHeight
TYPE(C_PTR), VALUE :: font
END FUNCTION glutStrokeHeight
END INTERFACE

!  void    glutBitmapString( void* font, const unsigned char *string )
PUBLIC glutBitmapString
INTERFACE
SUBROUTINE glutBitmapString(font, string) BIND(C,NAME="glutBitmapString")
IMPORT
! CHARACTER(C_CHAR), INTENT(IN) :: string
CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: string
TYPE(C_PTR), VALUE :: font
END SUBROUTINE glutBitmapString
END INTERFACE

!  void    glutStrokeString( void* font, const unsigned char *string )
PUBLIC glutStrokeString
INTERFACE
SUBROUTINE glutStrokeString(font, string) BIND(C,NAME="glutStrokeString")
IMPORT
! CHARACTER(C_CHAR), INTENT(IN) :: string
CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: string
TYPE(C_PTR), VALUE :: font
END SUBROUTINE glutStrokeString
END INTERFACE


!  Geometry functions, see freeglut_geometry.c
!  void    glutWireRhombicDodecahedron( void )
PUBLIC glutWireRhombicDodecahedron
INTERFACE
SUBROUTINE glutWireRhombicDodecahedron() BIND(C,NAME="glutWireRhombicDodecahedron")
IMPORT
END SUBROUTINE glutWireRhombicDodecahedron
END INTERFACE

!  void    glutSolidRhombicDodecahedron( void )
PUBLIC glutSolidRhombicDodecahedron
INTERFACE
SUBROUTINE glutSolidRhombicDodecahedron() BIND(C,NAME="glutSolidRhombicDodecahedron")
IMPORT
END SUBROUTINE glutSolidRhombicDodecahedron
END INTERFACE

!  void    glutWireSierpinskiSponge ( int num_levels, GLdouble offset[3], GLdouble scale )
PUBLIC glutWireSierpinskiSponge
INTERFACE
SUBROUTINE glutWireSierpinskiSponge(num_levels, offset, scale) BIND(C,NAME="glutWireSierpinskiSponge")
IMPORT
INTEGER(GLint), VALUE :: num_levels
REAL(GLdouble), DIMENSION(3) :: offset
REAL(GLdouble), VALUE :: scale
END SUBROUTINE glutWireSierpinskiSponge
END INTERFACE

!  void    glutSolidSierpinskiSponge ( int num_levels, GLdouble offset[3], GLdouble scale )
PUBLIC glutSolidSierpinskiSponge
INTERFACE
SUBROUTINE glutSolidSierpinskiSponge(num_levels, offset, scale) BIND(C,NAME="glutSolidSierpinskiSponge")
IMPORT
INTEGER(GLint), VALUE :: num_levels
REAL(GLdouble), DIMENSION(3) :: offset
REAL(GLdouble), VALUE :: scale
END SUBROUTINE glutSolidSierpinskiSponge
END INTERFACE

!  void    glutWireCylinder( GLdouble radius, GLdouble height, GLint slices, GLint stacks)
PUBLIC glutWireCylinder
INTERFACE
SUBROUTINE glutWireCylinder(radius, height, slices, stacks) BIND(C,NAME="glutWireCylinder")
IMPORT
INTEGER(GLint), VALUE :: slices, stacks
REAL(GLdouble), VALUE :: radius, height
END SUBROUTINE glutWireCylinder
END INTERFACE

!  void    glutSolidCylinder( GLdouble radius, GLdouble height, GLint slices, GLint stacks)
PUBLIC glutSolidCylinder
INTERFACE
SUBROUTINE glutSolidCylinder(radius, height, slices, stacks) BIND(C,NAME="glutSolidCylinder")
IMPORT
INTEGER(GLint), VALUE :: slices, stacks
REAL(GLdouble), VALUE :: radius, height
END SUBROUTINE glutSolidCylinder
END INTERFACE


!  Extension functions, see freeglut_ext.c
!  GLUTproc glutGetProcAddress( const char *procName )
PUBLIC glutGetProcAddress
INTERFACE
FUNCTION glutGetProcAddress(procName) BIND(C,NAME="glutGetProcAddress")
IMPORT
TYPE(C_FUNPTR) :: glutGetProcAddress
! CHARACTER(C_CHAR), INTENT(IN) :: procName
CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: procName
END FUNCTION glutGetProcAddress
END INTERFACE


!  Joystick functions, see freeglut_joystick.c
!  USE OF THESE FUNCTIONS IS DEPRECATED !!!!! 
!  contact the "freeglut" developer community at freeglut-developer@lists.sourceforge.net,
!  switch to the OpenGLUT library, or else port your joystick functionality over to PLIB's
!  "js" library.

!  Initialization functions, see freeglut_init.c
!  void    glutInitContextVersion( int majorVersion, int minorVersion )
PUBLIC glutInitContextVersion
INTERFACE
SUBROUTINE glutInitContextVersion(majorVersion, minorVersion) BIND(C,NAME="glutInitContextVersion")
IMPORT
INTEGER(GLint), VALUE :: majorVersion, minorVersion
END SUBROUTINE glutInitContextVersion
END INTERFACE

!  void    glutInitContextFlags( int flags )
PUBLIC glutInitContextFlags
INTERFACE
SUBROUTINE glutInitContextFlags(flags) BIND(C,NAME="glutInitContextFlags")
IMPORT
INTEGER(GLint), VALUE :: flags
END SUBROUTINE glutInitContextFlags
END INTERFACE

!  void    glutInitContextProfile( int profile )
PUBLIC glutInitContextProfile
INTERFACE
SUBROUTINE glutInitContextProfile(profile) BIND(C,NAME="glutInitContextProfile")
IMPORT
INTEGER(GLint), VALUE :: profile
END SUBROUTINE glutInitContextProfile
END INTERFACE


!  GLUT API macro definitions -- the display mode definitions
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CAPTIONLESS = int(z'0400', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_BORDERLESS = int(z'0800', kind=GLenum)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_SRGB = int(z'1000', kind=GLenum)


!  ** END OF FILE **


!  ** End of freeglut_ext.h **

!  ** END OF FILE **


! Font variables in GLUT_fonts.c
TYPE(C_PTR), BIND(C), PUBLIC, PROTECTED :: GLUT_STROKE_ROMAN,         &
    GLUT_STROKE_MONO_ROMAN, GLUT_BITMAP_9_BY_15, GLUT_BITMAP_8_BY_13, &
    GLUT_BITMAP_TIMES_ROMAN_10, GLUT_BITMAP_TIMES_ROMAN_24,           &
    GLUT_BITMAP_HELVETICA_10, GLUT_BITMAP_HELVETICA_12,               &
    GLUT_BITMAP_HELVETICA_18

! A special callback function
EXTERNAL GLUT_NULL_FUNC
PUBLIC GLUT_NULL_FUNC

CONTAINS

SUBROUTINE glutInit_f03()
  INTEGER(C_INT), DIMENSION(1) :: argcp=1
  TYPE(C_PTR), DIMENSION(1), TARGET :: argv=C_NULL_PTR
  CHARACTER(C_CHAR), DIMENSION(1), TARGET :: empty_string=C_NULL_CHAR

  argv(1)=C_LOC(empty_string)
  CALL glutInit_gl(argcp, C_LOC(argv))

END SUBROUTINE
END MODULE OpenGL_glut

SUBROUTINE GLUT_NULL_FUNC() ! Dummy callback, not really legal
END SUBROUTINE
