! FreeIPC - A library to allow Fortran programs that use MPI
!           to easily access shared memory within multicore nodes
! Date - 23rd Oct Version 0.0
! Copyright(C) 2009 Ian Bush
!
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation; either
!    version 3.0 of the License, or (at your option) any later version.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
! This library is released under LGPL. The full text of LGPL may be found
! at http://www.gnu.org/licenses/lgpl-3.0.txt
!
! In particular, the terms of the GNU LGPL imply that you CAN use this library with
! proprietary code with some restrictions (you can link, but not use the code directly). 
! Please consult the above URL for exact details.

    Type( c_ptr ) :: shmaddr

    ! Check things are initialized
    If( .Not. Associated( FIPC_ctxt_world%initialized ) ) Then
       error = FIPC_not_initialized
       Return
    End If

    If( Size( n ) < rank ) Then
       error = FIPC_insufficient_dims
    End If

    Call segment_create( type, rank, Int( n( 1:rank ), c_long ), ctxt, shmaddr, error )
    If( error /= 0 ) Then
       Return
    End If

    Call c_f_pointer( shmaddr, a, Int( n( 1:rank ), c_long ) )

    error = FIPC_success

