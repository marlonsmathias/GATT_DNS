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

    Type( FIPC_ctxt ) :: ctxt

    Type( segment_list_type ), Pointer :: p

    ! Check things are initialized
    If( .Not. Associated( FIPC_ctxt_world%initialized ) ) Then
       error = FIPC_not_initialized
       Return
    End If

    If( .Not. Associated( a ) ) Then
       error = FIPC_free_null_pointer
       Return
    End If

    ! First find the segment
    ! Be careful to only look at segments that could possible be associated
    ! with this pointer
    p => seg_list
    Do While( Associated( p ) )
       If( p%data%type == type ) Then
          If( Size( p%data%sizes ) == Size( Shape( a ) ) ) Then
             If( All( p%data%sizes == Shape( a ) ) ) Then
                Call c_f_pointer( p%data%shmaddr, fp, p%data%sizes )
                If( Associated( a, fp ) ) Then
                   Exit
                End If
             End If
          End If
       End If
       p => p%next
    End Do

    If( .Not. Associated( p ) ) Then
       error = FIPC_seg_not_found
       Return
    End If

    If( debug ) Then
       ctxt = p%data%ctxt
    End If

    Call segment_free( p%data%shmid, error )
    If( error /= FIPC_success ) Then
       Return
    End If

    a => Null()

    If( debug ) Then
       Call print_seg_list( ctxt )
    End If

    error = FIPC_success

