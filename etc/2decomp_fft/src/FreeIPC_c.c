/*

FreeIPC - A library to allow Fortran programs that use MPI
          to easily access shared memory within multicore nodes
Date - 23rd Oct Version 0.0
Copyright(C) 2009 Ian Bush

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 3.0 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


This library is released under LGPL. The full text of LGPL may be found
at http://www.gnu.org/licenses/lgpl-3.0.txt

In particular, the terms of the GNU LGPL imply that you CAN use this library with
proprietary code with some restrictions (you can link, but not use the code directly). 
Please consult the above URL for exact details.

*/

/*
  The shared memory facilities are provide by use of the standard System V
  IPC facilities. See 
  
  http://www.opengroup.org/onlinepubs/009695399/functions/xsh_chap02_07.html#tag_02_07
  
  for details
*/


/* These routines are mostly very thin wrappers for System V IPC routines.
The main work occurs in the Fortran. The reason these wrappers are being used
( instead we might call the routines direct from Fortran ) is that Fortran 
interop with C does define what you need to use to be interoperable with such
things like size_t, key_t, pid_t etc. Therefore call the wrappers using
basic int and long and then let the prototypes in scope weave their magic
to make sure we get the correct types */

/* Apparently this is required to stop some warning from the gnu compiler */
#if !defined __USE_SVID && !defined __USE_XOPEN && __GNUC__ >= 2
#define _SVID_SOURCE
#define _XOPEN_SOURCE 600
#endif

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include <sys/types.h>
#include <sys/shm.h>
#include <sys/sem.h>

/* Stuff for setting up segs and sems */
#define N_BASE_FILES 8
static char *base_files[ N_BASE_FILES ] = { 
     "/tmp", "/usr", "/bin", "/lib", "/etc", "/dev" , "/sbin", "/mnt" };
static int next_proj = -1;
static int next_file =  0;


void  fipc_get_errval      ( int *eacces, int *eexist, int *einval, int *enfile, 
                             int *enoent, int *enomem, int *enospc );
int   fipc_get_seg         ( long size, int create, int exclusive, int perms );
void *fipc_attach_seg      ( int shmid, int flag );
int   fipc_detach_seg      ( void *shmaddr );
int   fipc_remove_seg      ( int shmid );
int   fipc_get_sem         ( int create, int exclusive, int perms );
int   fipc_crit_start      ( int semid );
int   fipc_crit_end        ( int semid );
int   fipc_remove_sem      ( int semid );
int   fipc_inquire_seg     ( int shmid, int n, long shm_data[] );
int   fipc_sizeof_c_int    ( void );
int   fipc_sizeof_c_long   ( void );
int   fipc_sizeof_c_double ( void );
int   fipc_sizeof_c_complex( void );

void  fipc_get_errval ( int *eacces, int *eexist, int *einval, int *enfile, 
                        int *enoent, int *enomem, int *enospc ) {

  /* Get the values of the macros used for marking Sys V errors */

  *eacces = - EACCES;
  *eexist = - EEXIST;
  *einval = - EINVAL;
  *enfile = - ENFILE;
  *enoent = - ENOENT;
  *enomem = - ENOMEM;
  *enospc = - ENOSPC;

}

int fipc_get_seg( long size, int create, int exclusive, int perms ) {

  /* 
     Attempt to create a shared memory segment of size size bytes. If create
     is set create a new one. If exclusive is set create it in exclusive mode.
     The segment will have permissions given by perms

     On a positive value means success, a negative one an error. On error the value
     corresponds to the negative of the errno set by shmget
  */
  
  /* List of files to be tried in ftok. This set of files are chosen are some of those
     that the Filesystem Hierarchy Standard 2.3 says must exist - 
     see http://proton.pathname.com/fhs. If needed can extend
     this, not all the required files are used */

  int key;
  int retval;

  /* Generate a key from the next file in the list, the project number and ftok.
  Only least significant 8 bits of the project are used as that's how ftok seems to work */
  if( next_proj == 256 ) {
    next_proj = 0;
    next_file++;
    if( next_file == N_BASE_FILES ) {
      next_file = 0;
    }
  }
  else {
    next_proj++;
  }
  key = ftok( base_files[ next_file ], next_proj );

  if( create ) {
    perms |= IPC_CREAT;
  }

  if( exclusive ) {
    perms |= IPC_EXCL;
  }

  retval = shmget( key, size, perms );

  if( retval == -1 ) {
    retval = -errno;
  }

  return retval;

}

void *fipc_attach_seg( int shmid, int read_only ) {

  /*
    Attach to segment shmid. If read_only is set attach in read only mode

    On success return a pointer to the segment. On error return a NULL pointer
  */

  void *shmaddr;
  int flags;

  flags = 0;

  if( read_only ) {
    flags |= SHM_RDONLY;
  }

  shmaddr = shmat( shmid, NULL, flags );

  if( shmaddr == ( void * ) -1 ) {
    shmaddr = NULL;
  }

  return shmaddr;
  
}

int fipc_detach_seg( void *shmaddr ) {

  /*
    Detach from a seg at address shmaddr.

    Return zero on success, -1 on error
  */
  
  return shmdt( shmaddr );

}

int fipc_remove_seg( int shmid ) {

  /*
    Mark segment shmid for removal
    
    Return zero on success, -1 on error
  */

  struct shmid_ds *buf = NULL;

  return shmctl( shmid, IPC_RMID, buf );

}

int fipc_get_sem( int create, int exclusive, int perms ) {

  union semun {
    int val;
    struct semid_ds *buf;
    unsigned short *array;
  } sem_union;

  key_t key;

  int semid;

  /* Generate a key from the next file in the list, the project number and ftok.
  Only least significant 8 bits of the project are used as that's how ftok seems to work */
  if( next_proj == 256 ) {
    next_proj = 0;
    next_file++;
    if( next_file == N_BASE_FILES ) {
      next_file = 0;
    }
  }
  else {
    next_proj++;
  }
  key = ftok( base_files[ next_file ], next_proj );

  if( create ) {
    perms |= IPC_CREAT;
  }

  if( exclusive ) {
    perms |= IPC_EXCL;
  }

  if( ( semid = semget( key, 1, perms ) ) == -1) { 
    return -errno;
  } 

  sem_union.val = 1;
  if( semctl( semid, 0, SETVAL, sem_union ) == -1 ) {
    return -errno;
  }

  return semid;

}

int fipc_crit_start( int semid ) {

  struct sembuf sem_op;

  sem_op.sem_num = 0;
  sem_op.sem_op  = -1;
  sem_op.sem_flg = SEM_UNDO;

  return semop( semid, &sem_op, 1 );

}

int fipc_crit_end( int semid ) {

  struct sembuf sem_op;

  sem_op.sem_num = 0;
  sem_op.sem_op  = 1;
  sem_op.sem_flg = SEM_UNDO;

  return semop( semid, &sem_op, 1 );

}

int fipc_remove_sem( int semid ) {

  union semun {
    int val;
    struct semid_ds *buf;
    unsigned short *array;
  } sem_union;

  return semctl( semid, 0, IPC_RMID, sem_union );

}

int fipc_inquire_seg( int shmid, int n, long shm_data[] ){

  /*
    Inquire the properties of segment shmid and return the data
    in an array shm_data of length n. n must be at least 7

    On success return zero and shm_data will hold the following data:
    shm_data[ 0 ] : The size of the segment in bytes
    shm_data[ 1 ] : process ID of last shared memory operation
    shm_data[ 2 ] : process ID of creator
    shm_data[ 3 ] : number of current attaches
    shm_data[ 4 ] : time of last shmat()
    shm_data[ 5 ] : time of last shmdt()
    shm_data[ 6 ] : time of last change by shmctl()
    These are chosen as they correspond to members of the shmid_ds
    structure returned by shmctl that the opengroup "Single Unix
    Specification, Version 2" says must exist. See  
    http://opengroup.org/onlinepubs/007908775/xsh/sysshm.h.html
    Note the last three are returned in the structure as time_t
    objects which the above specification says may be either
    integers or floating poin types. Here we convert to long 
    for uniformity.

    On error return a negative value.
    -1 : The shmctl failed 
    -2 : The supplied array is not big enough, i.e. n < 7
    -3 : Failed to allocate the memory for the shid_ds strcture needed by semctl
  */

  struct shmid_ds *buf;

  int retval;

  buf = malloc( sizeof( *buf ) ); 

  if( buf != NULL ) {

    if( n >= 7 ) {
      retval = shmctl( shmid, IPC_STAT, buf );
      if( retval == 0 ) {
	shm_data[ 0 ] = buf->shm_segsz;
	shm_data[ 1 ] = buf->shm_lpid;
	shm_data[ 2 ] = buf->shm_cpid;
	shm_data[ 3 ] = buf->shm_nattch;
	shm_data[ 4 ] = buf->shm_atime;
	shm_data[ 5 ] = buf->shm_dtime;
	shm_data[ 6 ] = buf->shm_ctime;
	retval = 0;
      }
      else {
	retval = -1;
      }
    }

    else {
      retval = -2;
    }

    free( buf );

  }
  else {
    retval = -3;
  }
    
  return retval;

}

/* Stupid functions to work around problems on stupid 
   broken Cray XT4 series */

int fipc_sizeof_c_int( void ) {

  return sizeof( int );

}

int fipc_sizeof_c_long( void ) {

  return sizeof( long );

}

int fipc_sizeof_c_double( void ) {

  return sizeof( double );

}


int fipc_sizeof_c_complex( void ) {

  return 2 * sizeof( double );

}


