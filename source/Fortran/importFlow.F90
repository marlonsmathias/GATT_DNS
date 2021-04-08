    if(nrank.eq.0) then! Run only in the root process
        allocate(sliceSizes(nproc-1))
        allocate(slicesJStarts(nproc-1))
        allocate(slicesJEnds(nproc-1))
        allocate(slicesKStarts(nproc-1))
        allocate(slicesKEnds(nproc-1))
        
        allocate(Ug(nx,ny,nz))
        allocate(Vg(nx,ny,nz))
        allocate(Wg(nx,ny,nz))
        allocate(Rg(nx,ny,nz))
        allocate(Eg(nx,ny,nz))
        
        allocate(insideWall(nx,ny,nz))
        
        call readFlow(nSave,t,Ug,Vg,Wg,Rg,Eg)
        
        do i = 1,nx
            do j = 1,ny
                do k = 1,nz
                    insideWall(i,j,k) = isnan(Ug(i,j,k))
                    if (insideWall(i,j,k)) then
                        NaN = Ug(i,j,k)
                        Ug(i,j,k) = 0
                        Vg(i,j,k) = 0
                        Wg(i,j,k) = 0
                        Rg(i,j,k) = 1
                        Eg(i,j,k) = 1
                    endif
                enddo
            enddo
        enddo
        
        ! Get local slice
        U = Ug(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))
        V = Vg(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))
        W = Wg(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))
        R = Rg(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))
        E = Eg(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))
    
        
        do i = 1,nproc-1
        
            !First, get the sizes of all slices
            call MPI_RECV(sliceSizes(i),    1, MPI_INT, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
            call MPI_RECV(slicesJStarts(i), 1, MPI_INT, i, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
            call MPI_RECV(slicesJEnds(i),   1, MPI_INT, i, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
            call MPI_RECV(slicesKStarts(i), 1, MPI_INT, i, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
            call MPI_RECV(slicesKEnds(i),   1, MPI_INT, i, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
            
            !Now send the data
            call MPI_SEND(tstep,                                                                                 1,             MPI_INT,   i, 6,  MPI_COMM_WORLD, ierror )
            call MPI_SEND(t,                                                                                     1,             MPI_REAL8, i, 7,  MPI_COMM_WORLD, ierror )
            call MPI_SEND(Ug(xstart(1):xend(1),slicesJStarts(i):slicesJEnds(i),slicesKStarts(i):slicesKEnds(i)), sliceSizes(i), MPI_REAL8, i, 8,  MPI_COMM_WORLD, ierror )
            call MPI_SEND(Vg(xstart(1):xend(1),slicesJStarts(i):slicesJEnds(i),slicesKStarts(i):slicesKEnds(i)), sliceSizes(i), MPI_REAL8, i, 9,  MPI_COMM_WORLD, ierror )
            call MPI_SEND(Wg(xstart(1):xend(1),slicesJStarts(i):slicesJEnds(i),slicesKStarts(i):slicesKEnds(i)), sliceSizes(i), MPI_REAL8, i, 10, MPI_COMM_WORLD, ierror )
            call MPI_SEND(Rg(xstart(1):xend(1),slicesJStarts(i):slicesJEnds(i),slicesKStarts(i):slicesKEnds(i)), sliceSizes(i), MPI_REAL8, i, 11, MPI_COMM_WORLD, ierror )
            call MPI_SEND(Eg(xstart(1):xend(1),slicesJStarts(i):slicesJEnds(i),slicesKStarts(i):slicesKEnds(i)), sliceSizes(i), MPI_REAL8, i, 12, MPI_COMM_WORLD, ierror )
            
        enddo
        
    else ! If not root, receive data from root
    
        !First, send data about the slice size
        call MPI_SEND( xsize(1)*xsize(2)*xsize(3), 1, MPI_INT, 0, 1, MPI_COMM_WORLD, ierror )
        call MPI_SEND( xstart(2),                  1, MPI_INT, 0, 2, MPI_COMM_WORLD, ierror )
        call MPI_SEND( xend(2),                    1, MPI_INT, 0, 3, MPI_COMM_WORLD, ierror ) 
        call MPI_SEND( xstart(3),                  1, MPI_INT, 0, 4, MPI_COMM_WORLD, ierror ) 
        call MPI_SEND( xend(3),                    1, MPI_INT, 0, 5, MPI_COMM_WORLD, ierror ) 
        
        !Then, receive flow data from root
        call MPI_RECV(tstep, 1,                      MPI_INT  , 0, 6,  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
        call MPI_RECV(t, 1,                          MPI_REAL8, 0, 7,  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
        call MPI_RECV(U, xsize(1)*xsize(2)*xsize(3), MPI_REAL8, 0, 8,  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
        call MPI_RECV(V, xsize(1)*xsize(2)*xsize(3), MPI_REAL8, 0, 9,  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
        call MPI_RECV(W, xsize(1)*xsize(2)*xsize(3), MPI_REAL8, 0, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
        call MPI_RECV(R, xsize(1)*xsize(2)*xsize(3), MPI_REAL8, 0, 11, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
        call MPI_RECV(E, xsize(1)*xsize(2)*xsize(3), MPI_REAL8, 0, 12, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
    
    endif
