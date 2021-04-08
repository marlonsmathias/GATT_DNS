    if(nrank.eq.0) then! Run only in the root process
        
        call readMeanFlow(Ug,Vg,Wg,Rg,Eg)
        
        do i = 1,nx
            do j = 1,ny
                do k = 1,nz
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
        Umean = Ug(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))
        Vmean = Vg(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))
        Wmean = Wg(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))
        Rmean = Rg(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))
        Emean = Eg(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))
    
        
        do i = 1,nproc-1
            !Send the data
            call MPI_SEND(Ug(xstart(1):xend(1),slicesJStarts(i):slicesJEnds(i),slicesKStarts(i):slicesKEnds(i)), sliceSizes(i), MPI_REAL8, i, 8,  MPI_COMM_WORLD, ierror )
            call MPI_SEND(Vg(xstart(1):xend(1),slicesJStarts(i):slicesJEnds(i),slicesKStarts(i):slicesKEnds(i)), sliceSizes(i), MPI_REAL8, i, 9,  MPI_COMM_WORLD, ierror )
            call MPI_SEND(Wg(xstart(1):xend(1),slicesJStarts(i):slicesJEnds(i),slicesKStarts(i):slicesKEnds(i)), sliceSizes(i), MPI_REAL8, i, 10, MPI_COMM_WORLD, ierror )
            call MPI_SEND(Rg(xstart(1):xend(1),slicesJStarts(i):slicesJEnds(i),slicesKStarts(i):slicesKEnds(i)), sliceSizes(i), MPI_REAL8, i, 11, MPI_COMM_WORLD, ierror )
            call MPI_SEND(Eg(xstart(1):xend(1),slicesJStarts(i):slicesJEnds(i),slicesKStarts(i):slicesKEnds(i)), sliceSizes(i), MPI_REAL8, i, 12, MPI_COMM_WORLD, ierror )
            
        enddo
        
    else ! If not root, receive data from root
        
        !Receive flow data from root
        call MPI_RECV(Umean, xsize(1)*xsize(2)*xsize(3), MPI_REAL8, 0, 8,  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
        call MPI_RECV(Vmean, xsize(1)*xsize(2)*xsize(3), MPI_REAL8, 0, 9,  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
        call MPI_RECV(Wmean, xsize(1)*xsize(2)*xsize(3), MPI_REAL8, 0, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
        call MPI_RECV(Rmean, xsize(1)*xsize(2)*xsize(3), MPI_REAL8, 0, 11, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
        call MPI_RECV(Emean, xsize(1)*xsize(2)*xsize(3), MPI_REAL8, 0, 12, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
    
    endif
