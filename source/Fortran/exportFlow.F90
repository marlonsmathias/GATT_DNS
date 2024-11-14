        ! Gather global flow from slices
        if(nrank.eq.0) then! Run only in the root process
        
        Ug(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) = U
        Vg(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) = V
        Wg(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) = W
        Rg(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) = R
        Eg(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) = E
        
        do i = 1,nproc-1
            call MPI_RECV(Ug(xstart(1):xend(1),slicesJStarts(i):slicesJEnds(i),slicesKStarts(i):slicesKEnds(i)), sliceSizes(i), MPI_REAL8, i, 1,  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
            call MPI_RECV(Vg(xstart(1):xend(1),slicesJStarts(i):slicesJEnds(i),slicesKStarts(i):slicesKEnds(i)), sliceSizes(i), MPI_REAL8, i, 2,  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
            call MPI_RECV(Wg(xstart(1):xend(1),slicesJStarts(i):slicesJEnds(i),slicesKStarts(i):slicesKEnds(i)), sliceSizes(i), MPI_REAL8, i, 3,  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
            call MPI_RECV(Rg(xstart(1):xend(1),slicesJStarts(i):slicesJEnds(i),slicesKStarts(i):slicesKEnds(i)), sliceSizes(i), MPI_REAL8, i, 4,  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
            call MPI_RECV(Eg(xstart(1):xend(1),slicesJStarts(i):slicesJEnds(i),slicesKStarts(i):slicesKEnds(i)), sliceSizes(i), MPI_REAL8, i, 5,  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
        enddo
        
        ! Add nans inside walls
        do i = 1,nx
            do j = 1,ny
                do k = 1,nz
                    if (insideWall(i,j,k)) then
                        Ug(i,j,k) = NaN
                        Vg(i,j,k) = NaN
                        Wg(i,j,k) = NaN
                        Rg(i,j,k) = NaN
                        Eg(i,j,k) = NaN
                    endif
                enddo
            enddo
        enddo

        if(saveDerivs) then

            Udotg(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) = Udot
            Vdotg(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) = Vdot
            Wdotg(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) = Wdot
            Rdotg(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) = Rdot
            Edotg(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) = Edot

            do i = 1,nproc-1
                call MPI_RECV(Udotg(xstart(1):xend(1),slicesJStarts(i):slicesJEnds(i),slicesKStarts(i):slicesKEnds(i)), sliceSizes(i), MPI_REAL8, i, 6,  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
                call MPI_RECV(Vdotg(xstart(1):xend(1),slicesJStarts(i):slicesJEnds(i),slicesKStarts(i):slicesKEnds(i)), sliceSizes(i), MPI_REAL8, i, 7,  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
                call MPI_RECV(Wdotg(xstart(1):xend(1),slicesJStarts(i):slicesJEnds(i),slicesKStarts(i):slicesKEnds(i)), sliceSizes(i), MPI_REAL8, i, 8,  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
                call MPI_RECV(Rdotg(xstart(1):xend(1),slicesJStarts(i):slicesJEnds(i),slicesKStarts(i):slicesKEnds(i)), sliceSizes(i), MPI_REAL8, i, 9,  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
                call MPI_RECV(Edotg(xstart(1):xend(1),slicesJStarts(i):slicesJEnds(i),slicesKStarts(i):slicesKEnds(i)), sliceSizes(i), MPI_REAL8, i, 10,  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
            enddo
            
            ! Add nans inside walls
            do i = 1,nx
                do j = 1,ny
                    do k = 1,nz
                        if (insideWall(i,j,k)) then
                            Udotg(i,j,k) = NaN
                            Vdotg(i,j,k) = NaN
                            Wdotg(i,j,k) = NaN
                            Rdotg(i,j,k) = NaN
                            Edotg(i,j,k) = NaN
                        endif
                    enddo
                enddo
            enddo

            ! Save the flow to file
            call writeFlowDeriv(nSave,t,Ug,Vg,Wg,Rg,Eg,Udotg,Vdotg,Wdotg,Rdotg,Edotg)

        else
            ! Save the flow to file
            call writeFlow(nSave,t,Ug,Vg,Wg,Rg,Eg)
        endif
        
        else ! If not the root, send data
            call MPI_SEND(U, xsize(1)*xsize(2)*xsize(3), MPI_REAL8, 0, 1, MPI_COMM_WORLD, ierror)
            call MPI_SEND(V, xsize(1)*xsize(2)*xsize(3), MPI_REAL8, 0, 2, MPI_COMM_WORLD, ierror)
            call MPI_SEND(W, xsize(1)*xsize(2)*xsize(3), MPI_REAL8, 0, 3, MPI_COMM_WORLD, ierror)
            call MPI_SEND(R, xsize(1)*xsize(2)*xsize(3), MPI_REAL8, 0, 4, MPI_COMM_WORLD, ierror)
            call MPI_SEND(E, xsize(1)*xsize(2)*xsize(3), MPI_REAL8, 0, 5, MPI_COMM_WORLD, ierror)

            if(saveDerivs) then
                call MPI_SEND(Udot, xsize(1)*xsize(2)*xsize(3), MPI_REAL8, 0, 6, MPI_COMM_WORLD, ierror)
                call MPI_SEND(Vdot, xsize(1)*xsize(2)*xsize(3), MPI_REAL8, 0, 7, MPI_COMM_WORLD, ierror)
                call MPI_SEND(Wdot, xsize(1)*xsize(2)*xsize(3), MPI_REAL8, 0, 8, MPI_COMM_WORLD, ierror)
                call MPI_SEND(Rdot, xsize(1)*xsize(2)*xsize(3), MPI_REAL8, 0, 9, MPI_COMM_WORLD, ierror)
                call MPI_SEND(Edot, xsize(1)*xsize(2)*xsize(3), MPI_REAL8, 0, 10, MPI_COMM_WORLD, ierror)
            endif

        endif
