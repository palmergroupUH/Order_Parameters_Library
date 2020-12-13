SUBROUTINE calc_q3_cluster(n_atoms,n_q3_neigh,q3_cutoff,n_bonds,box,x,largest_cluster)

  IMPLICIT NONE

  ! Passed
  INTEGER,INTENT(IN) :: n_atoms
  INTEGER,INTENT(IN) :: n_q3_neigh
  INTEGER,INTENT(IN) :: n_bonds
  DOUBLE PRECISION,INTENT(IN) :: q3_cutoff
  DOUBLE PRECISION,INTENT(IN) :: box
  DOUBLE PRECISION,INTENT(IN),DIMENSION(0:3*n_atoms-1) :: x
  INTEGER,INTENT(OUT) :: largest_cluster

  ! Local
  INTEGER :: i,j,k
  INTEGER :: iatom,jatom, iindx,jindx
  DOUBLE PRECISION,DIMENSION(3) :: xi,xj,dxij
  DOUBLE PRECISION :: rij,rijsq

  ! Neighbor list
  INTEGER,DIMENSION(n_atoms) :: n_nneigh
  INTEGER,DIMENSION(75,n_atoms) :: nneigh_lst
  DOUBLE PRECISION,DIMENSION(4,75,n_atoms) :: rij_lst
  INTEGER :: nneigh_temp
  DOUBLE PRECISION,DIMENSION(4) :: rij_temp
  DOUBLE PRECISION :: q3_cutoffsq
 
  ! Spherical Harmonics
  DOUBLE PRECISION :: Pi,factor, fcij
  DOUBLE PRECISION :: pre0,pre1,pre2,pre3
  DOUBLE PRECISION :: pre(-3:3)
  DOUBLE PRECISION :: cost(3),sint(3)
  DOUBLE PRECISION :: cosp(3), sinp(3)
  COMPLEX*16 :: Y(-3:3),comaux
  COMPLEX*16 :: sum_Y(-3:3),sum_q3(-3:3,n_atoms)
  DOUBLE PRECISION,DIMENSION(75,n_atoms) :: dij
  DOUBLE PRECISION,DIMENSION(n_atoms) :: q3iatom

  ! Clustering
  INTEGER,DIMENSION(n_atoms) :: n_dij_bonds
  INTEGER :: n_clusters,n_ice
  INTEGER :: clusterID,clusterID_old, d6ijcounter
  INTEGER,DIMENSION(n_atoms) :: cluster_size, cluster_ID
  INTEGER,DIMENSION(n_atoms, n_atoms) :: cluster_lst
  LOGICAL,DIMENSION(n_atoms) :: InCluster, visited

  q3_cutoffsq =  q3_cutoff*q3_cutoff

  Pi = 3.1415926535897930d0
  factor = 2.0d0*dSQRT(Pi/7.0d0)
  fcij = factor*factor
  pre(-3) = (1.0d0/8.0d0)*dSQRT(35.0d0/Pi)
  pre(-2) = (1.0d0/4.0d0)*dSQRT(105.0d0/(2.*Pi))
  pre(-1) = (1.0d0/8.0d0)*dSQRT(21.0d0/Pi)
  pre(0) =  (1.0d0/4.0d0)*dSQRT(7.0d0/Pi)
  pre(1) = -pre(-1)
  pre(2) =  pre(-2)
  pre(3) = -pre(-3)

  ! Initialize arrays
  n_nneigh = 0
  nneigh_lst = 0
  rij_lst = 0.0d0
  dij = 0.0d0
  q3iatom = 0.0d0
  sum_q3 = 0.0d0
  dij = 0.0d0

  ! Loop over iatom
  DO iatom = 1,n_atoms-1
    iindx = 3*(iatom-1)

    ! Store iatom's coordinates
    xi(1) = x(iindx)
    xi(2) = x(iindx+1)
    xi(3) = x(iindx+2)

    ! Loop over jatom
    DO jatom = iatom+1,n_atoms
      jindx = 3*(jatom-1)

      ! Store jatom's coordinates
      xj(1) = x(jindx)
      xj(2) = x(jindx+1)
      xj(3) = x(jindx+2)

      ! Compute ij vector & separation distance
      dxij(:) = xj(:) - xi(:)
      dxij(:) = dxij(:) - dNINT(dxij(:)/box)*box
      rijsq = dxij(1)*dxij(1) + dxij(2)*dxij(2) + dxij(3)*dxij(3)

      ! If they are within the cutoff distance
      IF(rijsq .GT. q3_cutoffsq) CYCLE

        rij = dSQRT(rijsq)

        ! Update information for iatom

        ! Update the number of neighbors for iatom
        n_nneigh(iatom) = n_nneigh(iatom) + 1

        ! Include jatom as a neighbor of iatom
        nneigh_lst(n_nneigh(iatom),iatom)  = jatom

        ! Store the distance components for iatom
        rij_lst(1, n_nneigh(iatom), iatom) = dxij(1)
        rij_lst(2, n_nneigh(iatom), iatom) = dxij(2)
        rij_lst(3, n_nneigh(iatom), iatom) = dxij(3)
        rij_lst(4, n_nneigh(iatom), iatom) = rij

        ! Update information for jatom

        ! Update the number of neighbors for jatom
        n_nneigh(jatom) = n_nneigh(jatom) + 1

        ! Include iatom as a neighbor of jatom
        nneigh_lst(n_nneigh(jatom),jatom) = iatom

        ! Store the distance components for jatom
        rij_lst(1, n_nneigh(jatom), jatom) = -dxij(1)
        rij_lst(2, n_nneigh(jatom), jatom) = -dxij(2)
        rij_lst(3, n_nneigh(jatom), jatom) = -dxij(3)
        rij_lst(4, n_nneigh(jatom), jatom) = rij


      !ENDIF

    ENDDO
  ENDDO 

  ! Check the minimum number of neighbors are present
  IF(MINVAL(n_nneigh(:)) .LT. n_q3_neigh .AND. n_q3_neigh .GT. 0) THEN
    WRITE(*,'(A,I5,A,I3,A)') 'FATAL ERROR: ', MINLOC(n_nneigh(:)), ' HAS LESS THAN ', n_q3_neigh,' NEAREST NEIGHBORS IN SUBROUTINE Q3, ADJUST Q3_CUTOFF'
    STOP
  ENDIF

  ! Check that the maximum number of neighbors is not exceeded
  IF(MAXVAL(n_nneigh(:)) .GT. 75) THEN
    WRITE(*,'(A,I5,A)') 'FATAL ERROR: ', MAXLOC(n_nneigh(:)), ' HAS MORE THAN 75 NEAREST NEIGHBORS IN SUBROUTINE Q3'
    STOP
  ENDIF

  ! Sort the neighbors using insertion sort if using only the closest n_q3_neigh
  ! neighbors
  IF(n_q3_neigh .GT. 0 ) THEN
    DO iatom = 1,n_atoms
      DO k = 2,n_nneigh(iatom)
        nneigh_temp = nneigh_lst(k,iatom)
        rij_temp(:) = rij_lst(:,k,iatom)
        DO j = k-1,1,-1
          IF(rij_lst(4,j,iatom) .LE. rij_temp(4)) GOTO 10
          nneigh_lst(j + 1, iatom) = nneigh_lst(j,iatom)
          rij_lst(:,j + 1,iatom) = rij_lst(:, j, iatom)
        ENDDO
        jatom = 0
        10 CONTINUE
        nneigh_lst(j + 1, iatom) = nneigh_temp
        rij_lst(:,j + 1, iatom) = rij_temp(:)
      ENDDO
      n_nneigh(iatom) = n_q3_neigh
    ENDDO
  ENDIF


  ! Calculate q3 and dij

  ! Loop over the atoms
  DO iatom = 1,n_atoms
    iindx = 3*(iatom-1)

    ! Zero harmonic sum
    sum_Y = 0.0d0

    ! Loop over the nearest neighbors of iatom
    DO j = 1, n_nneigh(iatom)
      jatom = nneigh_lst(j,iatom)
      jindx = 3*(jatom-1)

      ! Get separation vector for list
      dxij(1) = rij_lst(1,j,iatom)
      dxij(2) = rij_lst(2,j,iatom)
      dxij(3) = rij_lst(3,j,iatom)
      rij = rij_lst(4,j,iatom)

      ! Calculate cos(theta) and sin(theta)
      cost(1) = dxij(3)/rij
      sint(1) = (1.0d0-cost(1)*cost(1))**0.5d0

      ! Calculate cos(phi) and sin(phi)
      IF (dABS(sint(1)) .LT. 1.0d-9) THEN
        cosp(1) = 1.0d0
        sinp(1) = 0.0d0
      ELSE
        cosp(1) = dxij(1)/(rij*sint(1))
        sinp(1) = dxij(2)/(rij*sint(1))
      ENDIF

      ! Calculate powers of sin(theta)
      sint(2) = sint(1)*sint(1)
      sint(3) = sint(2)*sint(1)

      ! Calculate powers of cos(theta)
      cost(2) = cost(1)*cost(1)
      cost(3) = cost(2)*cost(1)

      ! Calculate sin/cos (2phi)
      sinp(2) = 2.0d0*sinp(1)*cosp(1)
      cosp(2) = 2.0d0*cosp(1)**2 - 1.0d0

      ! Calculate sin/cos (3phi)
      sinp(3) = sinp(2)*cosp(1) + cosp(2)*sinp(1)
      cosp(3) = cosp(2)*cosp(1) - sinp(2)*sinp(1)

      ! Update the harmonic sum
      Y(-3) = pre(-3)*sint(3)*dCMPLX(cosp(3),-sinp(3))
      Y(-2) = pre(-2)*sint(2)*cost(1)*dCMPLX(cosp(2),-sinp(2))
      Y(-1) = pre(-1)*sint(1)*(5.0d0*cost(2)-1.0)*dCMPLX(cosp(1),-sinp(1))
      Y(0) = pre(0)*(5.*cost(3)-3.0d0*cost(1))
      Y(1) = pre(1)*sint(1)*(5.0d0*cost(2)-1.0)*dCMPLX(cosp(1),sinp(1))
      Y(2) = pre(2)*sint(2)*cost(1)*dCMPLX(cosp(2),sinp(2))
      Y(3) = pre(3)*sint(3)*dCMPLX(cosp(3),sinp(3))

      sum_Y = sum_Y + Y

    ENDDO

    sum_q3(:,iatom) = sum_Y(:)

    ! Calculate q3 for each atom
    DO i = -3,3
      q3iatom(iatom) = q3iatom(iatom) + sum_q3(i,iatom)*dCONJG(sum_q3(i,iatom))
    ENDDO
    q3iatom(iatom) = factor*dSQRT(q3iatom(iatom))/DBLE(n_nneigh(iatom))
  ENDDO

  ! Loop over the atoms
  n_dij_bonds = 0
  n_ice = 0
  DO iatom = 1,n_atoms
    DO j = 1, n_nneigh(iatom)
      jatom = nneigh_lst(j,iatom)
      DO i = -3,3
        dij(j,iatom) = dij(j,iatom) + sum_q3(i,iatom)*CONJG(sum_q3(i,jatom))
      ENDDO
      dij(j,iatom) = fcij*dij(j,iatom)/(q3iatom(iatom)*q3iatom(jatom))/DBLE(n_nneigh(iatom)*n_nneigh(jatom))
      IF((dij(j,iatom) .LE. -0.82d0) .OR. ((dij(j,iatom) .GE. -0.145d0) .AND. (dij(j,iatom) .LE. -0.065d0))) THEN
         n_dij_bonds(iatom) = n_dij_bonds(iatom) + 1
      ENDIF
    ENDDO
  ENDDO


  ! Initialize variables
  InCluster = .FALSE. !FALSE if atom i not in a cluster 
  n_clusters = 0
  cluster_size = 0
  cluster_lst = 0
  cluster_ID = 0

  !Loop over the atoms
  DO iatom = 1,n_atoms

    clusterID = 0

    IF (n_dij_bonds(iatom) .LT. n_bonds) CYCLE

      InCluster(iatom) = .TRUE.
  
      !Search iatoms's n.n. to see if they are in cluster(s)
      DO j = 1, n_nneigh(iatom)
        jatom = nneigh_lst(j,iatom)
 
        IF (InCluster(jatom) .EQV. .TRUE.) THEN !one of iatom's n.n. is in a cluster

          IF (clusterID .NE. 0) THEN !another one of iatoms's n.n. is also in a cluster

            IF (clusterID .NE. cluster_ID(jatom)) THEN !two of iatoms's n.n. belong to "different" clusters...combine

              !keep lower cluster ID number
              IF (clusterID .LT. cluster_ID(jatom)) THEN
                clusterID_old = cluster_ID(jatom)
              ELSE
                clusterID_old = clusterID
                clusterID = cluster_ID(jatom)
              ENDIF

              ! Merge clusters, keeping lower cluster ID number
              DO k = 1,cluster_size(clusterID_old)
                cluster_lst(cluster_size(clusterID) + k,clusterID) = cluster_lst(k,clusterID_old)
                cluster_ID(cluster_lst(k,clusterID_old)) = clusterID
              ENDDO
              cluster_size(clusterID) = cluster_size(clusterID) + cluster_size(clusterID_old)
              cluster_lst(:,clusterID_old) = 0 ! Clear old cluster data, reduce count
              cluster_size(clusterID_old) = 0     
              n_clusters = n_clusters - 1      

              !Place LAST cluster in cluster list into newly empty slot
              IF (clusterID_old .NE. (n_clusters+1)) THEN ! Must do move
                DO k = 1,cluster_size(n_clusters+1)
                 cluster_lst(k,clusterID_old) = cluster_lst(k,(n_clusters+1))
                 cluster_ID(cluster_lst(k,clusterID_old)) = clusterID_old
                ENDDO
                cluster_size(clusterID_old) = cluster_size(n_clusters+1)  ! Clear data for moved cluster 
                cluster_lst(:,(n_clusters+1)) = 0 
                cluster_size(n_clusters+1) = 0    
              ENDIF
            ENDIF
       
          ELSE
            clusterID = cluster_ID(jatom)
          ENDIF

        ENDIF
      ENDDO  ! j neighbor

      IF (clusterID .EQ. 0) then !if none of iatom's n.n. are in a cluster, imol begins new cluster
        n_clusters = n_clusters + 1
        clusterID = n_clusters
      ENDIF

      !Add iatom to the cluster
      cluster_ID(iatom) = clusterID
      cluster_size(clusterID) = cluster_size(clusterID) + 1
      cluster_lst(cluster_size(clusterID),clusterID) = iatom !cluster_lst(cluster number, list of atoms in cluster)

  ENDDO ! iatom
  
  largest_cluster = MAXVAL(cluster_size)
  
END SUBROUTINE

SUBROUTINE calc_Q6_cg(n_atoms,n_q6_neigh,q6_cutoff,n_bonds,box,x, Q6iatom) 

  IMPLICIT NONE
    
  ! Passed
  INTEGER,INTENT(IN) :: n_atoms
  INTEGER,INTENT(IN) :: n_q6_neigh
  INTEGER,INTENT(IN) :: n_bonds
  DOUBLE PRECISION,INTENT(IN) :: q6_cutoff
  DOUBLE PRECISION,INTENT(IN) :: box
  DOUBLE PRECISION,INTENT(IN),DIMENSION(0:3*n_atoms-1) :: x
  DOUBLE PRECISION,INTENT(OUT),DIMENSION(n_atoms) :: Q6iatom

  ! Local
  INTEGER :: i,j,k,m 
  INTEGER :: iatom,jatom, iindx,jindx, count_nb, num_NB 
  DOUBLE PRECISION,DIMENSION(3) :: xi,xj,dxij
  DOUBLE PRECISION :: rij,rijsq, dist 

  ! Neighbor list
  INTEGER,DIMENSION(n_atoms) :: n_nneigh
  INTEGER,DIMENSION(75,n_atoms) :: nneigh_lst
  INTEGER, DIMENSION(75) :: temp_nb_lst 
  DOUBLE PRECISION,DIMENSION(4,75,n_atoms) :: rij_lst
  INTEGER :: nneigh_temp
  DOUBLE PRECISION,DIMENSION(4) :: rij_temp
  DOUBLE PRECISION :: q6_cutoffsq
  

  ! Spherical Harmonics
  DOUBLE PRECISION :: Pi,factor, fcij, sum_Q6_norm 
  DOUBLE PRECISION :: pre_and_post(-6:6)
  DOUBLE PRECISION :: cost(6),sint(6)
  DOUBLE PRECISION :: cosp(6), sinp(6)
  COMPLEX*16 :: sum_Y(-6:6) , Y(-6:6), Q6_i_norm(-6:6), Q6_j_norm(-6:6), sum_neighbor_qlm(-6:6)  
  
  complex*16 :: sum_q6(-6:6,1:n_atoms) , Q6_cg(-6:6,1:n_atoms) 
  DOUBLE PRECISION :: scalar_Q6i_Q6j
  !INTEGER, DIMENSION(n_atoms) :: n_dij_bonds
  !DOUBLE PRECISION, DIMENSION(75,n_atoms) :: dij

  Q6iatom = 0.0d0

  Q6_cutoffsq = Q6_cutoff*Q6_cutoff

  Pi = 3.1415926535897930d0 

  factor = 2.0d0*dSQRT(Pi/13.0d0) 

  ! Initialize arrays
  n_nneigh = 0
  nneigh_lst = 0
  rij  = 0.0d0

  ! Loop over iatom
  DO iatom = 1,n_atoms-1
    iindx = 3*(iatom-1)

    ! Store iatom's coordinates
    xi(1) = x(iindx)
    xi(2) = x(iindx+1)
    xi(3) = x(iindx+2)  
 
    ! Loop over jatom
    DO jatom = iatom+1,n_atoms
      jindx = 3*(jatom-1)

      ! Store jatom's coordinates
      xj(1) = x(jindx)
      xj(2) = x(jindx+1)
      xj(3) = x(jindx+2)

      ! Compute ij vector & separation distance
      dxij(:) = xj(:) - xi(:)
      dxij(:) = dxij(:) - dNINT(dxij(:)/box)*box
      rijsq = dxij(1)*dxij(1) + dxij(2)*dxij(2) + dxij(3)*dxij(3)
      
       ! If they are within the cutoff distance
       IF(rijsq .GT. Q6_cutoffsq) CYCLE

         rij = dSQRT(rijsq)
   
         ! Update information for iatom

         ! Update the number of neighbors for iatom
         n_nneigh(iatom) = n_nneigh(iatom) + 1

         ! Include jatom as a neighbor of iatom
         nneigh_lst(n_nneigh(iatom),iatom)  = jatom

         ! Store the distance components for iatom
         rij_lst(1, n_nneigh(iatom), iatom) = dxij(1)
         rij_lst(2, n_nneigh(iatom), iatom) = dxij(2)
         rij_lst(3, n_nneigh(iatom), iatom) = dxij(3)
         rij_lst(4, n_nneigh(iatom), iatom) = rij

         ! Update informtion for jatom

         ! Update the number of neighbors for jatom
         n_nneigh(jatom) = n_nneigh(jatom) + 1

         ! Include iatom as a neighbor of jatom
         nneigh_lst(n_nneigh(jatom),jatom) = iatom

         ! Store the distance components for jatom
         rij_lst(:,n_nneigh(jatom),jatom) = rij_lst(:,n_nneigh(iatom),iatom)

     ENDDO
   ENDDO
  
  ! Check the minimum number of neighbors are present
  IF(MINVAL(n_nneigh(:)) .LT. n_Q6_neigh .AND. n_Q6_neigh .GT. 0) THEN
    WRITE(*,'(A,I5,A,I3,A)') 'FATAL ERROR: ', MINLOC(n_nneigh(:)), ' HAS LESS THAN ', n_Q6_neigh,' NEAREST NEIGHBORS IN SUBROUTINE Q6, ADJUST Q6_CUTOFF'
    STOP
  ENDIF

  ! Check that the maximum number of neighbors is not exceeded
  IF(MAXVAL(n_nneigh(:)) .GT. 75) THEN
    WRITE(*,'(A,I5,A)') 'FATAL ERROR: ', MAXLOC(n_nneigh(:)), ' HAS MORE THAN 75 NEAREST NEIGHBORS IN SUBROUTINE Q6'
    STOP
  ENDIF

  ! Sort the neighbors using insertion sort if using only the closest n_Q6_neigh neighbors
  IF(n_Q6_neigh .GT. 0 ) THEN
    DO iatom = 1,n_atoms
      DO k = 2,n_nneigh(iatom)
        nneigh_temp = nneigh_lst(k,iatom)
        rij_temp(:) = rij_lst(:,k,iatom)
        DO j = k-1,1,-1
          IF(rij_lst(4,j,iatom) .LE. rij_temp(4)) GOTO 10
          nneigh_lst(j + 1, iatom) = nneigh_lst(j,iatom)
          rij_lst(:,j + 1,iatom) = rij_lst(:, j, iatom)
        ENDDO
        jatom = 0
        10 CONTINUE
        nneigh_lst(j + 1, iatom) = nneigh_temp
        rij_lst(:,j + 1, iatom) = rij_temp(:)
      ENDDO
      n_nneigh(iatom) = n_Q6_neigh
    ENDDO
  ENDIF


  ! Calculate Q6
  
  ! Initialize the sum
  
  sum_Q6 = 0.0d0

  ! Loop over the atomecules
  DO iatom = 1,n_atoms

    ! Initialize the harmonic sum
    sum_Y = 0.0d0
  
    ! Loop over the 4 nearest neighbors of iatom
    DO j = 1,n_nneigh(iatom)

      ! Calculate cos(theta) and sin(theta)
      cost(1) = rij_lst(3,j,iatom)/rij_lst(4,j,iatom)
      sint(1) = (1.0d0-cost(1)*cost(1))**0.5d0

      ! Calculate cos(phi) and sin(phi)
      IF (dABS(sint(1)) .LT. 1.0d-9) THEN
        cosp(1) = 1.0d0
        sinp(1) = 0.0d0

      ELSE
        cosp(1) = rij_lst(1,j,iatom)/(rij_lst(4,j,iatom)*sint(1))
        sinp(1) = rij_lst(2,j,iatom)/(rij_lst(4,j,iatom)*sint(1))

      END IF

      ! calculate powers of sin(theta)
      sint(2) = sint(1)*sint(1)
      sint(3) = sint(2)*sint(1)
      sint(4) = sint(3)*sint(1) 
      sint(5) = sint(4)*sint(1) 
      sint(6) = sint(5)*sint(1) 

      ! calculate powers of cos(theta)
      cost(2) = cost(1)*cost(1)
      cost(3) = cost(2)*cost(1)
      cost(4) = cost(3)*cost(1)
      cost(5) = cost(4)*cost(1)
      cost(6) = cost(5)*cost(1)


      ! calculate sin/cos (2phi)
      sinp(2) = 2.0d0*sinp(1)*cosp(1)
      cosp(2) = 2.0d0*cosp(1)**2 - 1.0d0

      ! calculate sin/cos (3phi)
      sinp(3) = sinp(2)*cosp(1) + cosp(2)*sinp(1)
      cosp(3) = cosp(2)*cosp(1) - sinp(2)*sinp(1)

      ! calculate sin/cos (4phi) 

      sinp(4) = sinp(3)*cosp(1) + cosp(3)*sinp(1) 
      cosp(4) = cosp(3)*cosp(1) - sinp(3)*sinp(1) 

      ! calculate sin/cos (5phi) 

      sinp(5) = sinp(4)*cosp(1) + cosp(4)*sinp(1) 
      cosp(5) = cosp(4)*cosp(1) - sinp(4)*sinp(1) 

      ! calculate sin/cos (6phi) 

      sinp(6) = sinp(5)*cosp(1) + cosp(5)*sinp(1) 
      cosp(6) = cosp(5)*cosp(1) - sinp(5)*sinp(1) 

      ! Update the harmonic sum
      ! From Mathmatica function: "SphericalHarmonicsY[4,m,theta, phi]"	

      ! pre and post factor: 

      pre_and_post(-6) = 0.483084113580066d0 * sint(6)  
      pre_and_post(-5) = 1.67345245810010d0 * sint(5) * cost(1)  
      pre_and_post(-4) = 0.356781262853998d0 * (11.0d0 * cost(2)- 1.0d0) * sint(4)  
      pre_and_post(-3) = 0.651390485867716d0 * (-3.0d0*cost(1) + 11.0d0 * cost(3)) * sint(3)  
      pre_and_post(-2) = 0.325695242933858d0 * (1.0d0 - 18.0d0 * cost(2) + 33.0d0* cost(4)) * sint(2)  
      pre_and_post(-1) = 0.411975516301141d0 * (5.0d0*cost(1) - 30.0d0 * cost(3) + 33.0d0*cost(5)) * sint(1)  
      pre_and_post(0) = 0.06356920226762842d0 * (105.0d0 * cost(2) - 315.0d0 * cost(4) + 231.0d0 * cost(6) - 5.0d0)
      pre_and_post(1) = -pre_and_post(-1)  
      pre_and_post(2) = pre_and_post(-2) 
      pre_and_post(3) = -pre_and_post(-3)
      pre_and_post(4) = pre_and_post(-4) 
      pre_and_post(5) = -pre_and_post(-5)
      pre_and_post(6) = pre_and_post(-6) 

      Y(-6) = dcmplx(cosp(6),-sinp(6))*pre_and_post(-6)  
      Y(-5) = dcmplx(cosp(5),-sinp(5))*pre_and_post(-5) 
      Y(-4) = dcmplx(cosp(4),-sinp(4))*pre_and_post(-4)  
      Y(-3) = dcmplx(cosp(3),-sinp(3))*pre_and_post(-3) 
      Y(-2) = dcmplx(cosp(2),-sinp(2))*pre_and_post(-2) 
      Y(-1) = dcmplx(cosp(1),-sinp(1))*pre_and_post(-1)
      Y(0) = pre_and_post(0) 
      Y(1) = dcmplx(cosp(1),sinp(1))* pre_and_post(1) 
      Y(2) = dcmplx(cosp(2),sinp(2))* pre_and_post(2)  
      Y(3) = dcmplx(cosp(3),sinp(3))* pre_and_post(3) 
      Y(4) = dcmplx(cosp(4),sinp(4))* pre_and_post(4)
      Y(5) = dcmplx(cosp(5),sinp(5))* pre_and_post(5) 
      Y(6) = dcmplx(cosp(6),sinp(6))* pre_and_post(6) 

      sum_Y = sum_Y + Y 

    END DO 
 
    sum_Q6(:, iatom) = sum_Q6(:, iatom) + sum_Y/DBLE(n_nneigh(iatom)) 

  END DO

  ! compute coarse-grained Q6

  Q6_cg = 0.0d0 
  
  DO iatom = 1, n_atoms 

    ! sum of q6 bond order includes particle i itsself
    sum_neighbor_qlm = sum_q6(:,iatom) 

    DO j = 1, n_nneigh(iatom)  
 
        jatom = nneigh_lst(j,iatom) 
        
        sum_neighbor_qlm = sum_neighbor_qlm + sum_q6(:,jatom) 
   
    END DO 
    
    Q6_cg(:, iatom) = sum_neighbor_qlm / n_nneigh(iatom)  
      
    ! compute normalization 

    sum_Q6_norm = 0.0d0

    DO i = -6, 6

        sum_Q6_norm = sum_Q6_norm + Q6_cg(i,iatom)*dCONJG(Q6_cg(i,iatom))

    END DO 

    ! compute the second order invariants of coarse-grained Q6
    Q6iatom(iatom) = factor*dSQRT(sum_Q6_norm)
 
  END DO
  
END SUBROUTINE 

SUBROUTINE calc_q12_cluster(n_atoms,n_q12_neigh,q12_cutoff,crys_cutoff,n_bonds,box,x,nxtl,largest_cluster) 

  ! Last Modified by Jingxiang Guo (8/9/2020)

  ! Note: 

  ! 1. Add an new arugment "crys_cutoff" in function call
  ! This a cutoff determining if two particles are in a cluster, and maybe
  ! different from "ql_cutoff".


  ! 2. For spherical harmonics Y, the factorial constant invovling l,m 
  ! are hard-coded as double-precision a scalar value.

  IMPLICIT NONE

  ! Passed
  INTEGER,INTENT(IN) :: n_atoms
  INTEGER,INTENT(IN) :: n_q12_neigh
  INTEGER,INTENT(IN) :: n_bonds
  DOUBLE PRECISION,INTENT(IN) :: q12_cutoff, crys_cutoff
  DOUBLE PRECISION,INTENT(IN) :: box
  DOUBLE PRECISION,INTENT(IN),DIMENSION(0:3*n_atoms-1) :: x
  INTEGER,INTENT(OUT) :: nxtl,largest_cluster

  ! Local
  INTEGER :: i,j,k,m 
  INTEGER :: iatom,jatom, iindx,jindx, count_nb, num_NB 
  DOUBLE PRECISION,DIMENSION(3) :: xi,xj,dxij
  DOUBLE PRECISION :: rij,rijsq, dist 

  ! Neighbor list
  INTEGER,DIMENSION(n_atoms) :: n_nneigh
  INTEGER,DIMENSION(75,n_atoms) :: nneigh_lst
  INTEGER, DIMENSION(75) :: temp_nb_lst 
  DOUBLE PRECISION,DIMENSION(4,75,n_atoms) :: rij_lst
  INTEGER :: nneigh_temp
  DOUBLE PRECISION,DIMENSION(4) :: rij_temp
  DOUBLE PRECISION :: q12_cutoffsq
  

  ! Spherical Harmonics
  DOUBLE PRECISION :: Pi,factor, fcij, sum_Q12_norm 
  DOUBLE PRECISION :: pre_and_post(-12:12)
  DOUBLE PRECISION :: cost(12),sint(12)
  DOUBLE PRECISION :: cosp(12), sinp(12)
  COMPLEX*16 :: sum_Y(-12:12) , Y(-12:12), Q12_i_norm(-12:12), Q12_j_norm(-12:12), sum_neighbor_qlm(-12:12)  
  DOUBLE PRECISION,DIMENSION(n_atoms) :: Q12iatom
  complex*16 :: sum_q12(-12:12,1:n_atoms) , Q12_cg(-12:12,1:n_atoms) 
  DOUBLE PRECISION :: scalar_Q12i_Q12j
  INTEGER, DIMENSION(n_atoms) :: n_dij_bonds
  DOUBLE PRECISION, DIMENSION(75,n_atoms) :: dij

  ! Clustering
  INTEGER :: n_clusters
  INTEGER :: clusterID,clusterID_old, d6ijcounter
  INTEGER,DIMENSION(n_atoms) :: cluster_size, cluster_ID
  INTEGER,DIMENSION(n_atoms, n_atoms) :: cluster_lst
  LOGICAL,DIMENSION(n_atoms) :: InCluster, visited

  q12_cutoffsq =  q12_cutoff*q12_cutoff

  Pi = 3.1415926535897930d0 

  factor = 2.0d0*dSQRT(Pi/25.0d0) 

  fcij = factor*factor 

  ! Initialize arrays
  n_nneigh = 0
  nneigh_lst = 0
  rij_lst = 0.0d0
  rij  = 0.0d0
  
  ! Loop over iatom
  DO iatom = 1,n_atoms-1
    iindx = 3*(iatom-1)

    ! Store iatom's coordinates
    xi(1) = x(iindx)
    xi(2) = x(iindx+1)
    xi(3) = x(iindx+2)
    
    ! Loop over jatom
    DO jatom = iatom+1,n_atoms
      jindx = 3*(jatom-1)

      ! Store jatom's coordinates
      xj(1) = x(jindx)
      xj(2) = x(jindx+1)
      xj(3) = x(jindx+2)

      ! Compute ij vector & separation distance
      dxij(:) = xj(:) - xi(:)
      dxij(:) = dxij(:) - dNINT(dxij(:)/box)*box
      rijsq = dxij(1)*dxij(1) + dxij(2)*dxij(2) + dxij(3)*dxij(3)

       ! If they are within the cutoff distance
       IF(rijsq .GT. q12_cutoffsq) CYCLE

         rij = dSQRT(rijsq)         

         ! Update information for iatom

         ! Update the number of neighbors for iatom
         n_nneigh(iatom) = n_nneigh(iatom) + 1

         ! Include jatom as a neighbor of iatom
         nneigh_lst(n_nneigh(iatom),iatom)  = jatom

         ! Store the distance components for iatom
         rij_lst(1, n_nneigh(iatom), iatom) = dxij(1)
         rij_lst(2, n_nneigh(iatom), iatom) = dxij(2)
         rij_lst(3, n_nneigh(iatom), iatom) = dxij(3)
         rij_lst(4, n_nneigh(iatom), iatom) = rij

         ! Update informtion for jatom

         ! Update the number of neighbors for jatom
         n_nneigh(jatom) = n_nneigh(jatom) + 1

         ! Include iatom as a neighbor of jatom
         nneigh_lst(n_nneigh(jatom),jatom) = iatom
 
         ! Store the distance components for jatom
         rij_lst(:,n_nneigh(jatom),jatom) = rij_lst(:,n_nneigh(iatom),iatom)

     ENDDO
   ENDDO
   
  ! Check the minimum number of neighbors are present
  IF(MINVAL(n_nneigh(:)) .LT. n_q12_neigh .AND. n_q12_neigh .GT. 0) THEN
    WRITE(*,'(A,I5,A,I3,A)') 'FATAL ERROR: ', MINLOC(n_nneigh(:)), ' HAS LESS THAN ', n_q12_neigh,' NEAREST NEIGHBORS IN SUBROUTINE Q12, ADJUST Q12_CUTOFF'
    STOP
  ENDIF

  ! Check that the maximum number of neighbors is not exceeded
  IF(MAXVAL(n_nneigh(:)) .GT. 75) THEN
    WRITE(*,'(A,I5,A)') 'FATAL ERROR: ', MAXLOC(n_nneigh(:)), ' HAS MORE THAN 75 NEAREST NEIGHBORS IN SUBROUTINE Q12'
    STOP
  ENDIF

  ! Sort the neighbors using insertion sort if using only the closest n_Q12_neigh
  ! neighbors
  IF(n_q12_neigh .GT. 0 ) THEN
    DO iatom = 1,n_atoms
      DO k = 2,n_nneigh(iatom)
        nneigh_temp = nneigh_lst(k,iatom)
        rij_temp(:) = rij_lst(:,k,iatom)
        DO j = k-1,1,-1
          IF(rij_lst(4,j,iatom) .LE. rij_temp(4)) GOTO 10
          nneigh_lst(j + 1, iatom) = nneigh_lst(j,iatom)
          rij_lst(:,j + 1,iatom) = rij_lst(:, j, iatom)
        ENDDO
        jatom = 0
        10 CONTINUE
        nneigh_lst(j + 1, iatom) = nneigh_temp
        rij_lst(:,j + 1, iatom) = rij_temp(:)
      ENDDO
      n_nneigh(iatom) = n_q12_neigh
    ENDDO
  ENDIF

  ! Calculate q12 and dij 
  Q12iatom = 0.0d0 
  dij = 0.0d0  
  sum_q12 = 0.0d0 
  ! Loop over the atoms
  DO iatom = 1,n_atoms
    iindx = 3*(iatom-1)

    ! Zero harmonic sum
    sum_Y = 0.0d0

    ! Loop over the nearest neighbors of iatom
    DO j = 1, n_nneigh(iatom)
      jatom = nneigh_lst(j,iatom)
      jindx = 3*(jatom-1)

      ! Get separation vector for list
      dxij(1) = rij_lst(1,j,iatom)
      dxij(2) = rij_lst(2,j,iatom)
      dxij(3) = rij_lst(3,j,iatom)
      rij = rij_lst(4,j,iatom)

      ! Calculate cos(theta) and sin(theta)
      cost(1) = dxij(3)/rij
      sint(1) = (1.0d0-cost(1)*cost(1))**0.5d0

      ! Calculate cos(phi) and sin(phi)
      IF (dABS(sint(1)) .LT. 1.0d-9) THEN
        cosp(1) = 1.0d0
        sinp(1) = 0.0d0
      ELSE
        cosp(1) = dxij(1)/(rij*sint(1))
        sinp(1) = dxij(2)/(rij*sint(1))
      ENDIF

      ! Calculate powers of sin(theta)
      sint(2) = sint(1)*sint(1)
      sint(3) = sint(2)*sint(1)
      sint(4) = sint(3)*sint(1)
      sint(5) = sint(4)*sint(1) 
      sint(6) = sint(5)*sint(1) 
      sint(7) = sint(6)*sint(1) 
      sint(8) = sint(7)*sint(1) 
      sint(9) = sint(8)*sint(1)
      sint(10) = sint(9)*sint(1)
      sint(11) = sint(10)*sint(1)
      sint(12) = sint(11)*sint(1)

      ! calculate powers of cos(theta)
      cost(2) = cost(1)*cost(1)
      cost(3) = cost(2)*cost(1)
      cost(4) = cost(3)*cost(1)
      cost(5) = cost(4)*cost(1)
      cost(6) = cost(5)*cost(1)
      cost(7) = cost(6)*cost(1)
      cost(8) = cost(7)*cost(1)
      cost(9) = cost(8)*cost(1)
      cost(10) = cost(9)*cost(1)
      cost(11) = cost(10)*cost(1)
      cost(12) = cost(11)*cost(1)


      ! Calculate sin/cos (2phi)
      sinp(2) = 2.0d0*sinp(1)*cosp(1)
      cosp(2) = 2.0d0*cosp(1)**2 - 1.0d0

      ! Calculate sin/cos (3phi)
      sinp(3) = sinp(2)*cosp(1) + cosp(2)*sinp(1)
      cosp(3) = cosp(2)*cosp(1) - sinp(2)*sinp(1)

      ! calculate sin/cos (4phi)

      sinp(4) = sinp(3)*cosp(1) + cosp(3)*sinp(1)
      cosp(4) = cosp(3)*cosp(1) - sinp(3)*sinp(1)

      ! calculate sin/cos (5phi)

      sinp(5) = sinp(4)*cosp(1) + cosp(4)*sinp(1)
      cosp(5) = cosp(4)*cosp(1) - sinp(4)*sinp(1)

      ! calculate sin/cos (6phi)
      sinp(6) = sinp(5)*cosp(1) + cosp(5)*sinp(1)
      cosp(6) = cosp(5)*cosp(1) - sinp(5)*sinp(1)

      ! calculate sin/cos (7phi)

      sinp(7) = sinp(6)*cosp(1) + cosp(6)*sinp(1)
      cosp(7) = cosp(6)*cosp(1) - sinp(6)*sinp(1)

      ! calculate sin/cos (8phi)

      sinp(8) = sinp(7)*cosp(1) + cosp(7)*sinp(1)
      cosp(8) = cosp(7)*cosp(1) - sinp(7)*sinp(1)

      ! calculate sin/cos (9phi)

      sinp(9) = sinp(8)*cosp(1) + cosp(8)*sinp(1)
      cosp(9) = cosp(8)*cosp(1) - sinp(8)*sinp(1)

      ! Calculate sin/cos (10phi)
      sinp(10) = sinp(9)*cosp(1) + cosp(9)*sinp(1)
      cosp(10) = cosp(9)*cosp(1) - sinp(9)*sinp(1)

      ! Calculate sin/cos (11phi)

      sinp(11) = sinp(10)*cosp(1) + cosp(10)*sinp(1)
      cosp(11) = cosp(10)*cosp(1) - sinp(10)*sinp(1)

      ! Calculate sin/cos (12phi)

      sinp(12) = sinp(11)*cosp(1) + cosp(11)*sinp(1)
      cosp(12) = cosp(11)*cosp(1) - sinp(11)*sinp(1)


      ! pre and post factor: 

      ! Update the harmonic sum
      ! From Mathmatica function: "SphericalHarmonicsY[12,m,theta, phi]"	

      pre_and_post(-12) = 0.566266663742191d0 * sint(12) 
      pre_and_post(-11) = 2.77412876903310d0 * sint(11)*cost(1)  
      pre_and_post(-10) = 0.409022972331817d0 * sint(10)*(-1.0d0+23.0d0*cost(2))  
      pre_and_post(-9) = 1.10763944520068d0 * sint(9)*cost(1)*(-3.0d0+23.0d0*cost(2))  
      pre_and_post(-8) = 0.362560114310785d0 * sint(8)*(1.0d0-42.0d0*cost(2)+161.0d0*cost(4))  
      pre_and_post(-7) = 0.725120228621570d0 * sint(7) * cost(1)*(5.0d0-70.0d0*cost(2)+161.0d0*cost(4))
      pre_and_post(-6) = 0.06791373178178367d0 * sint(6)*(-5.0d0+285.0d0*cost(2)-1995.0d0*cost(4)+3059.0d0*cost(6))
      pre_and_post(-5) = 0.762329748554085d0 * sint(5) * cost(1)*(-5.0d0 + 95.0d0 *cost(2)-399.0d0 *cost(4)+437.0d0 *cost(6))
      pre_and_post(-4) = 0.06536923664453508d0 * sint(4) * (5.0d0 - 340.0d0 *cost(2)+3230.0d0 *cost(4)-9044.0d0 *cost(6)+7429.0d0 *cost(8))
      pre_and_post(-3) = 0.08715898219271344d0 * sint(3) * cost(1)*(45.0d0 -1020.0d0 *cost(2)+ 5814.0d0 *cost(4)-11628.0d0 *cost(6)+ 7429.0d0 *cost(8))
      pre_and_post(-2) = 0.106747516436237d0 * sint(2) *(-3.0d0 +225.0d0 *cost(2)-2550.0d0 *cost(4)+ 9690.0d0 *cost(6)-14535.0d0 *cost(8)+7429.0d0 *cost(10))
      pre_and_post(-1) = 0.01720392001939924d0 * sint(1) * cost(1)*(-231.0d0 +5775.0d0 *cost(2)-39270.0d0 *cost(4)+106590.0d0 *cost(6)-124355.0d0 *cost(8)+52003.0d0 *cost(10))
      pre_and_post(0) = 0.001377415975458390d0 * (231.0d0 -18018.0d0 *cost(2)+225225.0d0 *cost(4)-1021020.0d0 *cost(6)+2078505.0d0 *cost(8)-1939938.0d0 *cost(10)+676039.0d0 *cost(12))  
      pre_and_post(1) = -pre_and_post(-1)
      pre_and_post(2) = pre_and_post(-2)
      pre_and_post(3) = -pre_and_post(-3)
      pre_and_post(4) = pre_and_post(-4)
      pre_and_post(5) = -pre_and_post(-5)
      pre_and_post(6) = pre_and_post(-6)
      pre_and_post(7) = -pre_and_post(-7)
      pre_and_post(8) = pre_and_post(-8)
      pre_and_post(9) = -pre_and_post(-9)
      pre_and_post(10) = pre_and_post(-10)
      pre_and_post(11) = -pre_and_post(-11)
      pre_and_post(12) = pre_and_post(-12)
      

      Y(-12) = dCMPLX(cosp(12),-sinp(12))* pre_and_post(-12) 
      Y(-11) = dCMPLX(cosp(11),-sinp(11))* pre_and_post(-11)
      Y(-10) = dCMPLX(cosp(10),-sinp(10))* pre_and_post(-10)
      Y(-9) = dCMPLX(cosp(9),-sinp(9))* pre_and_post(-9)
      Y(-8) = dCMPLX(cosp(8),-sinp(8))* pre_and_post(-8) 
      Y(-7) = dCMPLX(cosp(7),-sinp(7))* pre_and_post(-7)
      Y(-6) = dCMPLX(cosp(6),-sinp(6))* pre_and_post(-6) 
      Y(-5) = dCMPLX(cosp(5),-sinp(5))* pre_and_post(-5)
      Y(-4) = dCMPLX(cosp(4),-sinp(4))* pre_and_post(-4)
      Y(-3) = dCMPLX(cosp(3),-sinp(3))* pre_and_post(-3) 
      Y(-2) = dCMPLX(cosp(2),-sinp(2))* pre_and_post(-2) 
      Y(-1) = dCMPLX(cosp(1),-sinp(1))* pre_and_post(-1) 
      Y(0) = pre_and_post(0)
      Y(1) = dCMPLX(cosp(1),sinp(1))* pre_and_post(1)
      Y(2) = dCMPLX(cosp(2),sinp(2))* pre_and_post(2)
      Y(3) = dCMPLX(cosp(3),sinp(3))* pre_and_post(3)
      Y(4) = dCMPLX(cosp(4),sinp(4))* pre_and_post(4)
      Y(5) = dCMPLX(cosp(5),sinp(5))* pre_and_post(5) 
      Y(6) = dCMPLX(cosp(6),sinp(6))* pre_and_post(6) 
      Y(7) = dCMPLX(cosp(7),sinp(7))* pre_and_post(7) 
      Y(8) = dCMPLX(cosp(8),sinp(8))* pre_and_post(8) 
      Y(9) = dCMPLX(cosp(9),sinp(9))* pre_and_post(9) 
      Y(10) = dCMPLX(cosp(10),sinp(10))* pre_and_post(10) 
      Y(11) = dCMPLX(cosp(11),sinp(11))* pre_and_post(11) 
      Y(12) = dCMPLX(cosp(12),sinp(12))* pre_and_post(12) 

      sum_Y = sum_Y + Y 

    END DO 

    sum_q12(:, iatom) = sum_Y(:) / DBLE(n_nneigh(iatom)) 

  END DO 
  
  ! compute coarse-grained Q12  

  Q12_cg = 0.0d0 
  
  DO iatom = 1, n_atoms 

    ! sum of q12 bond order includes particle i itsself
    sum_neighbor_qlm = sum_q12(:,iatom) 

    DO j = 1, n_nneigh(iatom) 
 
        jatom = nneigh_lst(j,iatom) 
        
        sum_neighbor_qlm = sum_neighbor_qlm + sum_q12(:,jatom) 
   
    END DO 
    
    Q12_cg(:, iatom) = sum_neighbor_qlm / n_nneigh(iatom)  
      
    ! compute normalization 

    sum_Q12_norm = 0.0d0

    DO i = -12, 12

        sum_Q12_norm = sum_Q12_norm + Q12_cg(i,iatom)*dCONJG(Q12_cg(i,iatom))

    END DO 
    
    ! compute the second order invariants of coarse-grained Q12
    Q12iatom(iatom) = factor*dSQRT(sum_Q12_norm)
 
  END DO 
 
  ! Compute scalar product between Q12(i) and Q12(j) 
  
  nxtl = 0 

  n_dij_bonds = 0 

  DO iatom = 1, n_atoms 

    Q12_i_norm = factor * Q12_cg(:, iatom) / (Q12iatom(iatom)) 

    DO j = 1, n_nneigh(iatom)  

      jatom = nneigh_lst(j,iatom) 

      Q12_j_norm = factor * Q12_cg(:, jatom) / (Q12iatom(jatom))
   
      scalar_Q12i_Q12j = 0.0d0 

      DO m = -12, 12 

        scalar_Q12i_Q12j = scalar_Q12i_Q12j + Q12_i_norm(m) * dconjg(Q12_j_norm(m))

      END DO 
      
      ! dij for checking distribution
      dij(j,iatom) = scalar_Q12i_Q12j 
 
      if (scalar_Q12i_Q12j > 0.75) then 
        
        n_dij_bonds(iatom) = n_dij_bonds(iatom) + 1
            
      end if  
    
    END DO

    ! if at least 12 bonds found in particle i, then its crystal-like 
    if (n_dij_bonds(iatom) >= n_bonds) then 

        nxtl = nxtl + 1 

    end if 

  END DO
  
  ! Smaller cutoff may be used to determine if two particles are in a clusters 
  ! Rebuild neighbor list using smaller cutoff to search particle in a cluster 
  DO iatom = 1, n_atoms

    count_nb = 0

    temp_nb_lst = 0

    ! Note added ------- Jingxiang Guo (8/11/2020)

    ! Scenario 1: If the number of nearest neighbors is specified,
    ! only the nearest neighbors will be considered in
    ! a cluster. "n_nneigh" contains number of nearest neighbors for
    ! each particle  

    ! Scenario 2: If a cutoff is given without nearest neighbors specified,
    ! loop over all neighbors to find the separation smaller than 
    ! the crystal cutoff. "n_nneigh" contains the number of neighbors 
    ! that is superset for clustering assignment

    DO j = 1, n_nneigh(iatom)

      dist = rij_lst(4, j, iatom)

      IF (dist <= crys_cutoff) THEN  

        count_nb = count_nb + 1

        temp_nb_lst(count_nb) = nneigh_lst(j, iatom) 

      END IF

    END DO 

    n_nneigh(iatom) = count_nb

    nneigh_lst(:, iatom) = temp_nb_lst 

  END DO 
   
  ! Initialize variables
  largest_cluster = 0
  InCluster = .FALSE. !FALSE if atom i not in a cluster 
  n_clusters = 0
  cluster_size = 0
  cluster_lst = 0
  cluster_ID = 0

  !Loop over the atoms
  
  DO iatom = 1,n_atoms

    clusterID = 0

    IF (n_dij_bonds(iatom) .LT. n_bonds) CYCLE

      InCluster(iatom) = .TRUE.
  
      !Search iatoms's n.n. to see if they are in cluster(s)
      DO j = 1, n_nneigh(iatom)
        jatom = nneigh_lst(j,iatom)
 
        IF (InCluster(jatom) .EQV. .TRUE.) THEN !one of iatom's n.n. is in a cluster

          IF (clusterID .NE. 0) THEN !another one of iatoms's n.n. is also in a cluster

            IF (clusterID .NE. cluster_ID(jatom)) THEN !two of iatoms's n.n. belong to "different" clusters...combine

              !keep lower cluster ID number
              IF (clusterID .LT. cluster_ID(jatom)) THEN
                clusterID_old = cluster_ID(jatom)
              ELSE
                clusterID_old = clusterID
                clusterID = cluster_ID(jatom)
              ENDIF

              ! Merge clusters, keeping lower cluster ID number
              DO k = 1,cluster_size(clusterID_old)
                cluster_lst(cluster_size(clusterID) + k,clusterID) = cluster_lst(k,clusterID_old)
                cluster_ID(cluster_lst(k,clusterID_old)) = clusterID
              ENDDO
              cluster_size(clusterID) = cluster_size(clusterID) + cluster_size(clusterID_old)
              cluster_lst(:,clusterID_old) = 0 ! Clear old cluster data, reduce count
              cluster_size(clusterID_old) = 0     
              n_clusters = n_clusters - 1      

              !Place LAST cluster in cluster list into newly empty slot
              IF (clusterID_old .NE. (n_clusters+1)) THEN ! Must do move
                DO k = 1,cluster_size(n_clusters+1)
                 cluster_lst(k,clusterID_old) = cluster_lst(k,(n_clusters+1))
                 cluster_ID(cluster_lst(k,clusterID_old)) = clusterID_old
                ENDDO
                cluster_size(clusterID_old) = cluster_size(n_clusters+1)  ! Clear data for moved cluster 
                cluster_lst(:,(n_clusters+1)) = 0 
                cluster_size(n_clusters+1) = 0    
              ENDIF
            ENDIF
       
          ELSE
            clusterID = cluster_ID(jatom)
          ENDIF

        ENDIF
      ENDDO  ! j neighbor

      IF (clusterID .EQ. 0) then !if none of iatom's n.n. are in a cluster, imol begins new cluster
        n_clusters = n_clusters + 1
        clusterID = n_clusters
      ENDIF

      !Add iatom to the cluster
      cluster_ID(iatom) = clusterID
      cluster_size(clusterID) = cluster_size(clusterID) + 1
      cluster_lst(cluster_size(clusterID),clusterID) = iatom !cluster_lst(cluster number, list of atoms in cluster)
  
  ENDDO ! iatom
  
  largest_cluster = MAXVAL(cluster_size)
  
END SUBROUTINE

SUBROUTINE calc_global_Q6(n_atoms,n_Q6_neigh,Q6_cutoff,box,x,Q6avg)
  implicit none 

  ! Passed
  INTEGER,INTENT(IN) :: n_atoms  
  INTEGER, INTENT(in) :: n_Q6_neigh
  DOUBLE PRECISION,INTENT(IN) :: Q6_cutoff 
  DOUBLE PRECISION,INTENT(IN) :: box 
  DOUBLE PRECISION,INTENT(OUT) :: Q6avg 
  DOUBLE PRECISION,INTENT(IN),DIMENSION(0:3*n_atoms-1) :: x 
 
  ! Local
  INTEGER :: i,j,k 
  INTEGER :: iatom,jatom, iindx,jindx 
  DOUBLE PRECISION,DIMENSION(3) :: xi,xj,dxij 
  DOUBLE PRECISION :: rij,rijsq

  ! Neighbor list 
  INTEGER,DIMENSION(n_atoms) :: n_nneigh 
  INTEGER,DIMENSION(75,n_atoms) :: nneigh_lst 
  DOUBLE PRECISION,DIMENSION(4,75,n_atoms) :: rij_lst 
  INTEGER :: nneigh_temp 
  DOUBLE PRECISION,DIMENSION(4) :: rij_temp 
  DOUBLE PRECISION :: Q6_cutoffsq 

  ! Spherical Harmonics 
  DOUBLE PRECISION :: pre_and_post(-6:6)  
  DOUBLE PRECISION :: cost(6),sint(6)
  DOUBLE PRECISION :: cosp(6), sinp(6)
  COMPLEX*16 :: Y(-6:6)
  COMPLEX*16 :: sum_Y(-6:6),sum_Q6(-6:6)

  
  Q6_cutoffsq = Q6_cutoff*Q6_cutoff

  ! Initialize arrays
  n_nneigh = 0
  nneigh_lst = 0
  rij  = 0.0d0

  ! Loop over iatom
  DO iatom = 1,n_atoms-1
    iindx = 3*(iatom-1)

    ! Store iatom's coordinates
    xi(1) = x(iindx)
    xi(2) = x(iindx+1)
    xi(3) = x(iindx+2)  
 
    ! Loop over jatom
    DO jatom = iatom+1,n_atoms
      jindx = 3*(jatom-1)

      ! Store jatom's coordinates
      xj(1) = x(jindx)
      xj(2) = x(jindx+1)
      xj(3) = x(jindx+2)

      ! Compute ij vector & separation distance
      dxij(:) = xj(:) - xi(:)
      dxij(:) = dxij(:) - dNINT(dxij(:)/box)*box
      rijsq = dxij(1)*dxij(1) + dxij(2)*dxij(2) + dxij(3)*dxij(3)
      
       ! If they are within the cutoff distance
       IF(rijsq .GT. Q6_cutoffsq) CYCLE

         rij = dSQRT(rijsq)
   
         ! Update information for iatom

         ! Update the number of neighbors for iatom
         n_nneigh(iatom) = n_nneigh(iatom) + 1

         ! Include jatom as a neighbor of iatom
         nneigh_lst(n_nneigh(iatom),iatom)  = jatom

         ! Store the distance components for iatom
         rij_lst(1, n_nneigh(iatom), iatom) = dxij(1)
         rij_lst(2, n_nneigh(iatom), iatom) = dxij(2)
         rij_lst(3, n_nneigh(iatom), iatom) = dxij(3)
         rij_lst(4, n_nneigh(iatom), iatom) = rij

         ! Update informtion for jatom

         ! Update the number of neighbors for jatom
         n_nneigh(jatom) = n_nneigh(jatom) + 1

         ! Include iatom as a neighbor of jatom
         nneigh_lst(n_nneigh(jatom),jatom) = iatom

         ! Store the distance components for jatom
         rij_lst(:,n_nneigh(jatom),jatom) = rij_lst(:,n_nneigh(iatom),iatom)

     ENDDO
   ENDDO
  
  ! Check the minimum number of neighbors are present
  IF(MINVAL(n_nneigh(:)) .LT. n_Q6_neigh .AND. n_Q6_neigh .GT. 0) THEN
    WRITE(*,'(A,I5,A,I3,A)') 'FATAL ERROR: ', MINLOC(n_nneigh(:)), ' HAS LESS THAN ', n_Q6_neigh,' NEAREST NEIGHBORS IN SUBROUTINE Q6, ADJUST Q6_CUTOFF'
    STOP
  ENDIF

  ! Check that the maximum number of neighbors is not exceeded
  IF(MAXVAL(n_nneigh(:)) .GT. 75) THEN
    WRITE(*,'(A,I5,A)') 'FATAL ERROR: ', MAXLOC(n_nneigh(:)), ' HAS MORE THAN 75 NEAREST NEIGHBORS IN SUBROUTINE Q6'
    STOP
  ENDIF

  ! Sort the neighbors using insertion sort if using only the closest n_Q6_neigh neighbors
  IF(n_Q6_neigh .GT. 0 ) THEN
    DO iatom = 1,n_atoms
      DO k = 2,n_nneigh(iatom)
        nneigh_temp = nneigh_lst(k,iatom)
        rij_temp(:) = rij_lst(:,k,iatom)
        DO j = k-1,1,-1
          IF(rij_lst(4,j,iatom) .LE. rij_temp(4)) GOTO 10
          nneigh_lst(j + 1, iatom) = nneigh_lst(j,iatom)
          rij_lst(:,j + 1,iatom) = rij_lst(:, j, iatom)
        ENDDO
        jatom = 0
        10 CONTINUE
        nneigh_lst(j + 1, iatom) = nneigh_temp
        rij_lst(:,j + 1, iatom) = rij_temp(:)
      ENDDO
      n_nneigh(iatom) = n_Q6_neigh
    ENDDO
  ENDIF


  ! Calculate Q6
  
  ! Initialize the sum
  
  sum_Q6 = 0.0d0

  ! Loop over the atomecules
  DO iatom = 1,n_atoms

    ! Initialize the harmonic sum
    sum_Y = 0.0d0
  
    ! Loop over the 4 nearest neighbors of iatom
    DO j = 1,n_nneigh(iatom)

      ! Calculate cos(theta) and sin(theta)
      cost(1) = rij_lst(3,j,iatom)/rij_lst(4,j,iatom)
      sint(1) = (1.0d0-cost(1)*cost(1))**0.5d0

      ! Calculate cos(phi) and sin(phi)
      IF (dABS(sint(1)) .LT. 1.0d-9) THEN
        cosp(1) = 1.0d0
        sinp(1) = 0.0d0

      ELSE
        cosp(1) = rij_lst(1,j,iatom)/(rij_lst(4,j,iatom)*sint(1))
        sinp(1) = rij_lst(2,j,iatom)/(rij_lst(4,j,iatom)*sint(1))

      END IF

      ! calculate powers of sin(theta)
      sint(2) = sint(1)*sint(1)
      sint(3) = sint(2)*sint(1)
      sint(4) = sint(3)*sint(1) 
      sint(5) = sint(4)*sint(1) 
      sint(6) = sint(5)*sint(1) 

      ! calculate powers of cos(theta)
      cost(2) = cost(1)*cost(1)
      cost(3) = cost(2)*cost(1)
      cost(4) = cost(3)*cost(1)
      cost(5) = cost(4)*cost(1)
      cost(6) = cost(5)*cost(1)


      ! calculate sin/cos (2phi)
      sinp(2) = 2.0d0*sinp(1)*cosp(1)
      cosp(2) = 2.0d0*cosp(1)**2 - 1.0d0

      ! calculate sin/cos (3phi)
      sinp(3) = sinp(2)*cosp(1) + cosp(2)*sinp(1)
      cosp(3) = cosp(2)*cosp(1) - sinp(2)*sinp(1)

      ! calculate sin/cos (4phi) 

      sinp(4) = sinp(3)*cosp(1) + cosp(3)*sinp(1) 
      cosp(4) = cosp(3)*cosp(1) - sinp(3)*sinp(1) 

      ! calculate sin/cos (5phi) 

      sinp(5) = sinp(4)*cosp(1) + cosp(4)*sinp(1) 
      cosp(5) = cosp(4)*cosp(1) - sinp(4)*sinp(1) 

      ! calculate sin/cos (6phi) 

      sinp(6) = sinp(5)*cosp(1) + cosp(5)*sinp(1) 
      cosp(6) = cosp(5)*cosp(1) - sinp(5)*sinp(1) 

      ! Update the harmonic sum
      ! From Mathmatica function: "SphericalHarmonicsY[4,m,theta, phi]"	

      ! pre and post factor: 

      pre_and_post(-6) = 0.483084113580066d0 * sint(6)  
      pre_and_post(-5) = 1.67345245810010d0 * sint(5) * cost(1)  
      pre_and_post(-4) = 0.356781262853998d0 * (11.0d0 * cost(2)- 1.0d0) * sint(4)  
      pre_and_post(-3) = 0.651390485867716d0 * (-3.0d0*cost(1) + 11.0d0 * cost(3)) * sint(3)  
      pre_and_post(-2) = 0.325695242933858d0 * (1.0d0 - 18.0d0 * cost(2) + 33.0d0* cost(4)) * sint(2)  
      pre_and_post(-1) = 0.411975516301141d0 * (5.0d0*cost(1) - 30.0d0 * cost(3) + 33.0d0*cost(5)) * sint(1)  
      pre_and_post(0) = 0.06356920226762842d0 * (105.0d0 * cost(2) - 315.0d0 * cost(4) + 231.0d0 * cost(6) - 5.0d0)
      pre_and_post(1) = -pre_and_post(-1)  
      pre_and_post(2) = pre_and_post(-2) 
      pre_and_post(3) = -pre_and_post(-3)
      pre_and_post(4) = pre_and_post(-4) 
      pre_and_post(5) = -pre_and_post(-5)
      pre_and_post(6) = pre_and_post(-6) 

      Y(-6) = dcmplx(cosp(6),-sinp(6))*pre_and_post(-6)  
      Y(-5) = dcmplx(cosp(5),-sinp(5))*pre_and_post(-5) 
      Y(-4) = dcmplx(cosp(4),-sinp(4))*pre_and_post(-4)  
      Y(-3) = dcmplx(cosp(3),-sinp(3))*pre_and_post(-3) 
      Y(-2) = dcmplx(cosp(2),-sinp(2))*pre_and_post(-2) 
      Y(-1) = dcmplx(cosp(1),-sinp(1))*pre_and_post(-1)
      Y(0) = pre_and_post(0) 
      Y(1) = dcmplx(cosp(1),sinp(1))* pre_and_post(1) 
      Y(2) = dcmplx(cosp(2),sinp(2))* pre_and_post(2)  
      Y(3) = dcmplx(cosp(3),sinp(3))* pre_and_post(3) 
      Y(4) = dcmplx(cosp(4),sinp(4))* pre_and_post(4)
      Y(5) = dcmplx(cosp(5),sinp(5))* pre_and_post(5) 
      Y(6) = dcmplx(cosp(6),sinp(6))* pre_and_post(6) 

      sum_Y = sum_Y + Y 

    END DO 
 
    sum_Q6 = sum_Q6 + sum_Y/DBLE(n_nneigh(iatom)) 

  END DO 

  ! Normalize the Q6 sum
  sum_Q6 = sum_Q6/DBLE(n_atoms)

  ! Calculate the average Q6
  Q6avg = 0.0d0
  DO i = -6,6

    Q6avg = Q6avg + sum_Q6(i)*dCONJG(sum_Q6(i))

  ENDDO
  Q6avg = dSQRT(Q6avg)
 
  END SUBROUTINE  


