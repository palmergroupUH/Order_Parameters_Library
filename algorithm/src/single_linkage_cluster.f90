!------------------------------------------- Program Descriptions ----------------------------------------

! This library contains 2 different implementations of single-linkage cluster algorithm:  
    ! Method 1: A less memory efficient method:
	! Method 2: Allen & Tildesley method based on https://github.com/Allen-Tildesley/examples 
 

! Note: 
! 1. For method 1, 
!    a. The space complexity scales as n^2. For example, if n_atoms ~ 10^5, then
!    40 GB memory is needed. 
!    b. "n_bonds" and "n_dij_bonds" are used as critertia to determine if a selected
!    particle be considered in a cluster
! 
! 2. For method 2, 
!    a. a linked-list neighbor list is rebuilt for a given cutoff, system size,
!    and particle's coordinates. The time complexity scales as n^2  


module cluster
    use system, only: c_int, c_double, dp, sp
    implicit none 

    private

    public :: largest_cluster_high_mem, &  
            & Allen_Tidesley_cluster_nb_list, & 
            & Allen_Tidesley_cluster_ID_list, & 
            & Allen_Tidesley_largest_cluster, & 
            & Allen_Tidesley_cluster, & 
            & rebuild_cluster_nblist, & 
            & cluster_search, & 
            & CLUSTER_BUILD, & 
            & CLUSTER_SEARCH_old

contains

!---------------------------------------------------------------------------------------------------
!                               Method 1:  
!---------------------------------------------------------------------------------------------------
	subroutine largest_cluster_high_mem(n_atoms,n_nneigh, nneigh_lst, n_bonds,n_dij_bonds, largest_cluster)
		implicit none 
        ! Passed
		integer,intent(in) :: n_atoms, n_bonds
		integer,intent(in),dimension(:) :: n_nneigh
		integer,intent(in),dimension(:,:) :: nneigh_lst
		integer,intent(in),dimension(:) :: n_dij_bonds	

        ! Local	
		integer :: n_clusters,n_ice
		integer :: clusterID,clusterID_old, d6ijcounter
		integer,dimension(n_atoms) :: cluster_size, cluster_ID
		integer,dimension(n_atoms,n_atoms) :: cluster_lst
		logical,dimension(n_atoms) :: InCluster, visited
		integer :: iatom,jatom,i,j,k  

        ! Return
		integer,intent(out) :: largest_cluster 
		
		InCluster = .FALSE. !FALSE if atom i not in a cluster 

		!allocate(cluster_lst(1:n_atoms,1:n_atoms)) 

		n_clusters = 0
		cluster_size = 0
		cluster_lst = 0
		cluster_ID = 0

		!Loop over the atoms
		DO iatom = 1,n_atoms

		clusterID = 0
	
        ! Modify this for crystalline criteria
        
		if (n_dij_bonds(iatom) < n_bonds) CYCLE
           
		  InCluster(iatom) = .TRUE.

		  !Search iatoms's n.n. to see if they are in cluster(s)
		  DO j = 1, n_nneigh(iatom)
			jatom = nneigh_lst(j,iatom)

			if (InCluster(jatom) .EQV. .TRUE.) THEN !one of iatom's n.n. is in a cluster

			  if (clusterID .NE. 0) THEN !another one of iatoms's n.n. is also in a cluster

				if (clusterID .NE. cluster_ID(jatom)) THEN !two of iatoms's n.n. belong to "different" clusters...combine

				  !keep lower cluster ID number
				  if (clusterID .LT. cluster_ID(jatom)) THEN
					clusterID_old = cluster_ID(jatom)
				  ELSE
					clusterID_old = clusterID
					clusterID = cluster_ID(jatom)
				  end if

				  ! Merge clusters, keeping lower cluster ID number
				  DO k = 1,cluster_size(clusterID_old)
					cluster_lst(cluster_size(clusterID) + k,clusterID) = cluster_lst(k,clusterID_old)
					cluster_ID(cluster_lst(k,clusterID_old)) = clusterID
				  end DO
				  cluster_size(clusterID) = cluster_size(clusterID) + cluster_size(clusterID_old)
				  cluster_lst(:,clusterID_old) = 0 ! Clear old cluster data, reduce count
				  cluster_size(clusterID_old) = 0     
				  n_clusters = n_clusters - 1      

				  !Place LAST cluster in cluster list into newly empty slot
				  if (clusterID_old .NE. (n_clusters+1)) THEN ! Must do move
					DO k = 1,cluster_size(n_clusters+1)
					 cluster_lst(k,clusterID_old) = cluster_lst(k,(n_clusters+1))
					 cluster_ID(cluster_lst(k,clusterID_old)) = clusterID_old
					end DO
					cluster_size(clusterID_old) = cluster_size(n_clusters+1)  ! Clear data for moved cluster 
					cluster_lst(:,(n_clusters+1)) = 0 
					cluster_size(n_clusters+1) = 0    
				  end if
				end if
		   
			  ELSE
				clusterID = cluster_ID(jatom)
			  end if

			end if

		  end do  ! j neighbor

		  if (clusterID .EQ. 0) then !if none of iatom's n.n. are in a cluster, imol begins new cluster
			n_clusters = n_clusters + 1
			clusterID = n_clusters
		  end if

		  !Add iatom to the cluster
		  cluster_ID(iatom) = clusterID
		  cluster_size(clusterID) = cluster_size(clusterID) + 1
		  cluster_lst(cluster_size(clusterID),clusterID) = iatom !cluster_lst(cluster number, list of atoms in cluster)

		end do ! iatom

		largest_cluster = maxval(cluster_size)

		end subroutine

!---------------------------------------------------------------------------------------------------
!                               Method 2:  
!---------------------------------------------------------------------------------------------------

    pure logical function two_atoms_are_bonded(r1,r2,cutoff_sqr,box)
        implicit none
        ! Passed
        real(dp),intent(in) :: cutoff_sqr
        real(dp),intent(in),dimension(:) :: r1,r2,box

        ! Local
        integer :: i,j
        real(dp), dimension(1:size(r1)) :: r12
        real(dp) :: r12_sqr

        r12 = r1 - r2

        ! Periodical boundary condition

        r12 = r12 - dnint(r12/box)*box

        r12_sqr = r12(1)*r12(1) + r12(2)*r12(2) + r12(3)*r12(3)

        ! Return
        two_atoms_are_bonded = r12_sqr <= cutoff_sqr

        end function

    pure subroutine Allen_Tidesley_cluster_nb_list(total_atoms, &
                                                  & coord, &
                                                  & box, &
                                                  & cutoff_sqr, &
                                                  & clinklist)
        implicit none
        integer, intent(in) :: total_atoms
        real(dp),intent(in),dimension(:) :: box
        real(dp),intent(in),dimension(:,:) :: coord
        real(dp),intent(in) :: cutoff_sqr
        integer :: i,j,k
        integer,intent(out),dimension(:), allocatable :: clinklist

        if ( .not. allocated(clinklist)) allocate(clinklist(1:total_atoms))
    
        clinklist = [(i,i=1,total_atoms)]

        do i = 1,total_atoms - 1

            ! if particle i points to its self
            if (i == clinklist(i)) then
    
                ! set the pointer to i 
                j = i

                do

                    do k = i+1, total_atoms
    
                        ! k is not in any cluster 
                        if (clinklist(k) ==k) then

                            ! see if j and k are in cluster 
                            if (two_atoms_are_bonded(coord(:,j), &
                                                   & coord(:,k), &
                                                   & cutoff_sqr, &
                                                   & box)) then 

                                clinklist([k,j]) = clinklist([j,k])

                            end if

                        end if

                    end do

                    j = clinklist(j)

                    if ( j == i ) exit

                end do

            end if

        end do

        end subroutine

    subroutine Allen_Tidesley_cluster_ID_list(cnb_list,total_atoms,cluster_ID_list,cluster_member_head)
        implicit none
        ! Passed
        integer, intent(in) :: total_atoms
        integer, intent(in), dimension(:) :: cnb_list

        ! Local
        integer,dimension(1:total_atoms) :: head_in_cluster,done
        integer :: i,j,cluster_ID

        ! Return
        integer, intent(out), dimension(:), allocatable :: cluster_ID_list
        integer, intent(out), dimension(:),allocatable  :: cluster_member_head

        if ( .not. allocated( cluster_ID_list)) allocate(cluster_ID_list(1:total_atoms))

        cluster_ID_list = 0

        cluster_ID = 0

        head_in_cluster = 0

        ! Loop over all clusters
        do

            if (all(cluster_ID_list > 0)) exit

            i = minloc(cluster_ID_list,dim=1 )

            cluster_ID = cluster_ID + 1

            j = i

            head_in_cluster(cluster_ID) = j

            cluster_ID_list(j) = cluster_ID

            ! Loop over all atoms in a cluster
            do

                j = cnb_list(j)

                if (j==i) exit

                cluster_ID_list(j) = cluster_ID

            end do

        end do

        if (.not. allocated(cluster_member_head)) allocate(cluster_member_head(1:cluster_ID))

        cluster_member_head = 0

        cluster_member_head = head_in_cluster(1:cluster_ID)

        end subroutine

    subroutine Allen_Tidesley_largest_cluster(cluster_ID_ary, largest_cluster) 
        implicit none 
        ! Passed
        integer, intent(in), dimension(:) :: cluster_ID_ary

        ! Local 
        integer :: i, num_clusters 
        integer, dimension(:), allocatable :: largest_cluster_lst

        ! Return
        integer, intent(out) :: largest_cluster 

        ! The maximum value is equal to number of clusters 
        num_clusters = maxval(cluster_ID_ary)
       
        ! Allocate cluster size array  
        allocate(largest_cluster_lst(1:num_clusters)) 

        ! For each cluster 
        do i = 1, num_clusters  

            ! get the size of cluster 
            largest_cluster_lst(i) = count(cluster_ID_ary == i)

        end do 

        ! The maximum is the largest cluster size
        largest_cluster = maxval(largest_cluster_lst)
       
        end subroutine 

   
    subroutine Allen_Tidesley_cluster(total_atoms, &
                                    & coord, &
                                    & box, &
                                    & cutoff_sqr, &
                                    & largest_cluster) 
        implicit none 
        ! Passed
        integer, intent(in) :: total_atoms
        real(dp), intent(in), dimension(:) :: box
        real(dp), intent(in), dimension(:,:) :: coord
        real(dp), intent(in) :: cutoff_sqr
        integer :: i,j,k
    
        ! Local
        integer, dimension(:), allocatable :: clinklist
        integer, dimension(:), allocatable :: cluster_ID_list
        integer, dimension(:), allocatable  :: cluster_member_head

        ! Return
        integer, intent(out) :: largest_cluster 

        if (total_atoms == 0 .or. total_atoms ==1 ) then 

            largest_cluster = total_atoms

        else

            call Allen_Tidesley_cluster_nb_list(total_atoms, &
                                              & coord, &
                                              & box, &
                                              & cutoff_sqr, &
                                              & clinklist)

            call Allen_Tidesley_cluster_ID_list(clinklist, &
                                              & total_atoms, &
                                              & cluster_ID_list, &
                                              & cluster_member_head) 

            
            call Allen_Tidesley_largest_cluster(cluster_ID_list, & 
                                              & largest_cluster)

        end if

        end subroutine 

    ! Method 3: a more memory efficient version
    subroutine rebuild_cluster_nblist(maxnb, num_cl, num_NB_list, cutoff_cl, Rij, NB_LIST, select_cl_id)
        implicit none

        ! Passed
        integer, intent(in) :: maxnb 
        integer, intent(in) :: num_cl
        real(dp), intent(in) :: cutoff_cl
        real(dp), intent(in), dimension(:,:,:) :: Rij 
        integer, intent(in), dimension(:) :: select_cl_id

        ! Local:
        integer :: i, j, count_nb, select_i, j_nb  
        integer, dimension(1:maxnb) :: temp_nb_list  

        ! Return
        integer, intent(inout), dimension(:) :: num_NB_list
        integer, intent(inout), dimension(:,:) :: NB_list

        do i = 1, num_cl

            select_i = select_cl_id(i)

            count_nb = 0

            temp_nb_list = 0

            do j = 1, maxnb  

                j_nb = NB_LIST(j, select_i)

                if (any(j_nb == select_cl_id)) then

                    if (Rij(4, j, select_i) <= cutoff_cl) then

                        count_nb = count_nb + 1

                        temp_nb_list(count_nb) = j_nb

                    end if 

                end if 

            end do

            NB_list(:, select_i) = temp_nb_list

            num_NB_list(select_i) = count_nb

        end do 

        end subroutine 

    subroutine cluster_search(maxnb, num_atom, NB_list, num_NB_list, atom_id, largest_cluster)
        implicit none 
        ! Passed
        integer, intent(in) :: maxnb 
        integer, intent(in) :: num_atom 
        integer, intent(in), dimension(:) :: num_NB_list  
        integer, intent(in), dimension(:,:) :: NB_list
        
        ! Local
        integer, dimension(1:maxnb) :: neighbor_of_neighbors
        integer :: picked_atom, this_cluster_size, num_nb
        integer :: pick_nb, num_cluster, num_atom_left, nb_counter 
        
        ! declare an array containing atom ids in a single cluster
        integer, dimension(1:num_atom) :: current_cluster 
        integer, dimension(1:num_atom) :: Cluster_atom_ID 
        integer, dimension(1:num_atom) :: num_atom_in_cluster 
        real(dp) :: start, finish 

        ! Return
        
        ! a list of atom ids to be considered in clusters
        integer, intent(inout), dimension(1:num_atom) :: atom_id
        integer, intent(out) :: largest_cluster 

        num_cluster = 0 
        num_atom_in_cluster = 0
        num_atom_left = 0 
        Cluster_atom_ID = 0     

        ! loop over all possible clusters (terminate if all particles are considered)  

        call cpu_time(start) 
        do while (num_atom_left < num_atom)
    
            ! clear up current cluster
            current_cluster = 0

            ! Pick an initial atom to begin
            picked_atom  = minval(atom_id, atom_id /=0) 
            
            ! initially picked atom should be in current cluster
            current_cluster(1) = picked_atom 

            ! if this atom has no neighbor 
            if (num_NB_list(picked_atom) == 0) then
    
                this_cluster_size = 1 

            else 

                ! all neighbors of initially picked atoms should be in a same cluster

                current_cluster(2:num_NB_list(picked_atom)+1) = NB_list(1:num_NB_list(picked_atom),picked_atom)

                nb_counter = 1

                ! loop over a particle neighbors' neighbors in a cluster
                do while (nb_counter < count((current_cluster /= 0)))

                    ! clear up neighbors' neighbor
                    neighbor_of_neighbors = 0

                    nb_counter = nb_counter + 1

                    ! pick a neighbor except the initially picked atom from the current cluster

                    pick_nb = current_cluster(nb_counter)

                    ! find its number of neighbors

                    num_nb = num_NB_list(pick_nb)

                    ! find its neighbor:

                    neighbor_of_neighbors = NB_list(1:num_nb, pick_nb)

                    call add_members_to_cluster(current_cluster, num_nb, neighbor_of_neighbors) 

                end do 

                
                this_cluster_size = count(current_cluster /= 0)
                
            end if 

            ! increase the number of cluster
            num_cluster = num_cluster + 1 

            ! update the number of atoms in a cluster
            num_atom_in_cluster(num_cluster) = this_cluster_size 

            ! put all atom id in cluster in a large array containing all cluster's atom id  
            !Cluster_atom_ID(num_atom_left + 1: num_atom_left + this_cluster_size) = current_cluster(1:this_cluster_size)  

            ! number of atoms left over  ( have not been considered in any cluster yet)
            num_atom_left = sum(num_atom_in_cluster) 
            
            ! in next iteration, you can pick another an particle
            ! not in any current clusters from "atom_id"
            call remove_atoms_already_in_cluster(num_atom, atom_id, current_cluster) 
           
        end do  
    
        largest_cluster = maxval(num_atom_in_cluster) 
        print*, "largest_cluster,", largest_cluster 
        end subroutine 

    subroutine add_members_to_cluster(current_cluster, num_nb, neighbor_of_neighbors) 
        implicit none 

        ! Passed
        integer, intent(in) :: num_nb 
        integer, intent(in), dimension(:) :: neighbor_of_neighbors

        ! Local
        integer :: i, j, num_memebers, counter 

        ! Return
        integer, intent(inout), dimension(:) :: current_cluster 

        num_memebers = count(current_cluster /= 0) 

        ! loop over each neighbors
        do i = 1, num_nb 

            counter = 0

            if (all(neighbor_of_neighbors(i) /= current_cluster)) then  

                current_cluster(num_memebers+1) = neighbor_of_neighbors(i)

                num_memebers = num_memebers + 1 
             
            end if 

        end do 

        end subroutine 

    subroutine remove_atoms_already_in_cluster(num_atoms, atom_id, current_cluster)
        implicit none 
        integer, intent(in) :: num_atoms
        integer, intent(in), dimension(:) :: current_cluster        
        integer, intent(inout), dimension(:) :: atom_id  
        integer :: i, j 

        do j = 1, num_atoms  
        
            ! if atom id appears in current_cluster
            ! remove it by setting to be 0
            if (any(current_cluster == atom_id(j))) then 
                
                atom_id(j) = 0 

            end if 

        end do

        end subroutine 

    SUBROUTINE CLUSTER_BUILD(num_NB_list,cutoff_cl,Rij,NB_LIST_CL,num_NB_list_CL,NB_LIST,cluster_ID)
            IMPLICIT NONE
            REAL(8),INTENT(IN) :: cutoff_cl
            INTEGER,INTENT(IN),DIMENSION(:,:) :: NB_LIST
            REAL(8),INTENT(IN),DIMENSION(:,:,:) :: Rij
            INTEGER,INTENT(IN),DIMENSION(:) :: num_NB_list
            INTEGER,INTENT(IN),DIMENSION(:),OPTIONAL :: cluster_ID
            INTEGER,INTENT(OUT),DIMENSION(1:SIZE(num_NB_list)) :: num_NB_list_CL
            INTEGER,INTENT(OUT),DIMENSION(1:SIZE(NB_list,DIM=1),1:SIZE(num_NB_list)) :: NB_LIST_CL
            INTEGER :: NB,CID,i,j,counter,num_cl

            ! Build cluster neighbor list only for certain molecuels
            IF (PRESENT(cluster_ID)) THEN
                num_NB_list_CL = 0

                NB_LIST_CL = 0

                counter = 0

                DO i = 1,COUNT(cluster_ID /=0)

                    CID = cluster_ID(i)

                    num_cl = 0
                    DO j = 1, 40 
                        NB = NB_LIST(j,CID)

                        !Find the position of this nb in the ice label list
                        !(should be non-zero value if exist)
                        ! Check if this NB is in the ice label list
                        IF ( ANY(NB-cluster_ID ==0) ) THEN

                            ! Then check if this NB is within the cluster cutoff
                            ! distance
                            IF (Rij(4,j,CID) <= cutoff_cl) THEN

                                ! If yes, then create new neighbor list for
                                ! cluster
                                num_cl = num_cl + 1
                                NB_LIST_CL(num_cl,CID) = NB
                            END IF
                        END IF
                    END DO

                    num_NB_list_CL(CID) = num_cl

                END DO
            
            ELSE
                ! Assume all molecules are considered in Cluster
                num_NB_list_CL = num_NB_list
                NB_LIST_CL = NB_LIST

            END IF

    END SUBROUTINE CLUSTER_BUILD

    SUBROUTINE CLUSTER_SEARCH_old(num_atom,NB_list,num_NB_list,AID)
            IMPLICIT NONE

            INTEGER,INTENT(IN) :: num_atom
            INTEGER,INTENT(IN),DIMENSION(:) :: num_NB_list
            INTEGER,INTENT(IN),DIMENSION(:,:) :: NB_list
            INTEGER,DIMENSION(1:SIZE(NB_LIST,DIM=1)) :: CB
            INTEGER :: i,j,k,l,f,ica,nb,num_nb,this_cluster_size,num_cluster,num_atom_counter,counter
            INTEGER,DIMENSION(1:SIZE(num_NB_list)) :: c_current
            INTEGER,INTENT(INOUT),DIMENSION(:) :: AID
            INTEGER,DIMENSION(1:SIZE(num_NB_list)) :: Cluster_atom_ID
            INTEGER,DIMENSION(1:SIZE(num_NB_list)) :: num_atom_cluster
                num_cluster = 0
                num_atom_cluster  = 0
                num_atom_counter = 0
                Cluster_atom_ID = 0

                DO WHILE (num_atom_counter < num_atom)

                        c_current = 0

                        ica  = MINVAL(AID,AID/=0) ! Pick initial  atom to begin

                        c_current(1) = ica

                        IF ( num_NB_list(ica) == 0 ) THEN ! some molecules have no neighbors. They are in its own cluster

                                this_cluster_size = 1
                        ELSE

                            c_current(2:num_NB_list(ica)+1) = NB_list(1:num_NB_list(ica),ica)

                            counter = 1

                            DO WHILE (counter < COUNT((c_current /= 0)))

                                    CB = 0

                                    counter = counter + 1

                                    nb = c_current(counter)

                                    num_nb = num_NB_list(nb)

                                    CB = NB_list(:,nb)

                                    CALL CLUSTER_COMBINE_ARRAY(c_current,num_nb,CB)

                            END DO

                        this_cluster_size = COUNT(c_current /= 0)
                        END IF
                        num_cluster = num_cluster + 1

                        num_atom_cluster(num_cluster) = num_atom_cluster(num_cluster) + this_cluster_size

                        !Cluster_atom_ID(num_atom_counter + 1: num_atom_counter + 1+ this_cluster_size) = c_current !c_current(1:this_cluster_size)

                        num_atom_counter = SUM(num_atom_cluster)

                        CALL CLUSTER_DIFF_ARRAY(AID,c_current)
                END DO
                print*, maxval(num_atom_cluster) 

                END SUBROUTINE
            SUBROUTINE CLUSTER_COMBINE_ARRAY(c_current,num_nb,CB)

                IMPLICIT NONE
                INTEGER,INTENT(IN) :: num_nb
                INTEGER,INTENT(INOUT),DIMENSION(:) :: c_current
                INTEGER :: k,l,counter,D,E
                INTEGER,INTENT(IN),DIMENSION(:) :: CB

                DO k = 1,num_nb
                        counter = 0
                        DO l = 1,COUNT(c_current /= 0)
                                D = CB(k)  - c_current(l)
                                IF ( D==0 ) THEN
                                        CYCLE
                                ELSE
                                        counter = counter + 1
                                END IF
                        END DO

                        IF (counter == COUNT((c_current /=0))) THEN

                                c_current(COUNT(c_current /=0)+1) = CB(k)
                        END IF
                END DO

                END SUBROUTINE CLUSTER_COMBINE_ARRAY

        SUBROUTINE CLUSTER_DIFF_ARRAY(AID,c_current)
                IMPLICIT NONE
                INTEGER,INTENT(INOUT),DIMENSION(:) :: AID
                INTEGER,INTENT(IN),DIMENSION(:) :: c_current
                INTEGER :: m,n
                DO m = 1,COUNT(c_current /= 0 )
                        DO n= 1,SIZE(AID )
                                IF (c_current(m) == AID(n)) THEN
                                        AID(n) = 0
                                END IF
                        END DO
                END DO

            END SUBROUTINE


 

end module  

