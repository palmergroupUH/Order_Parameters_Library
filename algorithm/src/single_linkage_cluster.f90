!------------------------------------------- Program Descriptions ----------------------------------------

! This program contains 2 different implementations of single-linkage cluster algorithm:  
    ! Method 1: Based on existing neighbor list 
	! Method 2: Based on https://github.com/Allen-Tildesley/examples 
	! Method 3:  



! Note: 
! 1. For method 1: large memory may be required. For example, if n_atoms ~ 30^4,
! then at least 4 GB memory is needed. "n_bonds" and "n_dij_bonds" are used as
! critertia to determine if a selected particle be considered in a cluster
! calculations.
! 
! 2. For method 2: coordinates,  box dimension, and cutoff are required to rebuild the
! neighbor list. it is less efficient but more general. Routines
! "gen_cluster_ID_list" generates an array "cluster_ID_list" with size of all
! particles considered, in which each particle is assigned a integer ID number as the cluster ID it belongs to.
! 3. 


module cluster
    use system, only: c_int, c_double, dp, sp
    implicit none 

    private

    public :: single_linkage_cluster_by_Jeremy, &  
            & gen_cluster_nb_list, & 
            & gen_cluster_ID_list, & 
            & get_largest_cluster

contains

	subroutine single_linkage_cluster_by_Jeremy(n_atoms,n_nneigh, nneigh_lst, n_bonds,n_dij_bonds, largest_cluster)
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
		
		if (n_dij_bonds(iatom) == n_bonds) CYCLE

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
!                               Single Linkage cluster routines
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
!                               Method 1:  
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

    pure subroutine gen_cluster_nb_list(total_atoms, &
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

    subroutine gen_cluster_ID_list(cnb_list,total_atoms,cluster_ID_list,cluster_member_head)
        implicit none
        ! Passed
        integer,intent(in) :: total_atoms
        integer,intent(in),dimension(:)  :: cnb_list

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

    subroutine get_largest_cluster(cluster_ID_ary, largest_cluster) 
        implicit none 
        ! Passed
        integer, intent(in), dimension(:) :: cluster_ID_ary

        ! Local 
        integer :: i, num_clusters 
        integer, dimension(:), allocatable :: largest_cluster_lst

        ! Passed
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

end module 
