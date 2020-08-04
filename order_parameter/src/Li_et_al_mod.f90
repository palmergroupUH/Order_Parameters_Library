module Li_et_al
    use system, only: c_int, c_double, dp, sp
    use order_parameter, only: build_homo_neighbor_list,&
                             & check_nb_list, & 
                             & apply_nearest_neighbor_crit, &  
                             & initialize_Ylm, & 
                             & compute_ave_Ylm

    use cluster, only: single_linkage_cluster_by_Jeremy, & 
                     & gen_cluster_nb_list, &
                     & gen_cluster_ID_list, & 
                     & get_largest_cluster 

    implicit none 
    private 


contains 

    subroutine initialize_Li_et_al(l, sph_const, Plm_const) bind(c, name="initialize_Li_et_al")
        implicit none
        ! Passed
        integer(c_int), intent(in) :: l
        
        ! Return
        real(c_double), intent(out), dimension(-l:l) :: sph_const
        real(c_double), intent(out), dimension(0:l,4) :: Plm_const

        call initialize_Ylm(l, sph_const, Plm_const)

        end subroutine 

    subroutine call_Li_et_al(sph_const, & 
                           & Plm_const, &
                           & total_atoms, &
                           & maxnb, & 
                           & nnb, &
                           & l, &
                           & cutoff_sqr, &
                           & box, &
                           & xyz, &
                           & Iij) bind(c, name="call_Li_et_al")
                         
        implicit none 

        ! Passed:
        real(c_double), intent(in), dimension(-l:l) :: sph_const 
        real(c_double), intent(in), dimension(0:l,4) :: Plm_const 
        integer(c_int), intent(in) :: total_atoms
        integer(c_int), intent(in) :: maxnb
        integer(c_int), intent(in) :: nnb
        integer(c_int), intent(in) :: l
        real(c_double), intent(in) :: cutoff_sqr
        real(c_double), intent(in), dimension(1:3) :: box
        real(c_double), intent(in), dimension(1:3, total_atoms) :: xyz

        ! Local:
        integer, dimension(1:total_atoms) :: num_NB_list 
        integer, dimension(1:maxnb, 1:total_atoms) :: NB_list 
        real(dp), dimension(1:4, 1:maxnb, 1:total_atoms) :: Rij
        complex(dp), dimension(-l:l, total_atoms) :: Ylm
        real(dp), dimension(1:2) :: bins_range
        real(dp) :: cluster_cut, start, finish  
        real(dp) :: interval, upper, lower
        integer :: largest_cluster, k, counter 
        integer, dimension(1:20) :: bin_index
        integer, dimension(:),allocatable :: ID_list, clinklist, cluster_member_head, largest_cluster_lst
        integer, dimension(:),allocatable, target :: cluster_ID  
        integer, dimension(:),pointer :: cluster_ID_ptr  

        ! Return:
        real(c_double), intent(out), dimension(1:total_atoms) :: Iij 

        ! get the neighbor list
        
        
        call build_homo_neighbor_list(total_atoms, maxnb, cutoff_sqr, xyz, box, num_NB_list, NB_list, Rij) 
        
        ! check the neighbor list: 
        call check_nb_list(num_NB_list,maxnb,nnb) 
        
        ! apply the nearest neighbor:
        call apply_nearest_neighbor_crit(nnb, total_atoms, num_NB_list, NB_list, Rij)
        
        ! compute_sum_Ylm  
        call compute_ave_Ylm(sph_const, Plm_const, total_atoms, l, num_NB_list, Rij, Ylm)
       
        call calc_Iql(total_atoms,l,num_NB_list,NB_list,Ylm,Iij)  
        
        !call cpu_time(start)
        !cluster_cut = 3.5d0
        !call rebuild_neighbor_list_nnb(total_atoms, &
                                     !& cluster_cut, &
                                     !& maxnb, &
                                     !& num_NB_list, &
                                     !& NB_list, &
                                     !& Rij)

        !call cpu_time(finish) 
    
        !call single_linkage_cluster_by_Jeremy(total_atoms, num_NB_list, NB_list, 4, chill_id_list, largest_cluster)
    
        !print*, largest_cluster 

        !print*, largest_cluster  

        !counter = 0

        !allocate(ID_list(1:count(chill_id_list/=4)))

        !do k = 1,total_atoms

            !if ( chill_id_list(k) /= 4 ) then

                !counter = counter + 1

                !ID_list(counter) = k

            !end if

        !end do
        
        !cluster_cut = 3.5d0*3.5d0

        !call gen_cluster_nb_list(size(ID_list),xyz(:,ID_list),box,cluster_cut,clinklist)
       
        !call gen_cluster_ID_list(clinklist,size(ID_list),cluster_ID,cluster_member_head)

        !call get_largest_cluster(cluster_ID, largest_cluster) 

        !print*, "Allen's largest cluster", largest_cluster  

        end subroutine

    subroutine calc_Iql(total_atoms,l,num_NB_list,NB_list,qlm,Iij)
        implicit none
        ! Passed
        integer, intent(in) :: total_atoms,l
        integer, intent(in), dimension(:) :: num_NB_list
        integer, intent(in), dimension(:,:) :: NB_list
        complex(dp), intent(in), dimension(:,:) :: qlm
    
        ! Local
        integer :: i,j,nb,num_comp,num_NB, m
        real(dp) :: sum_nb,normalized_qlm_i,normalized_qlm_j, sum_ylm, start, finish
        complex(dp),dimension(-l:l) :: Ylm_i, Ylm_j

        ! Return
        real(dp), intent(out), dimension(1:total_atoms) :: Iij

        Iij = 0.0d0
        
        do i = 1, total_atoms

            Ylm_i = qlm(:,i)

            num_NB = num_NB_list(i)

            normalized_qlm_i = ZSQRT(SUM(Ylm_i*DCONJG(Ylm_i)))

            sum_nb = 0

            do j = 1, num_NB

                nb = NB_list(j,i)

                Ylm_j = qlm(:,nb)

                normalized_qlm_j = ZSQRT(SUM(Ylm_j*DCONJG(Ylm_j)))

                sum_ylm = 0.0d0 

                do m = -l, l 

                    sum_ylm = sum_ylm + Ylm_i(m)* dconjg(Ylm_j(m))
       
                end do 
                
                sum_nb = sum_nb + sum_ylm/(normalized_qlm_i*normalized_qlm_j)
                
            end do
            
            Iij(i) = sum_nb*1.0d0/num_NB

        end do

        end subroutine

end module

