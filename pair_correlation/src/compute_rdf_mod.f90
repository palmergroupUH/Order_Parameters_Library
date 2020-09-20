module RDF
    use system
    implicit none 
    private 

    public :: build_homo_pair_dist_hist, & 
            & build_hetero_pair_dist_hist, &   
            & build_homo_pair_dist_hist_in_chunk,& 
            & normalize_hist 

contains

    subroutine build_homo_pair_dist_hist(total_atoms,&
                                     & cutoff,&
                                     & num_bins,&
                                     & XYZ,&
                                     & L_xyz,&
                                     & rdf_histogram) &
                                     & bind(c,name="call_homo_pair_dist_hist") 
        implicit none
        
        ! Passed 
        integer(c_int),intent(in) :: total_atoms,num_bins 
        real(c_double),intent(in) :: cutoff
        real(c_double),intent(in),dimension(1:3,1:total_atoms) :: XYZ
        real(c_double),intent(in),dimension(1:3) :: L_xyz

        ! Local 
        real(dp),dimension(1:size(L_xyz)) :: xyz_temp,xyz_separate
        integer :: i,j,bin_index
        real(dp) :: distance_sqr,distance,cutoff_sqr,r_interval 

        ! Output 
        real(c_double),intent(out),dimension(1:num_bins) :: rdf_histogram

        r_interval = cutoff/num_bins 

        cutoff_sqr = cutoff*cutoff

        rdf_histogram = 0.0d0 

        do i = 1, total_atoms-1

            xyz_temp = xyz(1:3,i)

            do j = i + 1,total_atoms

                ! compute the separation vector between i,j  

                xyz_separate = xyz_temp - xyz(1:3,j)

                ! apply minimum image convention 

                xyz_separate = xyz_separate - L_xyz*dnint(xyz_separate/L_xyz)

                ! sum of squared separation vector  

                distance_sqr = xyz_separate(1)*xyz_separate(1) & 
                           & + xyz_separate(2)*xyz_separate(2) & 
                           & + xyz_separate(3)*xyz_separate(3)

                if (distance_sqr < cutoff_sqr ) then 
                    
                    bin_index = int(dsqrt(distance_sqr)/r_interval) + 1
                    
                    rdf_histogram(bin_index)  =  rdf_histogram(bin_index) + 2

                end if

            end do 

        end do

        end subroutine 

    subroutine build_hetero_pair_dist_hist(total_atoms_A, &
                                         & total_atoms_B, &  
                                         & cutoff,&
                                         & num_bins,&
                                         & xyz_A,&
                                         & xyz_B, & 
                                         & L_xyz,&
                                         & rdf_histogram) &
                                         & bind(c,name="hetero_pair_distance_histogram") 
        implicit none
        
        ! Passed 
        integer(c_int),intent(in) :: total_atoms_A, total_atoms_B, num_bins 
        real(c_double),intent(in) :: cutoff
        real(c_double),intent(in),dimension(1:3,1:total_atoms_A) :: xyz_A
        real(c_double),intent(in),dimension(1:3,1:total_atoms_B) :: xyz_B
        real(c_double),intent(in),dimension(1:3) :: L_xyz

        ! Local 
        real(dp),dimension(1:size(L_xyz)) :: xyz_temp,xyz_separate
        integer :: i,j,bin_index
        real(dp) :: distance_sqr,distance,cutoff_sqr,r_interval 

        ! Output 
        real(c_double),intent(out),dimension(1:num_bins) :: rdf_histogram

        r_interval = cutoff/num_bins 

        cutoff_sqr = cutoff*cutoff

        rdf_histogram = 0.0d0 

        do i = 1, total_atoms_A

            xyz_temp = xyz_A(1:3,i)

            do j = 1, total_atoms_B

                ! compute the separation vector between i,j  

                xyz_separate = xyz_temp - xyz_B(1:3,j)

                ! apply minimum image convention 

                xyz_separate = xyz_separate - L_xyz*dnint(xyz_separate/L_xyz)

                ! sum of squared separation vector  

                distance_sqr = xyz_separate(1)*xyz_separate(1) & 
                           & + xyz_separate(2)*xyz_separate(2) & 
                           & + xyz_separate(3)*xyz_separate(3)

                if (distance_sqr < cutoff_sqr ) then 
                    
                    bin_index = int(dsqrt(distance_sqr)/r_interval) + 1
                    
                    rdf_histogram(bin_index)  =  rdf_histogram(bin_index) + 1 

                end if

            end do 

        end do

        end subroutine 

    subroutine build_homo_pair_dist_hist_in_chunk(total_atoms,&
                                                 & num_configs,& 
                                                 & cutoff,&
                                                 & num_bins,&
                                                 & XYZ,&
                                                 & L_xyz,&
                                                 & rdf_histogram, &
                                                 & sum_volume)&
                                                 & bind(c,name="call_homo_pair_dist_hist_in_chunk") 
        implicit none

        ! Passed 
        integer(c_int),intent(in) :: total_atoms,num_bins,num_configs 
        real(c_double),intent(in) :: cutoff
        real(c_double),intent(in),dimension(1:3,1:total_atoms,1:num_configs) :: XYZ
        real(c_double),intent(in),dimension(1:3,1:num_configs) :: L_xyz

        ! Local 
        real(dp),dimension(1:size(L_xyz)) :: xyz_temp,xyz_separate
        integer :: i,j,bin_index,iconfig
        real(dp) :: distance_sqr,distance,cutoff_sqr,r_interval 
        real(dp),dimension(1:3) :: box  

        ! Output 
        real(c_double),intent(out),dimension(1:num_bins) :: rdf_histogram
        real(c_double),intent(out) :: sum_volume
        
        r_interval = cutoff/num_bins 

        cutoff_sqr = cutoff*cutoff

        rdf_histogram = 0.0d0 

        sum_volume = 0.0

        do iconfig = 1, num_configs 
    
            box = L_xyz(:,iconfig) 

            sum_volume = sum_volume + product(box)

            do i = 1, total_atoms-1

                xyz_temp = xyz(1:3,i,iconfig)

                do j = i + 1,total_atoms

                    ! compute the separation vector between i,j  

                    xyz_separate = xyz_temp - xyz(1:3,j,iconfig)

                    ! apply minimum image convention 

                    xyz_separate = xyz_separate - box*dnint(xyz_separate/box)

                    ! sum of squared separation vector  

                    distance_sqr = xyz_separate(1)*xyz_separate(1) & 
                               & + xyz_separate(2)*xyz_separate(2) & 
                               & + xyz_separate(3)*xyz_separate(3)

                    if (distance_sqr < cutoff_sqr ) then 
                        
                        bin_index = int(dsqrt(distance_sqr)/r_interval) + 1
                        
                        rdf_histogram(bin_index)  =  rdf_histogram(bin_index) + 2

                    end if

                end do 

            end do

        end do 
        
        end subroutine 

    subroutine normalize_hist(rdf_histogram, &
                             & num_bins, &
                             & cutoff, & 
                             & natoms, & 
                             & num_configs, & 
                             & bulk_density, & 
                             & gr, & 
                             & r2hr) bind(c,name="call_normalize_hist")
        implicit none  

        ! Passed 
        integer(c_int),intent(in) :: num_bins,natoms,num_configs
        real(c_double),intent(in),dimension(1:num_bins) :: rdf_histogram 
        real(c_double),intent(in) :: cutoff,bulk_density  

        ! Local
        real(dp) :: preconst,vshell_i,half_interval,upper_r,lower_r,center_r,r_interval  
        integer :: i 

        ! Output
        real(c_double),dimension(1:num_bins),intent(out) :: gr  
        real(c_double),dimension(1:num_bins),intent(out) :: r2hr   

        preconst = 4.0d0/3.0d0*3.14159265358979d0

        r_interval = cutoff/num_bins

        half_interval = r_interval/2.0d0

        do i = 1, num_bins 

            upper_r = i*r_interval 

            lower_r = (i-1)*r_interval

            center_r = half_interval + r_interval*(i-1)

            vshell_i = preconst*( ( upper_r*upper_r*upper_r ) - (lower_r*lower_r*lower_r) )  

            gr(i) = rdf_histogram(i) / (vshell_i*bulk_density*natoms*num_configs)

            r2hr(i) = center_r * center_r * (gr(i) - 1)

        end do 				

        end subroutine 

end module 
