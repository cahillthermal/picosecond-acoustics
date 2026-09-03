
%probe_matrices

%     ******************************************************************
function[m_probe,m_probe_end,d_m_probe]= probe_matrices(k_probe,theta_probe,...
    i_polzn_probe,epsilon_probe,d_epsilon_d_strain,bin_size,n_films_total,...
    d_m_probe,m_probe)
 
%     calculates probe matrices and change in probe matrices with strain
%     k_z_probe= z component of wave vector in bin
%     k_z_probe_prime= z component of wave vector in bin after perturbation
%     t_in= the transmission coefficient into the bin 
%     t_out= the transmission coefficient out of the bin 
%     r_in= the reflection coefficient for light coming
%       from the inside of the bin 
%     r_in= the reflection coefficient for light coming
%       from the outside of the bin 
%     ******************************************************************
      epsilon_probe_prime=zeros(1,40);
      bin_size_prime=zeros(1,40) ;
      
      c=cos(theta_probe);
      s=sin(theta_probe);

      for n_film=1:n_films_total
        c_1=sqrt(epsilon_probe(n_film)-s^2);
        c_2=sqrt(epsilon_probe(n_film));
        c_3=epsilon_probe(n_film)*c;
        k_z_probe=k_probe*c_1;
        if (i_polzn_probe == 1) 
          t_in=2.0*c/(c+c_1);
          t_out=2.0*c_1/(c+c_1);
          r_in=(c-c_1)/(c+c_1);
          r_out=-r_in;
        end 
        if (i_polzn_probe == 2) 
          t_in=2.0*c_2*c/(epsilon_probe(n_film)*c+c_1);
          t_out=2.0*c_2*c_1/(epsilon_probe(n_film)*c+c_1);
          r_in=(c_3-c_1)/(c_3+c_1);
          r_out=-r_in;
        end
        
        c_4=exp(i*k_z_probe*bin_size(n_film));
        c_5=c_4*c_4;

        t=t_in*t_out*c_4/(1.0-r_out^2*c_5);
        r=r_in*(1-c_5)/(1.0-r_out^2*c_5);
             
        m_probe(n_film,1,1)=t-r^2/t;
        m_probe(n_film,1,2)=r/t;
        m_probe(n_film,2,1)=-r/t;
        m_probe(n_film,2,2)=1/t;

        if (n_film == n_films_total) 
          m_probe_end(1,1)=t_in-r_out*r_in/t_out;
          m_probe_end(1,2)=r_out/t_out;
          m_probe_end(2,1)=-r_in/t_out;
          m_probe_end(2,2)=1.0/t_out;
        end 

      end 
 
%      write(*,*)
%      write(*,*) 'TRANSFER MATRICES PER BIN FOR PROBE LIGHT'
%      do n_film=1,n_films_total
%        write(*,100) ((m_probe(n_film,i,j),j=1,2),i=1,2)
% 100     format(2x,4('(',f7.4,',',f7.4,')'))
%      end do
%      write(*,100) ((m_probe_end(i,j),j=1,2),i=1,2)
%      write(*,*)      

%     now recalculate for strained bin  ********************************
      strain=0.01;  % strain in bin, chosen to be small
      for n_film=1:n_films_total
        epsilon_probe_prime(n_film)=epsilon_probe(n_film)+ ...  
        d_epsilon_d_strain(n_film)*strain;        % perturbed dielectric constant
        bin_size_prime(n_film)=bin_size(n_film)*(1.0+strain);  % perturbed bin size 
      end 

      for n_film=1:n_films_total
        c_1=sqrt(epsilon_probe_prime(n_film)-s^2);
        c_2=sqrt(epsilon_probe_prime(n_film));
        c_3=epsilon_probe_prime(n_film)*c;
        k_z_probe_prime=k_probe*c_1;
        if (i_polzn_probe == 1)  
          t_in=2.0*c/(c+c_1);
          t_out=2.0*c_1/(c+c_1);
          r_in=(c-c_1)/(c+c_1);
          r_out=-r_in;
        end 
        if (i_polzn_probe == 2)  
          t_in=2.0*c_2*c/(epsilon_probe_prime(n_film)*c+c_1);
          t_out=2.0*c_2*c_1/(epsilon_probe_prime(n_film)*c+c_1);
          r_in=(c_3-c_1)/(c_3+c_1);
          r_out=-r_in;
        end 
        
        c_4=exp(i*k_z_probe_prime*bin_size_prime(n_film));
        c_5=c_4*c_4;

        t=t_in*t_out*c_4/(1.0-r_out^2*c_5);
        r=r_in*(1-c_5)/(1.0-r_out^2*c_5);
             
        d_m_probe(n_film,1,1)=(t-r^2/t-m_probe(n_film,1,1))/strain;
        d_m_probe(n_film,1,2)=(r/t-m_probe(n_film,1,2))/strain;
        d_m_probe(n_film,2,1)=(-r/t-m_probe(n_film,2,1))/strain;
        d_m_probe(n_film,2,2)=(1.0/t-m_probe(n_film,2,2))/strain;

      end 
 
%      write(*,*)
%      write(*,*) 'PERTBN OF TRANSFER MATRICES PER BIN FOR PROBE LIGHT'
%      do n_film=1,n_films_total
%        write(*,100) ((d_m_probe(n_film,i,j),j=1,2),i=1,2)
%      end do

     
      end