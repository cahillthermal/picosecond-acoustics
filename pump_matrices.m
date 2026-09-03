% pump_matrices
%     ******************************************************************
function[m_pump,m_pump_end]= pump_matrices(k_pump,theta_pump,i_polzn_pump,...
     epsilon_pump,bin_size,n_films_total)
%     t_in is the transmission coefficient into the bin 
%     t_out is the transmission coefficient out of the bin 
%     r_in is the reflection coefficient for light coming
%     from the inside of the bin 
%     r_in is the reflection coefficient for light coming
%     from the outside of the bin 
%     *****************************************************************
    

      c=cos(theta_pump);
      s=sin(theta_pump);

      for n_film=1:n_films_total
        c_1=sqrt(epsilon_pump(n_film)-s^2);
        c_2=sqrt(epsilon_pump(n_film));
        c_3=epsilon_pump(n_film)*c;
        k_z_pump=k_pump*c_1;

        if (i_polzn_pump == 1)  
          t_in=2.0*c/(c+c_1);
          t_out=2.0*c_1/(c+c_1);
          r_in=(c-c_1)/(c+c_1);
          r_out=-r_in;
        end

        if (i_polzn_pump == 2)
          t_in=2.0*c_2*c/(epsilon_pump(n_film)*c+c_1);
          t_out=2.0*c_2*c_1/(epsilon_pump(n_film)*c+c_1);
          r_in=(c_3-c_1)/(c_3+c_1);
          r_out=-r_in;
        end
        
        c_4=exp(i*k_z_pump*bin_size(n_film));
        c_5=c_4*c_4;

        t=t_in*t_out*c_4/(1-r_out^2*c_5);
        r=r_in*(1-c_5)/(1-r_out^2*c_5);
             
        m_pump(n_film,1,1)=t-r^2/t;
        m_pump(n_film,1,2)=r/t;
        m_pump(n_film,2,1)=-r/t;
        m_pump(n_film,2,2)=1/t;

        if (n_film == n_films_total)
          m_pump_end(1,1)=t_in-r_out*r_in/t_out;
          m_pump_end(1,2)=r_out/t_out;
          m_pump_end(2,1)=-r_in/t_out;
          m_pump_end(2,2)=1.0/t_out;
        end 

      end
 
%      write(*,*)
%      write(*,*) 'TRANSFER MATRICES PER BIN FOR PUMP LIGHT'
%      do n_film=1,n_films_total
%        write(*,100) ((m_pump(n_film,i,j),j=1,2),i=1,2)
%100     format(2x,4('(',f7.4,',',f7.4,')'))
%      end do
%      write(*,100) ((m_pump_end(i,j),j=1,2),i=1,2)
%      write(*,*)       
 
      end