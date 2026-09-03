% sensitivity_function
%     *******************************************************************
function[f_sens]= sensitivity_function(n_bins,n_bins_total,n_films_total,...
     z,m_probe,d_m_probe,m_probe_end,m_probe_total,r_probe,f_sens,...
     max_probe_bin,fid5)
%     calculates the sensitivity function *******************************
%     *******************************************************************
      
    m_a=zeros(2003,2);
    m_a(:,:,2)=zeros(2003,2);
    m_b=zeros(2003,2);
    m_b(:,:,2)=zeros(2003,2);
 
    d_m_probe_total=zeros(2,2);m_work=zeros(2,2);
      
    

%     calculation of the m_a matrices **********************************
%     set m_a(n_bins_total+2) equal to the unit matrix *****************
      m_a(n_bins_total+2,:,:)=[1 0;0 1];
%     continue calculating for decreasing n ****************************
%     if trap has been applied multiply by the unit matrix until *******
%       we reach n = max_probe_bin *************************************
%     otherwise first multiply by probe matrix for the end, and then a *
%       sequence of matrices for each bin ******************************
      if (max_probe_bin == n_bins_total) 
        for m=1:2
          for j=1:2
            m_a(n_bins_total+1,m,j)=m_probe_end(m,j);
          end 
        end       
      else 
        for i=1:2
          for j=1:2
            m_a(n_bins_total+1,i,j)=m_a(n_bins_total+2,i,j);
          end 
        end       
      end 

      n=n_bins_total+1;

      for n_film=n_films_total:-1:1     
        for n_bin=1:n_bins(n_film)
          n=n-1;
          if (n <= max_probe_bin)
            for m=1:2
              for j=1:2
                m_a(n,m,j)=m_a(n+1,m,1)*m_probe(n_film,1,j)+ ...
                m_a(n+1,m,2)*m_probe(n_film,2,j);
              end 
            end        
          else
            for m=1:2
              for j=1:2
                m_a(n,m,j)=m_a(n+1,m,j);
              end 
            end                   
          end 
        end 
      end 

%     calculation of m_b matrices **************************************
      m_b(1,:,:)=[1 0;0 1];

%     if trap has been applied multiply by the unit matrix after *******
%       we reach n = max_probe_bin *************************************
%     otherwise multiply by a sequence of matrices for each bin ********
      n=1 ;  
      for n_film=1:n_films_total      
        for n_bin=1:n_bins(n_film)
          n=n+1;
          if (n <= max_probe_bin)  
            for m=1:2
              for j=1:2
                m_b(n,m,j)=m_probe(n_film,m,1)*m_b(n-1,1,j)+...
                m_probe(n_film,m,2)*m_b(n-1,2,j);
              end 
            end        
          else 
            for m=1:2
              for j=1:2
                m_b(n,m,j)=m_b(n-1,1,j);
              end 
            end 
          end           
        end 
      end

%     calculation of the sensitivity function **************************
%     find change in the transfer matrix *******************************
      n=1;
      for n_film=1:n_films_total      
        for n_bin=1:n_bins(n_film)
          n=n+1;
          if (n <= max_probe_bin)
            for m=1:2
              for j=1:2
                m_work(m,j)=d_m_probe(n_film,m,1)*m_b(n-1,1,j)...
                 +d_m_probe(n_film,m,2)*m_b(n-1,2,j);
              end  
            end 
            for m=1:2
              for j=1:2
                d_m_probe_total(m,j)=m_a(n+1,m,1)*m_work(1,j)+ ... 
                 m_a(n+1,m,2)*m_work(2,j);
              end 
            end 
%     find change in reflectivity and sensitivity **********************
            d_r_probe=(m_probe_total(2,1)*d_m_probe_total(2,2)...
            -m_probe_total(2,2)*d_m_probe_total(2,1))...
             /m_probe_total(2,2)^2;
            f_sens(n)=2.0*real(r_probe*conj(d_r_probe));
          else
            f_sens(n)=0.0;
          end
        end
      end                
 
%      the sensitivity function is not saved so as to avoid the creation 
%      of a large file ************************************************* 
      fprintf(fid5, '%e  %e \r',0.0,f_sens(1));
      for n=1:n_bins_total
        fprintf(fid5, '%e  %e \r',z(n),f_sens(n));
        fprintf(fid5, '%e  %e \r',z(n),f_sens(n+1));
      end                           
                                 
      
      end