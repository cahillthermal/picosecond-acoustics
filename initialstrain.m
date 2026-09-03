%initialstrain
%     ******************************************************************
function[eta_left,eta_right]= initialstrain(bin_size,n_bins,eta_left,eta_right,...
    n_films_total,r_pump,alpha,diffusion,m_pump,n_bins_total,max_pump_bin)
%     calculates the initial strain distribution in the structure ****** 
%     the effect of electron and heat diffusion is included ************      
%     spread() array describes the spreading out of the strain *********
%     that occurs as a result of the diffusion *************************
%     *****************************************************************
 
      spread=zeros(1,2002); 
      eta_left_diffn=zeros(1,2002);
      eta_right_diffn=zeros(1,2002);
       
      n=1; % counter for bins

      a_1_left=1;  % initial amplitudes to left of structure 
      a_2_left=r_pump;     
      
      for n_film=1:n_films_total
        for n_bin=1:n_bins(n_film)
          n=n+1;
%     only find energy deposited in bins within pump penetration region
          if (n <= max_pump_bin) 
%     calculate amplitudes to right of bin *****************************
            a_1_right=m_pump(n_film,1,1)*a_1_left+m_pump(n_film,1,2)*a_2_left;       
            a_2_right=m_pump(n_film,2,1)*a_1_left+m_pump(n_film,2,2)*a_2_left;      

%     find energy deposited in the bin *********************************    
%           write(*,100) cabs(a_1_left)**2,cabs(a_1_right)**2,
%    &      cabs(a_2_right)**2,cabs(a_2_left)**2      

            energy=abs(a_1_left)^2-abs(a_1_right)^2 ...
         +abs(a_2_right)^2-abs(a_2_left)^2; 
%     calculate the strain *********************************************
%     notice that this is relative to the stress free state ************
%     thus heating gives a compression which corresponds to a negative strain
            eta_left(n)=-0.5*alpha(n_film)*energy/bin_size(n_film);
            eta_right(n)=eta_left(n);
            a_1_left=a_1_right;  % use these for the next bin 
            a_2_left=a_2_right;  %
          else
            eta_left(n)=0.0;
            eta_right(n)=0.0;
          end  
        end  
      end  
        eta_left(1)
        eta_right(1)
%     allow for the effect of electron and thermal diffusion ***********
%     in this program only the diffusion within a bin is considered ****
%     heat deposited in a bin is assumed to be spread out with an ******
%       exponential distribution as a result of diffusion **************
      %for n=1:n_bins_total
      %  eta_left_diffn(n)=0.0;  % strain distribution after diffusion 
      %  eta_right_diffn(n)=0.0; %already done by MATLAB (8/22/07)
      %end
      
      x_max=40.0; % largest argument used for the hyperbolic functions
      n=1; %  set counter for bin where light is absorbed
      for n_film=1:n_films_total  % consider each film in turn 
        if (diffusion(n_film) == 0.0)  
          x_0=x_max;
        else 
          x_0=bin_size(n_film)/diffusion(n_film);
        end 
        for n_bin=1:n_bins(n_film)
          n=n+1;  % index of bin for where light is absorbed 
          sum=0.0;
          for n_bin_p=1:n_bins(n_film)
            n_p=n+n_bin_p-n_bin; % index of bin to which energy diffuses 
            if (n_bin_p <= n_bin)   % diffusion to the left ******
              if (n_bin*x_0 < x_max)    
                spread(n_p)=cosh(n_bin_p*x_0)/cosh(n_bin*x_0); % evaluate cosh
              else                  
                x_1=min([x_max (n_bin-n_bin_p)*x_0]);
                spread(n_p)=exp(-x_1);  % approximate cosh by exponentials 
              end 
            else            % diffusion to the right *******************
              if ((n_bins(n_film)-n_bin)*x_0 < x_max)    
                spread(n_p)=cosh((n_bins(n_film)-n_bin_p)*x_0) ... % evaluate cosh
                /cosh((n_bins(n_film)-n_bin)*x_0);
              else
                x_1=min([x_max (n_bin_p-n_bin)*x_0]);
                spread(n_p)=exp(-x_1); % approximate cosh by exponentials 
              end 
            end 
            sum=sum+spread(n_p);
          end 
%     normalize distribution function **********************************
          for n_bin_p=1:n_bins(n_film)
            n_p=n+n_bin_p-n_bin; % index of bin to which energy diffuses             
            spread(n_p)=spread(n_p)/sum;
          end  
%     re-distribute strain ********************************************* 
          for n_bin_p=1:n_bins(n_film)
            n_p=n+n_bin_p-n_bin; % index of bin to which energy diffuses 
            eta_left_diffn(n_p)=eta_left_diffn(n_p)+ eta_left(n)*spread(n_p);
            eta_right_diffn(n_p)=eta_right_diffn(n_p)+eta_right(n)*spread(n_p);
          end 
        end 
      end

%     set strain equal to spread out strain in all bins ****************
      for n=2:n_bins_total+1
        eta_left(n)=eta_left_diffn(n);
        eta_right(n)=eta_right_diffn(n);
      end  
    

    
      end
%     ******************************************************************