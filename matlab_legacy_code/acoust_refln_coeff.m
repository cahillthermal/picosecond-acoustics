%  acoust_refln_coeff

%     ******************************************************************
function[refl_left,refl_right,tran_left,tran_right]=acoust_refln_coeff(ac_imp,...
    vel,refl_left,refl_right,tran_left,tran_right,n_bins_total,fid3)
%     ******************************************************************
%     refl_left(n) is reflection of strain going to the left at
%     interface between n and n-1
%     refl_right(n) is reflection of strain going to the right at
%     interface between n and n+1      
%     tran_right(n) is transmission of strain going to the right at
%     interface between n and n+1      
%     tran_left(n) is transmission of strain going to the left at
%     interface between n and n-1

      for n=2:n_bins_total+1
        if (n == 2)  
          refl_left(n)=-1.0;
          tran_left(n)=0.0;
        else 
          refl_left(n)=(ac_imp(n-1)-ac_imp(n))/(ac_imp(n-1)+ac_imp(n));
          tran_left(n)=(2*ac_imp(n)/(ac_imp(n-1)+ac_imp(n)))...
          *(vel(n)/vel(n-1));
        end

        refl_right(n)=(ac_imp(n+1)-ac_imp(n))/(ac_imp(n+1)+ac_imp(n));      
        tran_right(n)=(2*ac_imp(n)/(ac_imp(n+1)+ac_imp(n)))...
        *(vel(n)/vel(n+1));

%       the acoustic reflection coefficients are not written to a file *
%       to avoid the creation of a large file **************************
        
      fprintf(fid3,'%d  %e  %e  %e  %e \r',n,refl_left(n),refl_right(n),...
          tran_left(n),tran_right(n));
      end 
%      write(*,*) 'REFLECTION/TRANSMISSION COEFFICIENTS AT BOUNDARIES'  
% <<<NOT SURE WHAT HUMPHREY WAS TRYING TO DO HERE...POSSIBLY MAKE ABSORBING
% BOUNDARY CONDITIONS>>>
%      for n=1:n_bins_total
%        if ((abs(refl_left(n))> 0.0001) | (abs(refl_right(n))>0.0001))        
%        end if
%      end do 
 %     return
      end