
%   pump_refln_coeff
%     ******************************************************************
function[m_pump_total,r_pump,max_pump_bin]=pump_refln_coeff(n_bins,...
    n_bins_total,n_films_total,m_pump,m_pump_end,fid4)
%     calculation includes a trap that avoids the transfer matrix becoming 
%       very large for thick films of highly-absorbing materials.
%     if max_pump_bin is less than n_bins_total then the trap has been applied.
%     if the trap has been applied m_pump_total is the product of 
%       transfer matrices for the bins 1 through max_pump_bin.  
%     if the trap has not been applied m_pump_total is the product of 
%       transfer matrices for all bins times the end transfer matrix. 
%     ******************************************************************
      
      m_temp_2=zeros(2,2);m_pump_total=zeros(2,2);
                               
      max_pump_bin=n_bins_total;
      
      m_temp_1=[1 0;0 1];

      n=0;
      
      for n_film=1:n_films_total
          
           if (abs(m_temp_1(1,1)-m_temp_1(2,1)*m_temp_1(1,2)/m_temp_1(2,2))<= 1.0e-3)
         break
         end
        for n_bin=1:n_bins(n_film)
            

%                crappy=m_temp_1
            crap=abs(m_temp_1(1,1)-m_temp_1(2,1)*m_temp_1(1,2)/m_temp_1(2,2));
          n=n+1;
          for m=1:2
            for j=1:2
              m_temp_2(m,j)=m_pump(n_film,m,1)*m_temp_1(1,j)+...
             m_pump(n_film,m,2)*m_temp_1(2,j);
            end 
          end
          m_temp_2(1,1);
%          m_temp_2(1,2)
%          m_temp_2(2,1)
%     set matrix m_temp_1 equal to m_temp_2 ****************************
          m_temp_1=m_temp_2;
%     trap for very small transmission *********************************
%     prevents overflows for thick opaque films ************************
%     write(*,105) m_temp_1(1,1),m_temp_1(1,2),
%    & m_temp_1(2,1),m_temp_1(2,2)

     if (abs(m_temp_1(1,1)-m_temp_1(2,1)*m_temp_1(1,2)/m_temp_1(2,2))<= 1.0e-3)
         m_pump_total=m_temp_1;
         break
     end
%              max_pump_bin=n;
%              m_pump_total=m_temp_1;
%     if max_pump_bin is less than n_bins_total stop ******************
%     if max_pump_bin equals n_bins_total multiply by m_pump_end matrix 
              if (n >= n_bins_total)
                  m_pump_total=m_pump_end*m_temp_1;
                  sprintf('whatyousay???')
              end   
            end 
        
      end
      max_pump_bin=n;
      
      m_temp_2;
%      write(*,101) ((m_pump_total(i,j),j=1,2),i=1,2)
%101   format(' TOTAL PUMP TRANSFER MATRIX',/,
%     & 2x,'(',f12.5,',',f12.5,')',2x,'(',f12.5,',',f12.5,')',/,
%     & 2x,'(',f12.5,',',f12.5,')',2x,'(',f12.5,',',f12.5,')',/,)

       r_pump=-m_pump_total(2,1)/m_pump_total(2,2);
      sprintf('  PUMP AMPLITUDE REFLECTION COEFFICIENT IS')
      r_pump
      sprintf(' # PUMP REFLECTION MAGNITUDE IS  %e  ', abs(r_pump)^2)
      fprintf(fid4,'\r # PUMP REFLECTION MAGNITUDE IS  %e  \r',abs(r_pump)^2)


      end 