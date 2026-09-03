%probe_refln_coeff

%     ******************************************************************
function[m_probe,m_probe_end,m_probe_total,r_probe,max_probe_bin]= ...
  probe_refln_coeff(n_bins,fid4,n_bins_total,n_films_total,...
     m_probe,m_probe_end,m_probe_total)
%     calculates the amplitude reflection coefficient for the probe ****
%     calculation includes a trap that avoids the transfer matrix becoming 
%       very large for thick films of highly-absorbing materials.
%     if max_probe_bin is less than n_bins_total then the trap has been applied.
%     if the trap has been applied m_probe_total is the product of 
%       transfer matrices for the bins 1 through max_probe_bin.  
%     if the trap has not been applied m_probe_total is the product of 
%       transfer matrices for all bins times the end transfer matrix. 
%     ******************************************************************

      m_temp_1=zeros(2,2);m_temp_2=zeros(2,2);
      
      max_probe_bin=n_bins_total;
      
      m_temp_1=[1 0;0 1];
      n=0;
      for n_film=1:n_films_total
           if (abs(m_temp_1(1,1)-m_temp_1(2,1)*m_temp_1(1,2)/m_temp_1(2,2))<= 1.0e-3)
               break
           end
        for n_bin=1:n_bins(n_film)
          n=n+1;
          for m=1:2
            for j=1:2
              m_temp_2(m,j)=m_probe(n_film,m,1)*m_temp_1(1,j)+...
            m_probe(n_film,m,2)*m_temp_1(2,j);
            end 
          end 
%     set matrix m_temp_1 equal to m_temp_2 ****************************
          m_temp_1=m_temp_2;
%     trap for very small transmission *********************************
%     prevents overflows for thick opaque films ************************
%         write(*,105) m_temp_1(1,1),m_temp_1(1,2),
%    &    m_temp_1(2,1),m_temp_1(2,2)
%105       format(8e10.4)
          if (abs(m_temp_1(1,1)-m_temp_1(2,1)*m_temp_1(1,2)/m_temp_1(2,2))<= 1.0e-3) 
              max_probe_bin=n;
              m_probe_total=m_temp_1;
              break
          end
%     if max_probe_bin is less than n_bins_total stop ******************
%     if max_probe_bin equals n_bins_total multiply by m_probe_end matrix
              if (max_probe_bin >= n_bins_total)
                  m_probe_total=m_probe_end*m_temp_1;
              end
           
        end 
      end 
      



%      write(*,101) ((m_probe_total(i,j),j=1,2),i=1,2)
%c101   format(' TOTAL PROBE TRANSFER MATRIX',/,
%     & 2x,'(',f12.5,',',f12.5,')',2x,'(',f12.5,',',f12.5,')',/,
%     & 2x,'(',f12.5,',',f12.5,')',2x,'(',f12.5,',',f12.5,')',/,)

      r_probe=-m_probe_total(2,1)/m_probe_total(2,2);
      
      sprintf('  PROBE AMPLITUDE REFLECTION COEFFICIENT IS %e  ', r_probe)
      sprintf(' #PROBE REFLECTION MAGNITUDE IS  %e  ', abs(r_probe)^2)
      fprintf(fid4,' # PROBE REFLECTION MAGNITUDE IS  %e  \r',abs(r_probe)^2)
    

     
      end 