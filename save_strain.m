% save_strain

%     ******************************************************************
function[i_save]= save_strain(i_save,z,n_bins_total,eta_left,eta_right,...
    eta_scale,fid2)
%     ******************************************************************

      eta_offset=1.2*i_save*eta_scale;
      i_save=i_save+1;

%     write to cumulative file *****************************************
      for n=2:n_bins_total+1
        fprintf(fid2,'%e   %e \r', z(n-1),eta_right(n)+eta_left(n)-eta_offset);
        fprintf(fid2,'%e   %e \r', z(n),eta_right(n)+eta_left(n)-eta_offset);
      end 

      for  n=2:n_bins_total
        fprintf(fid2,'%e   %e \r',z(n),-eta_offset); 
      end  

      return
      end