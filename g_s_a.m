%  g_s_a  (generaized stepping algorithm)

%     ******************************************************************
function[g_left,g_right]= g_s_a(n_bins,n_bins_total,n_films_total,bin_size,...
    tau,attenuation,velocity,refl_left,refl_right,tran_left,tran_right,...
    g_0,g_1,g_2,g_left,g_right)
%     calculate coefficients for generalized stepping algorithm 
%     ******************************************************************  
      n_bin=0;
      for n_film=1:n_films_total
        x=attenuation(n_film)/tau;
        if (x > 0.25)  
          error('Attenuation in film # %d is too large \r Value of x is ... %e \r It is necessary to increase the value of the time step. \r Program will crash.',n_film,x) 
        end
        y=velocity(n_film)*tau/bin_size(n_film);
        for n=1:n_bins(n_film)
          n_bin=n_bin+1;
          g_0(n_bin)=x+0.5*(y-2)*(y-1);
          g_1(n_bin)=1.0-2.0*x-(y-1.0)^2;
          g_2(n_bin)=x+0.5*y*(y-1);
        end
      end 

%     calculate g_left and g_right coefficients ************************
%     g_left(*,n) is for strain arriving at bin n and going left *******
%     g_left(1,n) .. is from leftgoing strain from bin n *************** 
%     g_left(2,n) .. is from leftgoing strain from bin n+1 *************
%     g_left(3,n) .. is from leftgoing strain from bin n+2 *************
%     g_left(4,n) .. is from rightgoing strain from bin n-1 ************ 
%     g_left(5,n) .. is from rightgoing strain from bin n **************
%     g_left(6,n) .. is from rightgoing strain from bin n+1************* 
%     g_right(*,n) is for strain arriving at bin n and going right *****
%     g_right(1,n) .. is from rightgoing strain from bin n ************* 
%     g_right(2,n) .. is from rightgoing strain from bin n-1 ***********
%     g_right(3,n) .. is from rightgoing strain from bin n-2 ***********
%     g_right(4,n) .. is from leftgoing strain from bin n+1 ************ 
%     g_right(5,n) .. is from leftgoing strain from bin n **************
%     g_right(6,n) .. is from leftgoing strain from bin n-1*************
      for n=2:n_bins_total+1
        g_left(1,n)=g_0(n)+g_2(n)*refl_left(n)*refl_right(n);
        if (n+1 <= n_bins_total) 
          g_left(2,n)=g_1(n+1)*tran_left(n+1);
        else 
          g_left(2,n)=0.0;
        end
        if (n+2 <= n_bins_total)
          g_left(3,n)=g_2(n+2)*tran_left(n+2)*tran_left(n+1);
        else 
          g_left(3,n)=0.0;
        end 
        if (n-1 >= 1) 
          g_left(4,n)=g_2(n-1)*tran_right(n-1)*refl_right(n);
        else 
          g_left(4,n)=0.0;
        end  
        g_left(5,n)=g_1(n)*refl_right(n);
        if (n+1 <= n_bins_total) 
          g_left(6,n)=g_2(n+1)*refl_right(n+1)*tran_left(n+1);
        else 
          g_left(6,n)=0.0;
        end 
        g_right(1,n)=g_0(n)+g_2(n)*refl_right(n)*refl_left(n);
        if (n-1 >= 1)
          g_right(2,n)=g_1(n-1)*tran_right(n-1);
        else 
          g_right(2,n)=0.0;
        end 
        if (n-2 >= 1)
          g_right(3,n)=g_2(n-2)*tran_right(n-2)*tran_right(n-1);
        else 
          g_right(3,n)=0.0;
        end 
        if (n+1 <= n_bins_total) 
          g_right(4,n)=g_2(n+1)*tran_left(n+1)*refl_left(n);
        else 
          g_right(4,n)=0.0;
        end 
        g_right(5,n)=g_1(n)*refl_left(n);
        if (n-1 >= 1) 
          g_right(6,n)=g_2(n-1)*refl_left(n-1)*tran_right(n-1);
        else 
          g_right(6,n)=0.0;
        end   
      end 
      
      
      end