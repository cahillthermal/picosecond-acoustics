%     ******************************************************************
%     program 'r_t_6c.m'  first matlab version of r_t:
%       based on program written by H.M. - Brown Univ.
%     ******************************************************************
%     finds the time evolution of a strain set up in a multilayer system
%     program includes the calculation of
%       1) the distribution of energy deposited in the system
%       2) the sensitivity function 
%       3) the change in optical reflectivity as a function of time
%       4) the attenuation of sound 
%       5) the effect of diffusion of heat before the sound is generated
%       6) the effect of roughness
%     program reads in material parameters from a set
%       of input files for each material.
%     program is similar to 'r_t_5fa.for' except that roughness is included 
%     the generalized stepping algorithm is used. 
%     this program has fixed size arrays based on 2000 bins
                     
%     LAST EDIT DATE 9/15/96

%     ac_imp(n)= acoustic impedance of bin # n
%     alpha(n_film)= thermal expansion coefficient divided by specific 
%       heat per unit volume for film # n_film 
%     attenuation(n_film)= amplitude attenuation per unit time in film 
%       # n_film divided by omega**2
%     bin_size(n_film)= thickness of each bin in film # n_film
%     boundary(n_film)= right hand boundary of film # n_film
%     d_epsilon_d_strain(n_film)= strain derivative of dielectric constant
%     diffusion(n_film)= diffusion length for film n_film
%     d_m_probe(n_film,i,j)= change in transfer matrix for a bin in film # 
%       n_film for probe light
%     d_refln= change in intensity reflection coefficient for the probe
%     epsilon_probe(n_film)= dielectric constant for probe for film # n_film
%     epsilon_pump(n_film)= dielectric constant for pump for film # n_film
%     eta_left(n_bin)= strain propagating to the left in bin # n_bin
%     eta_right(n_bin)= strain propagating to the right in bin # n_bin
%     f_sens(n_bin)= sensitivity function in bin # n_bin
%     g_0(n_film), g_1(n_film), g_2(n_film)= generalized stepping 
%       algorithm coefficients for film # n_film
%     g_left(n_bin), g_right(n_bin)= generalized stepping algorithm 
%     coefficients for bin # n_bin 
%     i_polzn_probe= 1 for probe polzn perp to the plane of the interface (sigma)
%     i_polzn_probe= 2 for probe polzn in the plane of the interface (pi)
%     i_polzn_pump= 1 for pump polzn perp to the plane of the interface (sigma)
%     i_polzn_pump= 2 for pump polzn in the plane of the interface (pi)
%     i_t= number of time steps taken  
%     k_probe= 2 pi/lambda_probe
%     k_pump= 2 pi/lambda_pump
%     lambda_probe= wavelength of probe light in free space
%     lambda_pump= wavelength of pump light in free space
%     material$(n_film)= material for film # n_film
%     m_probe(n,i,j)= transfer matrix for bin # n for probe light
%     m_probe_end(i,j)= transfer matrix for end for probe light
%     m_probe_total= transfer matrix for structure for the probe light 
%     m_pump(n,i,j)= transfer matrix for bin # n for pump light
%     m_pump_end(i,j)= transfer matrix for end for pump light
%     m_pump_total= transfer matrix for structure for the pump light 
%     n_bin= # of an individual bin
%     n_bins(n_film)= # of elements in film # n_film
%     n_bins_total= total # of bins in the structure
%     n_film= # of a particular film in structure
%     n_films_total= total # of films in the structure
%     n_tau_divide= time step in simulation is time step in data divided
%       by n_tau_divide
%     refl_left(n)= reflection coefficient of a wave going to the left at
%       interface between film # n and film # n-1
%     refl_right(n)= reflection coefficient of a wave going to the right at
%       interface between film # n and film # n+1      
%     r_probe= the amplitude reflection coefficient for the probe light
%     r_pump= the amplitude reflection coefficient for the pump light
%     roughness(n_film) = rms fluctuation in thickness of film # n_film
%     r_t(n) = averaged change in reflectivity at data time step n
%     n_samples = # of configurations considered for film thicknesses 
%     tau= time step for strain development
%     tau_data= time step in data acquisition
%     theta_probe= angle of probe light in vacuum in radians
%     theta_probe_deg= angle of probe light in vacuum in degrees
%     theta_pump= angle of pump light in vacuum in radians
%     theta_pump_deg= angle of pump light in vacuum in degrees
%     thick(n_film)= thickness of film # n_film
%     thick_0(n_film)= mean thickness of film # n_film
%     tran_left(n)= transmission coefficient of a wave going to the left at
%       interface between film # n and film # n-1
%     tran_right(n)= transmission coefficient of a wave going to the right at
%       interface between film # n and film # n+1      
%     t_save= time interval between saves to disk
%     t_stop= time to stop
%     vel(n)= sound velocity in bin # n
%     velocity(n_film)= velocity in film # n_film
%     z(n)= right hand end of bin # n
                 
clear all
fclose('all')
%%%Declare Matrices%%%
    ac_imp=zeros(1,2002);
    vel=zeros(1,2002);
    z=zeros(1,2002);
    refl_left=zeros(1,2002);
    refl_right=zeros(1,2002);
    tran_left=zeros(1,2002);
    tran_right=zeros(1,2002);
    attenuation=zeros(1,41); 
    alpha=zeros(1,41);   
    diffusion=zeros(1,41); 
    thick=zeros(1,41); 
    thick_0=zeros(1,41); 
    bin_size=zeros(1,41); 
    boundary=zeros(1,41); 
    roughness=zeros(1,41); 
    g_0=zeros(1,2002);
    g_1=zeros(1,2002);
    g_2=zeros(1,2002);
    g_left=zeros(6,2002); g_right=zeros(6,2002);
    eta_left=zeros(1,2002);
    eta_right=zeros(1,2002);
    eta_left_new=zeros(1,2002);
    eta_right_new=zeros(1,2002); 
    f_sens=zeros(1,2002);
    velocity=zeros(1,41);
    density=zeros(1,41);
    r_t=zeros(1,5001);
    
    %complex numbers
    epsilon_pump=zeros(1,41);
    epsilon_probe=zeros(1,41);
    d_epsilon_d_strain=zeros(1,41);
    
    m_pump=zeros(40,2);
    m_pump(:,:,2)=zeros(40,2);
    
    m_pump_end=zeros(2,2);
    m_pump_total=zeros(2,2);
    
    m_probe=zeros(40,2);
    m_probe(:,:,2)=zeros(40,2);
      
    m_probe_end=zeros(2,2);
    m_probe_total=zeros(2,2);
    
    d_m_probe=zeros(40,2);
    d_m_probe(:,:,2)=zeros(40,2);
    
    
    n_bins=zeros(1,41);
    %material=zeros{1,41};
    
    
       
%%%OPEN SOME FILES%%%
%     file with data about the structure to be studied *****************
    fid1=fopen('r_t_6.inp','r');  
%     file for output strain distributions *****************************
    fid2=fopen('r_t_eta.dat','w');  
%     file for calculated acoustic reflection and transmission coefficients
    fid3=fopen('r_t_ac.dat','w');
%     file for reflection change as a function of time
    fid4=fopen('r_t.dat','w');
%     file for sensitivity function as a function of position
    fid5=fopen('r_t_sens.dat','w'); 
     
%     calculation of the film thicknesses and bin sizes  ***************

    str1=fgetl(fid1);
    tau_data=sscanf(str1,'%e');
    str1=fgetl(fid1);
    n_tau_divide=sscanf(str1,'%d');
    str1=fgetl(fid1);
    t_save=sscanf(str1,'%e');
    str1=fgetl(fid1);
    t_stop=sscanf(str1,'%e');
    str1=fgetl(fid1);
    info=sscanf(str1,'%e %e %e');
    lambda_pump=info(1);theta_pump_deg=info(2);i_polzn_pump=info(3);
    str1=fgetl(fid1);
    info=sscanf(str1,'%e %e %e');
    lambda_probe=info(1);theta_probe_deg=info(2);i_polzn_probe=info(3);
    str1=fgetl(fid1);
    n_films_total=sscanf(str1,'%d');
    str1=fgetl(fid1);
    n_samples=sscanf(str1,'%d');
    for j=1:n_films_total
        str1=fgetl(fid1);
        info=sscanf(str1,'%e %e %s');
        thick_0(j)=info(1); roughness(j)=info(2);
        material{j}=char(info(3:size(info)));
        material{j}=material{j}';
    end
    fclose(fid1);
      k_pump=2.0*pi/lambda_pump;
      k_probe=2.0*pi/lambda_probe;
      theta_pump=theta_pump_deg*pi/180.0;
      theta_probe=theta_probe_deg*pi/180.0;
      
      
                                                     
%     read in material properties **************************************
      for j=1:n_films_total
          fid6=fopen(material{j},'r');
          str2=fgetl(fid6);
          density(j)=sscanf(str2,'%e');
          str2=fgetl(fid6);
          velocity(j)=sscanf(str2,'%e');
          str2=fgetl(fid6);
          info=sscanf(str2,'%e %e');
          n_pump=info(1); kappa_pump=info(2);
          epsilon_pump(j)=(n_pump+i*kappa_pump)^2;
          str2=fgetl(fid6);
          info=sscanf(str2,'%e %e');
          n_probe=info(1);kappa_probe=info(2);
          epsilon_probe(j)=(n_probe+i*kappa_probe)^2;
          
          str2=fgetl(fid6);
          info=sscanf(str2,'%e %e');
          d_n_d_strain=info(1); d_kappa_d_strain=info(2);
          d_epsilon_d_strain(j)=2.0*(n_probe+i*kappa_probe)*...
              (d_n_d_strain+i*d_kappa_d_strain);
          
          str2=fgetl(fid6);
          alpha(j)=sscanf(str2,'%e');
          str2=fgetl(fid6);
          diffusion(j)=sscanf(str2,'%e');
          str2=fgetl(fid6);
          attenuation(j)=sscanf(str2,'%e');
          fclose(fid6);
      end
         
                                                                                                 
     tau=tau_data/n_tau_divide;  % time step for algorithm
                                                                 
%     ******************************************************************
%     this is the beginning of the loop over thickness fluctuations ****
%     ******************************************************************
      n_data=t_stop/tau_data;
      n_mother=floor(t_stop/tau_data);
      sprintf ('How are you? %f %f',n_data,n_mother)
      
      if (n_data > 5000) 
        sprintf('Too many data times for array')
      end
      
      for n_average=1:n_samples
         sprintf('%d',n_average)
         for n_film=1:n_films_total  %choose thicknesses
             x=gaussrand;
          thick(n_film)=thick_0(n_film)+x*roughness(n_film);
         end
         %save thicknesses to file
         fprintf(fid4,'%f ',thick(1:n_films_total));
         
   %     find the number and size of bins in each film and the total ******
   %     number of bins ***************************************************
         n_bins_total=0;
         for n_film=1:n_films_total
            bin_size(n_film)=tau*velocity(n_film);
            n_bins(n_film)=floor(thick(n_film)/bin_size(n_film)+0.5);
            if (n_bins(n_film)==0)  
          error(' The thickness of film # %d has been chosen  so small \r that the film contains zero bins. \r This will crash the program.',n_film)
            
            end 
          bin_size(n_film)=thick(n_film)/n_bins(n_film);
          n_bins_total=n_bins(n_film)+n_bins_total;
         end
       
         sprintf(' Total number of bins is %d ',n_bins_total)
         
         %     calculate geometry of structure **********************************
         z(1)=0.0;
         n_bin=1;
         for n_film=1:n_films_total
           for n=1:n_bins(n_film)
             n_bin=n_bin+1;
             z(n_bin)=z(n_bin-1)+bin_size(n_film);
             ac_imp(n_bin)=density(n_film)*velocity(n_film);
             vel(n_bin)=velocity(n_film);
           end
         end

         
         %     add bins to get correct boundary conditions
         ac_imp(1)=0.0;
         boundary(1)=0.0;
         vel(1)=vel(2);
         ac_imp(n_bins_total+2)=ac_imp(n_bins_total+1);
         vel(n_bins_total+2)=vel(n_bins_total+1);
         
         
         %     calculate the pump optical transfer matrices *********************
        [m_pump,m_pump_end]=pump_matrices(k_pump,theta_pump,i_polzn_pump,...
         epsilon_pump,bin_size,n_films_total);
     
         %     calculate the reflection coefficient for the pump ****************
        [m_pump_total,r_pump,max_pump_bin]=pump_refln_coeff(n_bins,...
            n_bins_total,n_films_total,m_pump,m_pump_end,fid4);
        r_pump;
        m_pump_total;
        
        %     calculate the initial strain distribution ************************
      [eta_left,eta_right]=initialstrain(bin_size,n_bins,eta_left,eta_right,...
          n_films_total,r_pump,alpha,diffusion,m_pump,n_bins_total,max_pump_bin);

      
      
      
        %     calculate largest value of strain for plotting *******************
      eta_scale=0.0;
      for n=1:n_bins_total
        eta_scale=max([eta_scale abs(eta_left(n)+eta_right(n))]);
      end
      
      %     calculate the probe optical transfer matrices ********************
      [m_probe,m_probe_end,d_m_probe]= probe_matrices(k_probe,theta_probe,...
          i_polzn_probe,epsilon_probe,d_epsilon_d_strain,bin_size,...
          n_films_total,d_m_probe,m_probe);
      
      %     calculate the reflection coefficient for the probe ***************
      [m_probe,m_probe_end,m_probe_total,r_probe,max_probe_bin]= ...
          probe_refln_coeff(n_bins,fid4,...
          n_bins_total,n_films_total,m_probe,m_probe_end,m_probe_total);
      
      %     calculate the sensitivity function *******************************
      [f_sens]= sensitivity_function(n_bins,n_bins_total,n_films_total, ...
          z,m_probe,d_m_probe,m_probe_end,m_probe_total,r_probe,f_sens,...
          max_probe_bin,fid5);
      
      %     calculation of the acoustic reflection and transmission coefficients 
      [refl_left,refl_right,tran_left,tran_right]= acoust_refln_coeff(ac_imp,vel,...
           refl_left,refl_right,tran_left,tran_right,n_bins_total,fid3);
       
       ref=refl_left(2)
       ref=refl_right(2)
       
      %     calculate coefficients for generalized stepping algorithm ********
      [g_left,g_right]= g_s_a(n_bins,n_bins_total,n_films_total,bin_size,...
          tau,attenuation,velocity,refl_left,refl_right,tran_left,...
          tran_right,g_0,g_1,g_2,g_left,g_right);
      

      
      i_save=0;
      
%     begin stepping the strain pulses *********************************
      for i_data=1:n_data+1   
        t=(i_data-1)*tau_data;
%       save strain distribution ***************************************
        if (mod(t+0.5*tau,t_save) < tau) 
         sprintf('Time is %e psecs. Strain distribution saved.',t) 
%        the call to the save_strain subroutine is disabled in order **
%        in order to avoid the creation of very large files *********** 
        [i_save]=save_strain(i_save,z,n_bins_total,eta_left,eta_right,...
            eta_scale,fid2);
        end
        
        %       calculate change in reflectivity of the probe ******************
        d_refln=0.0;
        for n=1:max_probe_bin
        d_refln=d_refln+(eta_left(n)+eta_right(n))*f_sens(n);
        end 
        
        %       record change in reflectivity **********************************
        r_t(i_data)=r_t(i_data)+d_refln ;

        for n_divide=1:n_tau_divide
%         calculate new strain distribution ****************************
        [eta_left,eta_right,eta_left_new,eta_right_new]= step(n_bins_total,...
            g_left,g_right,eta_left,eta_right,eta_left_new,eta_right_new);                            
        end 

      end
        
      end 
%     ******************************************************************
%     this is the end of the loop over thickness fluctuations **********
%     ******************************************************************      

%     save time and reflectivity change ********************************
for i_data=1:n_data
    fprintf(fid4,'%e   %e \r',(i_data-1)*tau_data,r_t(i_data)/n_samples);
end

fclose('all');

clf;

fh = figure(2);
clf;
set(fh,'Units','inches','Position', [0 2 4 3])
set(fh, 'color', 'white');
plot([1:n_data]*tau_data,r_t(1:3500))
hold on;
axis([0 700 -2.3e-4 2.3e-4])

time = 62;
plot([time 1.00001*time], [-100,100],'r','LineWidth',1);

time = 124
plot([time 1.00001*time], [-100,100],'r','LineWidth',1);
% time = 62*3
% plot([time 1.00001*time], [0,100*in(x1)],'r','LineWidth',1);


time = 211.7
plot([time 1.00001*time], [-100,100],'r','LineWidth',1);

time = 211.7+62
plot([time 1.00001*time], [-100,100],'r','LineWidth',1);

time = 335
plot([time 1.00001*time], [-100,100],'r','LineWidth',1);


% time = 26.5;
% plot([time 1.00001*time], [-100,100],'r','LineWidth',1);
% 

      
