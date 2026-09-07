function [S_f] = transmission_muller_matrix_and_stokes_param(n_air,n_medium,theta_skylight,S_i)
    
    n1=n_air; %Air refractive index
    n2=n_medium; %glass refractive index
    
    theta_in=theta_skylight;
    
    theta_t=asind(n1*sind(theta_in)/n2); %Snell's law
    
    t_s= 2*n1*cosd(theta_in)./(n1*cosd(theta_in)+n2*cosd(theta_t)); %Perpendicular component of transmitted light
    
    t_p= 2*n1*cosd(theta_in)./(n2*cosd(theta_in)+n1*cosd(theta_t)); %Parallel component of transmitted light
    
    M_r11=t_s.^2+t_p.^2;
    M_r12=t_s.^2-t_p.^2;
    
    M_r21=M_r12;
    M_r22=M_r11;
    
    M_r33=2.*t_s.*t_p;
    M_r34=0;
    
    M_r43=0;
    M_r44=M_r33;
    
    %%%Transmission Mueller matrix
    M_transmitted=(1/2)*[M_r11 M_r12 0 0; M_r21 M_r22 0 0; 0 0 M_r33 M_r34; 0 0 M_r43 M_r44];

    S_f=M_transmitted*S_i;
    S_f=S_f';

   end
