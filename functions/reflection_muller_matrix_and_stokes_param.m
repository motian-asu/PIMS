function [S_f] = reflection_muller_matrix_and_stokes_param(n_air,n_medium,theta_skylight,S_i)
            %Calculate M_refrac and Multiply with S_i
            %%Paper: The polarization patterns of skylight reflected off wave water surface
            n1=n_air; %Refractive index of the incident medium
            n2=n_medium; %Refractive index of the transmitted medium
            theta_inc=theta_skylight;
            theta_trans=asind(n1*sind(theta_inc)/n2); %Snell's law

            rs=(n1*cosd(theta_inc)-n2*cosd(theta_trans))./(n1*cosd(theta_inc)+n2*cosd(theta_trans));
            rp=(n2*cosd(theta_inc)-n1*cosd(theta_trans))./(n2*cosd(theta_inc)+n1*cosd(theta_trans));
            M_ref_11=rs^2+rp^2;
            M_ref_12=rp^2-rs^2;
            M_ref_21=M_ref_12;
            M_ref_22=M_ref_11;
            M_ref_33=2.*rs.*rp;
            M_ref_44=2.*rs.*rp;

            M_reflect=(1/2)*[
                            M_ref_11 M_ref_12 0 0
                            M_ref_21 M_ref_22 0 0
                            0 0 M_ref_33 0
                            0 0 0 M_ref_44];
                        
            S_f=M_reflect*S_i;
            S_f=S_f';
   end
