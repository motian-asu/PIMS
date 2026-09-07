function [cam_az,cam_zen,I_Pol,Q_Pol,U_Pol,V_Pol,cam_zen_Pol,cam_az_Pol,DoLP_Pol,AoP_Pol,DoCP_Pol,DoP_Pol,un_normalized_I_pol,un_normalized_Q_pol,un_normalized_U_pol,un_normalized_V_pol]=get_skylight_Rayleigh_Scat(theta,phi,cam_az_T,cam_zen_T)
    
%%  Condition for the simulation  
    error=0;
    lamda=480e-9; %%Wavelength
    Temp=5782; %Sun's blackbody temperature in K on TOA (top of Atmosphere)
    h=6.6261e-34; %plank's constant in Js
    c=3e8; %light velocity in ms-1;
    kb=1.3806e-23; %boltzmann's contant in JK-1
    %%% Solar irradiance measured at earth surface
    %%% from : http://www.oceanopticsbook.info/view/light_and_radiometry/level_2/blackbody_radiation#E:planckirradwavenm2
    Rsun=6.95e5; %suns radius in km
    Rearth=1.496e8; % radius of earth's orbit in km
    I=(Rsun/Rearth)^2*2*pi*h*c^2/lamda^5*(1/(exp(h*c/(lamda*kb*Temp))-1))*lamda; % in Wm-2
    % %particle detail
    r=100e-9; %Particle radius
    mr=1.53; %Real part of complex refractive index m = mr + i*mi
    mi=0.007;%Imaginary part of complex refractive index m = mr + i*mi
    m=mr+1i*mi;%Complex refractive index m = mr + i*mi
    % For unpolarized light I=I; For LP|| I,Q=I; For LP_per I=I,Q=-I; 
    I=I; Q=0; U=0;V=0;
    %%%Add Irridiance of SUN : I vs lamda -data or formula
    S=[I;  Q;  U;  V];
    %%
cam_zen_T=cam_zen_T;
cam_az_T=cam_az_T;

flag1=0;
for ii=1:length(cam_az_T)
    if cam_az_T(ii)<0
        cam_az_T(ii)=cam_az_T(ii)+360;
        flag1=1;
    end
end

%if there is a negative number in the cam_az_T vector then we need to sort it
if flag1==1
    cam_az_T=sort(cam_az_T);
end


        S_M_T=cell(length(cam_az_T),length(cam_zen_T));
        DoLP_M_T=cell(length(cam_az_T),length(cam_zen_T));
        DoP_M_T=cell(length(cam_az_T),length(cam_zen_T));
        DoCP_M_T=cell(length(cam_az_T),length(cam_zen_T));
        Q_M_T=cell(length(cam_az_T),length(cam_zen_T));
        U_M_T=cell(length(cam_az_T),length(cam_zen_T));
        V_M_T=zeros(length(cam_az_T),length(cam_zen_T));
        AoP_M_T=cell(length(cam_az_T),length(cam_zen_T));
        I_M_T=cell(length(cam_az_T),length(cam_zen_T));
        Il_M_T=cell(length(cam_az_T),length(cam_zen_T));
        Ir_M_T=cell(length(cam_az_T),length(cam_zen_T));
        rho_M_T=cell(length(cam_az_T),length(cam_zen_T));
        Iunp_M_T=cell(length(cam_az_T),length(cam_zen_T));
        Il_unp_M_T=cell(length(cam_az_T),length(cam_zen_T));
        Ir_unp_M_T=cell(length(cam_az_T),length(cam_zen_T));
        Il_lin_M_T=cell(length(cam_az_T),length(cam_zen_T));
        Ir_lin_M_T=cell(length(cam_az_T),length(cam_zen_T));
        mu_M_T=cell(length(cam_az_T),length(cam_zen_T));
        beta_M_T=cell(length(cam_az_T),length(cam_zen_T));
        kR_M_T=cell(length(cam_az_T),length(cam_zen_T));
        sigma_M_T=cell(length(cam_az_T),length(cam_zen_T));
        sigma_S_M_T=cell(length(cam_az_T),length(cam_zen_T));
        
        if error==1
    
    er_DoP_M_T=cell(length(cam_az_T),length(cam_zen_T));
    er_AoP_M_T=cell(length(cam_az_T),length(cam_zen_T));
    er_Q_M_T=cell(length(cam_az_T),length(cam_zen_T));
    er_U_M_T=cell(length(cam_az_T),length(cam_zen_T));

        end
       
for i=1:length(cam_az_T)
    
    phi_S=cam_az_T(i);
    
    for j=1:length(cam_zen_T)
        
        theta_S=cam_zen_T(j);
        
        %As some of the phi_S are problematic (produces invalid rotation angles) we try to calculate close to these values
        %but not exactly at those values of phi_S
        if phi_S==phi
            phi_S=phi_S+0.001;
            cam_az_T(i)=phi_S;
        elseif phi_S==phi-180
            phi_S=phi_S+0.001;
            cam_az_T(i)=phi_S;
        elseif theta_S==0
            theta_S=theta_S+0.001;
            cam_zen_T(j)=theta_S;
        elseif theta_S==90
            theta_S=theta_S-0.001;
            cam_zen_T(j)=theta_S;
        end

% mu: scattering angle
mu=(acosd(sind(theta)*sind(theta_S)*cosd(phi_S-phi)+cosd(theta)*cosd(theta_S)));

%Scattering matrix or Muller Matrix (w.r.t scattering plane: Which holds the incident and scattering light)
%Rayleigh Scattering (partilce smaller than incident wavelength) Matrix--for homogenous or radially inhomogeneous spherical
%particle
F11=(1+(cosd(mu)).^2)/2;
F12=-(1-(cosd(mu)).^2)/2;
F21=F12;
F22=F11;
F33=cosd(mu);
F34=0;
F43=0;
F44=F33;

F=[F11 F12 0 0; F21 F22 0 0; 0 0 F33 F34; 0 0 F43 F44];

% %particle detail
% r=Particle radius
% mr=Real part of complex refractive index m = mr + i*mi
% mi=Imaginary part of complex refractive index m = mr + i*mi
% m=mr+1i*mi=Complex refractive index m = mr + i*mi

a=2*pi*r/lamda; % Particle Size parameter

F_Rayl= a.^6*abs(((m.^2-1)/(m.^2+2)))^2.*F; %Rayley Scattering Matrix

%Stokes parameter after scattering (w.r.t scattering plane)
S_S=F_Rayl*S; %Matrix multiplication

%R is the radial distance from the scattering particle
%%From Equation: 2.74 (_Stokes parameters of skylight based on simulations and polarized radiometer measurements)
kR=sqrt(F11*S(1)/S_S(1)); %Will be used to normalize F matrix

%Converting Scattering matrix in scattering plane to Phase matrix in local Meridian Plane

% From: _Stokes parameters of skylight based on simulations and polarized
% radiometer measurements:
% The scattering matrix is defined with respect to the scattering plane. This plane transforms for different scattering
% events. Thus, the corresponding coordinate systems describing the incident
% and the scattered light are also not fixed. In order to represent interactions
% of light and particles in a uniform reference coordinate system, the scattering matrix should be transformed into the phase matrix that is defined
% relative to the local meridian plane.

%Calculation of ratation angle: incident (sigma) and scattering (sigma_S) rotation angle

%From: http://www.oceanopticsbook.info/view/light_and_radiometry/level_2/polarization_scattering_geometry#Eq:rotations1

phi_n=phi;

if phi>180
    phi=phi-180;
end

del_phi=phi_S-phi;

if ((mu ~= 0 || mu ~= 180) && ((0<=del_phi && del_phi<=180) && del_phi>=0))
    
    if (theta_S~=0 && theta~=0)
        sigma=acosd((cosd(theta_S)-cosd(theta)*cosd(mu))/(sind(mu)*sind(theta)));
        sigma_S=180+acosd((cosd(theta)-cosd(theta_S)*cosd(mu))/(sind(mu)*sind(theta_S)));       
    end
    
    if theta_S==0 % we don't have the condition when phi_S==0 or 360
        sigma=0;
        sigma_S=acosd(cosd(theta)*cosd(phi_S-phi));
    end
    
    if theta==0
        sigma=acosd(cosd(theta_S)*cosd(phi_S-phi));
        sigma_S=0;
    end
    
    if (theta_S==0 && theta==0)
    
        sigma=(phi_S-phi);
  
    end
%     if sind(theta_S)==0
%         sigma=(acosd(cosd(theta_S)));
%         sigma_S=acosd(-cosd(theta)*cosd(phi_S-phi));
%     elseif sind(theta)==0
%         sigma_S=(acosd(cosd(theta)));
%         sigma=acosd(-cosd(theta)*cosd(phi_S-phi));
%     end

elseif ((mu ~= 0 || mu ~= 180) && ((del_phi>0 &&(180<del_phi && del_phi<360)) || del_phi<0))
% elseif ((mu ~= 0 || mu ~= 180) && (del_phi>0 &&(180<del_phi && del_phi<360)))
    
    if (theta_S~=0 && theta~=0)
    sigma=180-acosd((cosd(theta_S)-cosd(theta)*cosd(mu))/(sind(mu)*sind(theta)));
    sigma_S=360-acosd((cosd(theta)-cosd(theta_S)*cosd(mu))/(sind(mu)*sind(theta_S)));
    end
    
    if ((phi_S==0 || phi_S==360) && theta_S==0)
        sigma=0;
        sigma_S=acosd(cosd(theta)*cosd(phi));
    end
    
    if ((phi_S~=0 || phi_S~=360) && theta_S==0)
        sigma=0;
        if phi_S>phi
            sigma_S=acosd(cosd(theta)*cosd(phi+(360-phi_S)));
        else
            sigma_S=acosd(cosd(theta)*cosd(phi-phi_S));
        end
    end
    
    if ((phi_S==0 || phi_S==360) && theta==0)
        sigma=acosd(cosd(theta_S)*cosd(phi));
        sigma_S=0;
    end
    
    if ((phi_S~=0 || phi_S~=360) && theta==0)
        sigma_S=0;
        if phi_S>phi
            sigma=acosd(cosd(theta_S)*cosd(phi+(360-phi_S)));
        else
            sigma=acosd(cosd(theta_S)*cosd(phi-phi_S));
        end
    end
    
    if (theta_S==0 && theta==0 && (phi_S==0 || phi_S==360))
    
        sigma=phi;
        
    elseif (theta_S==0 && theta==0 && (phi_S~=0 || phi_S~=360))
        if phi_S>phi
            sigma=phi+(360-phi_S);
        else
            sigma=phi-phi_S;
        end
    end
    
%     if sind(theta_S)==0
%         sigma=(acosd(cosd(theta_S)));
%         sigma_S=acosd(cosd(theta)*cosd(phi_S-phi)); %from my derivation
%             
%     elseif sind(theta)==0
%         sigma_S=(acosd(cosd(theta)));
%         sigma=acosd(-cosd(theta)*cosd(phi_S-phi));
%     end
%         
% elseif ((mu ~= 0 || mu ~= 180) && del_phi<0)
%     
% %     if (theta_S~=0 && theta~=0)
%     sigma=acosd((cosd(theta_S)-cosd(theta)*cosd(mu))/(sind(mu)*sind(theta)));
%     sigma_S=-180+acosd((cosd(theta)-cosd(theta_S)*cosd(mu))/(sind(mu)*sind(theta_S)));
% %     end

elseif (mu==0 || mu==180)
    
    sigma=0;
    sigma_S=0;
    
end
phi=phi_n;
%%% Rotation matrix calculated from sigma and sigma_S: we are looking to
%%% the beam (or we see that beam is coming towards us/detector)
L_sigma_S=[1 0 0 0; 0 cosd(2*(sigma_S)) sind(2*(sigma_S)) 0; 0 -sind(2*(sigma_S)) cosd(2*(sigma_S)) 0; 0 0 0 1];
L_sigma=[1 0 0 0; 0 cosd(2*(sigma)) sind(2*(sigma)) 0; 0 -sind(2*(sigma)) cosd(2*(sigma)) 0; 0 0 0 1];

if phi>180
    L_sigma_S=[1 0 0 0; 0 cosd(2*(180-sigma_S)) sind(2*(180-sigma_S)) 0; 0 -sind(2*(180-sigma_S)) cosd(2*(180-sigma_S)) 0; 0 0 0 1];
    L_sigma=[1 0 0 0; 0 cosd(2*(180-sigma)) sind(2*(180-sigma)) 0; 0 -sind(2*(180-sigma)) cosd(2*(180-sigma)) 0; 0 0 0 1];
end
%%% Rotation matrix calculated from sigma and sigma_S: we are looking into
%%% the beam (or we see that beam is going towards us/detector)
% if phi>180
% L_sigma_S=[1 0 0 0; 0 cosd(2*(sigma_S)) -sind(2*(sigma_S)) 0; 0 sind(2*(sigma_S)) cosd(2*(sigma_S)) 0; 0 0 0 1];
% L_sigma=[1 0 0 0; 0 cosd(2*sigma) -sind(2*sigma) 0; 0 sind(2*sigma) cosd(2*sigma) 0; 0 0 0 1];
% end

%Phase matrix from scattering matrix
% P=L_sigma_S*F_Rayl*L_sigma;
%Phase matrix from scattering matrix
if theta==0 && theta_S==0
    P=F_Rayl*L_sigma;
else
    P=L_sigma_S*F_Rayl*L_sigma;
end

%Final stokes parameter

S_S_C=P*S;

%%%Uncomment if you don't want to normalize
kR=1;

S_S_F=S_S_C./(kR)^2; %Nomrlaized Equation: 2.88 (_Stokes parameters of skylight based on simulations and polarized radiometer measurements)

non_normalized_S=S_S_F;

S_S_F=S_S_F./S_S_F(1); %This is to normalize the Stokes parameters

DoP_S= sqrt(S_S_F(2)^2+S_S_F(3)^2+S_S_F(4)^2)/S_S_F(1);
DoLP_S=sqrt(S_S_F(2)^2+S_S_F(3)^2)/S_S_F(1);
DoCP_S=S_S_F(4)/S_S_F(1);

Il_S=(S_S_F(1)+S_S_F(2))/2; %Parallel component of polarized radiance (I); Function of scattering angle; 
Ir_S=(S_S_F(1)-S_S_F(2))/2 ;%Perpendicular component of polarized radiance (I); Constant w.r.t scattering angle;

rho_S=Ir_S/Il_S; %Linear depolarization ratio - unpolarized: Il_S=Ir_S; Perpendicular: Il_S=0; Partially Pol: Il_S<Ir_S

%Scattered light is partially polarized so, I_S = Iunp_S + Ilin_S;

%Unpolarized portion of scattered light-> Iunp_S; 
Iunp_S=S(1)*(1-DoLP_S);%Eq 2.37 (_Stokes parameters of skylight based on simulations and polarized radiometer measurements)

Il_unp_S= Iunp_S/2;%Parralel portion of unpolarized light
Ir_unp_S= Iunp_S/2; %Perpendicular portion of unpolarized light

%Polarized portion of scattered light-> Ilin_S;
Il_lin_S=Il_S-Il_unp_S; %Parralel portion of polarized light
Ir_lin_S=Ir_S-Ir_unp_S ;%Perpendicular portion of polarized light

%angle of polarization
xhi_S=0.5*myatand_0to180(S_S_F(3),S_S_F(2)); %for 0 to 180 degree

% % %%uncomment his if statement-> For AoP in -90 to 90 degree
% if xhi_S>90
%     xhi_S=-(180-xhi_S);
% end

%Elipticity angle
beta_S= 0.5*((atand(S_S_F(4)./sqrt(S_S_F(2).^2+S_S_F(3).^2))));

S_M_T{i,j}=S_S_F;

if isnan(DoP_S)==1
    DoP_M_T{i,j}=DoP_S;
    DoLP_M_T{i,j}=DoLP_S;
    DoCP_M_T{i,j}=DoCP_S;
    AoP_M_T{i,j}=xhi_S;
    beta_M_T{i,j}=beta_S;
    Q_M_T{i,j}=S_S_F(2);
    U_M_T{i,j}=S_S_F(3);
    V_M_T(i,j)=S_S_F(4);
else
    
    error_percent=0.1;
    
    er_DoP=error*(error_percent*(-1+2 * rand));
    DoP_M_T{i,j}=DoP_S+er_DoP; %iff error is 1 then it will add 1% random error between -1 and 1

    DoLP_M_T{i,j}=DoLP_S+error*(error_percent*(-1+2 * rand));%iff error is 1 then it will add 1% random error between -1 and 1
    DoCP_M_T{i,j}=DoCP_S+error*(error_percent*(-1+2 * rand));%iff error is 1 then it will add 1% random error between -1 and 1

    er_AoP=error*(error_percent*(-180+360 * rand));
    xhi_S=xhi_S+er_AoP;%iff error is 1 then it will add 1% random error between -180 and 180

    %as we are adding/subtratcting 1% of 180, sometimes the AoP will become
    %more than 180 or less than 0. So, to bring them in the range of 0 to 180
    %we do the following
    if error==1 && xhi_S>180
        xhi_S=xhi_S-180;
    elseif error==1 && xhi_S<0
        xhi_S=xhi_S+180;
    end

    AoP_M_T{i,j}=xhi_S; %error added in the previous line

    % AoP_M{i,j}=xhi_S+error*(0.01*(-90+180 * rand));%iff error is 1 then it will add 1% random error between -90 and 90

    beta_M_T{i,j}=beta_S+error*(error_percent*(-1+2 * rand));%iff error is 1 then it will add 1% random error between -45 and 45

    er_Q=error*(error_percent*(-1+2 * rand));
    Q_M_T{i,j}=S_S_F(2)+er_Q;%if error is 1 then it will add 1% random error between -1 and 1

    er_U=error*(error_percent*(-1+2 * rand));
    U_M_T{i,j}=S_S_F(3)+er_U;%iff error is 1 then it will add 1% random error between -1 and 1

end

if error==1
    
    er_DoP_M_T{i,j}=er_DoP;
    er_AoP_M_T{i,j}=er_AoP;
    er_Q_M_T{i,j}=er_Q;
    er_U_M_T{i,j}=er_U;

end

    un_normalized_I(i,j)=non_normalized_S(1);
    un_normalized_Q(i,j)=non_normalized_S(2);
    un_normalized_U(i,j)=non_normalized_S(3);
    un_normalized_V(i,j)=non_normalized_S(4);
    
    I_M_T{i,j}=S_S_F(1);

    Il_M_T{i,j}=Il_S;
    Ir_M_T{i,j}=Ir_S;
    rho_M_T{i,j}=rho_S;
    Iunp_M_T{i,j}=Iunp_S;
    Il_unp_M_T{i,j}=Il_unp_S;
    Ir_unp_M_T{i,j}=Ir_unp_S;
    Il_lin_M_T{i,j}=Il_lin_S;
    Ir_lin_M_T{i,j}=Ir_lin_S;

    mu_M_T{i,j}=mu;
    sigma_M_T{i,j}=sigma;
    sigma_S_M_T{i,j}=sigma_S;
    kR_M_T{i,j}=kR;

   
    end
end

cam_zen1=cam_zen_T(1:length(cam_zen_T));
cam_az1=cam_az_T(1:length(cam_az_T));

%%

cam_zen=cam_zen1';
cam_az=cam_az1';

%%%%For polar plot (in origin)
AoP_Pol=reshape(cell2mat(AoP_M_T)',1,length(cam_az)*length(cam_zen))';
DoP_Pol=reshape(cell2mat(DoP_M_T)',1,length(cam_az)*length(cam_zen))';
DoLP_Pol=reshape(cell2mat(DoLP_M_T)',1,length(cam_az)*length(cam_zen))';
DoCP_Pol=reshape(cell2mat(DoCP_M_T)',1,length(cam_az)*length(cam_zen))';
Q_Pol=reshape(cell2mat(Q_M_T)',1,length(cam_az)*length(cam_zen))';
U_Pol=reshape(cell2mat(U_M_T)',1,length(cam_az)*length(cam_zen))';

if error==1
    er_DoP_Pol=reshape(cell2mat(er_DoP_M_T)',1,length(cam_az)*length(cam_zen))';
    er_AoP_Pol=reshape(cell2mat(er_AoP_M_T)',1,length(cam_az)*length(cam_zen))';
    er_Q_Pol=reshape(cell2mat(er_Q_M_T)',1,length(cam_az)*length(cam_zen))';
    er_U_Pol=reshape(cell2mat(er_U_M_T)',1,length(cam_az)*length(cam_zen))';
else
    er_DoP_Pol=zeros(length(cam_az)*length(cam_zen),1);
    er_AoP_Pol=zeros(length(cam_az)*length(cam_zen),1);
    er_Q_Pol=zeros(length(cam_az)*length(cam_zen),1);
    er_U_Pol=zeros(length(cam_az)*length(cam_zen),1);
end

mu_Pol=reshape(cell2mat(mu_M_T)',1,length(cam_az)*length(cam_zen))';
sigma_Pol=reshape(cell2mat(sigma_M_T)',1,length(cam_az)*length(cam_zen))';
sigma_S_Pol=reshape(cell2mat(sigma_S_M_T)',1,length(cam_az)*length(cam_zen))';
I_Pol=reshape(cell2mat(I_M_T)',1,length(cam_az)*length(cam_zen))';

for m1=1:length(cam_az)
    for m2=1:length(cam_zen)
    
        cam_zen_Pol(m2+length(cam_zen)*(m1-1))=cam_zen(m2);
        cam_az_Pol(m2+length(cam_zen)*(m1-1))=cam_az(m1);
        
    end
end

cam_zen_Pol=cam_zen_Pol';
cam_az_Pol=cam_az_Pol';

V_M_T=V_M_T';
V_Pol=V_M_T(:);

un_normalized_I=un_normalized_I';
un_normalized_I_pol=un_normalized_I(:);

un_normalized_Q=un_normalized_Q';
un_normalized_Q_pol=un_normalized_Q(:);

un_normalized_U=un_normalized_U';
un_normalized_U_pol=un_normalized_U(:);

un_normalized_V=un_normalized_V';
un_normalized_V_pol=un_normalized_V(:);

%%
% %DoP
% subplot(2,3,1)
% % figure
% contourf(cam_az,cam_zen,cell2mat(DoP_M_T)',500,'edgecolor','none');
% % contourf(cell2mat(DoP_M_T)',500,'edgecolor','none');
% colorbar
% colormap(parula(512))
% xlabel('camra__azimuth');
% ylabel('camera__zenith');
% title('DoP')
% set(gca, 'FontSize', 12);
% 
% %%
% %DoLP
% subplot(2,3,2)

%     figure
%     contourf(cam_az,cam_zen,cell2mat(DoLP_M_T)',500,'edgecolor','none');
%     % contourf(cell2mat(DoLP_M_T),500,'edgecolor','none');
%     colorbar
%     colormap(jet(512))
%     xlabel('camra__azimuth');
%     ylabel('camera__zenith');
%     title('DoLP')
%     set(gca, 'FontSize', 12);

% %%
% % %DoCP
% subplot(2,3,3)
% % figure
% contourf(cam_az,cam_zen,cell2mat(DoCP_M_T)',500,'edgecolor','none');
% % contourf(cell2mat(DoCP_M_T),500,'edgecolor','none');
% colorbar
% colormap(parula(512))
% xlabel('camra__azimuth');
% ylabel('camera__zenith');
% title('DoCP')
% set(gca, 'FontSize', 12);
% 
% %%
% %AoP: Angle of polarization 
% subplot(2,3,4)
    % figure
    % contourf(cam_az,cam_zen,cell2mat(AoP_M_T)',500,'edgecolor','none');
    % % contourf(cell2mat(AoP_M_T),500,'edgecolor','none');
    % colorbar
    % colormap(jet(512))
    % xlabel('camra__azimuth');
    % ylabel('camera__zenith');
    % title('AoP')
    % set(gca, 'FontSize', 12);

end

%%
