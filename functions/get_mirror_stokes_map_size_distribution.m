function [I_total,Q_total,U_total,V_total] = get_mirror_stokes_map_size_distribution(phi_mirror,theta_mirror,phi,theta,m,n_med,x,k0,S)
%% Get the rotation matrix
%Looking or camera direction
phi_S=phi_mirror;
theta_S=theta_mirror;

L_sigma_S=cell(length(phi_S),length(theta_S));
L_sigma=cell(length(phi_S),length(theta_S));


for i=1:length(phi_S)
    for j=1:length(theta_S)
        
        %As some of the phi_S are problematic (produces invalid rotation angles) we try to calculate close to these values
        %but not exactly at those values of phi_S
        if phi_S(i)==phi
            phi_S(i)=phi_S(i)+0.001;
        elseif phi_S(i)==phi-180
            phi_S(i)=phi_S(i)+0.001;
        elseif theta_S(j)==0
            theta_S(j)=theta_S(j)+0.001;
        elseif theta_S(j)==90
            theta_S(j)=theta_S(j)-0.001;
        end
        
        if theta==0
            theta=theta+1e-3;
        elseif theta==90
            theta=theta+1e-3;
        elseif theta==180
            theta=theta+1e-3;
        end
            
        mu(i,j)=(acosd(sind(theta)*sind(theta_S(j))*cosd(phi_S(i)-phi)+cosd(theta)*cosd(theta_S(j))));
        [L_sigma_S{i,j}, L_sigma{i,j}]=get_rotation_matrix(phi,theta,phi_S(i),theta_S(j), mu(i,j));
                
    end
end

scattering_angle=mu;
n_r=m/n_med;

% I_total=zeros(length(phi_S),length(theta_S));
% Q_total=zeros(length(phi_S),length(theta_S));
% U_total=zeros(length(phi_S),length(theta_S));
% V_total=zeros(length(phi_S),length(theta_S));

for p_count=1:length(x)
    x_select=x(p_count);
    
    %% Single Mie Scattering
    I=zeros(length(phi_S),length(theta_S));
    Q=zeros(length(phi_S),length(theta_S));
    U=zeros(length(phi_S),length(theta_S));
    V=zeros(length(phi_S),length(theta_S));
    for w=1:length(phi_S)
        for v=1:length(theta_S)
            if x_select~=0
                u=cos(pi/180*scattering_angle(w,v));

                %% Mie Scattering cross section
                % result=mie(m,x);
                result=mie(n_r,x_select);
                qsca=result(5);

                %% S1 and S2 
                % result_s1_s2=mie_S12(m, x, u);
                result_s1_s2=mie_S12(n_r, x_select, u);

                S1=result_s1_s2(1);
                S2=result_s1_s2(2);

                S1_T(w,v)=S1;
                S2_T(w,v)=S2;
                %% Mie Scattering Matrix
%                 M11=2*pi/k0^2/qsca*(abs(S1)^2+abs(S2)^2);
%                 M12=2*pi/k0^2/qsca*(abs(S2)^2-abs(S1)^2);
%                 M33=2*pi/k0^2/qsca*(S2*conj(S1)+S1*conj(S2));
%                 M34=2*pi/k0^2/qsca*(S2*conj(S1)-S1*conj(S2))*1i;

                M11=1/2*(abs(S1)^2+abs(S2)^2);
                M12=1/2*(abs(S2)^2-abs(S1)^2);
                M33=1/2*(S2*conj(S1)+S1*conj(S2));
                M34=1/2*(S2*conj(S1)-S1*conj(S2))*1i;

                M=[M11 M12 0 0
                    M12 M11 0 0
                    0 0 M33 M34
                    0 0 -M34 M33];

                %%
                S_S=(L_sigma_S{w,v}*M*L_sigma{w,v})*S; % scattering in air (skylight)

                I(w,v)=S_S(1);
                Q(w,v)=S_S(2);
                U(w,v)=S_S(3);
                V(w,v)=S_S(4);
            else    
                if phi_S(w)==phi && theta_S(v)==theta
                    I(w,v)=1;
                    Q(w,v)=0;
                    U(w,v)=0;
                    V(w,v)=0;
                end
            end
        end
    end
%     I_total=I_total+I*x(p_count,2);
%     Q_total=Q_total+Q*x(p_count,2);
%     U_total=U_total+U*x(p_count,2);
%     V_total=V_total+V*x(p_count,2);
    
    I_total{p_count,1}=I;
    Q_total{p_count,1}=Q;
    U_total{p_count,1}=U;
    V_total{p_count,1}=V;
       
end
end
