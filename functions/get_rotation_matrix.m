function [L_sigma_S, L_sigma]=get_rotation_matrix(phi,theta,phi_S,theta_S,mu)

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

elseif ((mu ~= 0 || mu ~= 180) && ((del_phi>0 &&(180<del_phi && del_phi<360)) || del_phi<0))

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

if phi>180
    sigma=180-sigma;
    sigma_S=180-sigma_S;
else
    sigma=sigma;
    sigma_S=sigma_S;
end
end
