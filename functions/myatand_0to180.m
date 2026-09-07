function v=myatand_0to180(y,x)
   
   %Specially for Angle of polarization calculation
   %https://en.wikipedia.org/wiki/Atan2
   %_Stokes parameters of skylight based on simulations and polarized radiometer measurements
   %Table 2.1: Ranges of the angle of polarization xhi determined by signs of Q and U for the linearly polarized light
   %---returns an angle in degree between 0 and 180 for atand
   
   v=zeros(size(x));
   
   v(x>0 & y>0) = atand( y(x>0 & y>=0) ./ x(x>0 & y>=0) );
   v(x>0 & y==0) = 0;
   v(x>0 & y<0) = 2*180+atand( y(x>0 & y<0) ./ x(x>0 & y<0) );
   
   v(x<0 & y>0) = 180+atand( y(x<0 & y>0) ./ x(x<0 & y>0) );
   v(x<0 & y==0) = 180;   
   v(x<0 & y<0)  = 180+atand( y(x<0 & y<0) ./ x(x<0 & y<0) );
   
   v(x==0 & y>0) = 180/2;
   v(x==0 & y==0) = NaN;
   v(x==0 & y<0) = 180+180/2;   
   
   v(isnan(x)==1 & isnan(y)==1) = NaN;
   
   end
   
   %%
