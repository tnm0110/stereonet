function [tr_p, pl_p] = plane2pole(tr,pl)

%  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
%  
%         Simple function to calculate pole of a given plane.
%  [tr_p pl_p] = plane2pole(tr,pl) returns trend (tr_p) and plunge (pl_p) 
%  of thepoles given strike and dip of a plane
%    
%          Input and output angles should be in radians

% ######################################################################

%Calulate the direction cosine of the pole
 cd = cos(pl);
 ce = -sin(pl) * cos(tr);
 cn = sin(pl) * sin(tr);

%convert NED to spherical coordinates to get trend and plunge of the pole
[tr_p,pl_p] = ned2sph(cn,ce,cd);

end