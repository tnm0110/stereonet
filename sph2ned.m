function [n,e,d] = sph2ned(tr,pl,varargin)

%   
%    Simple function to convert from spherical to cartesian coordinate 
%    system. 
%    [n,e,d] = sph2ned(tr,pl) computes direction cosines of the line given 
%    trend (tr) and plunge(pl)
%    'varargin' is the optional input if given ( = 'pole' ) will calculate 
%    the direction cosine of the pole of the plane given stirke and dip.   
%        
%             tr and pl must be in radians

if isempty(varargin)
    %disp(['Calculate the direction cosine a the line']);
    n = cos(pl) * cos(tr);
    e = cos(pl) * sin(tr);
    d = sin(pl);

else
    %disp(['Calculate the direction cosine of the pole of a plane']);
    
    n = sin(pl) * sin(tr);
    e = -sin(pl) * cos(tr);
    d = cos(pl);


end

end