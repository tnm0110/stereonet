function [tr,pl] = ned2sph(n,e,d)

%
%   simple function to convert from cartesian(NED) to spherical coordinates 
%   [tr,pl] = ned2sph(n,e,d) returns the trend (tr) and plunge (pl) of 
%   a line given NED direction cosines (n,e,d) == (North, East, Down)
%   
%             Trend and plunge are returned in radians
%
%   
%
% ######################################################################

 
% calculate plunge
pl = asin(d);  

% calculate trend

% handling E-W trending ambiguity when n = 0
if n == 0.0 
    if e < 0.0 
        tr = (3/2)*pi;  % trend ==> west
    else
        tr = pi/2;      % trend ==> east
    end
else
    tr = atan(e/n);          
    if n < 0.0        
        tr = tr+pi;     % solving negative trend
    end

%  making trend in between 0 and 2*pi    

if tr < 0.0
    tr = tr + 2*pi;
elseif tr >= 2*pi
    tr = tr - 2*pi;
end

end



end