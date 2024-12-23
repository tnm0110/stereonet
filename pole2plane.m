function [st, dip] = pole2plane(tr,pl)

%  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
%  
%         Simple function to calculate plane of a given pole.
%   [st dip] = pole2plane(tr,pl) returns strike (st) and dip 
%   of the plane given the trend(tr) and plunge(pl) of a pole
%    
%          Input and output angles should be in radians

% ######################################################################

%Calculate plane given its pole
if pl >= 0
     dip = (pi/2) - pl;
     az_dip = tr - pi;
    else
     dip = (pi/2) + pl;
     az_dip = tr;
end

st = az_dip - (pi/2);

%  making strike in between 0 and 2*pi    
if st < 0.0
    st = st + 2*pi;
elseif st >= 2*pi
    st = st - 2*pi;
end

end