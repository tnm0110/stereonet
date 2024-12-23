function [] = plotcircle(dip,strike,varargin)

%   
%    Simple function to plot a plane on the stereonet
%    [] = plotcircle(dip,strike,varargin) plots a great circle containing
%    the plane given stirke and dip. 'varargin' is the optional input to 
%    assign the color ( 'k' == black, 'b' == blue ) of the great circle.
%    
%    Stirke and dip must be in degree


% handling Mathemical error for dip 90
if dip == 90
    dip = 89.99; 
end
s = 0 : 1: 180;
theta = -90 - strike:1:90 -strike;
%theta = 90 - s;
%theta = theta.*(pi/180);
% theta2 = 90*(pi/180) - linspace(0,pi,181); 
pl = atand(sind(s).*tand(dip)); 
R = (2/sqrt(2))*sind(45 - (pl/2)); % equal area

% back to the carteasian domain 
x= R.*cosd(theta);
y= R.*sind(theta);
plotaxis      % draw the stereonet

% draw the plane
    
if isempty(varargin)
    plot(x,y,'r','LineWidth',1.5);
else
    c = varargin{1};
    plot(x,y,c,'LineWidth',1.5);

end