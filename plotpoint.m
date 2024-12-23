function [] = plotpoint(pl,tr,varargin)

%   
%    Simple function to plot a line on the stereonet
%    [] = plotpoint(pl,tr,varargin) plots the poit of a line given the
%    plunge(pl) and trend (tr) of the line. 'varargin' is the optional 
%    input to assign the color ( 'k' == black, 'b' == blue ) of the point.
%    
%    Stirke and dip must be in degree

plotaxis  % plot streonet axis
theta = 90 - tr;

R = sqrt(2)*sind(45 - (pl/2)); % equal area
x= R*cosd(theta);
y= R*sind(theta);

if isempty(varargin)
    scatter(x,y,35,'red','filled','o','MarkerEdgeColor','k');
else
    c = varargin{1};
    scatter(x,y,35,c,'filled','o','MarkerEdgeColor','k');
end
hold on
end