clear;
close all

%% Exercise 1
% Plotting planes of the joint sets
% Plot planes and poles of the following measurements of a joint set in an
% outcrop 

joints = [52,145;55,130;54,137;53,133;49,131]; % joint sets

figure(1)
for i=1:size(joints,1)
plotcircle(joints(i,1),joints(i,2)); % plot planes using funcion "plotcircle"
end

% plotting the poles of the joint sets
poles = zeros(size(joints));
for i=1:size(joints,1)
tr = deg2rad(joints(i,2));          % converting trend into radians
pl = deg2rad(joints(i,1));         
[tr_p, pl_p]=plane2pole(tr,pl); %calculate pole using "plane2pole" function
plotpoint(rad2deg(pl_p),rad2deg(tr_p))    % plotting the poles
poles(i,:) = [tr_p pl_p];                 % saving all poles
end

% %  calculating the mean vector from the joint set using direction cosinses 

% Calulate the direction cosine of the poles

dir_co=zeros(size(poles,1),3);
for i=1:size(poles,1)   
    pl = poles(i,2);
    tr = poles(i,1);        
    cd = sin(pl);
    ce = cos(pl) * sin(tr);
    cn = cos(pl) * cos(tr); 
    dir_co(i,:) = [cn ce cd]; % saving all direction cosines 
end

dir_co_m = mean(dir_co);     % mean vector of dc

% convert back to spherical to get tr and pl using "ned2sph" function
[tr_m, pl_m] = ned2sph(dir_co_m(1),dir_co_m(2),dir_co_m(3));
d=['Mean vector, Trend: ' num2str(rad2deg(tr_m)),' ; Plunge: ' ...
    num2str(rad2deg(pl_m))];
disp(d)

plotpoint(rad2deg(pl_m),rad2deg(tr_m),'y');
disp('Yellow point is the mean veotor of the pole of the joint set planes.')

%% Exercise 2
% You have two joint sets with orientations:
% 45°, 60
% 330°, 45
% Plot the planes (great circles) and the poles to the two planes
% Estimate the orientation of the line of intersection between
% the planes (trend & plunge).
%
% Using the basic vector operation, calculate the trend and
% plunge of the line of intersection of the planes (this is orthogonal to the poles of both
% planes, so calculate the vector cross product of the planes’ poles).

joints2 = [60, 45; 45, 330];
% plotting planes of the joint sets
for i=1:size(joints2,1)
    figure(2);
    plotcircle(joints2(i,1),joints2(i,2));
end

% plotting the poles of the joint sets
poles2 = zeros(size(joints2));
for i=1:size(joints2,1)
tr = deg2rad(joints2(i,2));          % converting trend into radians
pl = deg2rad(joints2(i,1));          % converting plunge into radians
[tr_p, pl_p] = plane2pole(tr,pl);         % calculate the poles
plotpoint(rad2deg(pl_p),rad2deg(tr_p))    % plotting the poles
poles2(i,:) = [tr_p pl_p];                % saving all poles
end

% Calulate the direction cosine of the poles
dir_co2=zeros(size(poles2,1),3);          
for i=1:size(poles2,1)   
    pl = poles2(i,2);
    tr = poles2(i,1);

    cd = sin(pl);
    ce = cos(pl) * sin(tr);
    cn = cos(pl) * cos(tr); 
    dir_co2(i,:) = [cn ce cd];  % saving the two direction cosines 
end

p1 = dir_co2(1,:); p3 = dir_co2(2,:);

% corss porduct to get the line of intersection of the planes
c=cross(p3,p1,2); 

% project the vector down to the lower hemisphere  
if c(3) < 0
    c = c.*(-1);
end
r = sqrt(c(1).^2 + c(2).^2 + c(3).^2);        % get unit vector 

% convert back to spherical to get trend and plunge of the line
[tr_l, pl_l] = ned2sph(c(1)/r,c(2)/r,c(3)/r);
plotpoint(rad2deg(pl_l),rad2deg(tr_l),'b')

fprintf(['Orientation of the intersection between the planes\n Trend :' ...
    '%03.3f', ' ' ...
    'Plunge : %02.3f\n'], rad2deg(tr_l), rad2deg(pl_l));

%% Exercise 3 



% For the following fault plane:
%  280°, 40
% The striation on the fault plane has a rake angle of 80° from the strike direction and the upper
% plate moves downward.
% On a stereonet, plot the fault plane and this slip vector.
% What type of fault is this? If there is any lateral slip, is it left- or right-lateral?
% 
% If there is any dip slip, is it normal or revers?
% Add the pole to the fault plane on the stereonet.
% Add two more lines (points) on the stereonet as follows: Both of these points
% are located on the plane (great circle) that includes the pole of the fault plane and the slip
% line. One point, labeled ‘P’ is located 45° between the pole and the slip line. The other,
% labeled ‘T’ is located 135° from the slip line (and 45° from the pole of the fault plane).




rake = 80;
strike_f = 280; dip_f = 40;

% (a)  Fault plane and slip vector
figure(3)
plotcircle(dip_f,strike_f,'r');

% Note: rotation attribute may not work for earlier Matlab versions
text(-0.40,0.62,'Fault','FontSize',10,'FontWeight','bold',Rotation=9);
text(0.15,0.63,'Plane','FontSize',10,'FontWeight','bold',Rotation=-15);
% calculate the rake attitude 

%direction cosine of a line having strike with zero plunge
[V1(1), V1(2), V1(3)] = sph2ned(deg2rad(strike_f),0);
%direction cosine of a line with true dip direcion and dip angle
dip_dir = strike_f-90;
[V2(1), V2(2), V2(3)] = sph2ned(deg2rad(10),deg2rad(dip_f));

% calulate the trend and plunge of the straitions line
% Note : this require symoblic math toolbox to solve the equation
syms a b c
eqn1 = a^2 + b^2 + c^2 == 1;   % a b c is the direction cosine of the line 
eqn2 = V1(1)*a + V1(2)*b + V1(3)*c == cosd(rake);
eqn3 = V2(1)*a + V2(2)*b + V2(3)*c == cosd(rake-90);

sol = solve([eqn1, eqn2, eqn3], [a, b, c]);  % solve the equation 
rake_dir = [ double(vpa(sol.a,6)), double(vpa(sol.b,6)), ...
    double(vpa(sol.c,6)) ];

rake_dir = real(rake_dir(1,:)); % two solution will be the same 

% rotate back to the spherical 
[tr_rake, pl_rake] = ned2sph(rake_dir(1),rake_dir(2),rake_dir(3));

% plot the straition line
plotpoint(rad2deg(pl_rake),rad2deg(tr_rake),'b');
text(0,0.5,'R/P2','FontSize',11,'FontWeight','bold');

% (c) fault type

disp(['Its primary a dip-slip normal fault with a little left latral' ...
    ' strike-slip component']);
% (c) pole of the fault plane
[p1(1), p1(2), p1(3)] = sph2ned(deg2rad(strike_f),deg2rad(dip_f),'pole');
[tr_p1, pl_p1] = ned2sph(p1(1),p1(2),p1(3));   % convert back to spherical 
plotpoint(rad2deg(pl_p1),rad2deg(tr_p1),'b');  % plot the pole
text(-0.25,-0.55,'P1','FontSize',11,'FontWeight','bold');

% (d) calculate the equitarial plane that contains P and T axis 

p3 = cross(p1, rake_dir);                  % pole of the ep 
r = sqrt(p3(1).^2 + p3(2).^2 + p3(3).^2);  % get unit vector
[tr_p3, pl_p3] = ned2sph(p3(1)/r,p3(2)/r,p3(3)/r); % convert back to spherical
plotpoint(rad2deg(pl_p3),rad2deg(tr_p3),'k');      % pole of the EP
text(0.80,0.0,'P3','FontSize',11,'FontWeight','bold');

% get the strke and dip of ep from from the pole
[st_ep, dip_ep] = pole2plane(tr_p3,pl_p3);
plotcircle(rad2deg(dip_ep),rad2deg(st_ep),'g--') % plot the plane EP 


% plot the auxiliary plane
% get the strke and dip of the AP from from the pole
[st_ap, dip_ap] = pole2plane(tr_rake,pl_rake);
plotcircle(rad2deg(dip_ap),rad2deg(st_ap),'r') % plot AP 


% find the P axis (creates 45 between pole and slip)
% it makes 45 from the slip vector and (90-45) from the P1

% solve the equation to get it the P axis
% Note : this require symoblic math toolbox to solve the equation
P_ang = 45;
syms a1 b1 c1                 % a1,b1,c1 => DC of the P axis line to find
eqn1 = a1^2 + b1^2 + c1^2 == 1;  
eqn2 = rake_dir(1)*a1 + rake_dir(2)*b1 + rake_dir(3)*c1 == cosd(P_ang);
eqn3 = p1(1)*a1 + p1(2)*b1 + p1(3)*c1 == cosd(90 - P_ang);

sol = solve([eqn1, eqn2, eqn3], [a1, b1, c1]);  % solve the equation 
P_axis_dc = [ double(vpa(sol.a1,6)), double(vpa(sol.b1,6)), ...
    double(vpa(sol.c1,6)) ];
P_axis_dc = real(P_axis_dc(1,:));    % two solution will be the same

% rotate back to the spherical 
[P_axis_tr, P_axis_pl] = ned2sph(P_axis_dc(1),P_axis_dc(2),P_axis_dc(3));
plotpoint(rad2deg(P_axis_pl),rad2deg(P_axis_tr),'y');
text(-0.2,0.1,'P','FontSize',11,'FontWeight','bold');

% Similarly find T axis (creates 135 between pole and slip)
T_ang = 135;
syms a2 b2 c2                 % a2,b2,c2 => DC of the T axis line to find
eqn1 = a2^2 + b2^2 + c2^2 == 1;    
eqn2 = rake_dir(1)*a2 + rake_dir(2)*b2 + rake_dir(3)*c2 == cosd(T_ang);
eqn3 = p1(1)*a2 + p1(2)*b2 + p1(3)*c2 == cosd(90 - T_ang);

sol = solve([eqn1, eqn2, eqn3], [a2, b2, c2]);  % solve the equation 
T_axis_dc = [ double(vpa(sol.a2,6)), double(vpa(sol.b2,6)), ...
    double(vpa(sol.c2,6)) ];
T_axis_dc = real(T_axis_dc(1,:)); % two solution will be the same

% rotate back to the spherical 
[T_axis_tr, T_axis_pl] = ned2sph(T_axis_dc(1),T_axis_dc(2),T_axis_dc(3));
plotpoint(rad2deg(T_axis_pl),rad2deg(T_axis_tr),'o');
text(0,-0.9,'T','FontSize',11,'FontWeight','bold');


%% Exercise 4

% The apparent dip of a planar feature (e.g., bedding, fault plane, joint) in an outcrop
% (including a road cut) is the plunge of the line defined by the intersection of the two planes (i.e.,
% the outcrop/road cut face and the bedding/fault plane/joint). The apparent dip value is the plunge
% of that line.


% On your stereonet, plot a plane representing a road cut that strikes 30° (NE)
% and dips 75° (Right-Hand Rule). Label this “road cut”.

% plot a plane representing the bedding with orientation 240°, 60°.



% Calculate the apparent dip. 

% Calculate the rake (i.e., the angle between the strike line and
% the slip direction you just calculated). This can be calculated from the vector dot
% product of the strike vector (a line with trend = strike, plunge = 0) and the trend &
% plunge of the slip vector (i.e., the apparent dip).



st_rc = 30; dip_rc = 75;
st_bed = 240; dip_bed = 60;

%  [a b] plot the planes
figure(4)
plotcircle(dip_rc,st_rc,'r');
text(0.13,-0.3,'Road cut','FontSize',11,'FontWeight','bold',Rotation=60);
plotcircle(dip_bed,st_bed,'k');
text(-0.1,0.45,'Bedding','FontSize',11,'FontWeight','bold',Rotation=20);

% calculation of apparent dip

%pole of the road cut
[P_rc(1), P_rc(2)] = plane2pole(deg2rad(st_rc),deg2rad(dip_rc));
plotpoint(rad2deg(P_rc(2)), rad2deg(P_rc(1)),'r');

%pole of the bedding
[P_bed(1), P_bed(2)] = plane2pole(deg2rad(st_bed),deg2rad(dip_bed));
plotpoint(rad2deg(P_bed(2)), rad2deg(P_bed(1)),'k');


% convert poles spherical to NED system
[P_rc_dc(1),P_rc_dc(2),P_rc_dc(3)] = sph2ned(P_rc(1),P_rc(2));
[P_bed_dc(1),P_bed_dc(2),P_bed_dc(3)] = sph2ned(P_bed(1),P_bed(2));

app_dc = cross(P_rc_dc,P_bed_dc);   % cross product to find apparent dip

% make sure vector pointing downward
if app_dc(3) < 0
    app_dc = app_dc.*(-1);
end

[app_plane_tr, app_plane_pl]  = ned2sph(app_dc(1),app_dc(2),app_dc(3));

fprintf(['Apparent dip: ' ...
    '%03.3f\n'], rad2deg(app_plane_pl))

%plotpoint(rad2deg(app_plane_pl),rad2deg(app_plane_tr), 'g')

% calculte rake
[strike_dc(1), strike_dc(2), strike_dc(3)] = sph2ned(st_bed,0);
rake = dot(app_dc,strike_dc);
fprintf('Rake: %03.3f\n',rad2deg(rake))
