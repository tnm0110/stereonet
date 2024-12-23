function [] = plotaxis()

% simple script to plot stereonet

X=-1:.01:1;
Y1=sqrt(1-(X.^2));
Y2=-sqrt(1-(X.^2));
X2 = [X flip(X)]; Y = [ Y1 flip(Y2)];
ax_mark=[0 0; 0 1; 0 -1 ; -1 0; 1 0];
c=[0 0.4470 0.7410];
plot(X2,Y,'Color',c,'LineWidth',1.5);
%plot(X,Y1,X,Y2,'b','LineWidth',1.5);
hold on
scatter(ax_mark(:,1),ax_mark(:,2),50,'+k','LineWidth',1.5);
hold on
% plot annotation
xt = [-0.03 -0.05  1.1 -1.2 ];
yt = [1.1 -1.1 0  0];
str = {'0','180','90','270'};
text(xt,yt,str,'FontSize',12);
axis equal
axis off
box on
set(gca,'XTick',[], 'YTick', [])
hold on