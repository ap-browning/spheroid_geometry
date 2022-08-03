%% SETUP

red = [255,59,48] / 255;
blue = [0,122,255] / 255;
orange = [255,159,10] / 255;
green = [52,199,89] / 255;
brighten = @(x,f) min(f*x,1);

clf;

%% GREENSPAN

% Load results
load("fig6abc.mat")

%% Figure 4(a) -- All
subplot(2,2,1); cla; hold on;
surf(x1,y1,z1','FaceColor',blue,'EdgeColor',brighten(blue,0.8),'FaceAlpha',0.6);
surf(x2,y2,z2','FaceColor',red,'EdgeColor',brighten(red,0.8),'FaceAlpha',0.6)

% Plot true value and constant likelihood trace
scatter3(p(1),p(2),p(3),200,'k','filled')
plot3(sol1(1,:),sol1(2,:),sol1(3,:),'Color',[0.2,0.2,0.2],'LineWidth',5)
plot3(sol2(1,:),sol2(2,:),sol2(3,:),'Color',[0.2,0.2,0.2],'LineWidth',5)

%% Figure 4(b) -- Linearisation

subplot(2,2,2); cla; hold on;
surf(x1,y1,z1','FaceColor',blue,'EdgeColor',brighten(blue,0.8),'FaceAlpha',0.2);
surf(x2,y2,z2','FaceColor',red,'EdgeColor',brighten(red,0.8),'FaceAlpha',0.2)

% Plot true value and constant likelihood trace
scatter3(p(1),p(2),p(3),200,'k','filled')

% Tangent planes
[x,y] = meshgrid(linspace(-1,1,20));
x = p(1) + 0.1*x;
y = p(2) + 50*y;
z = @(x,y,p,v) -(v(1) * (x - p(1)) + v(2) * (y - p(2))) / v(3) + p(3);

surf(x,y,z(x,y,p,J(2,:)),'FaceColor',blue,'EdgeColor','none');
surf(x,y,z(x,y,p,J(1,:)),'FaceColor',red,'EdgeColor','none');

% Constant likelihood direction
d = d(1:3) / norm(d(1:3));
x = p(1)+40*d(1)*[-1,1];
y = p(2)+40*d(2)*[-1,1];
z = p(3)+40*d(3)*[-1,1];
plot3(x,y,z,'Color',[1.0,0.2,0.2],'LineWidth',5)

% Constant K-lambda direction
n1 = J(1,1:3);
n2 = J(2,1:3);
d2 = cross(n1,n2); d2 = d2 / norm(d2);
x = p(1)+40*d(1)*[-1,1];
y = p(2)+40*d(2)*[-1,1];
z = p(3)+40*d(3)*[-1,1];
plot3(x,y,z,'Color',[0.2,0.2,1.0],'LineWidth',5)

%% Figure 4(c) -- necrotic core

subplot(2,2,3); cla; hold on;
surf(x1,y1,z1','FaceColor',blue,'EdgeColor',brighten(blue,0.8),'FaceAlpha',0.6);
surf(x2,y2,z2','FaceColor',red,'EdgeColor',brighten(red,0.8),'FaceAlpha',0.6)
surf(x3,y3,z3','FaceColor',orange,'EdgeColor',brighten(orange,0.8),'FaceAlpha',0.6)

% Plot true value and constant likelihood trace
scatter3(p(1),p(2),p(3),200,'k','filled')

%% STYLING

for i = 1:3
    subplot(2,2,i)

    % Axis limits and view      
    axis([0.2,1,50,300,0,10])
    view([35,20]);

    % Axis labels
    xlabel("Q");
    ylabel("R1");
    zlabel("gamma");

    % Plot styling
    ax = gca;
    ax.BoxStyle = 'full'; box on; grid on;
    pbaspect([1,1,1])

    % Axis labels
    str = 'abcdef';
    title(sprintf("(%s)",str(i)))
end

%% CUBE AND ZOOM...

cube_x = [0.7,1];
cube_y = [50,200];
cube_z = [0,3];

cube = [cube_x,cube_y,cube_z];

subplot(2,2,3);
axis(cube)

subplot(2,2,2)
plot3([cube_x(1),cube_x(1)],[cube_y(1),cube_y(1)],cube_z,'k')
plot3([cube_x(1),cube_x(1)],[cube_y(2),cube_y(2)],cube_z,'k')
plot3([cube_x(2),cube_x(2)],[cube_y(1),cube_y(1)],cube_z,'k')
plot3([cube_x(2),cube_x(2)],[cube_y(2),cube_y(2)],cube_z,'k')

plot3([cube_x(1),cube_x(1)],cube_y,[cube_z(1),cube_z(1)],'k')
plot3([cube_x(1),cube_x(1)],cube_y,[cube_z(2),cube_z(2)],'k')
plot3([cube_x(2),cube_x(2)],cube_y,[cube_z(1),cube_z(1)],'k')
plot3([cube_x(2),cube_x(2)],cube_y,[cube_z(2),cube_z(2)],'k')

plot3(cube_x,[cube_y(1),cube_y(1)],[cube_z(1),cube_z(1)],'k')
plot3(cube_x,[cube_y(1),cube_y(1)],[cube_z(2),cube_z(2)],'k')
plot3(cube_x,[cube_y(2),cube_y(2)],[cube_z(1),cube_z(1)],'k')
plot3(cube_x,[cube_y(2),cube_y(2)],[cube_z(2),cube_z(2)],'k')


%% RADIAL-DEATH

% Load results
load("fig6d.mat")

% x2,y2,z2 -> matrices for plotting
[Y2,Z2] = meshgrid(y2,z2); X2 = x2';

%% Figure 4(a) -- All
subplot(2,2,4); cla; hold on;
surf(x1(1:2:end),y1(1:2:end),z1(1:2:end,1:2:end)','FaceColor',blue,'EdgeColor',brighten(blue,0.8),'FaceAlpha',0.6);
surf(X2(1:2:end,1:2:end),Y2(1:2:end,1:2:end),Z2(1:2:end,1:2:end),'FaceColor',red,'EdgeColor',brighten(red,0.8),'FaceAlpha',0.6)

% Plot true value and constant likelihood trace
scatter3(p(1),p(2),p(3),200,'k','filled')
plot3(sol1(1,:),sol1(2,:),sol1(3,:),'Color',[0.2,0.2,0.2],'LineWidth',5)
plot3(sol2(1,:),sol2(2,:),sol2(3,:),'Color',[0.2,0.2,0.2],'LineWidth',5)

% Axis
axis([0.9,1.4,0,6,10,200])
view([35,20]);

% Axis labels
xlabel("lambda");
ylabel("zeta");
zlabel("R1");

% Plot styling
ax = gca;
ax.BoxStyle = 'full'; box on; grid on;
pbaspect([1,1,1])

% Title
title("(d)")

%% EXPORT FIGURE 6

% % Set size
% f = gcf;
% f.Position = [1,1,1200,1200];
% 
% % Export
% set(f,'renderer','Painters')
% saveas(f,'fig6','epsc')