clear
clc
close all

test_case = 'CubeEllipsoid';

switch test_case    
    case 'CubeEllipsoid'
        SQ1.a = [1 2 3];
        SQ1.e = [0.25 0.25];
        SQ1.R = axang2rotm([1 0 0 pi/3])*axang2rotm([0 1 0 2*pi/3]);
        SQ1.t = [-1;1;2];

        E2.a = [1 2 0.5];
        E2.e = [1 1];
        E2.R = axang2rotm([1 0 0 7*pi/6])*axang2rotm([0 1 0 -4*pi/6])*axang2rotm([0 0 1 3*pi/5]);
        E2.t =  [9;3;6];
end

figure(1)
PlotSuperquadric(SQ1,'b');
PlotSuperquadric(E2,'r');

[result] = Collide(SQ1, E2)

PlotSuperquadric(result.E2_c,'g');

line_eb = [SQ1.t E2.t];
plot3(line_eb(1,:), line_eb(2,:), line_eb(3,:),'LineWidth',2.0,'Color',[0.75 0 0.25]);
xlabel('x');
ylabel('y');
zlabel('z');

grid on;
view(-15,15)
%axis(10*[-1 1 -1 1 -1 1]);
