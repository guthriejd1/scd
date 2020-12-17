clear
clc
close all

test_case = 'cube_ellipsoid';
test_data_folder = '../test/data/';

S02 = @(theta) [cos(theta) -sin(theta); sin(theta) cos(theta)];
switch test_case    
    case 'box_ellipse_2d'
        SQ1.a = [1 2];
        SQ1.e = [0.25];
        SQ1.R = S02(-1*pi/5);
        SQ1.t = [1;-1];
        
        E2.a = [0.5 2];
        E2.e = [1];
        E2.R = S02(-1*pi/3);
        E2.t = [-3;+3];
    case 'cube_ellipsoid'
        SQ1.a = [1 2 3];
        SQ1.e = [0.25 0.25];
        SQ1.R = axang2rotm([1 0 0 pi/3])*axang2rotm([0 1 0 2*pi/3]);
        SQ1.t = [-1;1;2];

        E2.a = [1 2 0.5];
        E2.e = [1 1];
        E2.R = axang2rotm([1 0 0 7*pi/6])*axang2rotm([0 1 0 -4*pi/6])*axang2rotm([0 0 1 3*pi/5]);
        E2.t =  [9;3;6];
    otherwise
        error(['Test ' test_case ' is not defined']);
end

switch numel(SQ1.t)
    case 2 % 2D
        ColorList = colormap(lines(10));
        figure(1)
        h_sq1 = PlotSuperquadric2D(SQ1,ColorList(1,:));
        h_e2 = PlotSuperquadric2D(E2,ColorList(2,:));

        [result] = Collide2D(SQ1, E2);

        h_e2c = PlotSuperquadric2D(result.E2_c,ColorList(3,:));

        line_eb = [SQ1.t E2.t];
        plot(line_eb(1,:), line_eb(2,:),'LineWidth',2.0,'Color','k','LineStyle','-');
        xlabel('x');
        ylabel('y');
        legend([h_sq1 h_e2 h_e2c], {'SQ_1: Superquadric','E_2: Ellipse','E_2^c: Ellipse At Contact Point'});
    case 3 % 3D
        figure(1)
        h_sq1 = PlotSuperquadric(SQ1,'b');
        h_e2 = PlotSuperquadric(E2,'r');

        [result] = Collide(SQ1, E2);

        h_e2c = PlotSuperquadric(result.E2_c,'g');

        line_eb = [SQ1.t E2.t];
        plot3(line_eb(1,:), line_eb(2,:), line_eb(3,:),'LineWidth',2.0,'Color',[0.75 0 0.25]);
        xlabel('x');
        ylabel('y');
        zlabel('z');

        grid on;
end

% Write data to test file
filename = [test_data_folder test_case '.txt'];
writematrix([], filename, 'Delimiter', ',');

Shapes{1} = SQ1;
Shapes{2} = E2;
Shapes{3} = result.E2_c;
for i = 1:numel(Shapes)
    writematrix(Shapes{i}.a, filename, 'Delimiter', ',', 'WriteMode', 'append');
    writematrix(Shapes{i}.e, filename, 'Delimiter', ',', 'WriteMode', 'append');
    writematrix(Shapes{i}.R, filename, 'Delimiter', ',', 'WriteMode', 'append');
    writematrix(Shapes{i}.t, filename, 'Delimiter', ',', 'WriteMode', 'append');
end
switch numel(SQ1.t)
    case 2
        writematrix([result.theta result.collision], filename, 'Delimiter', ',', 'WriteMode', 'append');
    case 3
        writematrix([result.omega result.eta result.collision], filename, 'Delimiter', ',', 'WriteMode', 'append');
end


