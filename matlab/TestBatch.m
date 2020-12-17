clear
clc
close all

test_case = 'ellipsoid_ellipsoid_batch';
test_data_folder = '../test/data/';

make_plot = false;
S02 = @(theta) [cos(theta) -sin(theta); sin(theta) cos(theta)];

switch test_case    
    case 'ellipsoid_ellipsoid_batch'
        SQ1.a = [1 2 3];
        SQ1.e = [1 1];
        SQ1.R = axang2rotm([1 0 0 pi/3])*axang2rotm([0 1 0 2*pi/3]);
        SQ1.t = [-1;1;2];

        E2.a = [1 2 0.5];
        E2.e = [1 1];
        E2.R = axang2rotm([1 0 0 7*pi/6])*axang2rotm([0 1 0 -4*pi/6])*axang2rotm([0 0 1 3*pi/5]);
        E2.t =  [9;3;6];
        i = 0;
        for px = linspace(-5,5,10)
            for py = linspace(-5,5,10)
                for pz = linspace(-5,5,10)
                    i = i+1;
                    test_data{i}.SQ1 = SQ1;
                    test_data{i}.E2 = E2;
                    test_data{i}.E2.t = [px;py;pz];
                end
            end
        end

    otherwise
        error(['Test ' test_case ' is not defined']);
end

for ntest = 1:numel(test_data)
    close all
    SQ1 = test_data{ntest}.SQ1;
    E2 = test_data{ntest}.E2;
    
switch numel(SQ1.t)
    case 2 % 2D
        PlotSuperquadric2D(SQ1,'b');
        PlotSuperquadric2D(E2,'r');

        [result] = Collide2D(SQ1, E2)

        PlotSuperquadric2D(result.E2_c,'g');

        line_eb = [SQ1.t E2.t];
        plot(line_eb(1,:), line_eb(2,:),'LineWidth',2.0,'Color',[0.75 0 0.25]);
        xlabel('x');
        ylabel('y');
    case 3 % 3D
        if make_plot
            PlotSuperquadric(SQ1,'b');
            PlotSuperquadric(E2,'r');
        end

        [result] = Collide(SQ1, E2)

        if make_plot
            PlotSuperquadric(result.E2_c,'g');

            line_eb = [SQ1.t E2.t];
            plot3(line_eb(1,:), line_eb(2,:), line_eb(3,:),'LineWidth',2.0,'Color',[0.75 0 0.25]);
            xlabel('x');
            ylabel('y');
            zlabel('z');

            grid on;
            view(-15,15)
        end
end

% Write data to test file
filename = [test_data_folder test_case '.txt'];

if ntest == 1
    writematrix([], filename, 'Delimiter', ',');
end

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

end

