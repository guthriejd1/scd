function [h_plot pts] = PlotSuperquadric2D(S,Color)
pt_index = 0;

for theta = linspace(-pi,pi,1e3)
    pt_index = pt_index + 1;
    [theta_mod c_theta_sign s_theta_sign] = AdjustAngle(theta);

    pts(1,pt_index) = [S.a(1)*c_theta_sign*(cos(theta_mod)).^S.e(1)];
    pts(2,pt_index) = [S.a(2)*s_theta_sign*(sin(theta_mod)).^S.e(1)];
end
% Rotate and translate
pts = S.R*pts + repmat(S.t,1,size(pts,2));
h_plot = fill(pts(1,:),pts(2,:),Color);
hold on;
end

function [omega_mod cos_sign sin_sign] = AdjustAngle(omega)
    if omega < 0
        omega = omega+2*pi;
    end
    if (omega >= 0) && ( omega <= pi/2)
        omega_mod = omega;
        cos_sign = 1;
        sin_sign = 1;
    elseif (omega >= pi/2) && ( omega <= pi)
        omega_mod = pi - omega;
        cos_sign = -1;
        sin_sign = 1;
    elseif (omega >= pi) && ( omega <= 3*pi/2)
        omega_mod = omega - pi;
        cos_sign = -1;
        sin_sign = -1;
    elseif (omega >= 3*pi/2) && ( omega <= 2*pi)
        omega_mod = 2*pi - omega;
        cos_sign = 1;
        sin_sign = -1;
    else
        error(['Angle : ' num2str(180/pi*omega)]);
    end
end