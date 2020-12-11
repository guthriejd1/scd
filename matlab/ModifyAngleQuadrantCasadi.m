function [theta_mod cos_sign sin_sign] = ModifyAngleQuadrantCasadi(theta)
    theta = if_else(theta < 0, theta + 2*pi, theta);
    theta = mod(theta,2*pi);
    theta_mod = if_else(theta < pi/2, theta, if_else(theta <= pi, pi - theta, if_else(theta <= 3*pi/2, theta - pi, 2*pi - theta) ) ); 
    cos_sign = if_else(theta <= pi/2, 1, if_else(theta <= pi, -1, if_else(theta <= 3*pi/2, -1, 1) ) ); 
    sin_sign = if_else(theta <= pi/2, 1, if_else(theta <= pi,  1, if_else(theta <= 3*pi/2, -1, -1) ) );
end

