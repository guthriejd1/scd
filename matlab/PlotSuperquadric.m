function [h_plot pts] = PlotSuperquadric(S,Color)
pt_index = 0;

eta_pts = [0 logspace(-10,log10(pi/2),200)];
omega_pts = [0 logspace(-10,log10(pi),200)];
for eta = unique([-eta_pts eta_pts])
    for omega = unique([-omega_pts omega_pts])
        pt_index = pt_index + 1;
        [omega_mod c_omega_sign s_omega_sign] = AdjustAngle(omega);
        [eta_mod c_eta_sign s_eta_sign] = AdjustAngle(eta);

        pts(1,pt_index) = [S.a(1)*c_eta_sign*c_omega_sign*(cos(eta_mod)).^S.e(1)*(cos(omega_mod)).^S.e(2)];
        pts(2,pt_index) = [S.a(2)*c_eta_sign*s_omega_sign*(cos(eta_mod)).^S.e(1)*(sin(omega_mod)).^S.e(2)];
        pts(3,pt_index) = [S.a(3)*s_eta_sign*(sin(eta_mod)).^S.e(1)] + omega_mod*0; 
    end
end
% Rotate and translate
pts = S.R*pts + repmat(S.t,1,size(pts,2));
plot3(pts(1,:),pts(2,:),pts(3,:),'k');
hold on;

x = pts(1,:); y = pts(2,:); z = pts(3,:);

tri = delaunay(x,y);
h = trisurf(tri, x, y, z);
hold on;
view(30, 30)
% 
h.FaceColor = Color;
h.EdgeColor = 'none';
daspect([1 1 1])
view(3); 
axis tight
camlight 
lighting phong;

shading interp
light
lighting gouraud
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