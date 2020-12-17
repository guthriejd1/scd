function [result] = Collide(SQ1_W, E2_W)

% Pose of SQ1, E2 in world frame (W)
X_W_SQ1 = [SQ1_W.R SQ1_W.t; 0 0 0 1];
X_W_E2  = [E2_W.R  E2_W.t;  0 0 0 1];

% Pose of E2 in SQ1 frame
X_SQ1_E2 = (X_W_SQ1)^-1*X_W_E2;

SQ1.a = SQ1_W.a;
SQ1.e = SQ1_W.e;
SQ1.R = eye(3);
SQ1.t = zeros(3,1);

E2.a = E2_W.a;
E2.e = E2_W.e;
E2.R = X_SQ1_E2(1:3,1:3); 
E2.t = X_SQ1_E2(1:3,4);

opti = casadi.Opti();
eta = opti.variable(1,1);
omega = opti.variable(1,1);

[eta_mod c_eta_sign s_eta_sign] =  ModifyAngleQuadrantCasadi(eta);
[omega_mod c_omega_sign s_omega_sign] =  ModifyAngleQuadrantCasadi(omega);

opti.set_initial(eta, sign(E2.t(3))*45*pi/180);
opti.set_initial(omega, atan2(E2.t(2),E2.t(1)) );

r = min(E2.a);
T = E2.R * (diag(r./E2.a)) * E2.R';

x_ax = [SQ1.a(1)*c_eta_sign*c_omega_sign*(cos(eta_mod)).^SQ1.e(1)*(cos(omega_mod)).^SQ1.e(2)];
x_ay = [SQ1.a(2)*c_eta_sign*s_omega_sign*(cos(eta_mod)).^SQ1.e(1)*(sin(omega_mod)).^SQ1.e(2)];
x_az = [SQ1.a(3)*s_eta_sign*(sin(eta_mod)).^SQ1.e(1)]; 
x_a = [x_ax; x_ay; x_az];

grad_Phi = [1/SQ1.a(1)*c_eta_sign*c_omega_sign*(cos(eta_mod)).^(2-SQ1.e(1)) .* (cos(omega_mod)).^(2-SQ1.e(2));...
            1/SQ1.a(2)*c_eta_sign*s_omega_sign*(cos(eta_mod)).^(2-SQ1.e(1)) .* (sin(omega_mod)).^(2-SQ1.e(2));...
            1/SQ1.a(3)*s_eta_sign*(sin(eta_mod)).^(2-SQ1.e(1))];
        
x_eb = x_a + r*T^(-2)*grad_Phi/(norm(T^(-1)*grad_Phi));

f = cross(x_eb, E2.t);

opti.minimize(f'*f);

p_opts = struct('print_time',5);
s_opts = struct('print_level',5);
opti.solver('ipopt',p_opts,s_opts);
sol = opti.solve();

result.root_residuals = sol.value(f);
result.omega = sol.value(omega);
result.eta = sol.value(eta);

assert(norm(result.root_residuals) <= 1e-6, 'Tolerance not met');

% Find point that results in surface contact
x_c = sol.value(x_eb);
scale = norm(x_c)/norm(E2.t);
x_c = scale*E2.t;

if scale < 1
    result.collision = false;
else
    result.collision = true;
end

E2.t = x_c;
% Transform back to world frame
X_SQ1_E2c = [E2.R E2.t; 0 0 0 1];
X_W_E2c = X_W_SQ1*X_SQ1_E2c;

E2_c = E2;
E2_c.R = X_W_E2c(1:3,1:3);
E2_c.t = X_W_E2c(1:3,4);

result.E2_c = E2_c;
result.x_eb = X_W_E2c(1:3,4);


end

