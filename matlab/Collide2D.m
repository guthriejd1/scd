function [result] = Collide2D(SQ1_W, E2_W, theta)

switch nargin
    case 2
        % Solve root-finding problem to find theta
        opti = casadi.Opti();
        theta = opti.variable(1,1);   
    case 3
        % Return residual function f(theta) for a user-supplied theta
end

% Pose of SQ1, E2 in world frame (W)
X_W_SQ1 = [SQ1_W.R SQ1_W.t; 0 0 1];
X_W_E2  = [E2_W.R  E2_W.t;  0 0 1];

% Pose of E2 in SQ1 frame
X_SQ1_E2 = (X_W_SQ1)^-1*X_W_E2;

SQ1.a = SQ1_W.a;
SQ1.e = SQ1_W.e;
SQ1.R = eye(2);
SQ1.t = zeros(2,1);

E2.a = E2_W.a;
E2.e = E2_W.e;
E2.R = X_SQ1_E2(1:2,1:2); 
E2.t = X_SQ1_E2(1:2,3);


[theta_mod c_theta_sign s_theta_sign] =  ModifyAngleQuadrantCasadi(theta);

r = min(E2.a);
T = E2.R * (diag(r./E2.a)) * E2.R';

x_ax = [SQ1.a(1)*c_theta_sign*(cos(theta_mod)).^SQ1.e(1)];
x_ay = [SQ1.a(2)*s_theta_sign*(sin(theta_mod)).^SQ1.e(1)];
x_a = [x_ax; x_ay];

grad_Phi = 2/SQ1.e(1)*[1/SQ1.a(1)*c_theta_sign*(cos(theta_mod)).^(2-SQ1.e(1));...
                       1/SQ1.a(2)*s_theta_sign*(sin(theta_mod)).^(2-SQ1.e(1))];
        
x_eb = x_a + r*T^(-2)*grad_Phi/(norm(T^(-1)*grad_Phi));

f = E2.t(2)*x_eb(1) - E2.t(1)*x_eb(2);

switch nargin
    case 3
        result{1} = f;
        result{2} = x_eb;
        result{3} = E2.t;
        return;
    case 2
        opti.set_initial(theta, atan2(E2.t(2),E2.t(1)));
        opti.minimize(f'*f);

        p_opts = struct('print_time',5);
        s_opts = struct('print_level',5);
        opti.solver('ipopt',p_opts,s_opts);
        sol = opti.solve();

        result.root_residuals = sol.value(f);
        result.theta = sol.value(theta);

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
        X_SQ1_E2c = [E2.R E2.t; 0 0 1];
        X_W_E2c = X_W_SQ1*X_SQ1_E2c;

        E2_c = E2;
        E2_c.R = X_W_E2c(1:2,1:2);
        E2_c.t = X_W_E2c(1:2,3);

        result.E2_c = E2_c;
        result.x_eb = X_W_E2c(1:2,3);
end

end

