clc
clear all
close all

%% ===== INITIAL CONDITIONS =====
u0          = 170;          % forward velocity (m/s)  [Final Sizing: Cruise velocity phase 1]
v0          = 0;                % lateral velocity
alpha_trim  = 0.0857;           % trim AoA (rad)
w0          = u0*tan(alpha_trim);% vertical velocity
p0          = 0;                % roll rate
q0          = 0;                % pitch rate
r0          = 0;                % yaw rate
phi0        = 0;                % roll angle
theta0      = 0.0857;           % pitch angle = AoA at trim
psi0        = 0;                % yaw angle
x0          = 0;                % x position
y0          = 0;                % y position
z0          = 15000;            % altitude (ft)

initial_conditions = [u0 v0 w0 p0 q0 r0 phi0 theta0 psi0 x0 y0 z0];

%% ===== OPEN-LOOP ODE SOLVER =====
tspan   = [0 1000];
options = odeset('RelTol',1e-6,'AbsTol',1e-8);
[t_ol, states_ol] = ode45(@(t,y) equations_of_motion(t, y, [0;0;0]), ...
                           tspan, initial_conditions, options);

%% ===== LINEARISATION AT TRIM =====
dh       = 1e-4;
x_trim   = initial_conditions';
u_trim   = [0;0;0];   % [del_e; del_a; del_r] at trim
f0       = equations_of_motion(0, x_trim, u_trim);

% Full 12x12 A matrix (state Jacobian)
A_full = zeros(12,12);
for i = 1:12
    xp       = x_trim;
    xp(i)    = xp(i) + dh;
    fp       = equations_of_motion(0, xp, u_trim);
    A_full(:,i) = (fp - f0) / dh;
end

% Full 12x3 B matrix (control Jacobian) [del_e, del_a, del_r]
B_full = zeros(12,3);
for j = 1:3
    up    = u_trim;
    up(j) = up(j) + dh;
    fp    = equations_of_motion(0, x_trim, up);
    B_full(:,j) = (fp - f0) / dh;
end


%% ===== LONGITUDINAL SUBSYSTEM  (u, w, q, theta) -> indices 1,3,5,8 =====
idx_lon = [1,3,5,8];
A_lon   = A_full(idx_lon, idx_lon);
B_lon   = B_full(idx_lon, 1);       % elevator only

% Longitudinal LQR
%   State weights: u(m/s), w(m/s), q(rad/s), theta(rad)
%   Tune by increasing diagonal entries to tighten that channel
Q_lon = diag([0.01,  ...  % u  - loose (trim speed maintained by thrust)
              0.1,   ...  % w  - moderate
              10.0,  ...  % q  - tight  (damp short period)
              50.0]);     % theta - tight (attitude hold)
R_lon = 0.5;              % elevator authority cost

[K_lon, ~, eig_lon_cl] = lqr(A_lon, B_lon, Q_lon, R_lon);

fprintf('========================================\n');
fprintf('    LONGITUDINAL LQR GAIN  K_lon\n');
fprintf('========================================\n');
fprintf('  K = [Ku=%.4f  Kw=%.4f  Kq=%.4f  Ktheta=%.4f]\n', K_lon);

%% ===== LATERAL-DIRECTIONAL SUBSYSTEM  (v, p, r, phi) -> indices 2,4,6,7 =====
idx_lat = [2,4,6,7];
A_lat   = A_full(idx_lat, idx_lat);
B_lat   = B_full(idx_lat, [2,3]);   % aileron + rudder

% Lateral LQR
%   State weights: v(m/s), p(rad/s), r(rad/s), phi(rad)
Q_lat = diag([0.1,   ...  % v   - sideslip damping
              5.0,   ...  % p   - roll rate
              10.0,  ...  % r   - yaw damping
              80.0]);     % phi - bank angle hold
R_lat = diag([0.5, 1.0]); % aileron, rudder costs

[K_lat, ~, eig_lat_cl] = lqr(A_lat, B_lat, Q_lat, R_lat);

fprintf('\n========================================\n');
fprintf('  LATERAL-DIRECTIONAL LQR GAIN  K_lat\n');
fprintf('========================================\n');
fprintf('  Aileron  row: [Kv=%.4f  Kp=%.4f  Kr=%.4f  Kphi=%.4f]\n', K_lat(1,:));
fprintf('  Rudder   row: [Kv=%.4f  Kp=%.4f  Kr=%.4f  Kphi=%.4f]\n', K_lat(2,:));

%% ===== OPEN-LOOP EIGENVALUES =====
[~, D_lon_ol] = eig(A_lon);  eig_lon_ol = diag(D_lon_ol);
[~, D_lat_ol] = eig(A_lat);  eig_lat_ol = diag(D_lat_ol);

fprintf('\n========================================\n');
fprintf('   OPEN-LOOP LONGITUDINAL MODES\n');
fprintf('========================================\n');
for i = 1:4
    lam  = eig_lon_ol(i);
    wn   = abs(lam);
    zeta = -real(lam)/max(wn,1e-10);
    fprintf('  lambda = %+.4f %+.4fi  |  wn=%.4f rad/s  |  zeta=%.4f\n', ...
            real(lam), imag(lam), wn, zeta);
end
fprintf('\n   OPEN-LOOP LATERAL-DIR MODES\n');
fprintf('========================================\n');
for i = 1:4
    lam  = eig_lat_ol(i);
    wn   = abs(lam);
    zeta = -real(lam)/max(wn,1e-10);
    fprintf('  lambda = %+.4f %+.4fi  |  wn=%.4f rad/s  |  zeta=%.4f\n', ...
            real(lam), imag(lam), wn, zeta);
end

fprintf('\n========================================\n');
fprintf('   CLOSED-LOOP LONGITUDINAL MODES\n');
fprintf('========================================\n');
for i = 1:4
    lam  = eig_lon_cl(i);
    wn   = abs(lam);
    zeta = -real(lam)/max(wn,1e-10);
    fprintf('  lambda = %+.4f %+.4fi  |  wn=%.4f rad/s  |  zeta=%.4f\n', ...
            real(lam), imag(lam), wn, zeta);
end
fprintf('\n   CLOSED-LOOP LATERAL-DIR MODES\n');
fprintf('========================================\n');
for i = 1:4
    lam  = eig_lat_cl(i);
    wn   = abs(lam);
    zeta = -real(lam)/max(wn,1e-10);
    fprintf('  lambda = %+.4f %+.4fi  |  wn=%.4f rad/s  |  zeta=%.4f\n', ...
            real(lam), imag(lam), wn, zeta);
end
fprintf('========================================\n');

%% ===== CLOSED-LOOP SIMULATION (LQR ACTIVE) =====
% Disturb trim: +5 deg pitch, +3 deg roll to see LQR rejection
ic_disturbed        = initial_conditions;
ic_disturbed(8)     = theta0 + 5*pi/180;   % theta perturbed
ic_disturbed(7)     = phi0   + 3*pi/180;   % phi perturbed

tspan_cl = [0 200];   % shorter window — LQR should stabilise quickly

[t_cl, states_cl] = ode45( ...
    @(t,y) equations_of_motion_lqr(t, y, ...
           x_trim, K_lon, K_lat, idx_lon, idx_lat), ...
    tspan_cl, ic_disturbed, options);

%% ===== UNPACK OPEN-LOOP =====
u_ol     = states_ol(:,1);   v_ol = states_ol(:,2);   w_ol = states_ol(:,3);
p_ol     = states_ol(:,4);   q_ol = states_ol(:,5);   r_ol = states_ol(:,6);
phi_ol   = states_ol(:,7);   theta_ol = states_ol(:,8); psi_ol = states_ol(:,9);
x_ol     = states_ol(:,10);  y_ol = states_ol(:,11);  z_ol  = states_ol(:,12);
V_ol     = sqrt(u_ol.^2 + v_ol.^2 + w_ol.^2);
alpha_ol = atan2(w_ol, u_ol);
beta_ol  = asin(v_ol ./ V_ol);

%% ===== UNPACK CLOSED-LOOP =====
u_cl     = states_cl(:,1);   v_cl = states_cl(:,2);   w_cl = states_cl(:,3);
p_cl     = states_cl(:,4);   q_cl = states_cl(:,5);   r_cl = states_cl(:,6);
phi_cl   = states_cl(:,7);   theta_cl = states_cl(:,8); psi_cl = states_cl(:,9);
x_cl     = states_cl(:,10);  y_cl = states_cl(:,11);  z_cl  = states_cl(:,12);
V_cl     = sqrt(u_cl.^2 + v_cl.^2 + w_cl.^2);
alpha_cl = atan2(w_cl, u_cl);

% Recompute LQR control history for plotting
del_e_hist = zeros(length(t_cl),1);
del_a_hist = zeros(length(t_cl),1);
del_r_hist = zeros(length(t_cl),1);
for k = 1:length(t_cl)
    dx_lon = states_cl(k, idx_lon)' - x_trim(idx_lon);
    dx_lat = states_cl(k, idx_lat)' - x_trim(idx_lat);
    u_lqr_lon = -K_lon * dx_lon;
    u_lqr_lat = -K_lat * dx_lat;
    del_e_hist(k) = u_lqr_lon(1);
    del_a_hist(k) = u_lqr_lat(1);
    del_r_hist(k) = u_lqr_lat(2);
end

%% ===== PLOTS =====

% ── Figure 1: Linear Velocities (OL) ─────────────────────────────────────
figure(1);
sgtitle('Linear Velocities — Open Loop','FontWeight','bold','FontSize',13);
subplot(3,1,1); plot(t_ol,u_ol,'r','LineWidth',1.5); grid on;
xlabel('Time (s)'); ylabel('u (m/s)'); title('Forward Velocity');
subplot(3,1,2); plot(t_ol,v_ol,'b','LineWidth',1.5); grid on;
xlabel('Time (s)'); ylabel('v (m/s)'); title('Lateral Velocity');
subplot(3,1,3); plot(t_ol,w_ol,'g','LineWidth',1.5); grid on;
xlabel('Time (s)'); ylabel('w (m/s)'); title('Vertical Velocity');

% ── Figure 2: Total Velocity & Angles (OL) ───────────────────────────────
figure(2);
sgtitle('Velocity & Aerodynamic Angles — Open Loop','FontWeight','bold','FontSize',13);
subplot(3,1,1); plot(t_ol,V_ol,'m','LineWidth',1.5); grid on;
xlabel('Time (s)'); ylabel('V (m/s)'); title('Total Velocity');
subplot(3,1,2); plot(t_ol,alpha_ol*180/pi,'r','LineWidth',1.5); grid on;
xlabel('Time (s)'); ylabel('\alpha (deg)'); title('Angle of Attack');
subplot(3,1,3); plot(t_ol,beta_ol*180/pi,'b','LineWidth',1.5); grid on;
xlabel('Time (s)'); ylabel('\beta (deg)'); title('Sideslip Angle');

% ── Figure 3: Angular Rates (OL) ─────────────────────────────────────────
figure(3);
sgtitle('Angular Rates — Open Loop','FontWeight','bold','FontSize',13);
subplot(3,1,1); plot(t_ol,p_ol*180/pi,'b','LineWidth',1.5); grid on;
xlabel('Time (s)'); ylabel('p (deg/s)'); title('Roll Rate');
subplot(3,1,2); plot(t_ol,q_ol*180/pi,'g','LineWidth',1.5); grid on;
xlabel('Time (s)'); ylabel('q (deg/s)'); title('Pitch Rate');
subplot(3,1,3); plot(t_ol,r_ol*180/pi,'r','LineWidth',1.5); grid on;
xlabel('Time (s)'); ylabel('r (deg/s)'); title('Yaw Rate');

% ── Figure 4: Euler Angles (OL) ──────────────────────────────────────────
figure(4);
sgtitle('Euler Angles — Open Loop','FontWeight','bold','FontSize',13);
subplot(3,1,1); plot(t_ol,phi_ol*180/pi,'b','LineWidth',1.5); grid on;
xlabel('Time (s)'); ylabel('\phi (deg)'); title('Roll Angle');
subplot(3,1,2); plot(t_ol,theta_ol*180/pi,'g','LineWidth',1.5); grid on;
xlabel('Time (s)'); ylabel('\theta (deg)'); title('Pitch Angle');
subplot(3,1,3); plot(t_ol,psi_ol*180/pi,'r','LineWidth',1.5); grid on;
xlabel('Time (s)'); ylabel('\psi (deg)'); title('Yaw Angle');

% ── Figure 5: Position (OL) ───────────────────────────────────────────────
figure(5);
sgtitle('Position — Open Loop','FontWeight','bold','FontSize',13);
subplot(3,1,1); plot(t_ol,x_ol,'r','LineWidth',1.5); grid on;
xlabel('Time (s)'); ylabel('X (m)'); title('X Position');
subplot(3,1,2); plot(t_ol,y_ol,'b','LineWidth',1.5); grid on;
xlabel('Time (s)'); ylabel('Y (m)'); title('Y Position');
subplot(3,1,3); plot(t_ol,z_ol,'g','LineWidth',1.5); grid on;
xlabel('Time (s)'); ylabel('Altitude (ft)'); title('Altitude');

% ── Figure 6: 3D Trajectory (OL) ─────────────────────────────────────────
figure(6);
plot3(x_ol,y_ol,z_ol,'b','LineWidth',1.8); grid on;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Altitude (ft)');
title('3D Flight Trajectory — Open Loop','FontWeight','bold','FontSize',13);
hold on;
plot3(x_ol(1),y_ol(1),z_ol(1),'go','MarkerSize',10,'MarkerFaceColor','g');
plot3(x_ol(end),y_ol(end),z_ol(end),'rs','MarkerSize',10,'MarkerFaceColor','r');
legend('Trajectory','Start','End','Location','best');

% ── Figure 7: Eigenvalue Plot — OL vs CL ─────────────────────────────────
figure(7);
hold on; grid on;
xline(0,'k--','LineWidth',1);
yline(0,'k--','LineWidth',1);
scatter(real(eig_lon_ol), imag(eig_lon_ol), 120, 'b', 'o','filled','DisplayName','Lon OL');
scatter(real(eig_lat_ol), imag(eig_lat_ol), 120, 'r', '^','filled','DisplayName','Lat OL');
scatter(real(eig_lon_cl), imag(eig_lon_cl), 120, 'c', 's','filled','DisplayName','Lon CL (LQR)');
scatter(real(eig_lat_cl), imag(eig_lat_cl), 120, 'm', 'd','filled','DisplayName','Lat CL (LQR)');
xlabel('Real Part'); ylabel('Imaginary Part');
title('Eigenvalue Plot — Open Loop vs Closed Loop (LQR)','FontWeight','bold','FontSize',13);
legend('Location','best');
all_eigs = [eig_lon_ol; eig_lat_ol; eig_lon_cl; eig_lat_cl];
lbls     = {'L1_{OL}','L2_{OL}','L3_{OL}','L4_{OL}',...
            'LD1_{OL}','LD2_{OL}','LD3_{OL}','LD4_{OL}',...
            'L1_{CL}','L2_{CL}','L3_{CL}','L4_{CL}',...
            'LD1_{CL}','LD2_{CL}','LD3_{CL}','LD4_{CL}'};
for i = 1:16
    text(real(all_eigs(i))+0.002, imag(all_eigs(i))+0.002, lbls{i},'FontSize',8);
end

% ── Figure 8: LQR Disturbance Rejection — Euler Angles ───────────────────
figure(8);
sgtitle('LQR Disturbance Rejection — Euler Angles','FontWeight','bold','FontSize',13);
subplot(3,1,1);
plot(t_cl, phi_cl*180/pi, 'b','LineWidth',1.5); grid on;
yline(phi0*180/pi,'k--','Trim','LabelHorizontalAlignment','left');
xlabel('Time (s)'); ylabel('\phi (deg)'); title('Roll Angle (LQR Closed Loop)');

subplot(3,1,2);
plot(t_cl, theta_cl*180/pi, 'g','LineWidth',1.5); grid on;
yline(theta0*180/pi,'k--','Trim','LabelHorizontalAlignment','left');
xlabel('Time (s)'); ylabel('\theta (deg)'); title('Pitch Angle (LQR Closed Loop)');

subplot(3,1,3);
plot(t_cl, psi_cl*180/pi, 'r','LineWidth',1.5); grid on;
xlabel('Time (s)'); ylabel('\psi (deg)'); title('Yaw Angle (LQR Closed Loop)');

% ── Figure 9: LQR Control Surface Deflections ─────────────────────────────
figure(9);
sgtitle('LQR Control Surface Deflections','FontWeight','bold','FontSize',13);
subplot(3,1,1);
plot(t_cl, del_e_hist*180/pi, 'r','LineWidth',1.5); grid on;
xlabel('Time (s)'); ylabel('\delta_e (deg)'); title('Elevator Deflection');

subplot(3,1,2);
plot(t_cl, del_a_hist*180/pi, 'b','LineWidth',1.5); grid on;
xlabel('Time (s)'); ylabel('\delta_a (deg)'); title('Aileron Deflection');

subplot(3,1,3);
plot(t_cl, del_r_hist*180/pi, 'g','LineWidth',1.5); grid on;
xlabel('Time (s)'); ylabel('\delta_r (deg)'); title('Rudder Deflection');

% ── Figure 10: OL vs CL Comparison — Pitch & Roll ─────────────────────────
figure(10);
sgtitle('Open Loop vs LQR Closed Loop Comparison','FontWeight','bold','FontSize',13);

subplot(2,2,1);
t_ol_short = t_ol(t_ol <= 200);
plot(t_ol_short, theta_ol(t_ol<=200)*180/pi, 'r--','LineWidth',1.5); hold on;
plot(t_cl, theta_cl*180/pi, 'b','LineWidth',1.5); grid on;
yline(theta0*180/pi,'k:','Trim');
legend('Open Loop','LQR CL','Trim','Location','best');
xlabel('Time (s)'); ylabel('\theta (deg)'); title('Pitch Angle');

subplot(2,2,2);
plot(t_ol_short, phi_ol(t_ol<=200)*180/pi, 'r--','LineWidth',1.5); hold on;
plot(t_cl, phi_cl*180/pi, 'b','LineWidth',1.5); grid on;
yline(phi0*180/pi,'k:','Trim');
legend('Open Loop','LQR CL','Trim','Location','best');
xlabel('Time (s)'); ylabel('\phi (deg)'); title('Roll Angle');

subplot(2,2,3);
plot(t_ol_short, q_ol(t_ol<=200)*180/pi, 'r--','LineWidth',1.5); hold on;
plot(t_cl, q_cl*180/pi, 'b','LineWidth',1.5); grid on;
legend('Open Loop','LQR CL','Location','best');
xlabel('Time (s)'); ylabel('q (deg/s)'); title('Pitch Rate');

subplot(2,2,4);
plot(t_ol_short, p_ol(t_ol<=200)*180/pi, 'r--','LineWidth',1.5); hold on;
plot(t_cl, p_cl*180/pi, 'b','LineWidth',1.5); grid on;
legend('Open Loop','LQR CL','Location','best');
xlabel('Time (s)'); ylabel('p (deg/s)'); title('Roll Rate');

% ── Figure 11: 3D Trajectory — OL vs CL ──────────────────────────────────
figure(11);
plot3(x_ol, y_ol, z_ol, 'r--','LineWidth',1.5); hold on;
plot3(x_cl, y_cl, z_cl, 'b','LineWidth',1.8); grid on;
plot3(x_ol(1),y_ol(1),z_ol(1),'go','MarkerSize',10,'MarkerFaceColor','g');
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Altitude (ft)');
title('3D Trajectory — Open Loop vs LQR','FontWeight','bold','FontSize',13);
legend('Open Loop','LQR Closed Loop','Start','Location','best');


%% ===== Improved Eigenvalue Plot — OL vs CL =====
figure(12);
clf
hold on;
grid on;
box on;

% Coordinate axes
xline(0,'k--','LineWidth',1.5);
yline(0,'k--','LineWidth',1.5);

%% Plot eigenvalues

scatter(real(eig_lon_ol), imag(eig_lon_ol),...
    220,...
    'b','o',...
    'filled',...
    'MarkerEdgeColor','k',...
    'LineWidth',1.5,...
    'DisplayName','Longitudinal OL');

scatter(real(eig_lat_ol), imag(eig_lat_ol),...
    220,...
    'r','^',...
    'filled',...
    'MarkerEdgeColor','k',...
    'LineWidth',1.5,...
    'DisplayName','Lateral OL');

scatter(real(eig_lon_cl), imag(eig_lon_cl),...
    220,...
    'c','s',...
    'filled',...
    'MarkerEdgeColor','k',...
    'LineWidth',1.5,...
    'DisplayName','Longitudinal CL (LQR)');

scatter(real(eig_lat_cl), imag(eig_lat_cl),...
    220,...
    'm','d',...
    'filled',...
    'MarkerEdgeColor','k',...
    'LineWidth',1.5,...
    'DisplayName','Lateral CL (LQR)');

%% Labels

xlabel('Real Part','FontSize',13);
ylabel('Imaginary Part','FontSize',13);

title('Eigenvalue Plot — Open Loop vs Closed Loop (LQR)',...
      'FontWeight','bold',...
      'FontSize',15);

legend('Location','northwest',...
       'FontSize',11);

%% Appearance settings

set(gca,...
    'FontSize',12,...
    'LineWidth',1.5);

xlim([-4.8 0.5])
ylim([-5 5])

axis equal

%% ===== 6DOF EOM — OPEN LOOP (u_ctrl = external [del_e; del_a; del_r]) =====
function dydt = equations_of_motion(~, y, u_ctrl)

    del_e = u_ctrl(1);
    del_a = u_ctrl(2);
    del_r = u_ctrl(3);

    % ── UAV Parameters ──
    m   = 7800;
    Ix  = 54910.62746;
    Iy  = 17896.7971;
    Iz  = 68197.64348;
    Ixz = 0.1;
    g   = 9.79;
    b   = 9.8;
    rho = 0.6564;
    S   = 15.433;
    c   = 1.5748;
    T_max = 8000;

    % ── Aero Coefficients ──
    CL0 = 0.145404; CL_alpha = 4.97; CL_dele  =  0.017678;
    CD0 = 0.0130027; CD_alpha = 1.24065;
    CM0 = -0.032;   CM_alpha = -0.026798; CM_dele = -0.034968; CM_q = -0.30148;
    CY0 = 0; Cl0 = 0; Cn0 = 0;
    CY_beta = -0.01649;   Cl_beta = -0.24367; Cn_beta = 0.176;
    CY_p    = -0.229238; CY_r   =  0.1760;
    Cl_p    = -0.281678; Cl_r   =  0.118184;
    Cn_p    = -0.119354; Cn_r   = -0.506769;
    CY_dela = 0.001;    Cl_dela = 0.016131; Cn_dela = -0.013539;
    CY_delr = 0.197848; Cl_delr = -0.008512; Cn_delr =-0.007339;

    % ── Trim Thrust ──
    q_bar   = 8131.31;
    Va      = sqrt(q_bar / (0.5*rho));
    CL_trim = (2*m*g) / (rho*S*Va^2);
    del_e_trim = -((CM0*CL_alpha + CM_alpha*(CL_trim - CL0)) / ...
                   (CM_dele*CL_alpha - CM_alpha*CL_dele));
    CD_trim = CD0 + 0.065*CL_trim^2;
    D_trim  = 0.5*rho*Va^2*S*CD_trim;
    delta_t = D_trim / T_max;

    % ── Inertia Terms ──
    GAMA  = Ix*Iz - Ixz^2;
    GAMA1 = Ixz*(Ix - Iy + Iz) / GAMA;
    GAMA2 = (Iz*(Iz - Iy) + Ixz^2) / GAMA;
    GAMA3 = Iz / GAMA;
    GAMA4 = Ixz / GAMA;
    GAMA5 = (Iz - Ix) / Iy;
    GAMA6 = Ixz / Iy;
    GAMA7 = ((Ix - Iy)*Ix + Ixz^2) / GAMA;
    GAMA8 = Ix / GAMA;

    % ── Unpack States ──
    u = y(1); v = y(2); w = y(3);
    p = y(4); q = y(5); r = y(6);
    phi = y(7); theta = y(8); psi = y(9);

    % ── Derived ──
    V     = sqrt(u^2 + v^2 + w^2);
    alpha = atan2(w, u);
    beta  = asin(v / V);

    % ── Aero coefficients (state + control dependent) ──
    CL = CL0 + CL_alpha*alpha + CL_dele*del_e;
    CD = CD0 + CD_alpha*alpha + 0.065*CL^2;
    CY = CY0 + CY_beta*beta + CY_p*p*b/(2*V) + CY_r*r*b/(2*V) + CY_dela*del_a + CY_delr*del_r;
    Cl = Cl0 + Cl_beta*beta + Cl_p*p*b/(2*V) + Cl_r*r*b/(2*V) + Cl_dela*del_a + Cl_delr*del_r;
    Cn = Cn0 + Cn_beta*beta + Cn_p*p*b/(2*V) + Cn_r*r*b/(2*V) + Cn_dela*del_a + Cn_delr*del_r;
    Cm = CM0 + CM_alpha*alpha + CM_q*q*c/(2*V) + CM_dele*(del_e + del_e_trim);

    % ── Forces & Moments ──
    Lf = 0.5*rho*V^2*S*CL;
    D  = 0.5*rho*V^2*S*CD;
    Y  = 0.5*rho*V^2*S*CY;
    Mx = 0.5*rho*V^2*S*b*Cl;
    My = 0.5*rho*V^2*S*c*Cm;
    Mz = 0.5*rho*V^2*S*b*Cn;
    T  = T_max * delta_t;

    Fx = -D*cos(alpha) + Lf*sin(alpha) + T;
    Fz = -Lf*cos(alpha) - D*sin(alpha);
    Fy =  Y;

    Ax = Fx/m - g*sin(theta);
    Ay = Fy/m + g*sin(phi)*cos(theta);
    Az = Fz/m + g*cos(phi)*cos(theta);

    du = Ax + r*v - q*w;
    dv = Ay + p*w - r*u;
    dw = Az + q*u - p*v;

    dp = GAMA1*p*q - GAMA2*q*r + GAMA3*Mx + GAMA4*Mz;
    dq = GAMA5*p*r - GAMA6*(p^2 - r^2) + My/Iy;
    dr_dot = GAMA7*p*q - GAMA1*q*r + GAMA4*Mx + GAMA8*Mz;

    dphi   = p + q*sin(phi)*tan(theta) + r*cos(phi)*tan(theta);
    dtheta = q*cos(phi) - r*sin(phi);
    dpsi   = q*sin(phi)*sec(theta) + r*cos(phi)*sec(theta);

    dx_pos = u*cos(theta)*cos(psi) ...
           + v*(sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi)) ...
           + w*(cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi));
    dy_pos = u*cos(theta)*sin(psi) ...
           + v*(sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi)) ...
           + w*(cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi));
    dz_pos = -u*sin(theta) + v*sin(phi)*cos(theta) + w*cos(phi)*cos(theta);

    dydt = [du; dv; dw; dp; dq; dr_dot; dphi; dtheta; dpsi; dx_pos; dy_pos; dz_pos];
end

%% ===== 6DOF EOM — LQR CLOSED LOOP =====
function dydt = equations_of_motion_lqr(t, y, x_trim, K_lon, K_lat, idx_lon, idx_lat)
    % Compute LQR perturbation feedback
    dx_lon   = y(idx_lon) - x_trim(idx_lon);   % longitudinal state error
    dx_lat   = y(idx_lat) - x_trim(idx_lat);   % lateral state error

    u_lon    = -K_lon * dx_lon;   % scalar: elevator perturbation
    u_lat    = -K_lat * dx_lat;   % 2x1: [aileron; rudder] perturbation

    u_ctrl   = [u_lon(1); u_lat(1); u_lat(2)];

    % Saturate deflections to ±25 deg
    u_ctrl = max(min(u_ctrl, 25*pi/180), -25*pi/180);

    dydt = equations_of_motion(t, y, u_ctrl);
end