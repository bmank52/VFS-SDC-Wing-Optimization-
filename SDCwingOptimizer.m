% wing area optimizer
% Bryce Mankovsky
clc;
clear all;
close all;
%Constants
rho = 1.225; % kg/m^3
WminusWing = 582*9.8; % N
nu = 0.000014956; % m^2/s, kinematic
c = 339.14; % speed of sound, m/s

%Dimensions
l_fueslage = 3.5; %m
r_fueslage = 0.75; %m
t_half__fueslage = 0.625; %m
Sref_f = pi*r_fueslage*t_half__fueslage; %m^2

Sref_w = 12.78688525; %m^2
AR_vec = linspace(4,7.8,30);
v = linspace(20,75, 30); % m/s

[AR_grid, V_grid] = meshgrid(AR_vec, v);
LD_grid = zeros(size(AR_grid));   %for plotting later
C_Lgrid = zeros(size(AR_grid));
C_Dgrid = zeros(size(AR_grid));

b_wing = zeros(size(AR_vec));
c_wing = zeros(size(AR_vec));
Swet_w = 2.04*Sref_w;

sweep = 0;
t_c_wing = 0.11;
CL_max = 1.25;

for i = 1:length(AR_vec)

b_wing(i) = sqrt(AR_vec(i)*Sref_w);
c_wing(i) = Sref_w / b_wing(i);
end 

d_boom = 0.2;
t_boom = 0.3;
l_boom = 5;
Swet_boom = pi*d_boom*l_boom;

d_lg = 0.05;
l_lg = 1;
Sref_lg = (pi/4)*d_lg^2;

c_vtail = 1.3125;
b_vtail = 1.5;
sweep_vtail = 0;
t_c_vtail = 0.12;
Sref_vtail = b_vtail*c_vtail;
Swet_vtail = 2.04*Sref_vtail;

c_htail = 1.3125;
b_htail = 4;
sweep_htail = 0;
t_c_htail = 0.09;
Sref_htail = b_htail*c_htail;
Swet_htail = 2.04*Sref_htail;

d_tailboom = 0.3;
t_tailboom = 0.25;
l_tailboom = 1;
Stailboom = (pi/4)*d_tailboom*t_tailboom;
Swet_tailboom = 0.8657;
l_d_tailboom = l_tailboom/d_tailboom;

d_rotor = 1;
width_rotor= 0.1;
Sref_rotor = 8*d_rotor*width_rotor; % area swept by 6 rotors, flat plate
t_c_rotor = 0.1;

% 
Dtotal = zeros(length(v), length(AR_vec)); % Preallocate Dtotal
for i = 1:length(AR_vec)
    for j = 1:length(v)

M_0 = v(j)/c;
q = (1/2)*rho*v(j)^2;
W_estimate = WminusWing + 500;  % quick overestimate to be safe
C_L_check = (2*W_estimate) / (rho*Sref_w*v(j)^2);
if C_L_check > CL_max
    LD_grid(j,i) = NaN;
    continue
end
W = computeWeight_iterative(AR_vec(i), Sref_w, t_c_wing, q, 582) *9.8;
Re_length = v(j)/nu;

% Wing
Re = Re_length*c_wing(i);
cf = (204.8642*Re.^(-0.2805)+1.0409)/1000;
z = ((2-M_0^2)*cosd(sweep))/(sqrt(1-M_0^2*cosd(sweep)^2));
k = 1+z*t_c_wing+100*t_c_wing^4;
cDp_w = k.*cf.*Swet_w./Sref_w;
Dp_w =q*Sref_w*cDp_w;

% Fuselage
Re = Re_length*l_fueslage;
r_l = r_fueslage/l_fueslage;
cd0 = -13.353*r_l^9 + 98.818*r_l^8 + -313.37*r_l^7 + 556.41*r_l^6 + -606.66*r_l^5 + 419.07*r_l^4 + -182.56*r_l^3 + 48.371*r_l^2 + -7.1423*r_l+ 0.5127;
cd0_wing = cd0*Sref_f/Sref_w;
Dp_f = q*cd0_wing*Sref_w;

% Boom
% d_l = d_boom/l_boom;
% cd01(i) = -13.353*d_l^9 + 98.818*d_l^8 + -313.37*d_l^7 + 556.41*d_l^6 + -606.66*d_l^5 + 419.07*d_l^4 + -182.56*d_l^3 + 48.371*d_l^2 + -7.1423*d_l + 0.5127;
% Dp_b(i) = q*Sref_b*cd01(i);

% Landing Gear
% d_l = d_lg/l_lg;
% cd02(i) = -13.353*d_l^9 + 98.818*d_l^8 + -313.37*d_l^7 + 556.41*d_l^6 + -606.66*d_l^5 + 419.07*d_l^4 + -182.56*d_l^3 + 48.371*d_l^2 + -7.1423*d_l + 0.5127;
% Dp_lg(i) = q*Sref_lg*cd02(i);

% V Tail
Re = Re_length*b_vtail;
cf = (204.8642*Re^(-0.2805)+1.0409)/1000;
z = ((2-M_0^2)*cosd(sweep_vtail))/(sqrt(1-M_0^2*cosd(sweep_vtail)^2));
k = 1+z*t_c_vtail+100*t_c_vtail^4;
cDp_vtail = k*cf*Swet_vtail/Sref_w;
Dp_vtail = q*Sref_w*cDp_vtail;

%H Tail
Re = Re_length*b_htail;
cf = (204.8642*Re^(-0.2805)+1.0409)/1000;
z = ((2-M_0^2)*cosd(sweep_htail))/(sqrt(1-M_0^2*cosd(sweep_htail)^2));
k = 1+z*t_c_htail+100*t_c_htail^4;
cDp_htail= k*cf*Swet_htail/Sref_w;
Dp_htail = q*Sref_w*cDp_htail;

% Tail Boom
Re = Re_length*l_boom;
cf = (204.8642*Re^(-0.2805)+1.0409)/1000;
k = 1.01;
cDp_tailboom = k*cf*Swet_boom/Sref_w;
Dp_tailboom = q*Sref_w*cDp_tailboom;

% Rotors
Re = Re_length*d_rotor;
cf = (204.8642*Re^(-0.2805)+1.0409)/1000;
z = ((2-M_0^2)*cosd(0))/(sqrt(1-M_0^2*cosd(0)^2));
k = 1+z*t_c_rotor+100*t_c_rotor^4;
cDp_rotor = k*cf*Sref_rotor/Sref_w;
Dp_rotor = q*Sref_w*cDp_rotor;

% Induced Drag
C_L = (2*W)/(rho*Sref_w*v(j)*v(j));
if C_L > CL_max
    LD_grid(j,i) = 0;  % or 0 if you prefer
    continue
end
Cdi = (C_L^2)/(pi*AR_vec(i)*0.9);
Di = q*Sref_w*Cdi;

C_D = Cdi + cDp_rotor + cDp_tailboom + cDp_vtail + cd0_wing + cDp_w; 

% Prop Endurance Cond. (cl^3/2 / cd)
%enduranceCond(i) = C_L(i)^(3/2) / C_D(i);

Dp = Dp_tailboom+Dp_vtail+Dp_f+Dp_w+Dp_rotor;
Dtotal(j,i) = Dp + Di;
PowerRequired = ((Dtotal(j,i) + W * sind(9)) * v(j)) / 1000; % kW
% PR = TR*v = (D+Wsin(gamma))*v    ^^^ flight path angle = 9

LD_grid(j,i) = (C_L^(3/2)) / C_D;
C_Lgrid(j,i) = (C_L);
C_Dgrid(j,i) = (C_D);
    end
end

function W = computeWeight_iterative(AR, Sref_w_m2, t_c_wing, q_Pa, W_N_noWing)
    Sref_w_ft2 = Sref_w_m2 * 10.7639;
    q_psf = q_Pa * 0.0208854;
    W_lb_noWing = W_N_noWing / 4.44822;
    W_lb = W_lb_noWing;  % Initial guess: non-wing weight only
    tol = 1e-3;
    maxIter = 100;

    for iter = 1:maxIter
        W_wing_lb = 0.036 * Sref_w_ft2^0.578 * AR^0.6 * q_psf^0.006 * ...
                    (100*t_c_wing)^-0.3 * (3.75*W_lb)^0.49;

        W_lb_new = W_lb_noWing + W_wing_lb;

        if abs(W_lb_new - W_lb) < tol
            break
        end
        W_lb = W_lb_new;
    end

    W = W_lb * 4.44822;  % convert lb → N
end

[maxEnduranceCond, idx] = max(LD_grid(:));  % Find max L/D and its index in linear form
[row, col] = ind2sub(size(LD_grid), idx);  % Convert index to row (v) and column (AR)

best_AR = AR_vec(col);  % Use column index to extract AR from AR_vec
best_v = v(row);        % Use row index to extract velocity from v

[~, best_AR_idx] = min(abs(AR_vec - best_AR));

% Get corresponding span and chord
best_span = b_wing(best_AR_idx);
best_chord = c_wing(best_AR_idx);
best_w = computeWeight_iterative(best_AR, Sref_w, t_c_wing, .5*rho*best_v^2, WminusWing);

fprintf('Best Span = %.2f m\n', best_span);
fprintf('Best Chord = %.2f m\n', best_chord);
fprintf('Best Weight = %.2f N\n', best_w)

fprintf('Max CL^{3/2}/CD = %.2f at AR = %.2f and Velocity = %.2f m/s\n', maxEnduranceCond, best_AR, best_v);


figure
surf(AR_grid, V_grid, LD_grid)
xlabel('Aspect Ratio (AR)')
ylabel('Velocity (m/s)')
zlabel('CL^{1.5}/CD Ratio')
title('CL^{1.5}/CD vs AR vs Velocity')
colorbar
shading interp
grid on