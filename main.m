% Since multiple plots were obtained from the same base graph, I will
% end every comment with a number(s) which means that that section of the
% code must be run to get that numbered figure ie, if a comment ends with 
% (2, 5), it must be run so that the 2nd figure and 5th figure will be
% displayed correctly.

A_x = 0.5;
A_y = sin(pi/3);
C_x = 0;
C_y = 0;
S_x = 1;
S_y = 0;

plot([C_x S_x], [C_y, S_y], 'k', [S_x A_x], [S_y, A_y], 'k', [A_x C_x], [A_y C_y], 'k');
hold on;

m_AC = tan(pi/3);
m_CS = 0;
m_AS = tan(pi - pi/3);

for i=1:9
    x1 = i/10;
    y1 = 0;
    x2 = 0.5+sin(pi/6)*x1;
    y2 = sin(pi/3)*(1-x1);
    
    x3 = 1-x1;
    y3 = y1;
    x4 = 0.5-sin(pi/6)*x1;
    y4 = y2;
    
    plot([x1 x2], [y1 y2], 'Color', '#CDCDCD');
    plot([x3 x4], [y3 y4], 'Color', '#CDCDCD');
    plot([x2 x4], [y2 y4], 'Color', '#CDCDCD');
end

A_r = [0.55 0.5 0.4 0.3 0.2 0.1];
C_r = [0.35 0.43 0.57 0.68 0.79 0.895];
S_r = [0.1 0.07 0.03 0.02 0.01 0.005];

A_e = [0.6 0.5 0.4 0.3 0.2 0.1];
C_e = [0.13 0.04 0.03 0.02 0.015 0.01];
S_e = [0.27 0.46 0.57 0.68 0.785 0.89];

R = zeros(6, 2);
E = zeros(6, 2);
for i = 1:6
    R(i, 1) = (A_r(i)/tan(pi/3) + S_r(i)/sin(pi/3))*sin(pi/3);
    E(i, 1) = (A_e(i)/tan(pi/3) + S_e(i)/sin(pi/3))*sin(pi/3);
    R(i, 2) = A_r(i)*sin(pi/3);
    E(i, 2) = A_e(i)*sin(pi/3);
end

% Uncomment this section to get the first graph (1, 2, 3, 4, 5, 6, 7)
% plot(R(:, 1), R(:, 2), 'b*');
% plot(E(:, 1), E(:, 2), 'r*');
% xx_eq = R(6, 1):1e-2:E(6, 1);
% eq = spline([R(:, 1); E(:, 1)]', [R(:, 2); E(:, 2)]', xx_eq);
% plot(xx_eq, eq, 'g');

A_r_tie = [0.44 0.29 0.12];
A_r_tie_prime = A_r_tie*sin(pi/3);
A_e_tie = [0.56 0.4 0.18];
A_e_tie_prime = A_e_tie*sin(pi/3);
y_A_r = A_r_tie_prime;
y_A_e = A_e_tie_prime;

x_A_r = zeros(1, 3);
x_A_e = zeros(1, 3);

for i = 1:3
    fun_r = @(x) y_A_r(i) - spline([R(:, 1); E(:, 1)]', [R(:, 2); E(:, 2)]', x);
    fun_e = @(x) y_A_e(i) - spline([R(:, 1); E(:, 1)]', [R(:, 2); E(:, 2)]', x);
    x_A_r(i) = fzero(fun_r, 0);
    x_A_e(i) = fzero(fun_e, 0.8);
end

Rn_x = R(6, 1);
Rn_y = R(6, 2);

p_RnS = [(S_y-Rn_y)/(S_x-Rn_x) Rn_y-(S_y-Rn_y)*Rn_x/(S_x-Rn_x)];
p_tie1 = [(y_A_r(1)-y_A_e(1))/(x_A_r(1)-x_A_e(1)) y_A_e(1)-(y_A_r(1)-y_A_e(1))*x_A_e(1)/(x_A_r(1)-x_A_e(1))];
p_tie2 = [(y_A_r(2)-y_A_e(2))/(x_A_r(2)-x_A_e(2)) y_A_e(2)-(y_A_r(2)-y_A_e(2))*x_A_e(2)/(x_A_r(2)-x_A_e(2))];
p_tie3 = [(y_A_r(3)-y_A_e(3))/(x_A_r(3)-x_A_e(3)) y_A_e(3)-(y_A_r(3)-y_A_e(3))*x_A_e(3)/(x_A_r(3)-x_A_e(3))];

xx_tie = -1:1e-2:1;
yy_RnS = polyval(p_RnS, xx_tie);
yy_tie1 = polyval(p_tie1, xx_tie);
yy_tie2 = polyval(p_tie2, xx_tie);
yy_tie3 = polyval(p_tie3, xx_tie);

% Uncomment this section to see which tie line intersects the operating
% line (2)
% plot(xx_tie, yy_RnS, 'k');
% plot(xx_tie, yy_tie1, 'r');
% plot(xx_tie, yy_tie2, 'b');
% plot(xx_tie, yy_tie3, 'g');

delta_x = fzero(@(x) polyval(p_RnS, x) - polyval(p_tie1, x), -0.5);
delta_y = polyval(p_RnS, delta_x);

F_A = 0.35;
F_S = 0;
F_Aprime = F_A*sin(pi/3);
F_Sprime = 0;
F_x = F_Aprime/tan(pi/3);
F_y = F_Aprime;

p_deltaF = [(delta_y-F_y)/(delta_x-F_x) F_y-(delta_y-F_y)*F_x/(delta_x-F_x)];
E1_x = fzero(@(x) polyval(p_deltaF, x) - spline([R(:, 1); E(:, 1)]', [R(:, 2); E(:, 2)]', x), 0.8);
E1_y = polyval(p_deltaF, E1_x);

E1_Aprime = E1_y;
E1_Sprime = E1_x*sin(pi/3) - E1_y*cos(pi/3);
E1_A = E1_Aprime/sin(pi/3);
E1_S = E1_Sprime/sin(pi/3);
E1_C = 1 - (E1_A + E1_S);

% Uncomment this to find out E1 for minimum solvent case (3)
% plot([delta_x S_x], [delta_y S_y], 'k');
% plot([delta_x E1_x], [delta_y E1_y], '--b');

% Uncomment this to find out M (4)
% plot([F_x S_x], [F_y S_y], 'b');
% plot([Rn_x E1_x], [Rn_y E1_y], 'r');

p_FS = [(F_y-S_y)/(F_x-S_x) S_y-(F_y-S_y)*S_x/(F_x-S_x)];
p_RnE1 = [(Rn_y-E1_y)/(Rn_x-E1_x) E1_y-(Rn_y-E1_y)*E1_x/(Rn_x-E1_x)];
M_x = fzero(@(x) polyval(p_FS, x) - polyval(p_RnE1, x), 0.5);
M_y = polyval(p_FS, M_x);
M_Aprime = M_y;
M_Sprime = M_x*sin(pi/3) - M_y*cos(pi/3);
M_A = M_Aprime/sin(pi/3);
M_S = M_Sprime/sin(pi/3);
M_C = 1 - (M_A + M_S);

F = 1300;
S_min = F*(F_A-M_A)/M_A;

S = 1.5*S_min;
M_a_A = F*F_A/(1300+S);
M_a_S = S/(F+S);
M_a_C = 1 - (M_a_A + M_a_S);
M_a_Aprime = M_a_A*sin(pi/3);
M_a_Sprime = M_a_S*sin(pi/3);
M_a_x = M_a_Aprime/tan(pi/3)+M_a_Sprime/sin(pi/3);
M_a_y = M_a_Aprime;

pRnM_a = [(Rn_y-M_a_y)/(Rn_x-M_a_x) M_a_y-(Rn_y-M_a_y)*(M_a_x)/(Rn_x-M_a_x)];
E1_a_x = fzero(@(x) polyval(pRnM_a, x) - spline([R(:, 1); E(:, 1)]', [R(:, 2); E(:, 2)]', x), 0.8);
E1_a_y = polyval(pRnM_a, E1_a_x);

% Uncomment this to get the actual E1 (5, 7)
% plot([Rn_x E1_a_x], [Rn_y E1_a_y], 'r');

E1_a_Aprime = E1_a_y;
E1_a_Sprime = E1_a_x*sin(pi/3) - E1_a_y*cos(pi/3);
E1_a_A = E1_a_Aprime/sin(pi/3);
E1_a_S = E1_a_Sprime/sin(pi/3);
E1_a_C = 1 - (E1_a_A + E1_a_S);

p_FE1_a = [(F_y-E1_a_y)/(F_x-E1_a_x) F_y-(F_y-E1_a_y)*F_x/(F_x-E1_a_x)];
delta_a_x = fzero(@(x) polyval(p_FE1_a, x) - polyval(p_RnS, x), -1);
delta_a_y = polyval(p_FE1_a, delta_a_x);

% Uncomment this to get the actual difference point (6)
% plot([delta_a_x E1_a_x], [delta_a_y E1_a_y], '--b');
% plot([delta_a_x S_x], [delta_a_y S_y], '--b');

S_e_tie_prime = x_A_e*sin(pi/3) - y_A_e*cos(pi/3);
conjugate_Sprime = S_e_tie_prime;
conjugate_Aprime = A_r_tie_prime;
conjugate_x = conjugate_Aprime/tan(pi/3) + conjugate_Sprime/sin(pi/3);
conjugate_y = conjugate_Aprime;
p_conjugate = polyfit(conjugate_x, conjugate_y, 1);

new_tie1_Sprime = E1_a_Sprime;
new_tie1_x = fzero(@(x) new_tie1_Sprime - (x*sin(pi/3)-polyval(p_conjugate, x)*cos(pi/3)), 0.6);
new_tie1_y = polyval(p_conjugate, new_tie1_x);
R1_y = new_tie1_y;
R1_x = fzero(@(x) R1_y - spline([R(:, 1); E(:, 1)]', [R(:, 2); E(:, 2)]', x), 0.1);
% Uncomment this to get the first tie line (7)
% plot([E1_a_x R1_x], [E1_a_y R1_y], 'b');
R1_Aprime = R1_y;
R1_Sprime = R1_x*sin(pi/3) - R1_y*cos(pi/3);
R1_A = R1_Aprime/sin(pi/3);
R1_S = R1_Sprime/sin(pi/3);
R1_C = 1 - (R1_A + R1_S);


p_new_tie12 = [(delta_a_y-R1_y)/(delta_a_x-R1_x) R1_y-(delta_a_y-R1_y)*R1_x/(delta_a_x-R1_x)];
E2_x = fzero(@(x) polyval(p_new_tie12, x) - spline([R(:, 1); E(:, 1)]', [R(:, 2); E(:, 2)]', x), 0.8);
E2_y = polyval(p_new_tie12, E2_x);
E2_Aprime = E2_y;
E2_Sprime = E2_x*sin(pi/3) - E2_y*cos(pi/3);
E2_A = E2_Aprime/sin(pi/3);
E2_S = E2_Sprime/sin(pi/3);
E2_C = 1 - (E2_A + E2_S);
xx_new_tie12 = 0:1e-2:E2_x;

% Uncomment this to plot the line joining delta_a and R1 (7)
% plot(xx_new_tie12, polyval(p_new_tie12, xx_new_tie12), '--b');

E2_Sprime = E2_x*sin(pi/3)-E2_y*cos(pi/3);
new_tie2_Sprime = E2_Sprime;
new_tie2_x = fzero(@(x) new_tie2_Sprime - (x*sin(pi/3)-polyval(p_conjugate, x)*cos(pi/3)), 0.6);
new_tie2_y = polyval(p_conjugate, new_tie2_x);
R2_y = new_tie2_y;
R2_x = fzero(@(x) R2_y - spline([R(:, 1); E(:, 1)]', [R(:, 2); E(:, 2)]', x), 0.1);
% Uncomment this to plot the second tie line (7)
% plot([E2_x R2_x], [E2_y R2_y], 'b');
R2_Aprime = R2_y;
R2_Sprime = R2_x*sin(pi/3) - R2_y*cos(pi/3);
R2_A = R2_Aprime/sin(pi/3);
R2_S = R2_Sprime/sin(pi/3);
R2_C = 1 - (R2_A + R2_S);

p_new_tie23 = [(delta_a_y-R2_y)/(delta_a_x-R2_x) R2_y-(delta_a_y-R2_y)*R2_x/(delta_a_x-R2_x)];
E3_x = fzero(@(x) polyval(p_new_tie23, x) - spline([R(:, 1); E(:, 1)]', [R(:, 2); E(:, 2)]', x), 0.8);
E3_y = polyval(p_new_tie23, E3_x);
E3_Aprime = E3_y;
E3_Sprime = E3_x*sin(pi/3) - E3_y*cos(pi/3);
E3_A = E3_Aprime/sin(pi/3);
E3_S = E3_Sprime/sin(pi/3);
E3_C = 1 - (E3_A + E3_S);
xx_new_tie23 = 0:1e-2:E3_x;

% Uncomment this to plot the line joining delta_a and R2 (7)
% plot(xx_new_tie23, polyval(p_new_tie23, xx_new_tie23), '--b');

E3_Sprime = E3_x*sin(pi/3)-E3_y*cos(pi/3);
new_tie3_Sprime = E3_Sprime;
new_tie3_x = fzero(@(x) new_tie3_Sprime - (x*sin(pi/3)-polyval(p_conjugate, x)*cos(pi/3)), 0.6);
new_tie3_y = polyval(p_conjugate, new_tie3_x);
R3_y = new_tie3_y;
R3_x = fzero(@(x) R3_y - spline([R(:, 1); E(:, 1)]', [R(:, 2); E(:, 2)]', x), 0.1);
% Uncomment this to get the third tie line (7)
% plot([E3_x R3_x], [E3_y R3_y], 'b');
R3_Aprime = R3_y;
R3_Sprime = R3_x*sin(pi/3) - R3_y*cos(pi/3);
R3_A = R3_Aprime/sin(pi/3);
R3_S = R3_Sprime/sin(pi/3);
R3_C = 1 - (R3_A + R3_S);

xx_operating1 = 0:1e-2:E1_a_x;
xx_operating2 = 0:1e-2:S_x;

% Uncomment this to plot a portion of the  operating lines (7)
% plot(xx_operating1, polyval(p_FE1_a, xx_operating1), '--b');
% plot(xx_operating2, polyval(p_RnS, xx_operating2), '--b');