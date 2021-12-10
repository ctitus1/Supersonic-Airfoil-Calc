clear;%clc;

%% Parameters
% Design Parameters
M_inf = 1.45921505; % inflow mach number (1.7)
P_inf = 7.819281777; % kPa (atmospheric pressure at 60,000 ft) (7.2326)
% P_inf = 143.8527652; % kPa (atmospheric pressure at 0 ft) (101.325)
chord = 1; % m
wingspan = 18.287; % m

% Gas Parameters
gamma = 1.4; % gamma

% Plot Parameters
wait_delay = 0;
upper_colors = [0 1 0;
                0 0.4470 0.7410;
                0.8500 0.3250 0.0980;
                0.9290 0.6940 0.1250];
lower_colors = [1 0 1;
                0.4940 0.1840 0.5560;
                0.4660 0.6740 0.1880;
                0.3010 0.7450 0.9330];

% Test Parameters
% upper_phi_range = 1;
upper_phi_range = 2.866;
% upper_phi_range = [3:1:10]; % deg (angle upper surface makes with chord) should be about 0.5x the lower
% upper_phi_range = linspace(1,6,7); % deg (angle upper surface makes with chord) should be about 0.5x the lower
% lower_phi_range = 3;
lower_phi_range = 2.866;
% lower_phi_range = [3:1:10]; % deg (angle lower surface makes with chord) should be about 2x the upper
% lower_phi_range = linspace(1,6,7); % deg (angle lower surface makes with chord) should be about 2x the upper
% alpha_range = [0:.5:6]
alpha_range = linspace(0,1.35,51);
% alpha_range = 1.35;
% max_alpha_0 = 10.6
thickness_chord = 1/2*(tand(upper_phi_range)+tand(lower_phi_range))
% Data
data.upper_phi = upper_phi_range;
data.lower_phi = lower_phi_range;
data.cl = [];
data.cd = [];
data.LD = [];

%% Outer Loops
total_loops = length(upper_phi_range)*length(lower_phi_range)*length(alpha_range);
loop = 0;
u = 0;
for upper_phi = upper_phi_range
u = u + 1;
l = 0;
for lower_phi = lower_phi_range
l = l + 1;
%% Inner Loop
% alpha = max_alpha_0 - lower_phi; % for seeing max alpha before flow separation
for alpha = alpha_range
loop = loop + 1;
% Shock Calculations
upper_sol = calc(alpha,upper_phi,"upper",M_inf,P_inf,chord,gamma);
lower_sol = calc(alpha,lower_phi,"lower",M_inf,P_inf,chord,gamma);
% Visualization
if(true)
figure(1)
clf
hold on
draw_airfoil(upper_sol,upper_colors,chord)
draw_airfoil(lower_sol,lower_colors,chord)
hold off
title("Airfoil Geometry at with phi_{upper} = "+upper_phi+", phi_{lower} = "+lower_phi+", \alpha = "+alpha)
% legend(["upper surface","upper shock_{12}","upper expansion_{23_1}","upper expansion_{23_2}","upper shock_{34}","lower surface","lower shock_{12}","lower expansion_{23_1}","lower expansion_{23_2}","lower shock_{34}"])
xlim([-.1*chord,1.1*chord])
ylim([-.6*chord .6*chord])
axis square
end
cl = (1/cosd(lower_phi)*1/(gamma*M_inf^2)*(lower_sol.P2/P_inf*cosd(abs(lower_sol.flow2))+lower_sol.P3/P_inf*cosd(abs(lower_sol.flow3)))) - (1/cosd(upper_phi)*1/(gamma*M_inf^2)*(upper_sol.P2/P_inf*cosd(abs(upper_sol.flow2))+upper_sol.P3/P_inf*cosd(abs(upper_sol.flow3))));
cd = (1/cosd(lower_phi)*1/(gamma*M_inf^2)*(lower_sol.P2/P_inf*sind(abs(lower_sol.flow2))-lower_sol.P3/P_inf*sind(abs(lower_sol.flow3)))) + (1/cosd(upper_phi)*1/(gamma*M_inf^2)*(upper_sol.P2/P_inf*sind(abs(upper_sol.flow2))-upper_sol.P3/P_inf*sind(abs(upper_sol.flow3))));
LD = cl/cd;
% Print
fprintf("=== Values for loop = %d/%d (% 6.2f%s\n",loop,total_loops,100*loop/total_loops,"% done) ===")
fprintf("upper_phi = %f\n",upper_phi)
fprintf("lower_phi = %f\n",lower_phi)
fprintf("alpha = %f\n",alpha)
fprintf("cl = %f\n",cl)
fprintf("cd = %f\n",cd)
fprintf("L/D = %f\n",LD)
fprintf("==============================================\n\n")

%%
% if(cl < 0); cl = 0; end
% if(cd < 0); cd = 0; end
% if(LD < 0); LD = 0; end

% for varying phis
% data.cl(u,l) = cl;
% data.cd(u,l) = cd;
% data.LD(u,l) = LD;

% for plotting varying alpha
data.cl = [data.cl;cl];
data.cd = [data.cd;cd];
data.LD = [data.LD;LD];

data.upper_phi(u,l) = upper_phi;
data.lower_phi(u,l) = lower_phi;

pause(wait_delay);
end % end inner loop

end % end outer loops
end
%% Plotting

X = alpha_range;
X_label = "\alpha (deg)";

figure(2)
clf
plot(X,data.cl)
title("cl vs \alpha")
xlabel(X_label)
ylabel("cl")

figure(3)
clf
plot(X,data.cd)
title("cd vs \alpha")
ylabel("cd")

figure(4)
clf
plot(X,data.LD)
title("L/D vs \alpha")
ylabel("L/D")

% X = data.upper_phi;
% Y = data.lower_phi;
% X_label = "X = upper phi (deg)";
% Y_label = "Y = lower phi (deg)";
% 
% figure(2)
% clf
% surf(X,Y,data.cl)
% title("cl vs upper and lower half angle at \alpha = "+alpha)
% xlabel(X_label)
% ylabel(Y_label)
% zlabel("cl")
% 
% figure(3)
% clf
% surf(X,Y,data.cd)
% title("cd vs upper and lower half angle at \alpha = "+alpha)
% xlabel(X_label)
% ylabel(Y_label)
% zlabel("cd")
% 
% figure(4)
% clf
% surf(X,Y,data.LD)
% title("L/D vs upper and lower half angle at \alpha = "+alpha)
% xlabel(X_label)
% ylabel(Y_label)
% zlabel("L/D")

%% Functions
function sol = calc(alpha,phi,upper_or_lower,M_inf,P_inf,chord,gamma) % +theta if flow deflecetd CCW
% State 1 (inlet)
sol.alpha = alpha;
sol.upper_or_lower = upper_or_lower;
surface = 1;
if(strcmpi(upper_or_lower,"lower")) % if upper surface
    surface = -surface;
elseif(~strcmpi(upper_or_lower,"upper"))
    error("Invalid input to calc(). Input either upper or lower")
end
sol.surface = surface;
phi = surface*phi; % absolute angle of airfoil rel to chord
sol.phi = phi;

% geometry leading edge is (0,0)
xy12 = [0,0];
xy23 = [chord/2*secd(phi)*cosd(phi-alpha),chord/2*secd(phi)*sind(phi-alpha)];
xy34 = [chord*cosd(alpha),-chord*sind(alpha)];
sol.points.x = [xy12(1);xy23(1);xy34(1)];
sol.points.y = [xy12(2);xy23(2);xy34(2)];

% state 1
flow1 = 0;
sol.flow1 = flow1;
M1 = M_inf;
sol.M1 = M_inf;
P1 = P_inf;
sol.P1 = P_inf;

% state 2
dflowdir12 = phi - alpha; % absolute flow direction change
sol.dflowdir12 = dflowdir12;
flow2 = flow1 + dflowdir12;
sol.flow2 = flow2;
theta12 = surface*dflowdir12; % + theta for compression, - for expansion
sol.theta12 = theta12;
if (theta12 > 0) % if compression
    beta12 = Beqn(theta12,M1,gamma); % angle of shock relative to inflow
sol.beta12 = beta12;
    beta12_abs = surface*beta12 + flow1; % absolute angle of shock
sol.beta12_abs = beta12_abs;
    M2 = MOeqn(beta12,theta12,M1,gamma);
sol.M2 = M2;
    P2 = POeqn(P1,beta12,M1,gamma);
sol.P2 = P2;
else % both expansion fans will be on top of eachother if no deflection, check for this when plotting
    M2 = MEeqn(theta12,M1,gamma);
sol.M2 = M2;
    P2 = PReqn(P1,M1,M2,gamma);
sol.P2 = P2;
    u12_a = asind(1/M1);
sol.u12_a = u12_a;
    u12_a_abs = surface*u12_a + flow1;
sol.u12_a_abs = u12_a_abs;
    u12_b = asind(1/M2);
sol.u12_b = u12_b;
    u12_b_abs = surface*u12_b + flow2;
sol.u12_b_abs = u12_b_abs;
end

% state 3
dflowdir23 = -2*phi; % absolute flow direction change
sol.dflowdir23 = dflowdir23;
flow3 = flow2 + dflowdir23;
sol.flow3 = flow3;
theta23 = surface*dflowdir23; % + theta for compression, - for expansion
sol.theta23 = theta23;
if (theta23 > 0) % if compression
    beta23 = Beqn(theta23,M2,gamma); % angle of shock relative to inflow
sol.beta23 = beta23;
    beta23_abs = surface*beta23 + flow2; % absolute angle of shock
sol.beta23_abs = beta23_abs;
    M3 = MOeqn(beta23,theta23,M2,gamma);
sol.M3 = M3;
    P3 = POeqn(P2,beta23,M2,gamma);
sol.P3 = P3;
else % both expansion fans will be on top of eachother if no deflection, check for this when plotting
    M3 = MEeqn(theta23,M2,gamma);
sol.M3 = M3;
    P3 = PReqn(P2,M2,M3,gamma);
sol.P3 = P3;
    u23_a = asind(1/M2);
sol.u23_a = u23_a;
    u23_a_abs = surface*u23_a + flow2;
sol.u23_a_abs = u23_a_abs;
    u23_b = asind(1/M3);
sol.u23_b = u23_b;
    u23_b_abs = surface*u23_b + flow3;
sol.u23_b_abs = u23_b_abs;
end

% state 4
dflowdir34 = phi + alpha; % absolute flow direction change
sol.dflowdir34 = dflowdir34;
flow4 = flow3 + dflowdir34;
sol.flow4 = flow4;
theta34 = surface*dflowdir34; % + theta for compression, - for expansion
sol.theta34 = theta34;
if (theta34 > 0) % if compression
    beta34 = Beqn(theta34,M3,gamma); % angle of shock relative to inflow
sol.beta34 = beta34;
    beta34_abs = surface*beta34 + flow3; % absolute angle of shock
sol.beta34_abs = beta34_abs;
    M4 = MOeqn(beta34,theta34,M3,gamma);
sol.M4 = M4;
    P4 = POeqn(P3,beta34,M3,gamma);
sol.P4 = P4;
else % both expansion fans will be on top of eachother if no deflection, check for this when plotting
    M4 = MEeqn(theta34,M3,gamma);
sol.M4 = M4;
    P4 = PReqn(P3,M3,M4,gamma);
sol.P4 = P4;
    u34_a = asind(1/M3);
sol.u34_a = u34_a;
    u34_a_abs = surface*u34_a + flow3;
sol.u34_a_abs = u34_a_abs;
    u34_b = asind(1/M4);
sol.u34_b = u34_b;
    u34_b_abs = surface*u34_b + flow4;
sol.u34_b_abs = u34_b_abs;
end
end

function beta = Beqn(theta,M,gamma)
syms eqn(b)
eqn(b) = tand(theta) - 2*cotd(b)*(M^2*sind(b)^2-1)/(M^2*(gamma+cosd(2*b))+2);
beta_thetamax = asind(sqrt(1/(4*gamma*M^2)*((gamma+1)*M^2-(3-gamma)+sqrt((gamma+1)*((gamma+1)*M^4-2*(3-gamma)*M^2+gamma+9)))));
beta = double(vpasolve(eqn == 0,[0 beta_thetamax])); % assumes weak shock solution
if(isempty(beta))
    error("beta is empty")
end
end

function M2 = MOeqn(beta,theta,M1,gamma) % oblique shock relation
theta = abs(theta);
M2 = sqrt((2+(gamma-1)*M1^2*sind(beta)^2)/((2*gamma*M1^2*sind(beta)^2-(gamma-1))*sind(beta-theta)^2));
end

function M2 = MEeqn(theta12,M1,gamma) % expansion/compression wave relation
syms M2
eqn = theta12 + (v(M2,gamma)-v(M1,gamma)); % +theta for CCW, -theta for CW
v_max = rad2deg(pi/2*(sqrt((gamma+1)/(gamma-1))-1));
M2 = double(vpasolve(eqn == 0,[0 v_max]));
if(isempty(M2))
    error("Mf is empty")
end
end

function P2 = POeqn(P1,beta12,M1,gamma) % normal shock relation for P
Mn1 = M1*sind(beta12);
P2 = P1*(1+(2*gamma)/(gamma+1)*(Mn1^2-1));
end

function P2 = PReqn(P1,M1,M2,gamma) % isentropic pressure relation
P2 = P1*((2+(gamma-1)*M1^2)/(2+(gamma-1)*M2^2))^(gamma/(gamma-1));
end

function v = v(M,gamma) %Prandtl-Meyer function
v = sqrt((gamma+1)/(gamma-1))*atand(sqrt((gamma-1)/(gamma+1)*(M^2-1)))-atand(sqrt(M^2-1));
end

function line(xi,yi,angle,length,color)
xf = xi+length*cosd(angle);
yf = yi+length*sind(angle);
plot([xi xf],[yi yf],'Color',color,'LineWidth',4)
end

function draw_airfoil(sol,colors,chord)
shock_length = 1*chord;
plot(sol.points.x,sol.points.y,'Color',colors(1,:),'LineWidth',4)

if(sol.theta12 > 0)
    line(sol.points.x(1),sol.points.y(1),sol.beta12_abs+sol.flow1,shock_length,colors(2,:))
elseif(sol.flow1 ~= sol.flow2)
    line(sol.points.x(1),sol.points.y(1),sol.u12_a_abs+sol.flow1,shock_length,colors(2,:))
    line(sol.points.x(1),sol.points.y(1),sol.u12_b_abs+sol.flow2,shock_length,colors(2,:))
end
if(sol.theta23 > 0)
    line(sol.points.x(2),sol.points.y(2),sol.beta23_abs+sol.flow2,shock_length,colors(3,:))
elseif(sol.flow2 ~= sol.flow3)
    line(sol.points.x(2),sol.points.y(2),sol.u23_a_abs+sol.flow2,shock_length,colors(3,:))
    line(sol.points.x(2),sol.points.y(2),sol.u23_b_abs+sol.flow3,shock_length,colors(3,:))
end
if(sol.theta34 > 0)
    line(sol.points.x(3),sol.points.y(3),sol.beta34_abs+sol.flow3,shock_length,colors(4,:))
elseif(sol.flow3 ~= sol.flow4)
    line(sol.points.x(3),sol.points.y(3),sol.u34_a_abs+sol.flow3,shock_length,colors(4,:))
    line(sol.points.x(3),sol.points.y(3),sol.u34_b_abs+sol.flow4,shock_length,colors(4,:))
end
end
%% Unused Functions
function M2 = MNeqn(M1,gamma) % normal shock relation
M2 = sqrt((2+(gamma-1)*M1^2)/(2*gamma*M1^2-(gamma-1)));
end
function T2 = TNeqn(T1,Mn1,gamma) % normal shock relation for T
T2 = T1*(1+(2*gamma)/(gamma+1)*(Mn1^2-1))*(2+(gamma-1)*Mn1^2)/((gamma+1)*Mn1^2);
end
function T0 = T0eqn(T,M,gamma) % stagnation temperature relation
T0 = T*(1+(gamma-1)/2*M^2);
end
function P0 = P0eqn(P,M,gamma) % stagnation pressure relation
P0 = P*(1+(gamma-1)/2*M^2)^(gamma/(gamma-1));
end
function s2 = Seqn(s1,T2,T1,P2,P1,gamma,R) % computes s2 = s2 + ds
s2 = s1 + R*log((T2/T1)^(gamma/(gamma-1))*(P2/P1)^-1);
end
function T2 = TReqn(T1,M1,M2,gamma) % isentropic temperature relation
T2 = T1*(2+(gamma-1)*M1^2)/(2+(gamma-1)*M2^2);
end