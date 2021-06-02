clc
clear
close all
%% Defining constants 
z_k = 0.2; 
m = 0.1;
k = 10;
b = 10; 
z0 = 0.1; 
dot_z0 = 0; 
init = [z0, dot_z0];
g = 9.81; 
time_span = [0 10]; 
i = @(time_span) 10*sin(10*time_span);
p = 1.5e-5; 
e0 = 0.5e-3; 

%% Maintaining the stationary point

z_eq = 0.08; 
dot_z_eq = 0; 
ddot_z_eq = 0; 
% i_eq = (z_eq + e0) * sqrt(abs((k * (z_k - z_eq) - m * g) / p));

syms i_eq
i_eq = double(solve(m * ddot_z_eq == -p * (i_eq/(z_eq + e0)) ^ 2 - m*g - k * (z_eq - z_k) - b * dot_z_eq, i_eq)); 
disp('Equilibrium value of current is .. ')
disp(i_eq(2))
%% Derivation of displacement using ode45 

figure(1)
func = @(t,z) rhs_equation(t,z,m,k,b,e0,z_k,p,i,g);
[t,z] = ode45(func,time_span,init);
subplot(2,2,1)
plot(t,z(:,1),'r','DisplayName','z', 'LineWidth',2);
grid on
title('Response of the magnetic suspension system to sine input') 
ylabel('Displacement z (m)');
xlabel('Time (s)')
lgd = legend;
lgd.FontSize = 10;
lgd.Title.String = 'Legend';
subplot(2,2,3)
plot(t,z(:,2),'r','DisplayName','z_{dot}', 'LineWidth',2); 
grid on
ylabel('z_{dot} (m/s)');
xlabel('Time (s)')
lgd = legend;
lgd.FontSize = 10;
lgd.Title.String = 'Legend';

% For equilibrium plot
func = @(t,z) rhs_equation_eq(t,z,m,k,b,e0,z_k,p,i,g, i_eq);
[t,z] = ode45(func,time_span,init);
subplot(2,2,2)
plot(t,z(:,1),'b', 'DisplayName','z', 'LineWidth',2); 
grid on
title('Response of the magnetic suspension system to fixed equilibrium input') 
ylabel('Displacement z (m)');
xlabel('Time (s)')
ylim([0.075 0.101])
lgd = legend;
lgd.FontSize = 10;
lgd.Title.String = 'Legend';
subplot(2,2,4)
plot(t,z(:,2),'b','DisplayName','z_{dot}', 'LineWidth',2); 
grid on
ylabel('z_{dot} (m/s)');
xlabel('Time (s)')
lgd = legend;
lgd.FontSize = 10;
lgd.Title.String = 'Legend';

% definition of a function for ode45
function dzdt = rhs_equation(t,z,m,k,b,e0,z_k,p,i,g)
    F = p*(i(t)/(z(1)+e0))^2;
    dot_z = z(2);
    ddot_z = (-F + -m * g - k * (z(1) - z_k) - b * z(2)) / m;
    dzdt = [dot_z, ddot_z]'; 
end


% definition of a equilibrium function for ode45
function dzdt = rhs_equation_eq(t,z,m,k,b,e0,z_k,p,i,g, i_eq)
    F = p*(i_eq(2)/(z(1)+e0))^2;
    dot_z = z(2);
    ddot_z = (-F + -m * g - k * (z(1) - z_k) - b * z(2)) / m;
    dzdt = [dot_z, ddot_z]'; 
end
