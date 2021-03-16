% Finite Difference Analysis of Transient Heat Transfer in 1D
% The problem concerns solving the heat transfer in 1D using finite difference method.
% The problem is an extension of the topic covered in the lectures:
% finite difference method - lecture 5-6
% and heat transfer - lecture 10.
clear

X = 0.3;            %[m]
dx = 0.01;          %[m]
M = X / dx;         %num of setps for x

T = 48 * 60 * 60;   %[s]
dt = 1;             %[s]
N = T / dt;         %num of steps fot t

%cegla
end_of_brick = 2/3*M;
k(1:end_of_brick) = 0.2;     %heat conductivity [W / (m*K)]
Cp(1:end_of_brick) = 880;    %heat capacity [J / (kg*K)]
p(1:end_of_brick) = 1800;    %density [kg / m^3]

%styropian izolacyjny
k(end_of_brick-1:M) = 0.04;        
Cp(end_of_brick-1:M) = 1300;
p(end_of_brick-1:M) = 25;

% srating temperature
u = zeros(N, M);

u(:,1) = 273.15 + 20;
u(1,2:M-1) = 273.15 + -10;
u(:,M) = 273.15 + -15;

u(2,1:M) = u(1,1:M);

% heat source
f = zeros(N, M);
heat_source_x = 1/3 * M;
f(N/2:N,heat_source_x-2:heat_source_x+2) = 273.15 + 300;

% simulation
w = ((2 * dt) ./ (Cp(2:M-1) .* p(2:M-1) * dx * dx));
for n = 2 : N - 1    
    u(n+1,2:M-1) = w .*...
        ((k(2:M-1) .* (u(n,3:M) - 2 * u(n,2:M-1) + u(n,1:M-2))) + ...
        (0.25 * (k(3:M) - k(1:M-2)) .* (u(n,3:M) - u(n,1:M-2)))) + ...
        f(n+1,2:M-1) ./ (Cp(2:M-1) .* p(2:M-1)) + ...
        u(n,2:M-1);
end

% end result
u = u - 273.15; %K -> *C

figure(1)
surf(linspace(0,X,M), linspace(0,T,N), u(:,1:M), 'EdgeColor', "none")
title('u(x,t)')
xlabel('x[m]')
ylabel('t[s]')
zlabel('temperature [*C]')

figure(2)
plot(linspace(0,X,M), u(N,:).','Color','red')
title('Last iteration material temperature')
xlabel('x[m]')
ylabel('u[*C]')