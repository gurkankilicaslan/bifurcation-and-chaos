clear; 
clc; clf;
%% A SMALL REMINDER: WHEN YOU RUN THE CODE, PLEASE DO NOT PLAY WITH PLOTS OR CHANGE THEIR POSITIONS UNTIL RUNNING IS DONE
% SINCE IT CAN HAVE PROBLEM WHILE OPERATING !!!
a = 0.2; b= 0.2; c= 5; %defining parameters
x0 = [3 3 3]; %defining the initial values
tspan = [0:0.01:1000]; %defining the range and iteration

[t,x] = ode45(@(t,x) odefcn(x,a, b, c), tspan, x0); %solution with Runge-Kutta 4/5
figure(1) % to open the result in a seperate window 
plot3(x(:,1),x(:,2),x(:,3),'LineWidth',1.2); % to plot the solution "Map of the Rössler system"
xlabel("X"); %x axis name
ylabel("Y"); %y axis name
zlabel("Z"); %z axis name
axis equal;  % equating the axis lengths
set(gca,'FontSize',15); %defining the font size

d = exp(-3:0.25:2.5)';
C(1:size(d,1),1) = 0;
N = size(x,1);
Nc = 1000;
r = randperm(N);
for n = 1 : Nc
    for m = 1 : size(d,1)
        C(m,1) = C(m,1) + sum((x(:,1)-x(r(n),1)).^2 ...
            + (x(:,2)-x(r(n),2)).^2 + (x(:,3)-x(r(n),3)).^2 < d(m,1)^2);
    end
end
C = C / Nc;
figure(2); % RESULT FOR PART A
plot(log(d),log(C),'.-r','MarkerSize',50); %fitting the fracture
set(gca,'FontSize',20);



dz = (x(3:N,3) - x(2:N-1,3)).*(x(2:N-1,3) - x(1:N-2,3)); % finding extremum
ddz = x(3:N,3) + x(1:N-2,3) - 2*x(2:N-1,3); % finding maximum
zn = x((dz<0).*(ddz<0)==1,3); %defining z_n

figure(3) % RESULT FOR PART B 
plot(zn(1:end-1),zn(2:end),'.k'); % plotting z_n vs z_n_+_1
xlabel('z_n') %x axis name
ylabel('z_n_+_1') %y axis name
line([-5 20],[-5 20],'Color','r'); %drawing the middle line
axis equal;



a = 0.2; b= 0.2; %defining parameters a and b
tspan = [0:0.001:100]; %defining the range and iteration with increase 0.001

figure(4) % RESULT FOR PART C
hold on;
c_value = (2:0.01:6); % changing c from 2 to 6 with 0.01 increase at each step

for nc = 1 : size(c_value,2) % solution via Runge-Kutta 4/5 at each step
    c = c_value(nc);
    [t,x] = ode45(@(t,x) odefcn(x,a,b,c), tspan, x0);
    x0 = x(end,:);
    [t,x] = ode45(@(t,x) odefcn(x,a,b,c), tspan, x0);
    N = size(x,1);
    dz = (x(3:N,1) - x(2:N-1,1)).*(x(2:N-1,1) - x(1:N-2,1)); % find extremum
    ddz = x(3:N,1) + x(1:N-2,1) - 2*x(2:N-1,1); % find maximum
    zn = x((dz<0).*(ddz<0)==1,1);
    plot(c,zn(1:end),'.');
    drawnow; % adding each step on and on continuously
    x0 = x(end,:); %changing the initial value to the end value of result of the previous step
end

result = (exp(-0.04104)-exp(-0.02961))/(exp(-0.04343)-exp(-0.04104));
result % Feigenbaum constant is found

function dxdt = odefcn(x,a,b,c) %defining the equation system
dxdt = zeros(2,1);
dxdt(1) = -x(2) - x(3);
dxdt(2) = x(1) + a*x(2);
dxdt(3) = b + x(3)*(x(1) - c);
end

