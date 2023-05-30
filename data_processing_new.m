%% Initialization
clear all;clc;close all

%% Input parameters

g = 10.5;

%us acc 123 giro 123 magn 123 rotaz 1234 p h t 

%filename = 'dataLog_m_numerogeneratoaacaso.txt';
filename = 'vediamo.TXT';
data = readmatrix(filename);

time = data(:,1)/(10^6);
acc_raw = data(:,2:4);
acc = acc_raw;
giro = data(:,5:7);
rot_quat = data(:,8:11);
pres = data(:,12);
h = data(:,13);
temp = data(:,14);

%acc = (acc(:,1),acc(:,3),acc(:,2));  rivedi

%% Data processing

data_steps = size(data,1);

h_0 = h(1,1);
h = h - h_0;

x_in = [1,0,0];
z_in = [0,1,0];
y_in = [0,0,1];

% define absolute axis
rotm_in = quat2rotm(rot_quat(1,:));

x_abs = x_in * rotm_in;
y_abs = y_in * rotm_in;
z_abs = z_in * rotm_in;

acc_abs = [0 0 0];
vel_abs = [0 0 0];
pos_abs = [0 0.2 0];
plot(time,h,'-c')
title('Altitude in time')
xlabel('t [s]')
ylabel('h [m]')

% find relative directions, abs acc,abs vel e abs pos
for i = 2:data_steps-1
    rotm = quat2rotm(rot_quat(i,:));
    x_rel = x_in * rotm;
    y_rel = y_in * rotm;
    z_rel = z_in * rotm;
    
    ax_rel = acc(i,1) * x_rel;
    ay_rel = acc(i,2) * y_rel;
    az_rel = acc(i,3) * z_rel;
    
    acc_rel = ax_rel + ay_rel + az_rel;
    
    acc_abs(i,1) = dot(acc_rel,x_abs);
    acc_abs(i,2) = dot(acc_rel,y_abs);
    acc_abs(i,3) = dot(acc_rel,z_abs);
    acc_abs(i,2) = -acc_abs(i,2) - g;
    
    acc_mod(i) = 0;
    for k = 1:3
        acc_mod(i) = acc_mod(i) + acc_abs(i,k)^2;
    end
    acc_mod(i) = sqrt(acc_mod(i));
    
  
    vel_mod(i) = 0;
    delta_t = time(i+1) - time(i);

    for j = 1:3
        vel_abs(i,j) = vel_abs(i-1,j) + 1/2 * (acc_abs(i-1,j) + ...
            acc_abs(i,j)) * delta_t;
        vel_mod(i) = vel_mod(i) + vel_abs(i,j)^2;
        pos_abs(i,j) = pos_abs(i-1,j) + vel_abs(i-1,j) * delta_t + ...
            0.5 * acc_abs(i-1,j) * delta_t^2;
    end
    vel_mod(i) = sqrt(vel_mod(i));
end
    
%% Graphycs everywhere

% figure(1)
% plot(time,h,'-c')
% title('Altitude in time')
% xlabel('t [s]')
% ylabel('h [m]')

figure(2)
plot3(pos_abs(:,1),pos_abs(:,3),pos_abs(:,2),'-r');
grid
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
title('Trajectory')

figure(3)
plot(time(2:data_steps),acc_mod,'-g')
xlabel('t [s]')
ylabel('acc [m/s^2]')
title('Acceleration module in time')

figure(4)
plot(time,acc(:,2),'-k')
xlabel('t [s]')
ylabel('acc_z [m/s^2]')
title('Vertical acceleration in time')

figure(5)
plot(time(2:data_steps),vel_mod,'-y')
xlabel('t [s]')
ylabel('vel [m/s]')
title('Velocity module in time')

figure(6)
plot(time(2:data_steps),vel_abs(:,1),'-k')
title('vx')

figure(7)
plot(time(2:data_steps),vel_abs(:,3),'-k')
title('vy')

figure(8)
plot(time(2:data_steps),vel_abs(:,2),'-k')
title('vz')

figure(9)
plot(time(2:data_steps),pos_abs(:,2),'-k')
title('heigth in time')

figure(10)
plot(time,h,'-c')
title('Altitude in time without abisso')
xlabel('t [s]')
ylabel('h [m]')
ylim([-10 60])