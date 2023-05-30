%% Initialization
clear all;clc;close all

%% Input parameters

%us acc 123 giro 123 magn 123 rotaz 1234 p h t 

%filename = 'dataLog_m_numerogeneratoaacaso.txt';
test_2_static = 34047-500;
test_2_s = 34047;
test_2_f = 35905+100;
static_start = test_2_static;
flight_start = test_2_s;
flight_end = test_2_f;
offset_ln = flight_start - static_start;
flight_ln = flight_end - flight_start;


filename = 'test_2.TXT';
data_complete = readmatrix(filename);
data = data_complete(static_start:flight_end,:);

time = data(:,1)/(10^6);
acc_raw = data(:,2:4);
acc = acc_raw;
giro = data(:,5:7);
rot_quat_incomp = data(:,8:10);
pres = data(:,11);
h = data(:,12);
temp = data(:,13);

for i = 1:size(rot_quat_incomp,1)
    q_1 = rot_quat_incomp(i,1);
    q_2 = rot_quat_incomp(i,2);
    q_3 = rot_quat_incomp(i,3);
    q_0 = real(sqrt(1-(q_1^2+q_2^2+q_3^2)));
    rot_quat(i,:) = [q_0 q_1 q_2 q_3];
end

%acc = (acc(:,1),acc(:,3),acc(:,2));  rivedi

%% Data processing

%% define g, ax, ay, az while stopped

acc_raw_mod = zeros(1,offset_ln);

for i = 1:offset_ln-1
    for j = 1:3
        acc_raw_mod(1,i) = acc_raw_mod(1,i) + acc_raw(i,j)^2;
    end
    acc_raw_mod(1,i) = sqrt(acc_raw_mod(1,i));
        
end
g_calib = mean(acc_raw_mod);

ax_s = mean(acc_raw(1:offset_ln, 1));
ay_s = mean(acc_raw(1:offset_ln, 2));
az_s = mean(acc_raw(1:offset_ln, 3));

%% define absolute axis

x_eul = [1,0,0];
y_eul = [0,1,0];
z_eul = [0,0,1];

rotm_in = zeros(1,offset_ln);

for i = 1:offset_ln-1
    
    rotm = quat2rotm(rot_quat(i,:));
    x_abs(i,:) = x_eul * rotm;
    y_abs(i,:) = y_eul * rotm;
    z_abs(i,:) = z_eul * rotm;
    
    gx_abs(i,:) = acc(i,1) * x_abs(i,:);
    gy_abs(i,:) = acc(i,2) * y_abs(i,:);
    gz_abs(i,:) = acc(i,3) * z_abs(i,:);
    
end

gx_avg = mean(gx_abs,1,'omitnan');
gy_avg = mean(gy_abs,1,'omitnan');
gz_avg = mean(gz_abs,1,'omitnan');

g = gx_avg + gy_avg + gz_avg;
g_mod = norm(g);
g_dir = g/g_mod;

z_cart = g_dir;
x_cart = x_eul - dot(x_eul, z_cart) * z_cart;
x_cart = x_cart/norm(x_cart);
y_cart = y_eul - dot(y_eul, z_cart) * z_cart - dot(y_eul, x_cart) * x_cart;
y_cart = y_cart/norm(y_cart);

%% old
% 
% rotm_in = quat2rotm(rot_quat(1,:));
% 
% x_abs = x_cart * rotm_in;
% y_abs = y_in * rotm_in;
% z_abs = z_cart * rotm_in;

%% find relative directions, abs acc,abs vel e abs pos

data = data_complete(flight_start:flight_end,:);

time = data(:,1)/(10^6);
acc_raw = data(:,2:4);
acc = acc_raw;
giro = data(:,5:7);
rot_quat_incomp = data(:,8:10);
pres = data(:,11);
h = data(:,12);
temp = data(:,13);

for i = 1:size(rot_quat_incomp,1)
    q_1 = rot_quat_incomp(i,1);
    q_2 = rot_quat_incomp(i,2);
    q_3 = rot_quat_incomp(i,3);
    q_0 = real(sqrt(1-(q_1^2+q_2^2+q_3^2)));
    rot_quat(i,:) = [q_0 q_1 q_2 q_3];
end


data_steps = size(data,1);

%% pressure analysis

h_0 = h(1,1);
h = h - h_0;

%%

for i = 1:data_steps
    rotm = quat2rotm(rot_quat(i,:));
    x_ist = x_eul * rotm;
    y_ist = y_eul * rotm;
    z_ist = z_eul * rotm;
    
    ax_ist = acc(i,1) * x_ist;
    ay_ist = acc(i,2) * y_ist;
    az_ist = acc(i,3) * z_ist;
    
    acc_ist = ax_ist + ay_ist + az_ist;
    
    acc_cart(i,1) = dot(acc_ist,x_cart);
    acc_cart(i,2) = dot(acc_ist,y_cart);
    acc_cart(i,3) = dot(acc_ist,z_cart);
    
    %%test
    acc_mod_corrected(i) = 0;
    for k = 1:3
        acc_mod_corrected(i) = acc_mod_corrected(i) + acc_cart(i,k)^2;
    end
    acc_mod_corrected(i) = sqrt(acc_mod_corrected(i));
    
    acc_cart(i,3) = acc_cart(i,3) - g_calib;
    
    acc_mod(i)=0;
    for k = 1:3
        acc_mod(i) = acc_mod(i) + acc_cart(i,k)^2;
    end
    acc_mod(i) = sqrt(acc_mod(i));
    
end

vel_cart(i,:) = [0 0 0];
pos_cart(i,:) = [0 0 0.2];

for i = 1:data_steps-1
    
    delta_t = time(i+1) - time(i);
    
    vel_mod(i)=0;
    for j = 1:3
        vel_cart(i+1,j) = vel_cart(i,j) + 1/2 * (acc_cart(i,j) + ...
            acc_cart(i+1,j)) * delta_t;
        vel_mod(1,i) = vel_mod(1,i) + vel_cart(i,j)^2;
        pos_cart(i+1,j) = pos_cart(i,j) + vel_cart(i,j) * delta_t + ...
            0.5 * acc_cart(i,j) * delta_t^2;
    end
    vel_mod(1,i) = sqrt(vel_mod(1,i));
end
    
%% Graphycs everywhere

% figure(1)
% plot(time,h,'-c')
% title('Altitude in time')
% xlabel('t [s]')
% ylabel('h [m]')

figure(2)
plot3(pos_cart(:,1),pos_cart(:,2),pos_cart(:,3),'-r');
grid
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
title('Trajectory')

figure(3)
plot(time,acc_mod,'-g')
xlabel('t [s]')
ylabel('acc [m/s^2]')
title('Acceleration module in time')

figure(4)
plot(time,acc_cart(:,3),'-k')
xlabel('t [s]')
ylabel('acc_z [m/s^2]')
title('Vertical acceleration in time')

figure(5)
plot(time(1:data_steps-1),vel_mod,'-y')
xlabel('t [s]')
ylabel('vel [m/s]')
title('Velocity module in time')

figure(6)
plot(time,vel_cart(:,1),'-k')
title('vx')

figure(7)
plot(time,vel_cart(:,2),'-k')
title('vy')

figure(8)
plot(time,vel_cart(:,3),'-k')
title('vz')

figure(9)
plot(time,pos_cart(:,3),'-k')
title('heigth in time')

figure(10)
plot(time,h,'-c')
title('Pressure ltitude in time')
xlabel('t [s]')
ylabel('h [m]')

figure(11)
plot(time,acc_mod_corrected,'-r')
title('Accelerazione corretta?')

figure(12)
plot(time,acc_cart(:,1),'-k')
xlabel('t [s]')
ylabel('acc_x [m/s^2]')
title('ax in time')

figure(13)
plot(time,acc_cart(:,2),'-k')
xlabel('t [s]')
ylabel('acc_y [m/s^2]')
title('ay in time')