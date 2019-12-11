clear; clc;
'Setting inital conditions'
SERIALcom = false;
MYO = true;
PORT = "COM4";
TAU = .5;
K = 3;
r = 1; %.999; %incremental decrease of velocity
sampling = 1000;
FPS = 50;
TIME = 10;
accel_holder = zeros(3, 2);
s = tf('s');
F = K / (TAU * s + 1);
digital = c2d(F, 1/FPS);
h1 = cell2mat(digital.Numerator);
h2 = cell2mat(digital.Denominator);
filt1 = h1(2);
filt2 = - h2(2);
history = zeros(3, FPS * TIME);
V_history = zeros(3, FPS * TIME);
P_history = zeros(3, FPS * TIME);
accHunfil = zeros(3, FPS * TIME);
%linear postion, linear velocity, rotation
P = [0 0 0]; V = [0 0 0]; accel = [0 0 0]; feedback = [0 0 0];
% T = [0 0 0]; TV = [0 0 0]; A = 0;
I = [0 0 0]; t = 0:1/FPS:1;

% Constants
% Control system gains
Kv1 = 16;
Kv2 = 16;
Kv3 = 16;
Kp1 = 40;
Kp2 = 40;
Kp3 = 40;

%Inertia tensor 2:
I2xx = 0.00053933;
I2xy = -0.00000005;
I2xz = 0.00000397;
I2yx = -0.00000005;
I2yy = 0.00055977;
I2yz = -0.00000903;
I2zx = 0.00000397;
I2zy = -0.00000903;
I2zz = 0.00005433;

%Inertia tensor 3:
I3xx = 0.00078121;
I3xy = -0.00000052;
I3xz = 0.00000095;
I3yx = -0.00000052;
I3yy = 0.00027907;
I3yz = -0.00001532;
I3zx = 0.00000095;
I3zy = -0.00001532;
I3zz = 0.00094348;

d1 = 0.1859;
d2 = 0.0339;
d3 = 0.1422;

th1d = 0;
th2d = 0;
th3d = 0;

%% This opens up serial communication with the arduino
if SERIALcom
    s = serialport(PORT,57600)
end
%% Here the MyoMex is initializated and the gravity acceleration is sampled
if MYO
    mm = MyoMex();
    m = mm.myoData;
    pause(2)
    sample_holder = [0 0 0];
    for i = 0:sampling
        tic
        sample_holder = sample_holder + m.accel_fixed;
        while toc <= 1/FPS
            %waits until the time is the expected time of a frame
            continue;
        end
    end
    change = sample_holder / sampling
    G = change(3)
end
%% The main loop starts here
pause(.5);
'loop begin'

for i = 0 : TIME * FPS
    tic; %start counting the timer
    if MYO
        %% Lowpass filter
        
        P_history = circshift(P_history, -1, 2);
        P_history(:,1) = transpose(P);
        V_history = circshift(V_history, -1, 2);
        V_history(:,1) = transpose(V);
        history = circshift(history, -1, 2);
        history(:,1) = transpose(accel_holder(:,2));
        accel_holder(:,2) = filt1 * accel_holder(:,1) + filt2 * accel_holder(:,2);
        accel_holder(:,1) = transpose(m.accel_fixed);
        accHunfil = circshift(accHunfil, -1, 2);
        accHunfil(:,1) = accel_holder(:,1);
        V = V * r;
        for j = 1:3
            V(j) = V(j) + (accel_holder(j,2) - change(j) - ((j==3) * G))/ (FPS^2);
            P(j) = P(j) + V(j);
            %
            % %             TV(j) = TV(j) + tAccel(j) / FPS;    % TV -> thjd
            % %             T(j) = T(j) + TV(j) / FPS;          % T -> thj
        end
    end
    %% dynamics go here and updat I1 I2 I3 and A
    %% INVERSE KINEMATICS
    
    x = P(1,1);
    y = P(1,2);
    z = P(1,3);
    
    %Theta 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    th1 = atan2(y,x);
    
    %Theta 2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    L1 = 0.2198;
    L2 = 0.2778;
    L4 = sqrt(y^2 + x^2);
    L5 = z-0.0538;
    L3 = sqrt((L4^2)+(L5^2));
    
    phi2 = real(acos((L1^2+L3^2-L2^2)/(2*L1*L3)));
    phi1 = atan2(L4,L5);
    th2 = phi1-phi2;
    
    %Theta 3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    phi3 = real(acos((L1^2+L2^2-L3^2)/(2*L1*L2)));
    th3 = 180-phi3;
    
    %% TORQUES
    
    % Control System Variables
    e1 = th1-feedback(4,1)*0.088*(pi/180);
    e2 = th2-feedback(5,1)*0.088*(pi/180);
    e3 = th3-feedback(6,1)*0.088*(pi/180);
    ed1 = th1d-feedback(1,1);
    ed2 = th2d-feedback(2,1);
    ed3 = th3d-feedback(3,1);
    
    
    %Acceleration vector for control system
    qdd2 = [ed1*Kv1+e1*Kp1;
        ed2*Kv2+e2*Kp2;
        ed3*Kv3+e3*Kp3];
    
    th1dd = qdd2(1,1);
    th2dd = qdd2(2,1);
    th3dd = qdd2(3,1);
    
    tau1 = (-2*th1d*(d1^2*m1 - I2xx + I2yy)*cos(th2)*sin(th2) + (I2xy + I2yx)*cos(th2)^2*th1d - ((I2xy + I2yx)*sin(th2)^2*th1d)/2 + ((-th1d*(I2xy + I2yx)*sin(th2) + th2d*(I2zy + I2yz))*sin(th2))/2 - th2d*(I2xz + I2zx)*cos(th2)/2 - th1d*(d3^2*m2 - I3xx + I3yy)*sin(2*th2 + 2*th3) + th1d*(I3xy + I3yx)*cos(2*th2 + 2*th3) - 2*th1d*d3*m2*(d1 + d2)*sin(th3 + 2*th2) + ((th2d + th3d)*(I3yz + I3zy)*sin(th2 + th3))/2 - ((th2d + th3d)*(I3zx + I3xz)*cos(th2 + th3))/2 - m2*th1d*(d1 + d2)^2*sin(2*th2))*th2d + (-th1d*(d3^2*m2 - I3xx + I3yy)*sin(2*th2 + 2*th3) + th1d*(I3xy + I3yx)*cos(2*th2 + 2*th3) - th1d*d3*m2*(d1 + d2)*sin(th3 + 2*th2) + ((th2d + th3d)*(I3yz + I3zy)*sin(th2 + th3))/2 - ((th2d + th3d)*(I3zx + I3xz)*cos(th2 + th3))/2 - m2*th1d*d3*(d1 + d2)*sin(th3))*th3d + ((d1^2*m1 - I2xx + I2yy)*cos(th2)^2 + (I2xy + I2yx)*sin(th2)*cos(th2) + I2xx + ((d3^2*m2 - I3xx + I3yy)*cos(2*th2 + 2*th3))/2 + ((I3xy + I3yx)*sin(2*th2 + 2*th3))/2 + d3*m2*(d1 + d2)*cos(th3 + 2*th2) + m2*(d1 + d2)^2*cos(2*th2)/2 + m2*d3*(d1 + d2)*cos(th3) + ((d3^2 + (d1 + d2)^2)*m2)/2 + I3xx/2 + I3yy/2)*th1dd + (-((I2zy + I2yz)*cos(th2))/2 - ((I2xz + I2zx)*sin(th2))/2 - ((I3yz + I3zy)*cos(th2 + th3))/2 - ((I3zx + I3xz)*sin(th2 + th3))/2)*th2dd + (-((I3yz + I3zy)*cos(th2 + th3))/2 - ((I3zx + I3xz)*sin(th2 + th3))/2)*th3dd;
    tau2 = (((I2zy + I2yz)*th1d*sin(th2))/2 - th1d*(I2xz + I2zx)*cos(th2)/2 + th1d*(I3yz + I3zy)*sin(th2 + th3)/2 - th1d*(I3zx + I3xz)*cos(th2 + th3)/2)*th2d + (th1d*(I3yz + I3zy)*sin(th2 + th3)/2 - th1d*(I3zx + I3xz)*cos(th2 + th3)/2 - m2*(4*th2d + 2*th3d)*d3*(d1 + d2)*sin(th3)/2)*th3d + (-((I2zy + I2yz)*cos(th2))/2 - ((I2xz + I2zx)*sin(th2))/2 - ((I3yz + I3zy)*cos(th2 + th3))/2 - ((I3zx + I3xz)*sin(th2 + th3))/2)*th1dd + (d1^2*m1 + I2zz + 2*m2*d3*(d1 + d2)*cos(th3) + ((4*d3^2 + 4*(d1 + d2)^2)*m2)/4 + I3zz)*th2dd + (m2*d3*(d1 + d2)*cos(th3) + d3^2*m2 + I3zz)*th3dd + th1d^2*(d1^2*m1 - I2xx + I2yy)*cos(th2)*sin(th2) - th1d^2*(I2xy + I2yx)*cos(th2)^2/2 - ((-th1d*(I2xy + I2yx)*sin(th2) + th2d*(I2zy + I2yz))*th1d*sin(th2))/2 + th1d*th2d*(I2xz + I2zx)*cos(th2)/2 - m1*G2*cos(th2)*d1 + th1d^2*(d3^2*m2 - I3xx + I3yy)*sin(2*th2 + 2*th3)/2 - th1d^2*(I3xy + I3yx)*cos(2*th2 + 2*th3)/2 + th1d^2*d3*m2*(d1 + d2)*sin(th3 + 2*th2) - th1d*(th2d + th3d)*(I3yz + I3zy)*sin(th2 + th3)/2 + th1d*(th2d + th3d)*(I3zx + I3xz)*cos(th2 + th3)/2 + m2*th1d^2*(d1 + d2)^2*sin(2*th2)/2 + m2*G3*(-cos(th2)*(d1 + d2) - cos(th2 + th3)*d3);
    tau3 = (th1d*(I3yz + I3zy)*sin(th2 + th3)/2 - th1d*(I3zx + I3xz)*cos(th2 + th3)/2)*th2d + (th1d*(I3yz + I3zy)*sin(th2 + th3)/2 - th1d*(I3zx + I3xz)*cos(th2 + th3)/2 - m2*th2d*d3*(d1 + d2)*sin(th3))*th3d + (-((I3yz + I3zy)*cos(th2 + th3))/2 - ((I3zx + I3xz)*sin(th2 + th3))/2)*th1dd + (m2*d3*(d1 + d2)*cos(th3) + d3^2*m2 + I3zz)*th2dd + (d3^2*m2 + I3zz)*th3dd + th1d^2*(d3^2*m2 - I3xx + I3yy)*sin(2*th2 + 2*th3)/2 - th1d^2*(I3xy + I3yx)*cos(2*th2 + 2*th3)/2 + th1d^2*d3*m2*(d1 + d2)*sin(th3 + 2*th2)/2 - th1d*(th2d + th3d)*(I3yz + I3zy)*sin(th2 + th3)/2 + th1d*(th2d + th3d)*(I3zx + I3xz)*cos(th2 + th3)/2 + m2*(th1d^2 + 2*th2d*(th2d + th3d))*d3*(d1 + d2)*sin(th3)/2 - m2*G3*cos(th2 + th3)*d3;
    
    Torq = [tau1; tau2; tau3]
    
    %MX64 conversion: Joint 1
    amps1 = 0.871*Torq(1,1)+0.075;
    current1 = amps1/0.00336;
    
    %MX106 conversion: Joint 2
    amps2 = 0.006*tau2^4-0.067*tau2^3+0.33*tau2^2-0.192*tau2+0.726;
    current2 = amps2/0.00336;
    
    %MX64 conversion: Joint 3
    amps3 = 0.871*Torq(3,1)+0.075;
    current3 = amps3/0.00336;
    %% Verboise
    if ~rem(i, FPS) & MYO
        [i/FPS TIME]
    end
    %% Transmission
    if SERIALcom
        write(s, int16(5000), "int16");     %first WORD of the transmission
        write(s, int16(current1), "int16");       %Current of motor 1
        write(s, int16(current2), "int16");       %Current of motor 2
        write(s, int16(current3), "int16");       %Current of motor 3
        write(s, int16(A), "int16");        %Angle of griper(simetrical)
        feedback = read(s, 6, 'int16')      %theata and theatad of motor 1, 2 and 3
    end
    while toc <= 1/FPS
        %waits until the time is the expected time of a frame
        continue;
    end
end
%% Clear everything before ending the program.
if MYO
    mm.delete();
end
%% Plot everythin
a = true;
ua = false;
v = true;
p = true;
plot_t = 0 : (1 / FPS): TIME - (1 / FPS);
multiple = 10;
plot_accel = history(1,:) * multiple * a;
plot_uAccel = accHunfil(1,:) * multiple * ua;
plot_vel = V_history(1,:) * multiple * FPS * v;
plot_posi = P_history(1,:) * multiple * p;
figure
grid on
grid minor
plot(plot_t, plot_accel, plot_t, plot_uAccel, plot_t, plot_vel, plot_t, plot_posi);
% xlim([0 TIME])
% ylim([-15 15])
% clear;


