%% SMPC Rocket Landing - Cleaned Code
% When trying to run the matlab code, it will throw some errors because of our use of global variables.  
% 
% To run the code,  first run the whole document.  It will throw an error.  Then run subsection called “Step”.  This subsection has all of the global variables in it.  Run this section.   It will throw an error.  Then run the full document.  It will finally have the global variables registered and will be able to run.
% 
% This bug is a problem with matlabs global variables.  They recommend not using them but we found ourselves in the position where we needed to use them.  

%% Constants
clear;
clc;
m = 100000; %kg -> can be modified
g = [0;0;9.81];
v = [0;0;10]; %Initially v=0 m/s
p = 1.2; %kg/m^3
C_D = .75; %arbitrary coefficient of drag
D = 3.667; %m -> Diameter of Rocket
H = 70; %m -> Height of Rocket
C_M = 30; %Center of mass, for now it will be 30 m

%% Drag Induced Force and Moment Function

function [F_D,M_D] = Drag_Force_Moment(psi, theta, phi, v)
    C_M = 30; %Center of mass, for now it will be 30 m
    D = 3.667; %m -> Diameter of Rocket
    H = 70; %m -> Height of Rocket
    p = 1.2; %kg/m^3
    C_D = .75; %arbitrary coefficient of drag

    R_yaw = [cos(psi), -sin(psi), 0; 
        sin(psi), cos(psi), 0; 
        0,0,1];

    R_pitch = [cos(theta), 0, sin(theta);
        0,1,0;
        -sin(theta),0,cos(theta)];

    R_roll = [1,0,0;
        0, cos(phi), -sin(phi);
        0, sin(phi), cos(phi)];
    
    R = R_yaw*R_pitch*R_roll;
    
    l_0 = [0;0;1];
    l = R*l_0;
    
    alpha = acos(dot(v/norm(v),l));
    
    A = pi*D*H*abs(sin(alpha)) + pi/4*D^2*abs(cos(alpha));
    
    F_D = .5*p*C_D*A*v.^2;
    M_D = cross(l*(35-C_M)*sin(alpha),F_D);
end

%% Thrust Induced Force and Moment
%e1,e2,and F are control variables

function [F_T, M_T] = Thrust_Force(psi, theta, phi, e1, e2, F)
    C_M = 30; %Center of mass, for now it will be 30 m

    R_yaw = [cos(psi), -sin(psi), 0; 
        sin(psi), cos(psi), 0; 
        0,0,1];
    R_pitch = [cos(theta), 0, sin(theta);
        0,1,0;
        -sin(theta),0,cos(theta)];
    R_roll = [1,0,0;
        0, cos(phi), -sin(phi);
        0, sin(phi), cos(phi)];
    R_e1 = [cos(e1), 0, sin(e1);
        0,1,0;
        -sin(e1), 0, cos(e1)];
    R_e2 = [1,0,0;
        0, cos(e2), -sin(e2);
        0, sin(e2), cos(e2)];

    R_old = R_yaw*R_pitch*R_roll;
    l = R_old*[0;0;1];

    R = R_pitch*R_e1*R_yaw*R_roll*R_e2;
    
    l_Thrust = R*[0;0;1];
    
    F_T = F.*l_Thrust;
    M_T = cross(l*C_M*sin(acos(dot(l,l_Thrust))),F_T);
end
%% Updated Acceleration Function

function [sigma_dd, a] = Acceleration_Update(psi, theta, phi, e1, e2, F, x_d)
    m = 100000; %kg -> can be modified
    g = [0;0;9.81];
    D = 3.667; %m -> Diameter of Rocket
    H = 70; %m -> Height of Rocket

    r = D/2; %Radius of Rocket
    I = [1/12*m*(3*r^2 + H^2),0,0;
        0, 1/12*m*(3*r^2 + H^2), 0;
        0,0,1/2*m*r^2];

    [F_T, M_T] = Thrust_Force(psi, theta, phi, e1, e2, F);
    [F_D, M_D] = Drag_Force_Moment(psi, theta, phi, x_d);
    
    sigma_dd = inv(I)*(M_T - M_D); %[psi, theta, phi]
    
    ma = -m*g + F_D + F_T;
    a = ma/m;
end


%% Symbolic Acceleration Function
%Generating Functions using this method results in exponentially faster
%computation as subs() does not need to be used

syms psiv thetav phiv e1 e2 Fv x_dx x_dy x_dz

x_dv = [x_dx;x_dy;x_dz]; %Velocity Vector.  The v term means nothing

[sigma_dd, a] = Acceleration_Update(psiv, thetav, phiv, e1, e2, Fv, x_dv);

sigma1_dd_func = matlabFunction(sigma_dd(1));
sigma2_dd_func = matlabFunction(sigma_dd(2));
sigma3_dd_func = matlabFunction(sigma_dd(3));

a1_func = matlabFunction(a(1));
a2_func = matlabFunction(a(2));
a3_func = matlabFunction(a(3));

%% Expedited Acceleration Update (with Control Capabilities)
%Control variables are passed as symbolic values.  All other values are
%passed as values that they are at the current state

function [sigma_dd, a] = Controlled_Acceleration_Update(psiv, thetav, phiv, e1, e2, Fv, x_d, sigma1_dd_func, sigma2_dd_func, sigma3_dd_func, a1_func, a2_func, a3_func)
    
    x_dx = x_d(1); x_dy = x_d(2); x_dz = x_d(3);

    sigma_dd = [sigma1_dd_func(Fv,e1,e2,phiv,psiv,thetav,x_dx,x_dy,x_dz);
    sigma2_dd_func(Fv,e1,e2,phiv,psiv,thetav,x_dx,x_dy,x_dz);
    sigma3_dd_func(Fv,e1,e2,phiv,psiv,thetav,x_dx,x_dy,x_dz)];
    
    a=[a1_func(Fv,e1,e2,phiv,psiv,thetav,x_dx,x_dy,x_dz);
    a2_func(Fv,e2,phiv,psiv,thetav,x_dx,x_dy,x_dz);
    a3_func(Fv,e1,e2,phiv,psiv,thetav,x_dx,x_dy,x_dz)];
end

%% Next Control Step

function [u_F, u_e1, u_e2] = Control_Step_Value(xdes, x_d, x, sigma_d, sigma, F_prev, e1_prev, e2_prev, sigma1_dd_func, sigma2_dd_func, sigma3_dd_func, a1_func, a2_func, a3_func)
    
    F_Temp = sym('F_Temp');
    e1_Temp = sym('e1_Temp');
    e2_Temp = sym('e2_Temp');
   
    dt = .1;
    psi = sigma(1); theta = sigma(2); phi = sigma(3);

    sig_des = [0;0;0];

    [sigma_dd, a] = Controlled_Acceleration_Update(psi, theta, phi, e1_Temp(1), e2_Temp(1), F_Temp(:,1), x_d, sigma1_dd_func, sigma2_dd_func, sigma3_dd_func, a1_func, a2_func, a3_func);
    
    %this sigma_dd and a specified at the current values, will be used for
    %the next time steps.  The only value in these functions that will
    %change are the symbolic control variables.

    s_dd_func = matlabFunction(sigma_dd); %s_dd_func(F_Temp,e1_Temp,e2_Temp)
    acc_func = matlabFunction(a); %acc_func(F_Temp,e1_Temp,e2_Temp)

    F_Temp = sym('F_Temp', [1,5]);
    e1_Temp = sym('e1_Temp', [1,5]);
    e2_Temp = sym('e2_Temp', [1,5]);

    
    x_d1 = x_d + dt*acc_func(F_Temp(1), e1_Temp(1), e2_Temp(1));
    x1 = x + dt*x_d1; 
    sigma_d1 = sigma_d + dt*s_dd_func(F_Temp(1), e1_Temp(1), e2_Temp(1));
    sigma1 = sigma + dt*sigma_d1;

    x_d2 = x_d1 + dt*acc_func(F_Temp(2), e1_Temp(2), e2_Temp(2));
    x2 = x1 + dt*x_d2;
    sigma_d2 = sigma_d1 + dt*s_dd_func(F_Temp(2), e1_Temp(2), e2_Temp(2));
    sigma2 = sigma1 + dt*sigma_d2;
    
    x_d3 = x_d2 + dt*acc_func(F_Temp(3), e1_Temp(3), e2_Temp(3));
    x3 = x2 + dt*x_d3;
    sigma_d3 = sigma_d2 + dt*s_dd_func(F_Temp(3), e1_Temp(3), e2_Temp(3));
    sigma3 = sigma2 + dt*sigma_d3;

    x_d4 = x_d3 + dt*acc_func(F_Temp(4), e1_Temp(4), e2_Temp(4));
    x4 = x3 + dt*x_d4;
    sigma_d4 = sigma_d3 + dt*s_dd_func(F_Temp(4), e1_Temp(4), e2_Temp(4));
    sigma4 = sigma3 + dt*sigma_d4;
    
    x_d5 = x_d4 + dt*acc_func(F_Temp(5), e1_Temp(5), e2_Temp(5));
    x5 = x4 + dt*x_d5;
    sigma_d5 = sigma_d4 + dt*s_dd_func(F_Temp(5), e1_Temp(5), e2_Temp(5));
    sigma5 = sigma4 + dt*sigma_d5;

    Traj_Weight = 500;
    Angle_Weight = 1;
    F_var_weight = 10^-7;%10^-7;
    e_var_weight = .1;%0^8;
    
    J = Traj_Weight*sum((x1 - xdes(:,1)).^2 + (x2 - xdes(:,2)).^2 + (x3 - xdes(:,3)).^2 + (x4 - xdes(:,4)).^2 + (x5 - xdes(:,5)).^2); %+ F_var_weight*(sum((F_Temp(1:(end-1))-F_Temp(2:end)).^2)) + e_var_weight*(sum((e1_Temp(2:end) - e1_Temp(1:(end-1))).^2 + (e1_prev-e1_Temp(1))^2 + (e2_Temp(2:end) - e2_Temp(1:(end-1))).^2 + (e2_prev-e2_Temp(1))^2)) + sum(e1_Temp(:).^20) + sum(e2_Temp(:).^20); %... + Angle_Weight*sum((sigma1 - sig_des).^2 + (sigma2 - sig_des).^2 + (sigma3 - sig_des).^2 + (sigma4 - sig_des).^2 + (sigma5 - sig_des).^2)
    %J = J + Angle_Weight*sum((sigma1 - sig_des).^2 + (sigma2 - sig_des).^2 + (sigma3 - sig_des).^2 + (sigma4 - sig_des).^2 + (sigma5 - sig_des).^2);
    J = J + F_var_weight*(sum((F_Temp(1:(end-1))-F_Temp(2:end)).^2) + (F_prev - F_Temp(1)));
    J = J + e_var_weight*(sum((e1_Temp(2:end) - e1_Temp(1:(end-1))).^2 + (e1_prev-e1_Temp(1))^2 + (e2_Temp(2:end) - e2_Temp(1:(end-1))).^2 + (e2_prev-e2_Temp(1))^2));
    %J = J + sum(e1_Temp(:).^20) + sum(e2_Temp(:).^20); 
    
    J_func = matlabFunction(J);
    J_grad = gradient(J, [F_Temp(:); e1_Temp(:); e2_Temp(:)]);
   

    J_grad_Func = matlabFunction(J_grad); %J_grad_Func(F_Temp1,F_Temp2,F_Temp3,F_Temp4,F_Temp5,e1_Temp1,e1_Temp2,e1_Temp3,e1_Temp4,e1_Temp5,e2_Temp1,e2_Temp2,e2_Temp3,e2_Temp4,e2_Temp5)

    J_Func = matlabFunction(J);
    d = 1;
    u = [F_prev*ones(1,5), e1_prev*ones(1,5), e2_prev*ones(1,5)];
    
    % figure
    % 
    % x1_out = subs(x1, [F_Temp(1), F_Temp(2), F_Temp(3), F_Temp(4), F_Temp(5), e1_Temp(1), e1_Temp(2), e1_Temp(3), e1_Temp(4), e1_Temp(5), e2_Temp(1), e2_Temp(2), e2_Temp(3), e2_Temp(4), e2_Temp(5)], [u(1), u(2), u(3), u(4), u(5), u(6), u(7), u(8), u(9), u(10), u(11), u(12), u(13), u(14), u(15)]);
    % x2_out = subs(x2, [F_Temp(1), F_Temp(2), F_Temp(3), F_Temp(4), F_Temp(5), e1_Temp(1), e1_Temp(2), e1_Temp(3), e1_Temp(4), e1_Temp(5), e2_Temp(1), e2_Temp(2), e2_Temp(3), e2_Temp(4), e2_Temp(5)], [u(1), u(2), u(3), u(4), u(5), u(6), u(7), u(8), u(9), u(10), u(11), u(12), u(13), u(14), u(15)]);
    % x3_out = subs(x3, [F_Temp(1), F_Temp(2), F_Temp(3), F_Temp(4), F_Temp(5), e1_Temp(1), e1_Temp(2), e1_Temp(3), e1_Temp(4), e1_Temp(5), e2_Temp(1), e2_Temp(2), e2_Temp(3), e2_Temp(4), e2_Temp(5)], [u(1), u(2), u(3), u(4), u(5), u(6), u(7), u(8), u(9), u(10), u(11), u(12), u(13), u(14), u(15)]);
    % x4_out = subs(x4, [F_Temp(1), F_Temp(2), F_Temp(3), F_Temp(4), F_Temp(5), e1_Temp(1), e1_Temp(2), e1_Temp(3), e1_Temp(4), e1_Temp(5), e2_Temp(1), e2_Temp(2), e2_Temp(3), e2_Temp(4), e2_Temp(5)], [u(1), u(2), u(3), u(4), u(5), u(6), u(7), u(8), u(9), u(10), u(11), u(12), u(13), u(14), u(15)]);
    % x5_out = subs(x5, [F_Temp(1), F_Temp(2), F_Temp(3), F_Temp(4), F_Temp(5), e1_Temp(1), e1_Temp(2), e1_Temp(3), e1_Temp(4), e1_Temp(5), e2_Temp(1), e2_Temp(2), e2_Temp(3), e2_Temp(4), e2_Temp(5)], [u(1), u(2), u(3), u(4), u(5), u(6), u(7), u(8), u(9), u(10), u(11), u(12), u(13), u(14), u(15)]);
    % x_out_initial = double([x1_out, x2_out, x3_out, x4_out, x5_out]);
    % plot3(x_out_initial(1,:), x_out_initial(2,:), x_out_initial(3,:))
    % hold on
    % xlabel("X")
    % ylabel("Y")
    % zlabel("Z (Altitude)")
    % %title("Controlled Trajectory Regression vs Origional Trajectory")
    % grid on
    
    

    alpha = [10^6;10^6;10^6;10^6;10^6;1;1;1;1;1;1;1;1;1;1];
    new_gradient = zeros(1,10);
    for i = 1:10000
        if(norm(d,2)<1)
            break;
        end
        %fprintf("Iteration %f | Fx1 = %f | Fx2 = %f | Fx3 = %f | Fx4 = %f | Fx5 = %f| e1_1 = %f | e1_2 = %f | e1_3 = %f | e1_4 = %f | e1_5 = %f | e2_1 = %f | e2_2 = %f | e2_3 = %f | e2_4 = %f | e2_5 = %f \n", [i, u(1), u(2), u(3), u(4), u(5), u(6), u(7), u(8), u(9), u(10), u(11), u(12), u(13), u(14), u(15)])
         failed_count_e1 = [0,0,0,0,0];
         failed_count_e2 = [0,0,0,0,0];
         for i = 1:1000
             e1_plus_noise = u(6:10) + randn(1,5)*.1;
             e2_plus_noise = u(11:15) + randn(1,5)*.1;
             failed_count_e1 = failed_count_e1 + (abs(e1_plus_noise) > pi/2);
             failed_count_e2 = failed_count_e2 + (abs(e2_plus_noise) > pi/2);
         end
         Pf_e1 = failed_count_e1/1000;
         Pf_e2 = failed_count_e2/1000;

         new_gradient(1:5) = (Pf_e1 > .2).*u(6:10).^19;
         new_gradient(6:10) = (Pf_e2 > .2).*u(11:15).^19;


        grad = J_grad_Func(u(1), u(2), u(3), u(4), u(5), u(6), u(7), u(8), u(9), u(10), u(11), u(12), u(13), u(14), u(15));
        grad(6:15) = grad(6:15) + new_gradient';
        normalized_grad = [grad(1:5,1); grad(6:15,1)/norm(grad(6:15,1),2)];
        d = -normalized_grad.*alpha;
        u = u + d';
       
        %Problem is that I am using a force vector for control, when it
        %should just be force magnitude, and the force vector will be
        %determined from orientation and e values.
        %J_final = J_Func(u(1), u(2), u(3), u(4), u(5), u(6), u(7), u(8), u(9), u(10), u(11), u(12), u(13), u(14), u(15));
        

    end
     u_F = u(1);
     u_e1 = u(6);
     u_e2 = u(11);

     
    % %double(subs(x1, [F_Temp(1), e1_Temp(1), e2_Temp(1)], [, u_e1, u_e2]))
    % 
    % plot3(xdes(1,:), xdes(2,:), xdes(3,:))
    % x1_out = subs(x1, [F_Temp(1), F_Temp(2), F_Temp(3), F_Temp(4), F_Temp(5), e1_Temp(1), e1_Temp(2), e1_Temp(3), e1_Temp(4), e1_Temp(5), e2_Temp(1), e2_Temp(2), e2_Temp(3), e2_Temp(4), e2_Temp(5)], [u(1), u(2), u(3), u(4), u(5), u(6), u(7), u(8), u(9), u(10), u(11), u(12), u(13), u(14), u(15)]);
    % x2_out = subs(x2, [F_Temp(1), F_Temp(2), F_Temp(3), F_Temp(4), F_Temp(5), e1_Temp(1), e1_Temp(2), e1_Temp(3), e1_Temp(4), e1_Temp(5), e2_Temp(1), e2_Temp(2), e2_Temp(3), e2_Temp(4), e2_Temp(5)], [u(1), u(2), u(3), u(4), u(5), u(6), u(7), u(8), u(9), u(10), u(11), u(12), u(13), u(14), u(15)]);
    % x3_out = subs(x3, [F_Temp(1), F_Temp(2), F_Temp(3), F_Temp(4), F_Temp(5), e1_Temp(1), e1_Temp(2), e1_Temp(3), e1_Temp(4), e1_Temp(5), e2_Temp(1), e2_Temp(2), e2_Temp(3), e2_Temp(4), e2_Temp(5)], [u(1), u(2), u(3), u(4), u(5), u(6), u(7), u(8), u(9), u(10), u(11), u(12), u(13), u(14), u(15)]);
    % x4_out = subs(x4, [F_Temp(1), F_Temp(2), F_Temp(3), F_Temp(4), F_Temp(5), e1_Temp(1), e1_Temp(2), e1_Temp(3), e1_Temp(4), e1_Temp(5), e2_Temp(1), e2_Temp(2), e2_Temp(3), e2_Temp(4), e2_Temp(5)], [u(1), u(2), u(3), u(4), u(5), u(6), u(7), u(8), u(9), u(10), u(11), u(12), u(13), u(14), u(15)]);
    % x5_out = subs(x5, [F_Temp(1), F_Temp(2), F_Temp(3), F_Temp(4), F_Temp(5), e1_Temp(1), e1_Temp(2), e1_Temp(3), e1_Temp(4), e1_Temp(5), e2_Temp(1), e2_Temp(2), e2_Temp(3), e2_Temp(4), e2_Temp(5)], [u(1), u(2), u(3), u(4), u(5), u(6), u(7), u(8), u(9), u(10), u(11), u(12), u(13), u(14), u(15)]);
    % x_out = double([x1_out, x2_out, x3_out, x4_out, x5_out]);
    % plot3(x_out(1,:), x_out(2,:), x_out(3,:))
    % hold off
    % MSE_Initial = sum(sum((x_out_initial - xdes).^2)/5);
    % MSE_Final = sum(sum((x_out - xdes).^2)/5);
    % fprintf("MSE Initially: %f \n", MSE_Initial);
    % fprintf("MSE After Optimal u is found: %f \n", MSE_Final)
    % legend("Origional Controlled Trajectory", "Desired Trajectory", "Optimally Controlled Trajectory")
    % J_grad_Func(u(1), u(2), u(3), u(4), u(5), u(6), u(7), u(8), u(9), u(10), u(11), u(12), u(13), u(14), u(15));
end

%% Desired Trajectory

%Trajectory
iter = 1000;
t = linspace(0,100,iter);

x_Traj = linspace(0,100,iter);
y_Traj = linspace(0,100,iter);
z_Traj = 1000 + 10*x_Traj;%1020.*1./((1+exp((x_Traj/10-5))).*(1+exp((y_Traj/10-5))));
 figure;
 scatter3(x_Traj,y_Traj,z_Traj,20, t,'filled');
 colorbar;
%% Uncontrolled Dynamics
x_d = [1;0;10]; %Initial Velocity
x = [0;0;1000]; %Initial Positon w.r.t Ground
sigma_d = [0;0;0]; %intial angle positions [theta, psi, phi]
sigma = [0.01;.01;0.01]; %Initial change in angle positions


y_net = [];
t_net = [];
for i = 1:25
    y0 = zeros(15,1);
    y0(1:3) = x(:);
    y0(4:6) = x_d(:);
    y0(7:9) = sigma(:);
    y0(10:12) = sigma_d(:);
    y0(13:15) = [0;0;0];
    
    t_span = [0,.1];
    
    [t,y] = ode45(@system_of_eq, t_span, y0);
   
    %fprintf("Uncontrolled Iteration %f\n", i)
    %fprintf("x1 = %f | x2 = %f | x3 = %f | x_d1 = %f | x_d2 = %f | x_d3 = %f | s1 = %f | s2 = %f | s3 = %f | s_d1 = %f | sd_2 = %f | s_d3 = %f\n", [y(end,1),y(end,2),y(end,3), y(end,4),y(end,5),y(end,6), y(end,7),y(end,8),y(end,9), y(end,10),y(end,11),y(end,12)])
   
   
    x = [y(end,1);y(end,2);y(end,3)];
    x_d = [y(end,4);y(end,5);y(end,6)];
    sigma = [y(end,7);y(end,8);y(end,9)];
    sigma_d = [y(end,10);y(end,11);y(end,12)];

    y_net = [y_net; y];
    t_net = [t_net; t+.1*(i-1)];
end


%scatter3(y_net(:,1),y_net(:,2),y_net(:,3), 30, t_net, 'filled')
%colorbar;

 figure
 plot3(y_net(:,1),y_net(:,2),y_net(:,3))
% 
 hold on


%% Testing
%% Initial Condition
dt = .1;
x_d = [1;0;10]; %Initial Velocity
x = [0;0;1000]; %Initial Positon w.r.t Ground
sigma_d = [0;0;0]; %intial angle positions [theta, psi, phi]
sigma = [0.01;.01;0.01]; %Initial change in angle positions
F_prev = 1000000;
e1_prev = 0;
e2_prev = 0;

%% Control Step




%% Step
%Run Everything
%Run this Section
%Run Everything
e1_collected = [];
e2_collected = [];


%RK45
% 
% x_d1 = x_d + dt*acc_func(u_F, u_e1, u_e2);
% x1 = x + dt*x_d1; 
% sigma_d1 = sigma_d + dt*s_dd_func(u_F, u_e1, u_e2);
% sigma1 = sigma + dt*sigma_d1;


 dt = .1;
 
global s1_dd_func
 s1_dd_func = sigma1_dd_func;
global s2_dd_func
s2_dd_func = sigma2_dd_func;
global s3_dd_func
s3_dd_func = sigma3_dd_func;
global a1_function
a1_function = a1_func;
global a2_function
a2_function = a2_func;
global a3_function
a3_function = a3_func;


 function dydt = system_of_eq(t,y)
    global s1_dd_func
    global s2_dd_func
    global s3_dd_func
    global a1_function
    global a2_function
    global a3_function
    sigma1_dd_func = s1_dd_func;
    sigma2_dd_func = s2_dd_func;
    sigma3_dd_func = s3_dd_func;
    a1_func = a1_function;
    a2_func = a2_function;
    a3_func = a3_function;
    
    x = y(1:3);
    x_d = y(4:6);
    sigma = y(7:9);
    sigma_d = y(10:12);
    psi = sigma(1); theta = sigma(2); phi = sigma(3);
    u = y(13:15);
    u_F = u(1); u_e1 = u(2); u_e2 = u(3);
    [sigma_dd, a] = Controlled_Acceleration_Update(psi, theta, phi, u_e1, u_e2, u_F, x_d, sigma1_dd_func, sigma2_dd_func, sigma3_dd_func, a1_func, a2_func, a3_func);
    dydt = zeros(15,1);
    dydt(1:3) = x_d(:);
    dydt(4:6) = a(:);
    dydt(7:9) = sigma_d(:);
    dydt(10:12) = sigma_dd(:);

 end
y_net_ctrl = [];
t_net_ctrl = [];
for i = 1:25
    
    xdes = [x_Traj(i),y_Traj(i),z_Traj(i); x_Traj(i+1),y_Traj(i+1),z_Traj(i+1); x_Traj(i+2),y_Traj(i+2),z_Traj(i+2); x_Traj(i+3),y_Traj(i+3),z_Traj(i+3); x_Traj(i+4),y_Traj(i+4),z_Traj(i+4)]';
    [u_F, u_e1, u_e2] = Control_Step_Value(xdes, x_d, x, sigma_d, sigma, F_prev, e1_prev, e2_prev, sigma1_dd_func, sigma2_dd_func, sigma3_dd_func, a1_func, a2_func, a3_func);
    
    
    y0 = zeros(15,1);
    y0(1:3) = x(:);
    y0(4:6) = x_d(:);
    y0(7:9) = sigma(:);
    y0(10:12) = sigma_d(:);
    y0(13:15) = [u_F; u_e1; u_e2];
    
    t_span = [0,.1];
    
    [t,y] = ode45(@system_of_eq, t_span, y0);
    fprintf("Controlled Iteration %f\n", i)
    fprintf("u_F = %f | u_e1 = %f | u_e2 = %f | x1 = %f | x2 = %f | x3 = %f | x_d1 = %f | x_d2 = %f | x_d3 = %f | s1 = %f | s2 = %f | s3 = %f | s_d1 = %f | sd_2 = %f | s_d3 = %f \n\n", [u_F, u_e1, u_e2, y(end,1),y(end,2),y(end,3), y(end,4),y(end,5),y(end,6), y(end,7),y(end,8),y(end,9), y(end,10),y(end,11),y(end,12)])
   
   
    x = [y(end,1);y(end,2);y(end,3)];
    x_d = [y(end,4);y(end,5);y(end,6)];
    sigma = [y(end,7);y(end,8);y(end,9)];
    sigma_d = [y(end,10);y(end,11);y(end,12)];

    y_net_ctrl = [y_net_ctrl; y];
    t_net_ctrl = [t_net_ctrl; t+.1*(i-1)];
    F_prev = u_F;

    e1_collected = [e1_collected, u_e1];
    e2_collected = [e2_collected, u_e2];

    e1_prev = u_e1;
    e2_prev = u_e2;
end

%% 

%scatter3(y_net_ctrl(:,1),y_net_ctrl(:,2),y_net_ctrl(:,3), 30, t_net_ctrl, 'filled')
plot3(y_net_ctrl(:,1),y_net_ctrl(:,2),y_net_ctrl(:,3))
%colorbar;
%legend("Uncontrolled", "Controlled")
plot3(x_Traj(1:25), y_Traj(1:25), z_Traj(1:25));
legend("Uncontrolled Trajectory", "Controlled Trajectory", "Desired Trajectory")
ylabel("Y")
xlabel("X")
zlabel("Z (Altitude)")
hold off
grid on

%% 

t = linspace(0,2.5,25);
figure
     stairs(t,e1_collected, '-b')
     hold on
     stairs(t,e2_collected, '-k')
     yline(pi/2, '--r', 'pi/2 limit')
     yline(-pi/2, '--r', 'pi/2 limit')
     hold off
     xlabel("Time (s)")
     ylabel("Gimbal Angles (radians)")
     legend("Gimbal Angle e1", "Gimbal Angle e2")
     ylim([-3.5,3.5])
     


%% Results
% iter = 1093;
% t = linspace(0,2.5,iter);
% 
% x_Traj = linspace(0,2.5,iter);
% y_Traj = linspace(0,2.5,iter);
% z_Traj = 1000 + 10*x_Traj;%1020.*1./((1+exp((x_Traj/10-5))).*(1+exp((y_Traj/10-5))));
% 
% Desired_Traj = [x_Traj', y_Traj', z_Traj'];
% MSE_Controlled = sum(sum((Desired_Traj - y_net_ctrl(:,1:3)).^2))
% 
% iter = 1041;
% t = linspace(0,2.5,iter);
% 
% x_Traj = linspace(0,2.5,iter);
% y_Traj = linspace(0,2.5,iter);
% z_Traj = 1000 + 10*x_Traj;%1020.*1./((1+exp((x_Traj/10-5))).*(1+exp((y_Traj/10-5))));
% Desired_Traj = [x_Traj', y_Traj', z_Traj'];
% 
% MSE_Uncontrolled = sum(sum((Desired_Traj - y_net(:, 1:3)).^2))
