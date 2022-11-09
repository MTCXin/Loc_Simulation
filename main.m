clear all; clc;

%%%%%%%%%% Hyper param %%%%%%%%%%%
Simulation_Time=20;     % seconds of simulation (s)
Ts=0.1;                 % Sampling interval in (s)
BS_num=8;               % number of base station
F_clock=5e9;            % nominal frequency of clock (Hz)
BS_loc=[0,15,0.5;0,35,1.5;15,0,2.5;35,0,0.5;50,15,1.5;50,35,2.5;15,50,0.5;15,35,1.5];
MS_center=[25,25,1];      % assuming mobile station do a circular motion
MS_radius=6;    
MS_ang_velocity=0.5;
MS_height_ang_velocity=0.4;
MS_height_diff=0.5;
DELTA_TK_ORD = 1e-4;    % the order of delta_tk in (s)
FK_ORD = 1e-5;          % the order of fk in (s/s)
%%%%%%%%%%%% Constant %%%%%%%%%%%%
h_2=2e-23;  h_0=2e-20; 
St=2*h_0; Sf=8*pi*pi*h_2;
c=299792458;        %light speed
Qk=[(St*Ts+Sf*Ts^3/3),(Sf*Ts^2/2);(Sf*Ts^2/2),Sf*Ts];

%%%%%%%%%%% initialize %%%%%%%%%%%
time=0;
BSdelta_tk_1=rand(BS_num,1)*DELTA_TK_ORD;
BSf_k_1=rand(BS_num,1)*FK_ORD;
MSdelta_tk_1=rand()*DELTA_TK_ORD;
MSf_k_1=rand()*FK_ORD;
BSdelta_tk=zeros(BS_num,1);
BSf_k=zeros(BS_num,1);
MSdelta_tk=0;
MSf_k=0;

while time<Simulation_Time
    error=0;
    pseudodis=zeros(BS_num,1);
    %%%%%%%%%%% calcu the fk and deltk of time step
    noise=mvnrnd([0,0],Qk,BS_num+1);
    MSdelta_tk=MSdelta_tk_1+Ts*MSf_k_1+noise(BS_num+1,1);
    MSf_k=MSf_k_1+noise(BS_num+1,2);
    for i=1:BS_num
        BSdelta_tk(i)=BSdelta_tk_1(i)+Ts*BSf_k_1(i)+noise(i,1);
        BSf_k(i)=BSf_k_1(i)+noise(i,2);
    end

    %%%%%%%%%%% calcu the real location of MS
    real_x = MS_radius*cos(time*MS_ang_velocity)+MS_center(1);
    real_y = MS_radius*sin(time*MS_ang_velocity)+MS_center(2);
    real_z = MS_height_diff*cos(MS_height_ang_velocity*time)+MS_center(3);

    %%%%%%%%%%% calcu the noised pseodo distance
    for i=1:BS_num
        pseudodis(i)=((BS_loc(i,1)-real_x)^2+(BS_loc(i,2)-real_y)^2+(BS_loc(i,3)-real_z)^2)^(1/2)+(BSdelta_tk(i)-MSdelta_tk)*c;
    end
    [x,y,z]=EKF_direct(pseudodis,BS_num,BS_loc,Ts,Qk);

%     error=error+norm(MS-theta)^2;
%     RMSE(j)=(error/num)^(1/2);
end
semilogx(sigma_var,RMSE,'-O');
xlabel('std of measured noise');
ylabel('RMSE');
title('5 Base');
legend('TDOA-CHAN');