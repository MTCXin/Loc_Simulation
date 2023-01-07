clear all;  clc; close all;

%%%%%%%%%% Hyper param %%%%%%%%%%%
Simulation_Time=20;     % seconds of simulation (s)
Ts=0.05;                 % Sampling interval in (s)
BS_num=12;               % number of base station

Sigmaxy=5*Ts^2;              % indicate the xy_pos change variance of MS during TS 
Sigmaz=0.5*Ts^2;               % indicate the z_pos change variance of MS during TS 
variance_xy=0.5;            % indicate the ini xy_pos variance of MS 
variance_z=0.05;            % indicate the ini z_pos variance of MS 

% BS_loc=[0,15,0.5;0,35,1.5;15,0,2.5;35,0,0.5;50,15,1.5;50,35,2.5;15,50,0.5;35,50,1.5]; % 8 version
BS_loc=[0,12,0.5;0,25,1.5;0,38,2.5;12,0,2.5;25,0,0.5;38,0,1.5;50,12,1.5;50,25,0.5;50,38,2.5;12,50,0.5;25,50,1.5;38,50,2.5]; % 12 version
MS_center=[25,25,1];      % assuming mobile station do a circular motion
MS_radius=8;    
MS_ang_velocity=0.3;
MS_height_ang_velocity=0.04;
MS_height_diff=0.5;

SS_loc=[5,5,1.5]; %Stationary Station
DELTA_TK_ORD = 1e-5;    % the order of delta_tk in (s)
FK_ORD = 1e-6;          % the order of fk in (s/s)
RMSE_TAU = 0.04;
MEASURE_NOISE=10;
%%%%%%%%%%%% Constant %%%%%%%%%%%%
h_2=2e-23;  h_0=2e-20; 
St=2*h_0; Sf=8*pi*pi*h_2;
c=299792458;        %light speed
Qk=[(St*Ts+Sf*(Ts^3)/3),(Sf*(Ts^2)/2);(Sf*(Ts^2)/2),Sf*Ts];

%%%%%%%%%%% initialize %%%%%%%%%%%
time=0;
% BSdelta_tk_1=rand(BS_num,1)*DELTA_TK_ORD;  %random
% BSdelta_tk_1=ones(BS_num,1)*DELTA_TK_ORD;
BSdelta_tk_1=0.5*ones(BS_num,1)*DELTA_TK_ORD+0.5*rand(BS_num,1)*DELTA_TK_ORD;  %assigned
% BSf_k_1=rand(BS_num,1)*FK_ORD;   %random
% BSf_k_1=ones(BS_num,1)*FK_ORD;
BSf_k_1=0.5*ones(BS_num,1)*FK_ORD+0.5*rand(BS_num,1)*FK_ORD;
% MSdelta_tk_1=rand()*DELTA_TK_ORD;
MSdelta_tk_1=0.5*ones(1)*DELTA_TK_ORD+0.5*rand()*DELTA_TK_ORD;
SSdelta_tk_1=0.5*ones(1)*DELTA_TK_ORD+0.5*rand()*DELTA_TK_ORD;
% MSdelta_tk_1=ones(1)*DELTA_TK_ORD;
% SSdelta_tk_1=ones(1)*DELTA_TK_ORD;
% MSf_k_1=rand()*FK_ORD;
MSf_k_1=0.5*ones()*FK_ORD+0.5*rand()*FK_ORD;
SSf_k_1=0.5*ones()*FK_ORD+0.5*rand()*FK_ORD;
% MSf_k_1=ones()*FK_ORD;
% SSf_k_1=ones()*FK_ORD;
BSdelta_tk=zeros(BS_num,1);
BSf_k=zeros(BS_num,1);
MSdelta_tk=0;
MSf_k=0;
SSdelta_tk=0;
SSf_k=0;
real_track=zeros(Simulation_Time/Ts,3);
esti_track=zeros(Simulation_Time/Ts,3);
RMSE_xy=zeros(Simulation_Time/Ts,1);
RMSE_xy_real=zeros(Simulation_Time/Ts,1);
RMSE_z=zeros(Simulation_Time/Ts,1);
RMSE_z_real=zeros(Simulation_Time/Ts,1);
RMSE_delta=zeros(Simulation_Time/Ts,1);
RMSE_fk=zeros(Simulation_Time/Ts,1);
RMSE_delta_real=zeros(Simulation_Time/Ts,1);
RMSE_fk_real=zeros(Simulation_Time/Ts,1);
PErr_xy=zeros(Simulation_Time/Ts,1);
PErr_z=zeros(Simulation_Time/Ts,1);
time_track=zeros(Simulation_Time/Ts,1);
time_step=1;

while time<Simulation_Time
    error=0;
    pseudodis=zeros(BS_num,1);
    pseudodis_SS=zeros(BS_num,1);
    %%%%%%%%%%% calcu the fk and deltk of time step
    noise=mvnrnd([0,0],Qk,BS_num+2);
    pseudor_noise=mvnrnd(0,MEASURE_NOISE,BS_num*2);
    MSdelta_tk=MSdelta_tk_1+Ts*MSf_k_1+noise(BS_num+1,1);
    MSf_k=MSf_k_1+noise(BS_num+1,2);
    SSdelta_tk=SSdelta_tk_1+Ts*MSf_k_1+noise(BS_num+2,1);
    SSf_k=SSf_k_1+noise(BS_num+2,2);
    for i=1:BS_num
        BSdelta_tk(i)=BSdelta_tk_1(i)+Ts*BSf_k_1(i)+noise(i,1);
        BSf_k(i)=BSf_k_1(i)+noise(i,2);
    end
    BSdelta_tk_diff=(BSdelta_tk(1:end-1)-BSdelta_tk(end));
    BSf_k_diff=(BSf_k(1:end-1)-BSf_k(end));
    %%%%%%%%%%% calcu the real location of MS
    real_x = MS_radius*cos(time*MS_ang_velocity)+MS_center(1);
    real_y = MS_radius*sin(time*MS_ang_velocity)+MS_center(2);
    real_z = MS_height_diff*cos(MS_height_ang_velocity*time)+MS_center(3);

    %%%%%%%%%%% calcu the noised pseodo distance
    for i=1:BS_num
        pseudodis(i)=((BS_loc(i,1)-real_x)^2+(BS_loc(i,2)-real_y)^2+(BS_loc(i,3)-real_z)^2)^(1/2)+(MSdelta_tk-BSdelta_tk(i))*c+pseudor_noise(i);
        pseudodis_SS(i)=((BS_loc(i,1)-SS_loc(1))^2+(BS_loc(i,2)-SS_loc(2))^2+(BS_loc(i,3)-SS_loc(3))^2)^(1/2)+(SSdelta_tk-BSdelta_tk(i))*c+pseudor_noise(i+BS_num);
    end
%     [x,y,z,Gxx,Gyy,Gzz,BSf_k_pre,BSdelta_tk_pre]=EKF_direct(pseudodis,BS_num,BS_loc,Ts,Qk,Sigmaxy,Sigmaz,variance_xy,variance_z,DELTA_TK_ORD,FK_ORD,BSdelta_tk_1,BSf_k_1,MSdelta_tk_1,MSf_k_1,real_x,real_y,real_z,MEASURE_NOISE);
%     [x,y,z,Gxx,Gyy,Gzz]=EKF_direct_diff(pseudodis,BS_num,BS_loc,Ts,Qk,Sigmaxy,Sigmaz,variance_xy,variance_z,DELTA_TK_ORD,FK_ORD,BSdelta_tk_1,BSf_k_1,MSdelta_tk_1,MSf_k_1,real_x,real_y,real_z,MEASURE_NOISE);
%     [x,y,z,Gxx,Gyy,Gzz,BSf_k_pre,BSdelta_tk_pre]=doubleEKF_diff_clock_pos_ver2(pseudodis,pseudodis_SS,BS_num,BS_loc,SS_loc,Ts,Qk,Sigmaxy,Sigmaz,variance_xy,variance_z,DELTA_TK_ORD,FK_ORD,BSdelta_tk_1,BSf_k_1,MSdelta_tk_1,MSf_k_1,real_x,real_y,real_z,MEASURE_NOISE);
    [x,y,z,Gxx,Gyy,Gzz,BSf_k_pre,BSdelta_tk_pre]=EWMA_EKF_diff(pseudodis,pseudodis_SS,BS_num,BS_loc,SS_loc,Ts,Qk,Sigmaxy,Sigmaz,variance_xy,variance_z,DELTA_TK_ORD,FK_ORD,BSdelta_tk_1,BSf_k_1,MSdelta_tk_1,MSf_k_1,real_x,real_y,real_z,MEASURE_NOISE);
    real_track(time_step,:)=[real_x,real_y,real_z];
    esti_track(time_step,:)=[x,y,z];

    MSdelta_tk_1=MSdelta_tk;
    MSf_k_1=MSf_k;
    for i=1:BS_num
        BSdelta_tk_1(i)=BSdelta_tk(i);
        BSf_k_1(i)=BSf_k(i);
    end
    time_track(time_step)=time;
    
    RMSE_xy_real(time_step)=((real_x-x)^2+(real_y-y)^2);
    RMSE_z_real(time_step)=((real_z-z)^2);
    RMSE_delta_real(time_step)=sum((BSdelta_tk_diff-BSdelta_tk_pre).^2);
    RMSE_fk_real(time_step)=sum((BSf_k_diff-BSf_k_pre).^2);
    if time_step==1
    RMSE_xy(time_step)=((real_x-x)^2+(real_y-y)^2)^(1/2);
    RMSE_z(time_step)=((real_z-z)^2)^(1/2);
    RMSE_delta(time_step)=sum(((BSdelta_tk(1:end-1)-BSdelta_tk(end))-BSdelta_tk_pre).^2);
    RMSE_fk(time_step)=sum(((BSf_k(1:end-1)-BSf_k(end))-BSf_k_pre).^2);
    PErr_z(time_step)=(Gzz)^(1/2);
    PErr_xy(time_step)=(Gxx+Gyy)^(1/2);
    else
    RMSE_xy(time_step)=((real_x-x)^2+(real_y-y)^2)^(1/2)*RMSE_TAU+(1-RMSE_TAU)*RMSE_xy(time_step-1);
    RMSE_z(time_step)=((real_z-z)^2)^(1/2)*RMSE_TAU+(1-RMSE_TAU)*RMSE_z(time_step-1); 
    RMSE_delta(time_step)=sum(((BSdelta_tk(1:end-1)-BSdelta_tk(end))-BSdelta_tk_pre).^2)^(1/2)*RMSE_TAU+(1-RMSE_TAU)*RMSE_delta(time_step-1);
    RMSE_fk(time_step)=sum(((BSf_k(1:end-1)-BSf_k(end))-BSf_k_pre).^2)^(1/2)*RMSE_TAU+(1-RMSE_TAU)*RMSE_fk(time_step-1);
    PErr_z(time_step)=(Gzz)^(1/2)*RMSE_TAU+(1-RMSE_TAU)*PErr_z(time_step-1);
    PErr_xy(time_step)=(Gxx+Gyy)^(1/2)*RMSE_TAU+(1-RMSE_TAU)*PErr_xy(time_step-1);
    end
    time_step=time_step+1;
    time=time+Ts;
    
end

figure(1);
plot3(real_track(:,1),real_track(:,2),real_track(:,3),'blue');
hold on;
plot3(esti_track(:,1),esti_track(:,2),esti_track(:,3),'red');
hold on;
plot3(BS_loc(:,1),BS_loc(:,2),BS_loc(:,3),'*');
hold on;
plot3(MS_radius+MS_center(1), MS_center(2),MS_height_diff+MS_center(3),'x');
hold on;
plot3(SS_loc(1), SS_loc(2),SS_loc(3),'o');
legend('real','esti','BS location','start point','SS location');
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
title('track');
time_track=time_track(20:end);
RMSE_xy=RMSE_xy(20:end); 
RMSE_z=RMSE_z(20:end);
RMSE_delta=RMSE_delta(20:end); 
RMSE_fk=RMSE_fk(20:end); 
RMSExy=mean(RMSE_xy_real)^(1/2);
RMSEz=mean(RMSE_z_real)^(1/2);
RMSE_fk_real=mean(RMSE_fk_real)^(1/2);
RMSE_delta_real=mean(RMSE_delta_real)^(1/2);
% PErr_xy=PErr_xy(300:end);
% PErr_z=PErr_z(300:end);
figure(2);
plot(time_track,RMSE_xy,'red');
xlabel('time (s)');
ylabel('RMSE-xy (m)');
title('RMSE-xy')

figure(3);
plot(time_track,RMSE_z,'red');
xlabel('time (s)');
ylabel('RMSE-z (m)');
title('RMSE-z')

figure(4);
plot(time_track,RMSE_fk,'red');
xlabel('time (s)');
ylabel('RMSE-fk');
title('RMSE-fk')

figure(5);
plot(time_track,RMSE_delta,'red');
xlabel('time (s)');
ylabel('RMSE-delta (s)');
title('RMSE-delta')

% figure(4);
% plot(time_track,PErr_xy,'red');
% xlabel('time (s)');
% ylabel('Err-pre-xy (m)');
% title('Err-pre-xy')
% 
% figure(5);
% plot(time_track,PErr_z,'red');
% xlabel('time (s)');
% ylabel('Err-pre-z (m)');
% title('Err-pre-z')