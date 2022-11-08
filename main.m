clear all; clc;

%%%%%%%%%% Hyper param %%%%%%%%%%%
Simulation_Time=20;     % seconds of simulation (s)
Ts=0.1;                 % Sampling interval in (s)
BS_num=8;               % number of base station
BS_loc=[0,15;0,35;15,0;35,0;50,15;50,35;15,50;15,35];
MS_center=[25,25];      % assuming mobile station do a circular motion
MS_radius=6;    
MS_ang_velocity=0.1;


%%%%%%%%%%%% Constant %%%%%%%%%%%%
h_2=2e-23;  h_0=2e-20; 
St=2*h_0; Sf=8*pi*pi*h_2;
c=299792458;        %light speed


[m,~]=size(S);
num=1000;
RMSE=zeros(length(sigma_var),1);

for j=1:length(sigma_var)
    error=0;
    sigma_var1=sigma_var(j);
    for i=1:num
    r1=S-ones(m,1)*MS;
    r2=(sum(r1.^2,2)).^(1/2);
    r=r2(2:end,:)-ones(m-1,1)*r2(1,:)+sigma_var1*randn(m-1,1); %noised distance to MS
    sigma=sigma_var1^2;
    theta=TDOA_chan(S,r,sigma);
    error=error+norm(MS-theta)^2;
    end
    RMSE(j)=(error/num)^(1/2);
end
semilogx(sigma_var,RMSE,'-O');
xlabel('std of measured noise');
ylabel('RMSE');
title('5 Base');
legend('TDOA-CHAN');