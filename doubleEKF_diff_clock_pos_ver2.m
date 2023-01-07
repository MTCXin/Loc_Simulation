function [x0,y0,z0,Gxx,Gyy,Gzz,BSf_k_pre,BSdelta_tk_pre] = doubleEKF_diff_clock_pos_ver2(pseudodis,pseudodis_SS,BS_num,BS_loc,SS_loc,Ts,Qk,Sigmaxy,Sigmaz,variance_xy,variance_z,DELTA_TK_ORD,FK_ORD,BSdelta_tk_1,BSf_k_1,MSdelta_tk_1,MSf_k_1,real_x,real_y,real_z,MEASURE_NOISE)
%DOUBLEEKF_DIFF_CLOCK_POS_VER2 此处显示有关此函数的摘要
%   use EKF1 to estimate fk_diff, use EKF2 to estimate
%   clock_delta_diff and position of BS
%   position of SS is unknown
    

factor_time=1e3;
Ts=Ts.*factor_time; %convert to ?s
DELTA_TK_ORD=DELTA_TK_ORD.*factor_time;
BSdelta_tk_1=BSdelta_tk_1.*factor_time;
MSdelta_tk_1=MSdelta_tk_1.*factor_time;
Qk(1,1)=Qk(1,1)*factor_time^2;
Qk(1,2)=Qk(1,2)*factor_time;
Qk(2,1)=Qk(2,1)*factor_time;

factor_lenth=1e5;
pseudodis=pseudodis./factor_lenth;
pseudodis_SS=pseudodis_SS./factor_lenth;
BS_loc=BS_loc./factor_lenth;
SS_loc=SS_loc./factor_lenth;
Sigmaxy=Sigmaxy./(factor_lenth^2);
Sigmaz=Sigmaz./(factor_lenth^2);
variance_xy=variance_xy./factor_lenth^2;
variance_z=variance_z./factor_lenth^2;
real_x=real_x./factor_lenth;
real_y=real_y./factor_lenth;
real_z=real_z./factor_lenth;
MEASURE_NOISE=MEASURE_NOISE./factor_lenth^2;

c=299792458/(factor_time*factor_lenth); %(?m/?s)

BSdelta_tk_1_diff=BSdelta_tk_1(1:end-1)'-BSdelta_tk_1(end);
%     BSdelta_tk_1_diff=zeros(BS_num-1);
BSf_k_1_diff=BSf_k_1(1:end-1)'-BSf_k_1(end);
%     BSf_k_1_diff=zeros(BS_num-1);

persistent A1 Q1 R1 x1 P1 A2 Q2 R2 x2 P2 prev_pseudodis_SS firstRun

if isempty(firstRun)

    
    prev_pseudodis_SS=pseudodis_SS+c*BSf_k_1*Ts;
    x1=BSf_k_1_diff';   % hidden state:
    x2=[BSdelta_tk_1_diff,BSf_k_1_diff,real_x,real_y,real_z]';
    R1 = MEASURE_NOISE*eye(BS_num-1);
    R2 = MEASURE_NOISE*eye(BS_num*2-2);
    A1=eye(BS_num-1);
    A2=eye(BS_num*2+1);
    for i = 1:BS_num-1
        A2(i,i+BS_num-1)=Ts;
    end
    P1=zeros(BS_num-1);
    P2=zeros(BS_num*2+1);
    for i = 1:BS_num-1
        P1(i,i)=FK_ORD^2;
        P2(i+BS_num-1,i+BS_num-1)=FK_ORD^2;
        P2(i,i)=DELTA_TK_ORD^2;
    end
    P2(BS_num*2-1,BS_num*2-1)=variance_xy;
    P2(BS_num*2,BS_num*2)=variance_xy; 
    P2(BS_num*2+1,1+BS_num*2)=variance_z;
    Q1=zeros(BS_num-1);
    Q2=zeros(BS_num*2+1);
    for i=1:BS_num-1
        Q1(i,i)=2*Qk(2,2); 
        Q2(i,i)=2*Qk(1,1);
        Q2(i+BS_num-1,i+BS_num-1)=2*Qk(2,2); 
        Q2(i,i+BS_num-1)=2*Qk(1,2);
        Q2(i+BS_num-1,i)=2*Qk(1,2);
    end
    Q2(BS_num*2-1,BS_num*2-1)=Sigmaxy; Q2(BS_num*2,BS_num*2)=Sigmaxy; Q2(BS_num*2+1,BS_num*2+1)=Sigmaz;
    firstRun=1;
end
H1=zeros(BS_num-1,BS_num-1);
for i = 1:BS_num-1
    H1(i,i)=-1*c*Ts;
end

xp1=A1*x1;
Pp1=A1*P1*A1' + Q1;
K1=Pp1*H1'/(H1*Pp1*H1' + R1);

hx1=zeros(BS_num-1,1);
for i = 1:BS_num-1
    hx1(i)=-1*x1(i)*Ts*c;
end

pseudodis_diff1=(pseudodis_SS(1:end-1)-pseudodis_SS(end))-(prev_pseudodis_SS(1:end-1)-prev_pseudodis_SS(end));
x1=xp1 + K1*(pseudodis_diff1-hx1);
P1=Pp1 - K1*H1*Pp1;

H2=zeros(BS_num*2-2,BS_num*2+1);
for i = 1:BS_num-1
    H2(i,i)=-1*c;
    H2(i,i+BS_num-1)=-1*c*Ts;
    H2(i+BS_num-1,i+BS_num-1)=1;
    H2(i,BS_num*2-1)=(x2(BS_num*2-1)-BS_loc(i,1))/((x2(BS_num*2-1)-BS_loc(i,1))^2+(x2(BS_num*2)-BS_loc(i,2))^2+(x2(BS_num*2+1)-BS_loc(i,3))^2)^(1/2) - (x2(BS_num*2-1)-BS_loc(end,1))/((x2(BS_num*2-1)-BS_loc(end,1))^2+(x2(BS_num*2)-BS_loc(end,2))^2+(x2(BS_num*2+1)-BS_loc(end,3))^2)^(1/2);
    H2(i,BS_num*2)  =(x2(BS_num*2  )-BS_loc(i,2))/((x2(BS_num*2-1)-BS_loc(i,1))^2+(x2(BS_num*2)-BS_loc(i,2))^2+(x2(BS_num*2+1)-BS_loc(i,3))^2)^(1/2) - (x2(BS_num*2  )-BS_loc(end,2))/((x2(BS_num*2-1)-BS_loc(end,1))^2+(x2(BS_num*2)-BS_loc(end,2))^2+(x2(BS_num*2+1)-BS_loc(end,3))^2)^(1/2);
    H2(i,BS_num*2+1)=(x2(BS_num*2+1)-BS_loc(i,3))/((x2(BS_num*2-1)-BS_loc(i,1))^2+(x2(BS_num*2)-BS_loc(i,2))^2+(x2(BS_num*2+1)-BS_loc(i,3))^2)^(1/2) - (x2(BS_num*2+1)-BS_loc(end,3))/((x2(BS_num*2-1)-BS_loc(end,1))^2+(x2(BS_num*2)-BS_loc(end,2))^2+(x2(BS_num*2+1)-BS_loc(end,3))^2)^(1/2);
end

for i = 1:BS_num-1
    R2(i+BS_num-1,i+BS_num-1)=1e-9;   %P1(i,i); %FK_ORD^2  
end
xp2=A2*x2;
Pp2=A2*P2*A2' + Q2;
K2=Pp2*H2'/(H2*Pp2*H2' + R2);
hx2=zeros(BS_num-1,1);
for i = 1:BS_num-1
    hx2(i+BS_num-1)=x2(i+BS_num-1);
    hx2(i)= ((x2(BS_num*2-1)-BS_loc(i,1))^2+(x2(BS_num*2)-BS_loc(i,2))^2+(x2(BS_num*2+1)-BS_loc(i,3))^2)^(1/2) - ((x2(BS_num*2-1)-BS_loc(end,1))^2+(x2(BS_num*2)-BS_loc(end,2))^2+(x2(BS_num*2+1)-BS_loc(end,3))^2)^(1/2) - c*x2(i);
end

pseudodis_diff2=[(pseudodis(1:end-1)-pseudodis(end));x1];
x2=xp2 + K2*(pseudodis_diff2-hx2);
P2=Pp2 - K2*H2*Pp2;
prev_pseudodis_SS=pseudodis_SS;
x0=x2(BS_num*2-1)*factor_lenth;
y0=x2(BS_num*2)*factor_lenth;
z0=x2(BS_num*2+1)*factor_lenth;
Gxx=P2(BS_num*2-1,BS_num*2-1)*factor_lenth^2;
Gyy=P2(BS_num*2,BS_num*2)*factor_lenth^2;
Gzz=P2(BS_num*2+1,BS_num*2+1)*factor_lenth^2;
BSf_k_pre=x2(BS_num:BS_num*2-2);
BSdelta_tk_pre=x2(1:BS_num-1)/factor_time;
end



