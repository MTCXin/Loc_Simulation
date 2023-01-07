function [x0,y0,z0,Gxx,Gyy,Gzz] = doubleEKF_diff_clock_pos(pseudodis,pseudodis_SS,BS_num,BS_loc,SS_loc,Ts,Qk,Sigmaxy,Sigmaz,variance_xy,variance_z,DELTA_TK_ORD,FK_ORD,BSdelta_tk_1,BSf_k_1,MSdelta_tk_1,MSf_k_1,real_x,real_y,real_z,MEASURE_NOISE)
%DOUBLEEKF_DIFF_CLOCK_POS 此处显示有关此函数的摘要
%   use EKF1 to estimate clock_delta_diff and fk_diff, use EKF2 to estimate
%   position of BS
    
c=299792458;

persistent A1 Q1 R1 x1 P1 A2 Q2 R2 x2 P2 firstRun

if isempty(firstRun)

    BSdelta_tk_1_diff=BSdelta_tk_1(1:end-1)'-BSdelta_tk_1(end);
%     BSdelta_tk_1_diff=zeros(BS_num-1);
    BSf_k_1_diff=BSf_k_1(1:end-1)'-BSf_k_1(end);
%     BSf_k_1_diff=zeros(BS_num-1);
    x1=[BSdelta_tk_1_diff,BSf_k_1_diff]';   % hidden state:[deltaTk*n,fk*n, deltaTk, fk, x, y, z]( 2n+5 dims)
    x2=[real_x,real_y,real_z]';
    R1 = MEASURE_NOISE*eye(BS_num-1);
    R2 = MEASURE_NOISE*eye(BS_num-1);
    A1=eye(BS_num*2-2);
    A2=eye(3);
    for i = 1:BS_num-1
        A1(i,i+BS_num-1)=Ts;
    end

    P1=zeros(2*BS_num-2);
    P2=zeros(3);
    for i = 1:BS_num-1
        P1(i,i)=DELTA_TK_ORD^2;
        P1(i+BS_num-1,i+BS_num-1)=FK_ORD^2;
    end
    P2(1,1)=variance_xy;
    P2(2,2)=variance_xy; 
    P2(3,3)=variance_z;
    Q1=zeros(BS_num*2-2);
    Q2=zeros(3);
    for i=1:BS_num-1
        Q1(i,i)=2*Qk(1,1);
        Q1(i+BS_num-1,i+BS_num-1)=2*Qk(2,2); %TODO
        Q1(i,i+BS_num-1)=2*Qk(1,2);
        Q1(i+BS_num-1,i)=2*Qk(1,2);
    end
    Q2(1,1)=Sigmaxy; Q2(2,2)=Sigmaxy; Q2(3,3)=Sigmaz;
    firstRun=1;
end

H1=zeros(BS_num-1,BS_num*2-2);
for i = 1:BS_num-1
    H1(i,i)=-1*c;
end

xp1=A1*x1;
Pp1=A1*P1*A1' + Q1;
K1=Pp1*H1'/(H1*Pp1*H1' + R1);

hx1=zeros(BS_num-1,1);
for i = 1:BS_num-1
    hx1(i)= ((SS_loc(1)-BS_loc(i,1))^2+(SS_loc(2)-BS_loc(i,2))^2+(SS_loc(3)-BS_loc(i,3))^2)^(1/2) - ((SS_loc(1)-BS_loc(end,1))^2+(SS_loc(2)-BS_loc(end,2))^2+(SS_loc(3)-BS_loc(end,3))^2)^(1/2) - c*x1(i);
end

pseudodis_diff1=pseudodis_SS(1:end-1)-pseudodis_SS(end);
x1=xp1 + K1*(pseudodis_diff1-hx1);
P1=Pp1 - K1*H1*Pp1;

H2=zeros(BS_num-1,3);
for i = 1:BS_num-1
    H2(i,1)=(x2(1)-BS_loc(i,1))/((x2(1)-BS_loc(i,1))^2+(x2(2)-BS_loc(i,2))^2+(x2(3)-BS_loc(i,3))^2)^(1/2) - (x2(1)-BS_loc(end,1))/((x2(1)-BS_loc(end,1))^2+(x2(2)-BS_loc(end,2))^2+(x2(3)-BS_loc(end,3))^2)^(1/2);
    H2(i,2)=(x2(2)-BS_loc(i,2))/((x2(1)-BS_loc(i,1))^2+(x2(2)-BS_loc(i,2))^2+(x2(3)-BS_loc(i,3))^2)^(1/2) - (x2(2)-BS_loc(end,2))/((x2(1)-BS_loc(end,1))^2+(x2(2)-BS_loc(end,2))^2+(x2(3)-BS_loc(end,3))^2)^(1/2);
    H2(i,3)=(x2(3)-BS_loc(i,3))/((x2(1)-BS_loc(i,1))^2+(x2(2)-BS_loc(i,2))^2+(x2(3)-BS_loc(i,3))^2)^(1/2) - (x2(3)-BS_loc(end,3))/((x2(1)-BS_loc(end,1))^2+(x2(2)-BS_loc(end,2))^2+(x2(3)-BS_loc(end,3))^2)^(1/2);
end
xp2=A2*x2;
Pp2=A2*P2*A2' + Q2;
K2=Pp2*H2'/(H2*Pp2*H2' + R2);
hx2=zeros(BS_num-1,1);
for i = 1:BS_num-1
    hx2(i)= ((x2(1)-BS_loc(i,1))^2+(x2(2)-BS_loc(i,2))^2+(x2(3)-BS_loc(i,3))^2)^(1/2) - ((x2(1)-BS_loc(end,1))^2+(x2(2)-BS_loc(end,2))^2+(x2(3)-BS_loc(end,3))^2)^(1/2) - c*x1(i);
end

pseudodis_diff2=pseudodis(1:end-1)-pseudodis(end);
x2=xp2 + K2*(pseudodis_diff2-hx2);
P2=Pp2 - K2*H2*Pp2;
x0=x2(1);
y0=x2(2);
z0=x2(3);
Gxx=P2(1);
Gyy=P2(2);
Gzz=P2(3);
end

