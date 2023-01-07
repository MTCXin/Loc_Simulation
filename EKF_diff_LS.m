function [x0,y0,z0,Gxx,Gyy,Gzz] = EKF_diff_LS(pseudodis,BS_num,BS_loc,Ts,Qk,Sigmaxy,Sigmaz,variance_xy,variance_z,DELTA_TK_ORD,FK_ORD,BSdelta_tk_1,BSf_k_1,MSdelta_tk_1,MSf_k_1,real_x,real_y,real_z,MEASURE_NOISE)
%EKF_direct_diff
c=299792458;

persistent A Q R x P firstRun

if isempty(firstRun)

    BSdelta_tk_1_diff=BSdelta_tk_1(1:end-1)'-BSdelta_tk_1(end);
    BSf_k_1_diff=BSf_k_1(1:end-1)'-BSf_k_1(end);
    x=[BSdelta_tk_1_diff,BSf_k_1_diff,real_x,real_y,real_z]';   % hidden state:[deltaTk*n,fk*n, deltaTk, fk, x, y, z]( 2n+5 dims)
    R = MEASURE_NOISE*eye(BS_num-1);
    A=eye(BS_num*2+1);
    for i = 1:BS_num-1
        A(i,i+BS_num-1)=Ts;
    end

    P=zeros(2*BS_num+1);
    for i = 1:BS_num-1
        P(i,i)=DELTA_TK_ORD^3;
        P(i+BS_num-1,i+BS_num-1)=FK_ORD^3;
    end

    P(BS_num*2-1,BS_num*2-1)=variance_xy;
    P(BS_num*2,BS_num*2)=variance_xy; 
    P(BS_num*2+1,1+BS_num*2)=variance_z;
    Q=zeros(BS_num*2+1);
    for i=1:BS_num-1
        Q(i,i)=2*Qk(1,1);
        Q(i+BS_num-1,i+BS_num-1)=2*Qk(2,2); %TODO
        Q(i,i+BS_num-1)=2*Qk(1,2);
        Q(i+BS_num-1,i)=2*Qk(1,2);
    end
    Q(BS_num*2-1,BS_num*2-1)=Sigmaxy; Q(BS_num*2,BS_num*2)=Sigmaxy; Q(BS_num*2+1,BS_num*2+1)=Sigmaz;
    firstRun=1;
end

H=zeros(BS_num-1,BS_num*2+1);
for i = 1:BS_num-1
    H(i,i)=-1*c;
    H(i,BS_num*2-1)=(x(BS_num*2-1)-BS_loc(i,1))/((x(BS_num*2-1)-BS_loc(i,1))^2+(x(BS_num*2)-BS_loc(i,2))^2+(x(BS_num*2+1)-BS_loc(i,3))^2)^(1/2) - (x(BS_num*2-1)-BS_loc(end,1))/((x(BS_num*2-1)-BS_loc(end,1))^2+(x(BS_num*2)-BS_loc(end,2))^2+(x(BS_num*2+1)-BS_loc(end,3))^2)^(1/2);
    H(i,BS_num*2)  =(x(BS_num*2  )-BS_loc(i,2))/((x(BS_num*2-1)-BS_loc(i,1))^2+(x(BS_num*2)-BS_loc(i,2))^2+(x(BS_num*2+1)-BS_loc(i,3))^2)^(1/2) - (x(BS_num*2  )-BS_loc(end,2))/((x(BS_num*2-1)-BS_loc(end,1))^2+(x(BS_num*2)-BS_loc(end,2))^2+(x(BS_num*2+1)-BS_loc(end,3))^2)^(1/2);
    H(i,BS_num*2+1)=(x(BS_num*2+1)-BS_loc(i,3))/((x(BS_num*2-1)-BS_loc(i,1))^2+(x(BS_num*2)-BS_loc(i,2))^2+(x(BS_num*2+1)-BS_loc(i,3))^2)^(1/2) - (x(BS_num*2+1)-BS_loc(end,3))/((x(BS_num*2-1)-BS_loc(end,1))^2+(x(BS_num*2)-BS_loc(end,2))^2+(x(BS_num*2+1)-BS_loc(end,3))^2)^(1/2);
end

xp=A*x;
Pp=A*P*A' + Q;
K=Pp*H'/(H*Pp*H' + R);

hx=zeros(BS_num-1,1);
for i = 1:BS_num-1
    hx(i)= ((x(BS_num*2-1)-BS_loc(i,1))^2+(x(BS_num*2)-BS_loc(i,2))^2+(x(BS_num*2+1)-BS_loc(i,3))^2)^(1/2) - ((x(BS_num*2-1)-BS_loc(end,1))^2+(x(BS_num*2)-BS_loc(end,2))^2+(x(BS_num*2+1)-BS_loc(end,3))^2)^(1/2) - c*x(i);
end

pseudodis_diff=pseudodis(1:end-1)-pseudodis(end);
x=xp + K*(pseudodis_diff-hx);
P=Pp - K*H*Pp;

x0=x(BS_num*2-1);
y0=x(BS_num*2);
z0=x(BS_num*2+1);
Gxx=P(BS_num*2-1,BS_num*2-1);
Gyy=P(BS_num*2,BS_num*2);
Gzz=P(BS_num*2+1,BS_num*2+1);
end
end

