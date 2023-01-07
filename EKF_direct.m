function [x0,y0,z0,Gxx,Gyy,Gzz,BSf_k_pre,BSdelta_tk_pre] = EKF_direct(pseudodis,BS_num,BS_loc,Ts,Qk,Sigmaxy,Sigmaz,variance_xy,variance_z,DELTA_TK_ORD,FK_ORD,BSdelta_tk_1,BSf_k_1,MSdelta_tk_1,MSf_k_1,real_x,real_y,real_z,MEASURE_NOISE)
%EKF_direct 
c=299792458;

persistent A Q R x P firstRun

if isempty(firstRun)

    x=[BSdelta_tk_1',BSf_k_1',MSdelta_tk_1,MSf_k_1,real_x,real_y,real_z]';   % hidden state:[deltaTk*n,fk*n, deltaTk, fk, x, y, z]( 2n+5 dims)
    R = MEASURE_NOISE*eye(BS_num);
    A=eye(BS_num*2+5);
    for i = 1:BS_num
        A(i,i+BS_num)=Ts;
    end
    A(2*BS_num+1,2*BS_num+2)=Ts;

    P=zeros(2*BS_num+5);
    for i = 1:BS_num
        P(i,i)=DELTA_TK_ORD^3;
        P(i+BS_num,i+BS_num)=FK_ORD^3;
    end
    P(1+BS_num*2,1+BS_num*2)=DELTA_TK_ORD^3;
    P(2+BS_num*2,2+BS_num*2)=FK_ORD^3;
    P(3+BS_num*2,3+BS_num*2)=variance_xy;
    P(4+BS_num*2,4+BS_num*2)=variance_xy; 
    P(5+BS_num*2,5+BS_num*2)=variance_z;
    Q=zeros(BS_num*2+5);
    for i=1:BS_num
        Q(i,i)=Qk(1,1);
        Q(i+BS_num,i+BS_num)=Qk(2,2);
        Q(i,i+BS_num)=Qk(1,2);
        Q(i+BS_num,i)=Qk(1,2);
    end
    Q(BS_num*2+3,BS_num*2+3)=Sigmaxy; Q(BS_num*2+4,BS_num*2+4)=Sigmaxy; Q(BS_num*2+5,BS_num*2+5)=Sigmaz;
    Q(BS_num*2+1,BS_num*2+1)=Qk(1,1); Q(BS_num*2+2,BS_num*2+2)=Qk(2,2);
    Q(BS_num*2+1,BS_num*2+2)=Qk(1,2); Q(BS_num*2+2,BS_num*2+1)=Qk(1,2);
    firstRun=1;
end

H=zeros(BS_num,BS_num*2+5);
for i = 1:BS_num
    H(i,i)=-1*c;
    H(i,BS_num*2+1)=c;
    H(i,BS_num*2+3)=(x(BS_num*2+3)-BS_loc(i,1))/((x(BS_num*2+3)-BS_loc(i,1))^2+(x(BS_num*2+4)-BS_loc(i,2))^2+(x(BS_num*2+5)-BS_loc(i,3))^2)^(1/2);
    H(i,BS_num*2+4)=(x(BS_num*2+4)-BS_loc(i,2))/((x(BS_num*2+3)-BS_loc(i,1))^2+(x(BS_num*2+4)-BS_loc(i,2))^2+(x(BS_num*2+5)-BS_loc(i,3))^2)^(1/2);
    H(i,BS_num*2+5)=(x(BS_num*2+5)-BS_loc(i,3))/((x(BS_num*2+3)-BS_loc(i,1))^2+(x(BS_num*2+4)-BS_loc(i,2))^2+(x(BS_num*2+5)-BS_loc(i,3))^2)^(1/2);
end

xp=A*x;
Pp=A*P*A' + Q;
K=Pp*H'/(H*Pp*H' + R);

hx=zeros(BS_num,1);
for i = 1:BS_num
    hx(i)= ((x(BS_num*2+3)-BS_loc(i,1))^2+(x(BS_num*2+4)-BS_loc(i,2))^2+(x(BS_num*2+5)-BS_loc(i,3))^2)^(1/2) + c*(x(BS_num*2+1)-x(i));
end

x=xp + K*(pseudodis-hx);
P=Pp - K*H*Pp;

x0=x(BS_num*2+3);
y0=x(BS_num*2+4);
z0=x(BS_num*2+5);
Gxx=P(BS_num*2+3,BS_num*2+3);
Gyy=P(BS_num*2+4,BS_num*2+4);
Gzz=P(BS_num*2+5,BS_num*2+5);
BSf_k_pre=(x(BS_num+1:BS_num*2-1)-x(BS_num*2));
BSdelta_tk_pre=(x(1:BS_num-1)-x(BS_num));
end

