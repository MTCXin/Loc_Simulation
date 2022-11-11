function [x0,y0,z0] = EKF_direct(pseudodis,BS_num,BS_loc,Ts,Qk,Sigmaxy,Sigmaz)
%EKF_direct 
c=299792458;

persistent A Q R x P firstRun

if isempty(firstRun)

    x=[zeros(1,BS_num),zeros(1,BS_num),0,0,mean(BS_loc)];   % hidden state:[deltaTk*n,fk*n, deltaTk, fk, x, y, z]( 2n+5 dims)
    R = zeros(BS_num);
    A=eye(BS_num*2+5);
    for i = 1:BS_num
        A(i,i+BS_num)=Ts;
    end
    A(2*BS_num+1,2*BS_num+2)=Ts;
    P=(c*1e-4)^2*eye(BS_num);
    Q=zeros(BS_num*2+5);
    for i=1:BS_num
        Q(i,i)=Qk(1,1);
        Q(i+BS_num,i+BS_num)=Qk(2,2);
        Q(i,i+BS_num)=Qk(1,2);
        Q(i+BS_num,i)=Qk(1,2);
    end
    Q(BS_num+3,BS_num+3)=Sigmaxy; Q(BS_num+4,BS_num+4)=Sigmaxy; Q(BS_num+5,BS_num+5)=Sigmaz;
    Q(BS_num+1,BS_num+1)=Qk(1,1); Q(BS_num+2,BS_num+2)=Qk(2,2);
    Q(BS_num+1,BS_num+2)=Qk(1,2); Q(BS_num+2,BS_num+1)=Qk(1,2);
    firstRun=1;
end

H=zeros(BS_num,BS_num*2+5);
for i = 1:BS_num
    H(i,i)=c;
    H(BS_num*2+1,i)=-1*c;
    H(BS_num*2+3,i)=(x(BS_num*2+3)-BS_loc(i,1))/((x(BS_num*2+3)-BS_loc(i,1))^2+(x(BS_num*2+4)-BS_loc(i,2))^2+(x(BS_num*2+5)-BS_loc(i,3))^2)^(1/2);
    H(BS_num*2+4,i)=(x(BS_num*2+4)-BS_loc(i,2))/((x(BS_num*2+3)-BS_loc(i,1))^2+(x(BS_num*2+4)-BS_loc(i,2))^2+(x(BS_num*2+5)-BS_loc(i,3))^2)^(1/2);
    H(BS_num*2+5,i)=(x(BS_num*2+5)-BS_loc(i,3))/((x(BS_num*2+3)-BS_loc(i,1))^2+(x(BS_num*2+4)-BS_loc(i,2))^2+(x(BS_num*2+5)-BS_loc(i,3))^2)^(1/2);
end

xp=A*x;
Pp=A*P*A' + Q;
K=Pp*H'*inv(H*Pp*H' + R);

hx=zeros(BS_num,1);
for i = 1:BS_num
    hx(i)= ((x(BS_num*2+3)-BS_loc(i,1))^2+(x(BS_num*2+4)-BS_loc(i,2))^2+(x(BS_num*2+5)-BS_loc(i,3))^2)^(1/2) + c*(x(i)-x(BS_num*2+1));
end

x=xp + K*(pseudodis-hx);
P=Pp - K*H*Pp;

x0=x(BS_num*2+3);
y0=x(BS_num*2+4);
z0=x(BS_num*2+5);
end

