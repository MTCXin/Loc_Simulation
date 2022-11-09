function [x0,y0,z0] = EKF_direct(pseudodis,BS_num,BS_loc,Ts,Qk)
%EKF_direct 

persistent A Q R x P firstRun

if isempty(firstRun)

    x=[zeros(1,BS_num),zeros(1,BS_num),0,0,mean(BS_loc)];   % hidden state:[deltaTk*n,fk*n, deltaTk, fk, x, y, z]( 2n+5 dims)

    
%     |1 ... Ts
%     |0 1 ...0 Ts
%     |
%     |
%     |
    A=eye(BS_num*2+5);
    for i = 1:BS_num
        A(i,i+BS_num)=Ts;
    end
    A(2*BS_num+1,2*BS_num+2)=Ts;
    
    firstRun=1;
end

x0=0;
y0=0;
z0=0;
end

