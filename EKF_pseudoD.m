function [x0,y0,z0] = EKF_pseudoD(pseudodis,BS_num,BS_loc,Ts,Qk)
%EKF_pseudoD

persistent A Q R x P firstRun

if isempty(firstRun)

    x=[zeros(1,BS_num),zeros(1,BS_num),0,0,mean(BS_loc)];   % hidden state:[deltaTk*n,fk*n, deltaTk, fk, x, y, z]( 2n+5 dims)

    A=eye(BS_num*2+5);
    for i = 1:BS_num
        A(i,i+BS_num)=Ts;
    end
    A(2*BS_num+1,2*BS_num+2)=Ts;

    Q=
    
    firstRun=1;
end
