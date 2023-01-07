clear all; close all;
disp('Three points Measurement in 3d space.')
BS_num=12;
% BS_loc=[0,15,0.5;0,35,1.5;15,0,2.5;35,0,0.5;50,15,1.5;50,35,2.5;15,50,0.5;35,50,1.5];
BS_loc=[0,12,0.5;0,25,1.5;0,38,2.5;12,0,2.5;25,0,0.5;38,0,1.5;50,12,1.5;59,25,0.5;50,38,2.5;12,50,0.5;25,50,1.5;38,50,2.5];
a = 50;
x = BS_loc(:,1);
y = BS_loc(:,2);
z = BS_loc(:,3);


target = [25,25];
x0=linspace(0,50,100);
y0=linspace(0,50,100);
z0=ones(1,100);

HDOP =zeros(100,100);
h=zeros(BS_num,3);

for xx=1:1:length(x0)
    for yy=1:1:length(y0)

        for m=1:BS_num
            r=sqrt((x0(xx)-x(m))^2+(y0(yy)-y(m))^2+(z0(1)-z(m))^2);
            h(m,1)=(x(m)-x0(xx))/r;
            h(m,2)=(y(m)-y0(yy))/r;
            h(m,3)=(z(m)-z0(1))/r;
        end
       
        G = inv(h'*h);
        HDOP(xx,yy) = sqrt((G(1,1)+G(2,2)));

    end
end




figure(1);
mean_DOP=mean(HDOP,2)*ones(1,100);
pcolor(linspace(0,50,100),linspace(0,50,100),HDOP);
xlabel('x坐标/m');ylabel('y坐标/m');
shading interp;
colorbar;colormap('jet');
hold on;
plot(x,y,'r*');
title('HDOP Value');

