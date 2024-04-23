clear all
clc

theta1=0;
theta2=pi/2;
E1=111000000000;
E2=8000000000;
nu12=0.33;
G12=3000000000;
nu21=nu12*E2/E1;
q0=-1000;
a=0.5;
b=0.25;

%x=linspace(0,0.5,100);
%y=linspace(0,25,100);
z=[-0.001, -0.00051, -0.00049, -0.0001, 0, 0.0001, 0.00049, 0.00051, 0.001];
x=a/2;
y=b/2;

R1=[cos(theta1).^2 sin(theta1).^2 cos(theta1).*sin(theta1); sin(theta1).^2 cos(theta1).^2 -cos(theta1).*sin(theta1); -2.*sin(theta1).*cos(theta1) 2.*cos(theta1).*sin(theta1) cos(theta1).^2-sin(theta1).^2];

Qstar1=[E1/(1-nu12*nu21) nu12*E2/(1-nu12*nu21) 0; nu12*E2/(1-nu12*nu21) E2/(1-nu12*nu21) 0; 0 0 G12];

Q1=transpose(R1)*Qstar1*R1;


R2=[cos(theta2).^2 sin(theta2).^2 cos(theta2).*sin(theta2); sin(theta2).^2 cos(theta2).^2 -cos(theta2).*sin(theta2); -2.*sin(theta2).*cos(theta2) 2.*cos(theta2).*sin(theta2) cos(theta2).^2-sin(theta2).^2];

Qstar2=[E1/(1-nu12*nu21) nu12*E2/(1-nu12*nu21) 0; nu12*E2/(1-nu12*nu21) E2/(1-nu12*nu21) 0; 0 0 G12];

Q2=transpose(R2)*Qstar2*R2;

Apiastra1=Q1*0.0005+Q2*0.0005+Q2*0.0005+Q1*0.0005;
Apiastra2=Q2*0.0005+Q1*0.0005+Q1*0.0005+Q2*0.0005;

Bpiastra1=1/2.*(Q1.*(0.0005^2-0.001^2)+Q2.*(0-0.0005^2)+Q2.*(0.0005^2-0)+Q1.*(0.001^2-0.0005^2));

D1=1/3.*(Q1*(-0.0005^3-(-0.001)^3)+Q2*(0-(-0.0005)^3)+Q2*(0.0005^3-0)+Q1*(0.001^3-0.0005^3));
D2=1/3.*(Q2*(-0.0005^3-(-0.001)^3)+Q1*(0-(-0.0005)^3)+Q1*(0.0005^3-0)+Q2*(0.001^3-0.0005^3));

Spostamentomax1=q0/(D1(1,1)*(pi/a)^4+2*(D1(1,2)+2*D1(3,3))*(pi/a)^2*(pi/b)^2+D1(2,2)*(pi/b)^4);
Spostamentomax2=q0/(D2(1,1)*(pi/a)^4+2*(D2(1,2)+2*D2(3,3))*(pi/a)^2*(pi/b)^2+D2(2,2)*(pi/b)^4);

epsilonxx1=z.*((pi/a)^2*Spostamentomax1.*sin(pi* x./a).*sin(pi* y./b)); %piastra a
epsilonyy1=z.*((pi/b)^2*Spostamentomax1.*sin(pi* x./a).*sin(pi* y./b));
gammaxy1=z.*(-2*(pi/a)*(pi/b)*Spostamentomax1.*cos(pi* x./a).*cos(pi* y./b));

% figure 
% plot(epsilonxx1,z)
% grid on
% xlabel('epsilonxx')
% ylabel('z')
% 
% figure 
% plot(epsilonyy1,z)
% grid on
% xlabel('epsilonyy')
% ylabel('z')
% 
% figure 
% plot(gammaxy1,z)
% grid on
% xlabel('gammaxy')
% ylabel('z')

epsilonxx2=z.*((pi/a)^2*Spostamentomax2.*sin(pi* x./a).*sin(pi* y./b)); %piastra b
epsilonyy2=z.*((pi/b)^2*Spostamentomax2.*sin(pi* x./a).*sin(pi* y./b));
gammaxy2=z.*(-2*(pi/a)*(pi/b)*Spostamentomax2.*cos(pi* x./a).*cos(pi* y./b));

% figure 
% plot(epsilonxx2,z)
% grid on
% xlabel('epsilonxx')
% ylabel('z')
% 
% figure 
% plot(epsilonyy2,z)
% grid on
% xlabel('epsilonyy')
% ylabel('z')
% 
% figure 
% plot(gammaxy2,z)
% grid on
% xlabel('gammaxy')
% ylabel('z')

%sigmaxx1=[Q1(1,1); Q1(1,2); Q1(1,3)].*[epsilonxx1; epsilonyy1; gammaxy1];  %per prendere la prima riga di un vettore Q1(1,:)
%sigmayy1=[Q1(2,:)]*[epsilonxx1 ;epsilonyy1; gammaxy1];  % i ; non sono facoltativi quando separo i vettori
%tauxy1=[Q1(3,:)]*[epsilonxx1 ;epsilonyy1 ;gammaxy1];

%calcolo le sigmaxx lungo i quattro strati della prima piastra 
sigmaxx1piastra1=Q1(1,1)*epsilonxx1+Q1(1,2)*epsilonyy1+Q1(1,3)*gammaxy1;
sigmaxx2piastra1=Q2(1,1)*epsilonxx1+Q2(1,2)*epsilonyy1+Q2(1,3)*gammaxy1;
sigmaxx3piastra1=Q2(1,1)*epsilonxx1+Q2(1,2)*epsilonyy1+Q2(1,3)*gammaxy1;
sigmaxx4piastra1=Q1(1,1)*epsilonxx1+Q1(1,2)*epsilonyy1+Q1(1,3)*gammaxy1;

% figure 
% plot(sigmaxx1piastra1(1:2),z(1:2))
% xlabel('sigmaxx')
% ylabel('z')
% grid on
% hold on 
% plot(sigmaxx2piastra1(3:5),z(3:5))
% hold on
% plot(sigmaxx3piastra1(5:7),z(5:7))
% hold on
% plot(sigmaxx4piastra1(8:9),z(8:9))
% hold on
% yline(-0.0005,'LineWidth',1.5)
% yline(-0.001,'LineWidth',1.5)
% yline(0,'LineWidth',1.5)
% yline(0.0005,'LineWidth',1.5)
% yline(0.001,'LineWidth',1.5)

%calcolo le sigmayy lungo i 4 strati della prima piastra
% sigmayy1piastra1=Q1(1,2)*epsilonxx1+Q1(2,2)*epsilonyy1+Q1(2,3)*gammaxy1;
% sigmayy2piastra1=Q2(1,2)*epsilonxx1+Q2(2,2)*epsilonyy1+Q2(2,3)*gammaxy1;
% sigmayy3piastra1=Q2(1,2)*epsilonxx1+Q2(2,2)*epsilonyy1+Q2(2,3)*gammaxy1;
% sigmayy4piastra1=Q1(1,2)*epsilonxx1+Q1(2,2)*epsilonyy1+Q1(2,3)*gammaxy1;
% 
% figure 
% plot(sigmayy1piastra1(1:2),z(1:2))
% xlabel('sigmayy')
% ylabel('z')
% grid on
% hold on 
% plot(sigmayy2piastra1(3:5),z(3:5))
% hold on
% plot(sigmayy3piastra1(5:7),z(5:7))
% hold on
% plot(sigmayy4piastra1(8:9),z(8:9))
% hold on
% yline(-0.0005,'LineWidth',1.5)
% yline(-0.001,'LineWidth',1.5)
% yline(0,'LineWidth',1.5)
% yline(0.0005,'LineWidth',1.5)
% yline(0.001,'LineWidth',1.5)

tauxy1piastra1=Q1(3,1)*epsilonxx1+Q1(3,2)*epsilonyy1+Q1(3,3)*gammaxy1;
tauxy2piastra1=Q2(3,1)*epsilonxx1+Q2(3,2)*epsilonyy1+Q2(3,3)*gammaxy1;
tauxy3piastra1=Q2(3,1)*epsilonxx1+Q2(3,2)*epsilonyy1+Q2(3,3)*gammaxy1;
tauxy4piastra1=Q1(3,1)*epsilonxx1+Q1(3,2)*epsilonyy1+Q1(3,3)*gammaxy1;

% figure 
% plot(tauxy1piastra1(1:2),z(1:2))
% xlabel('tauxy')
% ylabel('z')
% grid on
% hold on 
% plot(tauxy2piastra1(3:5),z(3:5))
% hold on
% plot(tauxy3piastra1(5:7),z(5:7))
% hold on
% plot(tauxy4piastra1(8:9),z(8:9))
% hold on
% yline(-0.0005,'LineWidth',1.5)
% yline(-0.001,'LineWidth',1.5)
% yline(0,'LineWidth',1.5)
% yline(0.0005,'LineWidth',1.5)
% yline(0.001,'LineWidth',1.5)

%calcolo analogamente le sigma e le tau sulla piastra 2
sigmaxx1piastra2=Q2(1,1)*epsilonxx2+Q2(1,2)*epsilonyy2+Q2(1,3)*gammaxy2;
sigmaxx2piastra2=Q1(1,1)*epsilonxx2+Q1(1,2)*epsilonyy2+Q1(1,3)*gammaxy2;
sigmaxx3piastra2=Q1(1,1)*epsilonxx2+Q1(1,2)*epsilonyy2+Q1(1,3)*gammaxy2;
sigmaxx4piastra2=Q2(1,1)*epsilonxx2+Q2(1,2)*epsilonyy2+Q2(1,3)*gammaxy2;

% figure 
% plot(sigmaxx1piastra2(1:2),z(1:2))
% xlabel('sigmaxx')
% ylabel('z')
% grid on
% hold on 
% plot(sigmaxx2piastra2(3:5),z(3:5))
% hold on
% plot(sigmaxx3piastra2(5:7),z(5:7))
% hold on
% plot(sigmaxx4piastra2(8:9),z(8:9))
% hold on
% yline(-0.0005,'LineWidth',1.5)
% yline(-0.001,'LineWidth',1.5)
% yline(0,'LineWidth',1.5)
% yline(0.0005,'LineWidth',1.5)
% yline(0.001,'LineWidth',1.5)


sigmayy1piastra2=Q2(1,2)*epsilonxx2+Q2(2,2)*epsilonyy2+Q2(2,3)*gammaxy2;
sigmayy2piastra2=Q1(1,2)*epsilonxx2+Q1(2,2)*epsilonyy2+Q1(2,3)*gammaxy2;
sigmayy3piastra2=Q1(1,2)*epsilonxx2+Q1(2,2)*epsilonyy2+Q1(2,3)*gammaxy2;
sigmayy4piastra2=Q2(1,2)*epsilonxx2+Q2(2,2)*epsilonyy2+Q2(2,3)*gammaxy2;

% figure 
% plot(sigmayy1piastra2(1:2),z(1:2))
% xlabel('sigmayy')
% ylabel('z')
% grid on
% hold on 
% plot(sigmayy2piastra2(3:5),z(3:5))
% hold on
% plot(sigmayy3piastra2(5:7),z(5:7))
% hold on
% plot(sigmayy4piastra2(8:9),z(8:9))
% hold on
% yline(-0.0005,'LineWidth',1.5)
% yline(-0.001,'LineWidth',1.5)
% yline(0,'LineWidth',1.5)
% yline(0.0005,'LineWidth',1.5)
% yline(0.001,'LineWidth',1.5)

tauxy1piastra2=Q2(3,1)*epsilonxx2+Q2(3,2)*epsilonyy2+Q2(3,3)*gammaxy2;
tauxy2piastra2=Q1(3,1)*epsilonxx2+Q1(3,2)*epsilonyy2+Q1(3,3)*gammaxy2;
tauxy3piastra2=Q1(3,1)*epsilonxx2+Q1(3,2)*epsilonyy2+Q1(3,3)*gammaxy2;
tauxy4piastra2=Q2(3,1)*epsilonxx2+Q2(3,2)*epsilonyy2+Q2(3,3)*gammaxy2;

% figure 
% plot(tauxy1piastra2(1:2),z(1:2))
% xlabel('tauxy')
% ylabel('z')
% grid on
% hold on 
% plot(tauxy2piastra2(3:5),z(3:5))
% hold on
% plot(tauxy3piastra2(5:7),z(5:7))
% hold on
% plot(tauxy4piastra2(8:9),z(8:9))
% hold on
% yline(-0.0005,'LineWidth',1.5)
% yline(-0.001,'LineWidth',1.5)
% yline(0,'LineWidth',1.5)
% yline(0.0005,'LineWidth',1.5)
% yline(0.001,'LineWidth',1.5)