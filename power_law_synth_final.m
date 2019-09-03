close all
clear all
clc

Gamma = 2.675987E8;
Del = 50.9*10^-3;
del = 35.1*10^-3;

bvec = [0.049	-0.919	-0.391;0.726	0.301	-0.618;-0.683	0.255	-0.684;...
    0.845	-0.502	-0.186;-0.730	-0.619	-0.288;-0.051	0.039	0.998;...
    -0.018	0.871	-0.491;-0.444	0.494	0.747;-0.989	-0.086	-0.116;...
    -0.470	-0.855	0.221;0.412	0.400	0.819;-0.552	0.790	-0.267;...
    -0.123	-0.477	0.871;-0.848	0.141	0.510;-0.341	-0.788	-0.512;...
    0.361	-0.529	0.768;-0.472	0.850	0.234;-0.856	-0.481	0.189;...
    0.797	0.162	0.582;0.467	-0.009	-0.884;0.013	0.998	-0.056;...
    0.882	-0.387	0.267;0.017	-0.536	-0.844;-0.442	-0.651	0.617;...
    0.365	-0.058	0.929;0.977	-0.004	-0.213;-0.406	-0.902	-0.145;...
    -0.627	0.614	0.479;-0.354	0.772	-0.528;-0.658	-0.472	-0.586;...
    0.423	0.322	-0.847;0.212	-0.754	-0.622;0.912	-0.104	0.398;...
    -0.311	0.947	-0.077;0.679	0.632	-0.374;0.135	-0.286	0.949;...
    -0.647	0.230	0.727;0.904	0.397	0.158;-0.757	0.647	-0.087;...
    0.143	0.284	0.948;0.233	0.894	-0.382;0.664	-0.531	0.527;...
    0.157	0.710	0.686;-0.895	-0.214	0.392;0.594	0.080	0.801;...
    -0.006	0.688	-0.726;0.665	0.746	0.024;-0.277	0.276	0.920;...
    -0.962	0.268	0.046;-0.133	-0.970	-0.202;0.790	-0.405	-0.461;...
    -0.194	-0.193	0.962;-0.236	0.952	0.195;-0.884	-0.272	-0.379;...
    0.463	-0.307	0.831;0.700	0.066	-0.711;-0.200	0.928	-0.314;...
    0.550	0.705	0.449;-0.670	0.727	0.153;0.237	0.722	-0.650;...
    0.960	0.260	-0.100];
bvec60 = bvec(1:30,:);
bvec60 = [bvec60; -bvec60];
bvec30 = bvec(1:15,:);
bvec30 = [bvec30; -bvec30];

for i = 1:21
    b(1+(i-1)*length(bvec30):i*length(bvec30),1) = 0.5*(i-1)*10^-3;
end
G_mag = sqrt(b./Gamma^2/del.^2./(Del-del/3));

g = repmat(bvec30,21,1);

xdata(:,1:3) = g;
xdata(:,4) = G_mag*10^6;
xdata(:,5) = Del;
xdata(:,6) = del;

sigma = 1/250;
% sigma = 0;

% De_par = randn*0.2+2;
De_par = 2*1000;
% f1 = randn*0.1+0.65;
f1 = 0.65;
f2 = 1-f1;
% Da = randn*0.2+2;
Da = 2;
% kappa = randn*19 + 11;
kappa = 11;
% De_per = [randn*0.1+0.25, randn*0.1+0.5, randn*0.1+0.75]*1000;
De_per = [0.25, 0.5, 0.75]*1000;
eta = [0, 1, 1.5];
% eta = [0, 1, 1.5, 2, 2.5, 3];
OD = 2/pi*atan(1./kappa);

g=[xdata(:,1) xdata(:,2) xdata(:,3)];
G_mag= (xdata(:,4)./1000000); %(Tesla/micro meter)
Del=xdata(:,5);          %(in seconds)
del=xdata(:,6);          %(in seconds)
Gamma = 2.675987E8;
td=(Del-(del/3));        %(in seconds)
b=(td).*((Gamma.*del.*G_mag).^2); %(sec/mircometer^2)

D_intra = Da*1000;
% D_par = De_par*1000;

sinT = sin(0);
cosT = cos(0);
sinP = sin(0);
cosP = cos(0);
n = [cosP * sinT sinP * sinT  cosT]; % fiber direction

% the roots of the derivatives of the first kind Bessel function
am = [1.84118307861360, 5.33144196877749, 8.53631578218074, 11.7060038949077, ...
    14.8635881488839,  18.0155278304879, 21.1643671187891, 24.3113254834588,  ...
    27.4570501848623, 30.6019229722078, 33.7461812269726,  36.8899866873805, ...
    40.0334439409610, 43.1766274212415,  46.3195966792621, 49.4623908440429, ...
    52.6050411092602,  55.7475709551533, 58.8900018651876, 62.0323477967829,  ...
    65.1746202084584, 68.3168306640438, 71.4589869258787,  74.6010956133729, ...
    77.7431620631416, 80.8851921057280,  84.0271895462953, 87.1691575709855, ...
    90.3110993488875,  93.4530179063458, 96.5949155953313, 99.7367932203820, ...
    102.878653768715, 106.020498619541, 109.162329055405,  112.304145672561,...
    115.445950418834, 118.587744574512,  121.729527118091, 124.871300497614, ...
    128.013065217171,  131.154821965250, 134.296570328107, 137.438311926144, ...
    140.580047659913,143.721775748727, 146.863498476739,  150.005215971725, ...
    153.146928691331, 156.288635801966,  159.430338769213, 162.572038308643,...
    165.713732347338,  168.855423073845, 171.997111729391, 175.138794734935, ...
    178.280475036977, 181.422152668422, 184.563828222242,  187.705499575101];

% from Genu, Aboitiz 1992
bin = [2, 13, 26, 17, 13, 7, 5, 7, 2, 1, 1, 1, 1, 1, 1, 1]';
R = [0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.6]';

histR = [0.2,0.2, 0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,...
0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,...
0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,...
0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,...
1,1,1,1,1,1,1,1,1,1,1,1,1,...
1.2,1.2,1.2,1.2,1.2,1.2,1.2,...
1.4,1.4,1.4,1.4,1.4,...
1.6,1.6,1.6,1.6,1.6,1.6,1.6,1.8,1.8,2,2.2,2.4,2.6,2.8,3,3.6];

% edges = 0:0.1:4;
% figure,histogram(histR,edges)
% title('histogram of the axon diameters')
% set(gca, 'FontSize', 15)
 
b_delta = -0.5;

% the signal attenuation from Ven Gelderen equation for
% different axon diamters weighted by r^2

for l = 1:length(De_per)
    %%%%%%%%%%%%%%%%%%%%%%%% Extra Cellular ec %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    signal_ec = slow_exchange_b_tensor_SynthMeasWatsonHinderedDiffusion_PGSE([De_par De_per(l) kappa], g, b, n',b_delta);
    nG4=signal_ec(:);
    
    for j = 1:length(eta)
        for i = 1:length(R)
            %%%%%%%%%%%%%%%  IntraCellular ic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            signal_ic = b_tensor_R_SynthMeasWatsonSHCylNeuman_PGSE([D_intra R(i)*eta(j) kappa], g,G_mag,Del, del, n', am, b_delta);
            nG3=signal_ic(:);
            
            S_ic(i,:) = nG3 * R(i)^2 * bin(i);
        end
        S_ic_t(:,j) = sum(S_ic)'/sum(bin.*R.^2);
        S(:,j,l) = [S_ic_t(:,j) nG4]*[f1 f2]';
    end
end

alpha = 1;
beta = 0.2;
gama = 0;
opts = optimset('Display','off');
lb=[0 0];
ub=[2 2];
x2 = [alpha, beta];
[C, ia, ic] = unique(b);
p = length(C);

for l = 1:length(De_per)
    for j = 1:length(eta)
        S_noise = sqrt((squeeze(S(:,j,l)) + normrnd(0,sigma,length(S(:,j,l)),1)).^2 ...
            + normrnd(0,sigma,length(S(:,j,l)),1).^2);
        % we have to find the average signal value per shell.
        for i = 1:p-1
            S_shell(i,j,l) = mean(S_noise(ia(i):ia(i+1)-1));
        end
        S_shell(p,j,l) = mean(S_noise(ia(p):end));
        [xf(:,j,l),resnorm,~,exitflag,output] = ...
            lsqcurvefit(@power_law2,x2,C(15:p)*1000,S_shell(15:p,j,l),lb,ub,opts);
    end
end

figure
ha = tight_subplot(3,3,[.06 .04],[.08 .05],[.06 .02]);
axes(ha(1));
plot(1./C(15:p)/1000,S_shell(15:p,1,1),'bo-', 'LineWidth', 1)
hold on
signal = power_law2(xf(:,1,1),C*1000);
plot(1./C(15:p)/1000,signal(15:p),'r-.', 'LineWidth', 1)
txt = ['\alpha = ' sprintf('%.2f',xf(1,1,1))];
text(0.105,0.063 ,txt,'FontSize', 20)
txt = ['\beta = ' sprintf('%.2f',xf(2,1,1))];
text(0.105,0.053 ,txt,'FontSize', 20)
ylim([0.02,0.07])
set(gca, 'FontSize', 20)
ylabel('S(b)/S_0   (\eta = 0)')
title('D_e^\perp = 0.25 (\mum^2/ms)')
set(gca,'xtick',[])
legend('direction averaged signal','power-law fit')

axes(ha(2));
plot(1./C(15:p)/1000,S_shell(15:p,1,2),'bo-', 'LineWidth', 1)
hold on
signal = power_law2(xf(:,1,2),C*1000);
plot(1./C(15:p)/1000,signal(15:p),'r-.', 'LineWidth', 1)
txt = ['\alpha = ' sprintf('%.2f',xf(1,1,2))];
text(0.105,0.063 ,txt,'FontSize', 20)
txt = ['\beta = ' sprintf('%.2f',xf(2,1,2))];
text(0.105,0.053 ,txt,'FontSize', 20)
ylim([0.02,0.07])
set(gca, 'FontSize', 20)
title('D_e^\perp = 0.5 (\mum^2/ms)')
set(gca,'xtick',[])
set(gca,'ytick',[])

axes(ha(3));
plot(1./C(15:p)/1000,S_shell(15:p,1,3),'bo-', 'LineWidth', 1)
hold on
signal = power_law2(xf(:,1,3),C*1000);
plot(1./C(15:p)/1000,signal(15:p),'r-.', 'LineWidth', 1)
txt = ['\alpha = ' sprintf('%.2f',xf(1,1,3))];
text(0.105,0.063 ,txt,'FontSize', 20)
txt = ['\beta = ' sprintf('%.2f',xf(2,1,3))];
text(0.105,0.053 ,txt,'FontSize', 20)
ylim([0.02,0.07])
set(gca, 'FontSize', 20)
title('D_e^\perp = 0.75 (\mum^2/ms)')
set(gca,'xtick',[])
set(gca,'ytick',[])

% axes(ha(4));
% plot(1./C(15:p)/1000,sum(S_shell(15:p,1,:),3)/3,'bo-', 'LineWidth', 1)
% [xf3,~,~,~,~] = lsqcurvefit(@power_law2,x2,C(15:p)*1000,sum(S_shell(15:p,1,:),3)/3,lb,ub,opts);
% hold on
% signal = power_law2(xf3,C*1000);
% plot(1./C(15:p)/1000,signal(15:p),'r-.', 'LineWidth', 1)
% txt = ['\alpha = ' sprintf('%.2f',xf3(1))];
% text(0.105,0.063 ,txt,'FontSize', 20)
% txt = ['\beta = ' sprintf('%.2f',xf3(2))];
% text(0.105,0.053 ,txt,'FontSize', 20)
% ylim([0.02,0.07])
% set(gca, 'FontSize', 20)
% title('Mean signal of three D_e^\perp')
% set(gca,'xtick',[])
% set(gca,'ytick',[])

axes(ha(4));
plot(1./C(15:p)/1000,S_shell(15:p,2,1),'bo-', 'LineWidth', 1)
hold on
signal = power_law2(xf(:,2,1),C*1000);
plot(1./C(15:p)/1000,signal(15:p),'r-.', 'LineWidth', 1)
txt = ['\alpha = ' sprintf('%.2f',xf(1,2,1))];
text(0.105,0.063 ,txt,'FontSize', 20)
txt = ['\beta = ' sprintf('%.2f',xf(2,2,1))];
text(0.105,0.053 ,txt,'FontSize', 20)
ylim([0.02,0.07])
set(gca, 'FontSize', 20)
ylabel('S(b)/S_0   (\eta = 1)')
set(gca,'xtick',[])

axes(ha(5));
plot(1./C(15:p)/1000,S_shell(15:p,2,2),'bo-', 'LineWidth', 1)
hold on
signal = power_law2(xf(:,2,2),C*1000);
plot(1./C(15:p)/1000,signal(15:p),'r-.', 'LineWidth', 1)
txt = ['\alpha = ' sprintf('%.2f',xf(1,2,2))];
text(0.105,0.063 ,txt,'FontSize', 20)
txt = ['\beta = ' sprintf('%.2f',xf(2,2,2))];
text(0.105,0.053 ,txt,'FontSize', 20)
ylim([0.02,0.07])
set(gca, 'FontSize', 20)
set(gca,'xtick',[])
set(gca,'ytick',[])

axes(ha(6));
plot(1./C(15:p)/1000,S_shell(15:p,2,3),'bo-', 'LineWidth', 1)
hold on
signal = power_law2(xf(:,2,3),C*1000);
plot(1./C(15:p)/1000,signal(15:p),'r-.', 'LineWidth', 1)
txt = ['\alpha = ' sprintf('%.2f',xf(1,2,3))];
text(0.105,0.063 ,txt,'FontSize', 20)
txt = ['\beta = ' sprintf('%.2f',xf(2,2,3))];
text(0.105,0.053 ,txt,'FontSize', 20)
ylim([0.02,0.07])
set(gca, 'FontSize', 20)
set(gca,'xtick',[])
set(gca,'ytick',[])

% axes(ha(8));
% plot(1./C(15:p)/1000,sum(S_shell(15:p,2,:),3)/3,'bo-', 'LineWidth', 1)
% [xf3,~,~,~,~] = lsqcurvefit(@power_law2,x2,C(15:p)*1000,sum(S_shell(15:p,2,:),3)/3,lb,ub,opts);
% hold on
% signal = power_law2(xf3,C*1000);
% plot(1./C(15:p)/1000,signal(15:p),'r-.', 'LineWidth', 1)
% txt = ['\alpha = ' sprintf('%.2f',xf3(1))];
% text(0.105,0.063 ,txt,'FontSize', 20)
% txt = ['\beta = ' sprintf('%.2f',xf3(2))];
% text(0.105,0.053 ,txt,'FontSize', 20)
% ylim([0.02,0.07])
% set(gca, 'FontSize', 20)
% set(gca,'xtick',[])
% set(gca,'ytick',[])

axes(ha(7));
plot(1./C(15:p)/1000,S_shell(15:p,3,1),'bo-', 'LineWidth', 1)
hold on
signal = power_law2(xf(:,3,1),C*1000);
plot(1./C(15:p)/1000,signal(15:p),'r-.', 'LineWidth', 1)
txt = ['\alpha = ' sprintf('%.2f',xf(1,3,1))];
text(0.105,0.063 ,txt,'FontSize', 20)
txt = ['\beta = ' sprintf('%.2f',xf(2,3,1))];
text(0.105,0.053 ,txt,'FontSize', 20)
ylim([0.02,0.07])
set(gca, 'FontSize', 20)
ylabel('S(b)/S_0   (\eta = 1.5)')
xlabel('1/b (\mum^2/ms)')

axes(ha(8));
plot(1./C(15:p)/1000,S_shell(15:p,3,2),'bo-', 'LineWidth', 1)
hold on
signal = power_law2(xf(:,3,2),C*1000);
plot(1./C(15:p)/1000,signal(15:p),'r-.', 'LineWidth', 1)
txt = ['\alpha = ' sprintf('%.2f',xf(1,3,2))];
text(0.105,0.063 ,txt,'FontSize', 20)
txt = ['\beta = ' sprintf('%.2f',xf(2,3,2))];
text(0.105,0.053 ,txt,'FontSize', 20)
ylim([0.02,0.07])
set(gca, 'FontSize', 20)
xlabel('1/b (\mum^2/ms)')
set(gca,'ytick',[])

axes(ha(9));
plot(1./C(15:p)/1000,S_shell(15:p,3,3),'bo-', 'LineWidth', 1)
hold on
signal = power_law2(xf(:,3,3),C*1000);
plot(1./C(15:p)/1000,signal(15:p),'r-.', 'LineWidth', 1)
txt = ['\alpha = ' sprintf('%.2f',xf(1,3,3))];
text(0.105,0.063 ,txt,'FontSize', 20)
txt = ['\beta = ' sprintf('%.2f',xf(2,3,3))];
text(0.105,0.053 ,txt,'FontSize', 20)
ylim([0.02,0.07])
set(gca, 'FontSize', 20)
xlabel('1/b (\mum^2/ms)')
set(gca,'ytick',[])

% axes(ha(12));
% plot(1./C(15:p)/1000,sum(S_shell(15:p,3,:),3)/3,'bo-', 'LineWidth', 1)
% [xf3,~,~,~,~] = lsqcurvefit(@power_law2,x2,C(15:p)*1000,sum(S_shell(15:p,3,:),3)/3,lb,ub,opts);
% hold on
% signal = power_law2(xf3,C*1000);
% plot(1./C(15:p)/1000,signal(15:p),'r-.', 'LineWidth', 1)
% txt = ['\alpha = ' sprintf('%.2f',xf3(1))];
% text(0.105,0.063 ,txt,'FontSize', 20)
% txt = ['\beta = ' sprintf('%.2f',xf3(2))];
% text(0.105,0.053 ,txt,'FontSize', 20)
% ylim([0.02,0.07])
% set(gca, 'FontSize', 20)
% xlabel('1/b (\mum^2/ms)')
% set(gca,'ytick',[])

