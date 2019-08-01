% % TODO:
% %    0. Add cuts per antenna height
% %    1. Add integration over multiple heights - for complex target with
% %      multiple contributers.
% %    2.  

close all
clear
clc

%% Settings
freqs  =  (77:0.02:81)*1e9;   % (62:0.02:65)*1e9;  % [Hz]

h1_0 = 1.5;  % [m]
h2 = 1.5;    % [m]
D  = 10;  % [m]

N_ants = 20;
d_ants = 1.6e-3;  % 1.6e-3;  % [m]

Pol_ang = 0;  % To add treatment!

Ant_BW = 180;  % [deg]
eps_surface = 10;

theta_vec = (-80:1:80)*pi/180;
useApodization = 0;

%% Calc
k  = 2*pi*freqs/3e8;

Smat_tx_prll = zeros(N_ants,length(freqs));
Smat_tx_ppdc = zeros(N_ants,length(freqs));

for n = 1:N_ants
    h1 = h1_0 + n*d_ants;
    L0 = sqrt(D^2 + (h2-h1)^2);

    theta0 = atan((h1-h2)/D)*180/pi;

    d1 = D*h1/(h1 + h2);
    d2 = D - d1;
    L1 = sqrt(h1^2 + d1^2);
    L2 = sqrt(h2^2 + d2^2);

    theta2 = atan(d1/h1)*180/pi;
    BW_loss = -3*(abs((90-theta2)-theta0)/Ant_BW).^2;    % [dB]


    theta0 = theta0*pi/180;
    theta2 = theta2*pi/180;

    n_srf = sqrt(eps_surface);
%     R_prll = abs((cos(theta2) - n_srf*sqrt(1 - (1/n_srf*sin(theta2))^2)) / ...
%                     (cos(theta2) + n_srf*sqrt(1 - (1/n_srf*sin(theta2))^2)));
%     R_ppdc = abs((sqrt(1 - (1/n_srf*sin(theta2))^2) - n_srf*cos(theta2)) / ...
%                     (sqrt(1 - (1/n_srf*sin(theta2))^2) + n_srf*cos(theta2))); 

    R_prll = ((cos(theta2) - n_srf*sqrt(1 - (1/n_srf*sin(theta2))^2)) / ...
                    (cos(theta2) + n_srf*sqrt(1 - (1/n_srf*sin(theta2))^2)));
    R_ppdc = ((sqrt(1 - (1/n_srf*sin(theta2))^2) - n_srf*cos(theta2)) / ...
                    (sqrt(1 - (1/n_srf*sin(theta2))^2) + n_srf*cos(theta2))); 


    E__prll = exp(1i*k*L0) + R_prll*10^(BW_loss/10)*exp(1i*k*(L1+L2));
    E__ppdc = exp(1i*k*L0) + R_ppdc*10^(BW_loss/10)*exp(1i*k*(L1+L2));
    
    Smat_tx_prll(n,:) = E__prll;
    Smat_tx_ppdc(n,:) = E__ppdc;

end


%% "2D Imaging"
% win = CreateFreqWindow(freqs,{'nt4',[0 1]},1);
win = ones(1,length(freqs));

N_fft= 2^11;
Smat_tx_prll_t = fftshift(fft(win.*Smat_tx_prll,N_fft,2)/(N_fft)); 
Smat_tx_ppdc_t = fftshift(fft(win.*Smat_tx_ppdc,N_fft,2)/(N_fft)); 

t_end = 1/(freqs(2)-freqs(1));  
dt = t_end/N_fft;
t = -t_end/2:dt:t_end/2-dt;


Rx_pos_vec  = (d_ants*(0:N_ants-1));
if useApodization
    margins_y = 0.5;
    rx_vec_forwin = (Rx_pos_vec - mean(Rx_pos_vec)); 
    rx_vec_forwin = rx_vec_forwin / max(rx_vec_forwin);
    [w,msg] = gencoswin('hamming',50);
    win_rx  = interp1(linspace(-1/margins_y,1/margins_y,50),w,rx_vec_forwin);
else
    win_rx = ones(1,length(Rx_pos_vec));
end

AzRecMat = [];
for ll = 1:length(theta_vec)
    AzRecMat(ll,:) = win_rx.*exp(2i*pi*mean(freqs)/3e8*Rx_pos_vec*sin(theta_vec(ll)));   
end

Smat_t_theta_prll = AzRecMat*Smat_tx_prll_t;
Smat_t_theta_ppdc = AzRecMat*Smat_tx_ppdc_t;


%% Plotting
figure; subplot(3,1,1);
stem(0,h1,'bv','BaseValue',0,'MarkerSize',15); grid on; hold all
stem(D,h2,'ro','BaseValue',0,'MarkerSize',15); 
hold all; plot([0 , D],[h1 , h2],'c--','LineWidth',0.2);
plot([0,d1],[h1,0],'m--',[d1,D],[0,h2],'m--','LineWidth',0.2);
xlim([-0.5 D*1.5]);  ylim([0 , max([h1,h2])*1.5]); 
legend('Antenna','Target');
xlabel('range [m]'); ylabel('height [m]');

subplot(3,1,2);
plot(freqs*1e-9,db(Smat_tx_prll));  grid on;
xlabel('f[GHz]'); ylabel('Pr_norm [dB]');
title('Single-reflection interference, H-pol.');

subplot(3,1,3);
plot(freqs*1e-9,db(Smat_tx_ppdc));  grid on;
xlabel('f[GHz]'); ylabel('Pr_norm [dB]');
title('Single-reflection interference, V-pol.');


figure; subplot(2,1,1);
stem(0,h1,'bv','BaseValue',0,'MarkerSize',15); grid on; hold all
stem(D,h2,'ro','BaseValue',0,'MarkerSize',15); 
hold all; plot([0 , D],[h1 , h2],'c--','LineWidth',0.2);
plot([0,d1],[h1,0],'m--',[d1,D],[0,h2],'m--','LineWidth',0.2);
xlim([-0.5 D*1.5]);  ylim([0 , max([h1,h2])*1.5]); 
legend('Antenna','Target');
xlabel('range [m]'); ylabel('height [m]');

subplot(2,1,2);
plot(t*3e8/2,db(Smat_tx_prll_t),t*3e8/2,db(Smat_tx_ppdc_t));  grid on;
xlabel('x [m]'); ylabel('Pr_norm [dB]');
title('Time domain');
legend('H-pol','V-pol');



figure; subplot(1,2,1);
surface(t*3e8/2,theta_vec*180/pi,db(Smat_t_theta_prll) - max(max(db(Smat_t_theta_prll))),'EdgeColor','none');
xlabel('x [m]'); ylabel('\theta [deg]');
title('H-pol');
set(gca,'CLim',[-30,0]);

subplot(1,2,2);
surface(t*3e8/2,theta_vec*180/pi,db(Smat_t_theta_ppdc) - max(max(db(Smat_t_theta_ppdc))),'EdgeColor','none');
xlabel('x [m]'); ylabel('\theta [deg]');
title('V-pol');
set(gca,'CLim',[-30,0]);


%% Adding Polarization
Smat_t_theta = cos(Pol_ang*pi/180)*Smat_t_theta_ppdc + sin(Pol_ang*pi/180)*Smat_t_theta_prll;


