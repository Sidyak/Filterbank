%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K. Radmacher, 12.06.15
%
% compare ISO filterbank and own derivated filterbank
% INFO: file "Window.m" includes ISO coef.
%       file "Coef_32band_optimal_prototype.mat" includes own calculated 
%       optimized filter coef.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
format long
format compact

load('Coef_32band_optimal_prototype.mat','-ascii'); % load prototype coef.
h=Coef_32band_optimal_prototype;
Nfft=512;
% calc. filter coef.
C=zeros(Nfft,1);    
n_ind=1;
for j_ind=0:7
    for k_ind=0:63
        C(n_ind)=2*(-1).^j_ind*h((k_ind+1)+64*(j_ind));
        n_ind=n_ind+1;
    end
end
xmin=1;
xmax=512;
Ciso=Window()';     % load ISO coefs
f1=figure(1);  
subplot(311), stem(Ciso,'k')
xlim([xmin xmax])
title('\it{a) MPEG-1 ISO}','Interpreter', 'Latex', 'FontSize', 16);
xlabel('$$n\rightarrow$$','Interpreter', 'Latex', 'FontSize', 14);
ylabel('$$C_{ISO}[n]$$','Interpreter', 'Latex', 'FontSize', 14);

subplot(312), stem(C,'k')
xlim([xmin xmax])
title('\it{b) Aus eigener Herleitung}','Interpreter', 'Latex', 'FontSize', 16);
xlabel('$$n\rightarrow$$','Interpreter', 'Latex', 'FontSize', 14);
ylabel('$$C_{O}[n]$$','Interpreter', 'Latex', 'FontSize', 14);

subplot(313), stem((Ciso-C),'k')
xlim([xmin xmax])
title('\it{c) Abweichung}','Interpreter', 'Latex', 'FontSize', 16);
xlabel('$$n\rightarrow$$','Interpreter', 'Latex', 'FontSize', 14);
ylabel('$$C_{ISO}[n]-C_{O}[n]$$','Interpreter', 'Latex', 'FontSize', 14);

disp('done!')