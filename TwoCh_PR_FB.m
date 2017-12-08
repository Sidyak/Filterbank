%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K.Radmacher, 03.09.14
%
% 2 channel QMF-Bank design with perfect reconstruuction
% reference: N. Fliege, Example 6.2
% INFO: this example using a minimal phase filter prozedur 
%       file "MinimalPhaseFIR.m" needed!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
format long
format compact

Fs=48e3;
delta_stop = 0.0114; 
delta_pass = delta_stop;
fpass = 0.4*Fs/2;
fstop = 0.6*Fs/2; 
frequencys= [fpass fstop];
[N_FIR,fo,mo,w] = remezord( frequencys, [1 0], [delta_pass delta_stop], Fs );
if mod(N_FIR,2)==1
    N_FIR=N_FIR+1; % order must be even
end
N=N_FIR+1;   % number of coef.
B_TP = remez(N_FIR,fo,mo,w);
% create minimal phase filter (all zeros inside the EK)
b  = MinimalPhaseFIR( B_TP, delta_stop, 0);
% create analysis- & synthesis-filter
N=length(b);
h0=b;
for n=0:N-1
    h1(n+1)=(-1)^(N-1-n).*b(N-n);
    g0(n+1)=b(N-n);
    g1(n+1)=(-1)^n*b(n+1);
end
[H0 wh0]=freqz(h0);
[H1 wh1]=freqz(h1);
[G0 wg0]=freqz(g0);
[G1 wg1]=freqz(g1);
f1=figure(1);
subplot(211), plot(wh0/(2*pi),20*log10(abs(H0)),'b'); 
hold on
plot(wh1/(2*pi),20*log10(abs(H1)),'r'); grid
hold off
leg=legend('$$|H_0(e^{j \omega})|$$','$$|H_1(e^{j \omega})|$$');
set(leg,'Interpreter', 'Latex', 'FontSize',14);
xlabel('$$\omega \rightarrow$$','Interpreter', 'Latex', 'FontSize', 14);
ylabel('\it{Amplitude/dB}','Interpreter', 'Latex', 'FontSize', 14);
title('\it{a) Analyse-Filter}','Interpreter', 'Latex', 'FontSize', 16);

subplot(212), plot(wg0/(pi),20*log10(abs(G0)),'b'); 
hold on
plot(wg1/(pi),20*log10(abs(G1)),'r'); grid
hold off
leg=legend('$$|G_0(e^{j \omega})|$$','$$|G_1(e^{j \omega})|$$');
set(leg,'Interpreter', 'Latex', 'FontSize',14);
xlabel('$$\omega \rightarrow$$','Interpreter', 'Latex', 'FontSize', 14);
ylabel('\it{Amplitude/dB}','Interpreter', 'Latex', 'FontSize', 14);
title('\it{b) Synthese-Filter}','Interpreter', 'Latex', 'FontSize', 16);
f11=figure(11);
plot(wh0/(pi),20*log10(abs(H0)),'b'); 
hold on
plot(wh1/(pi),20*log10(abs(H1)),'r'); grid
hold off
leg=legend('$$|H_0(e^{j \omega})|$$','$$|H_1(e^{j \omega})|$$');
set(leg,'Interpreter', 'Latex', 'FontSize',14);
xlabel('$$\omega/\pi \rightarrow$$','Interpreter', 'Latex', 'FontSize', 14);
ylabel('\it{Amplitude/dB}','Interpreter', 'Latex', 'FontSize', 14);

f12=figure(12);
plot(wg0/(pi),20*log10(abs(G0)),'b'); 
hold on
plot(wg1/(pi),20*log10(abs(G1)),'r'); grid
hold off
leg=legend('$$|G_0(e^{j \omega})|$$','$$|G_1(e^{j \omega})|$$');
set(leg,'Interpreter', 'Latex', 'FontSize',14);
xlabel('$$\omega/\pi \rightarrow$$','Interpreter', 'Latex', 'FontSize', 14);
ylabel('\it{Amplitude/dB}','Interpreter', 'Latex', 'FontSize', 14);

f2=figure(2);
subplot(221), [hz1, hp1, ht1] = zplane(h0); title('a) h0 - Minimalphasig'); grid
set(findobj(hz1, 'Type', 'line'), 'Color', 'k'); 
set(findobj(hp1, 'Type', 'line'), 'Color', 'k');
set(findobj(ht1, 'Type', 'line'), 'Color', 'k');
subplot(222), [hz1, hp1, ht1] = zplane(h1); title('b) h1 - Maximalphasig'); grid
set(findobj(hz1, 'Type', 'line'), 'Color', 'k'); 
set(findobj(hp1, 'Type', 'line'), 'Color', 'k');
set(findobj(ht1, 'Type', 'line'), 'Color', 'k');
subplot(223), [hz1, hp1, ht1] = zplane(g0); title('c) g0 - Maximalphasig'); grid
set(findobj(hz1, 'Type', 'line'), 'Color', 'k'); 
set(findobj(hp1, 'Type', 'line'), 'Color', 'k');
set(findobj(ht1, 'Type', 'line'), 'Color', 'k');
subplot(224), [hz1, hp1, ht1] = zplane(g1); title('d) g1 - Minimalphasig'); grid
set(findobj(hz1, 'Type', 'line'), 'Color', 'k'); 
set(findobj(hp1, 'Type', 'line'), 'Color', 'k');
set(findobj(ht1, 'Type', 'line'), 'Color', 'k');

f3=figure(3);
subplot(221), stem(0:length(h0)-1,h0,'k','LineWidth',1); title('a) h0 - Minimalphasig'); grid
xlabel('n \rightarrow')
subplot(222), stem(0:length(h1)-1,h1,'k','LineWidth',1); title('b) h1 - Maximalphasig'); grid
xlabel('n \rightarrow')
subplot(223), stem(0:length(g0)-1,g0,'k','LineWidth',1); title('c) g0 - Maximalphasig'); grid
xlabel('n \rightarrow')
subplot(224), stem(0:length(g1)-1,g1,'k','LineWidth',1); title('d) g1 - Minimalphasig'); grid
xlabel('n \rightarrow')
% distortion- & alias-function
F_dist=2/2*(G0.*H0+G1.*H1);
j=1i; M=2;
for i=1:length(h0)
   W(i)=exp(j*2*pi/M*(-(i-1)));
end
[H0_ wh0_]=freqz(h0.*W);
[H1_ wh1_]=freqz(h1.*W);
F_alias=abs(1/2*(G0.*H0_+G1.*H1_));
f5=figure(5);
subplot(211), plot(wh0/(pi),20*log10(abs(F_dist)),'k'); grid
xlabel('$$\omega/\pi \rightarrow$$','Interpreter', 'Latex', 'FontSize', 14);
ylabel('$$F_{dist}(e^{j \omega})/dB $$','Interpreter', 'Latex', 'FontSize', 14);
title('\it{a) Verzerrungsfunktion }','Interpreter', 'Latex', 'FontSize', 16);
subplot(212), plot(wh0/(pi),(20*log10(F_alias)),'k'); grid
xlabel('$$\omega/\pi \rightarrow $$','Interpreter', 'Latex', 'FontSize', 14);
ylabel('$$F_{alias}(e^{j \omega})/dB $$','Interpreter', 'Latex', 'FontSize', 14);
title('\it{b) Aliasfunktion}','Interpreter', 'Latex', 'FontSize', 16);

disp('done!')
