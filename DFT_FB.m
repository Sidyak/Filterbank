%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K.Radmacher, 30.09.14
%
% M-th-band polyphase DFT Filterbank
% M=4
% including distortion- & alias-function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
format long
format compact

M=4;    % M-th bandfilter (M=2 -> Halfbandfilter)    
Fs=48e3;
f_mth_band=(Fs/2)*(1/(1*M));    % stopbandfreq. for halbbandfilter
x=1.2e3;
f=[f_mth_band-x f_mth_band+x];  % frequencys for pass- and stopband, must be symetic
a=[1 0];                        % for lowpass
rippledB=60;
dev=[10^-(rippledB/20) 10^-(rippledB/20)];	% deviation for pass- and stopband must be equal
[n,fo,ao,w] = firpmord(f,a,dev,Fs);
if mod(n,2)==1                  % order must be even / number of koeff. must be odd
    n=n+1;
end
b_lp=firpm(n, fo, ao);
B_LP=fftshift(db(abs(fft(b_lp,4096))));
freq=-2048:2047;
f1=figure(1);
subplot(311),
plot(freq/(2048), B_LP,'k');   
grid
xlabel('$$\omega/\pi \rightarrow$$','Interpreter', 'Latex', 'FontSize', 14);
title('\it{a) Amplitudengang}','Interpreter', 'Latex', 'FontSize', 14);
ylabel('$$|H_0(e^{j\omega})|/dB$$','Interpreter', 'Latex', 'FontSize', 14);
ylim([-80 5])
subplot(312), stem(0:n,b_lp,'k');   
grid
title('\it{b) Impulsantwort}','Interpreter', 'Latex', 'FontSize', 14);
ylabel('$$h_0[n]$$','Interpreter', 'Latex', 'FontSize', 14);
xlabel('$$n \rightarrow$$','Interpreter', 'Latex', 'FontSize', 14);
subplot(313),
[hz1, hp1, ht1] = zplane(b_lp);  
set(findobj(hz1, 'Type', 'line'), 'Color', 'k'); 
set(findobj(hp1, 'Type', 'line'), 'Color', 'k');
set(findobj(ht1, 'Type', 'line'), 'Color', 'k');
title('\it{c) Pol- Nullstellenplan}','Interpreter', 'Latex', 'FontSize', 14);
axis([-2 2 -2 2]); grid
% create polyphases
N_samp=length(b_lp);
E0=b_lp;
E0(2:M:N_samp)=0; E0(3:M:N_samp)=0; E0(4:M:N_samp)=0;
E1=b_lp;
E1(1:M:N_samp)=0; E1(3:M:N_samp)=0; E1(4:M:N_samp)=0;
E2=b_lp;
E2(1:M:N_samp)=0; E2(2:M:N_samp)=0; E2(4:M:N_samp)=0;
E3=b_lp;
E3(1:M:N_samp)=0; E3(2:M:N_samp)=0; E3(3:M:N_samp)=0;

j=1i;   % imaginary number
% create DFT Matrix
W=zeros(M,M)*j;
for i=1:M
    for k=1:M
        W(i,k)=exp(j*2*pi/M*(i-1)*(k-1));
    end
end
e=[E0' E1' E2' E3'];    % Polyphase vector
% do the DFT
h=zeros(M,length(e(:,1)));
for k=1:M
    h(k,:)=(W(k,1)*e(:,1)')+(W(k,2)*e(:,2)')+(W(k,3)*e(:,3)')+(W(k,4)*e(:,4)');    
end
% transform 
l=1024;
H0=fft(h(1,:),l);
H1=fft(h(2,:),l);
H2=fft(h(3,:),l);
H3=fft(h(4,:),l);
% Distortion function
F_dist=((H0)+(H1)+(H2)+(H3));
% Alias components
for i=1:length(h(1,:))
   W_(i)=exp(2i*pi/M*(i-1)*1); 
   W__(i)=exp(2i*pi/M*(i-1)*2); 
   W___(i)=exp(2i*pi/M*(i-1)*3); 
end
H0_     =(((fft(h(1,:).*W_,1024))));
H0__    =(((fft(h(1,:).*W__,1024))));
H0___   =(((fft(h(1,:).*W___,1024))));
H1_     =(((fft(h(2,:).*W_,1024))));
H1__    =(((fft(h(2,:).*W__,1024))));
H1___   =(((fft(h(2,:).*W___,1024))));
H2_     =(((fft(h(3,:).*W_,1024))));
H2__    =(((fft(h(3,:).*W__,1024))));
H2___   =(((fft(h(3,:).*W___,1024))));
H3_     =(((fft(h(4,:).*W_,1024))));
H3__    =(((fft(h(4,:).*W__,1024))));
H3___   =(((fft(h(4,:).*W___,1024))));
% Synthesefilter coincides with Analysefilter
G0=fft(h(1,:),1024); G1=fft(h(2,:),1024); G2=fft(h(3,:),1024); G3=fft(h(4,:),1024);
% alias function
F_alias=sqrt(1/M*(...
    abs(G0.*H0_  +G1.*H1_  +G2.*H2_  +G3.*H3_  )+...
    abs(G0.*H0__ +G1.*H1__ +G2.*H2__ +G3.*H3__ )+...
    abs(G0.*H0___+G1.*H1___+G2.*H2___+G3.*H3___)));
f2=figure(2);
subplot(312), plot((-512:511)/(512),db(fftshift(abs(F_dist))),'k');   
grid
xlabel('$$\omega/\pi \rightarrow$$','Interpreter', 'Latex', 'FontSize', 14);
title('\it{b) Verzerrungsfunktion}','Interpreter', 'Latex', 'FontSize', 14);
ylabel('$$F_{dist}(e^{j \omega})/dB$$','Interpreter', 'Latex', 'FontSize', 14);
subplot(313), plot((-512:511)/(512),db(fftshift(F_alias)),'k');   
grid
xlabel('$$\omega/\pi \rightarrow$$','Interpreter', 'Latex', 'FontSize', 14);
title('\it{c) Aliasfunktion}','Interpreter', 'Latex', 'FontSize', 14);
ylabel('$$F_{alias}(e^{j \omega})/dB$$','Interpreter', 'Latex', 'FontSize', 14);
subplot(311), p0=plot((-512:511)/(512), fftshift(db(abs(H0))),'k');   
hold on
title('\it{a) Teilb\"ander}','Interpreter', 'Latex', 'FontSize', 14);
p1=plot((-512:511)/(512), db(fftshift(abs(H1))),'b');   
p2=plot((-512:511)/(512), db(fftshift(abs(H2))),'r');   
p3=plot((-512:511)/(512), db(fftshift(abs(H3))),'g');   
grid
hold off
ylim([-80 5])
xlabel('$$\omega/\pi \rightarrow$$','Interpreter', 'Latex', 'FontSize', 14);
ylabel('$$|H_i(e^{j\omega})|/dB, i=0...3$$','Interpreter', 'Latex', 'FontSize', 14);
leg=legend('Location', [0.95 0.8296 0.005 0.05],'|H_0(e^{j \omega})|','|H_1(e^{j \omega})|','|H_2(e^{j \omega})|','|H_3(e^{j \omega})|','Location','Best');
set(leg,'FontAngle','italic','TextColor',[.3,.2,.1]);

disp('done')