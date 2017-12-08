%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K.Radmacher, 29.09.14
%
% Cosinus-modulated pseudo QMF-Filterbank 
% including distortion- & alias-function
% 
% INFO: for optimized prototype/filterbank 
%       uncommend row no. 30
%       files "opt_filter.m" and "ovlp_ripple.m" needed!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
format long
format compact

Fs=48e3;
M=4;
dev_dB=60;
dev=10^(-dev_dB/20);
tw=1.2e3;   % transition width
f=[(Fs/(2*M*2)-tw) (Fs/(2*M*2)+tw)];
[N_FIR,fo,ao,w] = firpmord(f,[1 0],[dev dev],Fs);
if mod(N_FIR,2)==1  	% order must be even / number of coef. must be odd
    N_FIR=N_FIR+1;
end
b=firpm(N_FIR, fo, ao);
%%%%%% optimized prototype filter %%%%%%
%%%%%%  UNCOMMEND FOLLOWING ROW   %%%%%%
%[b, passedge]=opt_filter(N_FIR,M); % using opt_filter.m and ovlp_ripple.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[B,w]=freqz(b,1024);
f1=figure(1);
subplot(311),plot(2*(-512:511)/1024, db(fftshift(abs(fft(b,1024)))),'k'); 
grid
xlabel('$$\omega/\pi \rightarrow$$','Interpreter', 'Latex', 'FontSize', 14);
title('\it{a) Amplitudengang}','Interpreter', 'Latex', 'FontSize', 14);
ylim([-140 5])
ylabel('$$|H_0(e^{j\omega})|/dB$$','Interpreter', 'Latex', 'FontSize', 14);
subplot(312), stem(0:N_FIR,b,'k'); 
grid
ylabel('$$h_0[n]$$','Interpreter', 'Latex', 'FontSize', 14);
title('\it{b) Impulsantwort}','Interpreter', 'Latex', 'FontSize', 14);
xlabel('$$n \rightarrow$$','Interpreter', 'Latex', 'FontSize', 14);
subplot(313),
[hz1, hp1, ht1] = zplane(b); %set(graph_tp,'Color','k', 'LineWidth',1); 
set(findobj(hz1, 'Type', 'line'), 'Color', 'k'); 
set(findobj(hp1, 'Type', 'line'), 'Color', 'k');
set(findobj(ht1, 'Type', 'line'), 'Color', 'k');
title('\it{c) Pol- Nullstellenplan}','Interpreter', 'Latex', 'FontSize', 14);
axis([-2 2 -2 2]); grid
% create polyphase filtercoefficients - without dezimation!
N_samp=length(b);
MM=2*M;
E0=b;
E0(2:MM:N_samp)=0; E0(3:MM:N_samp)=0; E0(4:MM:N_samp)=0;
E0(5:MM:N_samp)=0; E0(6:MM:N_samp)=0; E0(7:MM:N_samp)=0; E0(8:MM:N_samp)=0;
for i=MM+1:2*MM:N_samp
   E0(i)=-E0(i); 
end
E1=b;
E1(1:MM:N_samp)=0; E1(3:MM:N_samp)=0; E1(4:MM:N_samp)=0;
E1(5:MM:N_samp)=0; E1(6:MM:N_samp)=0; E1(7:MM:N_samp)=0; E1(8:MM:N_samp)=0;
for i=MM+2:2*MM:N_samp
   E1(i)=-E1(i); 
end
E2=b;
E2(1:MM:N_samp)=0; E2(2:MM:N_samp)=0; E2(4:MM:N_samp)=0;
E2(5:MM:N_samp)=0; E2(6:MM:N_samp)=0; E2(7:MM:N_samp)=0; E2(8:MM:N_samp)=0;
for i=MM+3:2*MM:N_samp
   E2(i)=-E2(i); 
end
E3=b;
E3(1:MM:N_samp)=0; E3(2:MM:N_samp)=0; E3(3:MM:N_samp)=0;
E3(5:MM:N_samp)=0; E3(6:MM:N_samp)=0; E3(7:MM:N_samp)=0; E3(8:MM:N_samp)=0;
for i=MM+4:2*MM:N_samp
   E3(i)=-E3(i); 
end
E4=b;
E4(1:MM:N_samp)=0; E4(2:MM:N_samp)=0; E4(3:MM:N_samp)=0;
E4(4:MM:N_samp)=0; E4(6:MM:N_samp)=0; E4(7:MM:N_samp)=0; E4(8:MM:N_samp)=0;
for i=MM+5:2*MM:N_samp
   E4(i)=-E4(i); 
end
E5=b;
E5(1:MM:N_samp)=0; E5(2:MM:N_samp)=0; E5(3:MM:N_samp)=0;
E5(4:MM:N_samp)=0; E5(5:MM:N_samp)=0; E5(7:MM:N_samp)=0; E5(8:MM:N_samp)=0;
for i=MM+6:2*MM:N_samp
   E5(i)=-E5(i); 
end
E6=b;
E6(1:MM:N_samp)=0; E6(2:MM:N_samp)=0; E6(3:MM:N_samp)=0;
E6(4:MM:N_samp)=0; E6(5:MM:N_samp)=0; E6(6:MM:N_samp)=0; E6(8:MM:N_samp)=0; 
for i=MM+7:2*MM:N_samp
   E6(i)=-E6(i); 
end
E7=b;
E7(1:MM:N_samp)=0; E7(2:MM:N_samp)=0; E7(3:MM:N_samp)=0;
E7(4:MM:N_samp)=0; E7(5:MM:N_samp)=0; E7(6:MM:N_samp)=0; E7(7:MM:N_samp)=0;
for i=MM+8:2*MM:N_samp
   E7(i)=-E7(i); 
end
% create cos.-mod.-matrix
T_a=zeros(M,2*M);
T_s=zeros(M,2*M);
for k=1:M   % number of subbands
    for n=1:2*M   % number of coef.
        teta=(-1)^(k-1)*pi/4;
        coef=(N_FIR/(1));
        T_a(k,n)=2*cos(pi/M*((k-1)+0.5)*((n-1)-(coef)/2)+teta);  % analyze
        T_s(k,n)=2*cos(pi/M*((k-1)+0.5)*((n-1)-(coef)/2)-teta);  % synthesis
    end
end
% create polyphase vector
e=[E0' E1' E2' E3' E4' E5' E6' E7']; 
% do the cos.-modulation
f=zeros(M,length(e(:,1)));
p=zeros(M,length(e(:,1)));
for k=1:M
    p(k,:)=(T_a(k,1)*e(:,1)')+(T_a(k,2)*e(:,2)')+(T_a(k,3)*e(:,3)')+(T_a(k,4)*e(:,4)')+(T_a(k,5)*e(:,5)')+(T_a(k,6)*e(:,6)')+(T_a(k,7)*e(:,7)')+(T_a(k,8)*e(:,8)');    
    f(k,:)=(T_s(k,1)*e(:,1)')+(T_s(k,2)*e(:,2)')+(T_s(k,3)*e(:,3)')+(T_s(k,4)*e(:,4)')+(T_s(k,5)*e(:,5)')+(T_s(k,6)*e(:,6)')+(T_s(k,7)*e(:,7)')+(T_s(k,8)*e(:,8)');    
end
f2=figure(2);
P0=fft(p(1,:),1024);
plot((-512:511)/512, db(fftshift(abs(P0))),'k'); 
hold on
title('$$M=4 Cos.-Mod.-Filterbank (pseudo quadrature mirror FB)$$','Interpreter', 'Latex', 'FontSize', 14);
P1=fft(p(2,:),1024);
plot((-512:511)/512, db(fftshift(abs(P1))),'b'); 
P2=fft(p(3,:),1024);
plot((-512:511)/512, db(fftshift(abs(P2))),'r'); 
P3=fft(p(4,:),1024);
plot((-512:511)/512, db(fftshift(abs(P3))),'g'); 
grid
hold off
legend('$$H0','H1','H2','H3$$','Interpreter', 'Latex', 'FontSize', 14);
xlabel('$$f/Fs \rightarrow$$','Interpreter', 'Latex', 'FontSize', 14);
ylabel('$$magnitude / dB$$','Interpreter', 'Latex', 'FontSize', 14);
% distortion function
F_dist=((P0)+(P1)+(P2)+(P3));
% create alias components
for i=1:length(p(1,:))
   W_(i)=exp(2i*pi/M*(i-1)*1); 
   W__(i)=exp(2i*pi/M*(i-1)*2); 
   W___(i)=exp(2i*pi/M*(i-1)*3); 
end
H0_     =(((fft(p(1,:).*W_,1024))));
H0__    =(((fft(p(1,:).*W__,1024))));
H0___   =(((fft(p(1,:).*W___,1024))));
H1_     =(((fft(p(2,:).*W_,1024))));
H1__    =(((fft(p(2,:).*W__,1024))));
H1___   =(((fft(p(2,:).*W___,1024))));
H2_     =(((fft(p(3,:).*W_,1024))));
H2__    =(((fft(p(3,:).*W__,1024))));
H2___   =(((fft(p(3,:).*W___,1024))));
H3_     =(((fft(p(4,:).*W_,1024))));
H3__    =(((fft(p(4,:).*W__,1024))));
H3___   =(((fft(p(4,:).*W___,1024))));
% synthesis filter coincides with analysis filter
G0=fft(f(1,:),1024); G1=fft(f(2,:),1024); G2=fft(f(3,:),1024); G3=fft(f(4,:),1024);
% alias function
F_alias=sqrt(1/M*(...
    abs(G0.*H0_  +G1.*H1_  +G2.*H2_  +G3.*H3_  )+...
    abs(G0.*H0__ +G1.*H1__ +G2.*H2__ +G3.*H3__ )+...
    abs(G0.*H0___+G1.*H1___+G2.*H2___+G3.*H3___)));
f3=figure(3);
subplot(312), plot((-512:511)/512,db(fftshift(abs(F_dist))),'k'); 
grid
xlabel('$$\omega/\pi \rightarrow$$','Interpreter', 'Latex', 'FontSize', 14);
ylabel('$$F_{dist}(e^{j \omega})/dB$$','Interpreter', 'Latex', 'FontSize', 14);
title('\it{b) Verzerrungsfunktion }','Interpreter', 'Latex', 'FontSize', 14);
subplot(313), plot((-512:511)/512,db(fftshift(F_alias)),'k'); 
grid
xlabel('$$\omega/\pi \rightarrow$$','Interpreter', 'Latex', 'FontSize', 14);
ylabel('$$F_{alias}(e^{j \omega})/dB$$','Interpreter', 'Latex', 'FontSize', 14);
title('\it{c) Aliasfunktion }','Interpreter', 'Latex', 'FontSize', 14);
subplot(311), p0=plot((-512:511)/512, db(fftshift(abs(P0))),'k'); 
hold on
title('\it{a) Teilb\"ander }','Interpreter', 'Latex', 'FontSize', 14);
p1=plot((-512:511)/512, db(fftshift(abs(P1))),'b'); 
p2=plot((-512:511)/512, db(fftshift(abs(P2))),'r'); 
p3=plot((-512:511)/512, db(fftshift(abs(P3))),'g'); 
grid
hold off
ylim([-140 5])
xlabel('$$\omega/\pi \rightarrow$$','Interpreter', 'Latex', 'FontSize', 14);
ylabel('$$|H_i(e^{j\omega})|/dB, i=0...3 $$','Interpreter', 'Latex', 'FontSize', 14);
leg=legend('Location', [0.95 0.830 0.005 0.05],'|H_0(e^{j \omega})|','|H_1(e^{j \omega})|','|H_2(e^{j \omega})|','|H_3(e^{j \omega})|','Location','Best');
set(leg,'FontAngle','italic','TextColor',[.3,.2,.1]);

disp('done')