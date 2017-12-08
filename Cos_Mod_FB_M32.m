%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K.Radmacher, 05.11.14
%
% M=32 Cosinuns-modulated Filterbank using optimized prototype
% including distortion- & alias-function
% INFO: determination of prototype filter coefs takes a while to be done!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
format long 
format compact

M=32;
N_FIR=512-1;    % Filter order
%%%%%% optimized prototype filter %%%%%%
% need  a few seconds to determine coefs
[b, passedge]=opt_filter(N_FIR,M); % using opt_filter.m and ovlp_ripple.m
[B w]=freqz(b,512);
f1=figure(1);
subplot(211), stem(0:511,b, 'k');  
grid
xlim([0 512])
xlabel('$$n \rightarrow$$','Interpreter', 'Latex', 'FontSize', 14);
ylabel('$$h_0[n]$$','Interpreter', 'Latex', 'FontSize', 14);
title('\it{a) Impulsantwort}','Interpreter', 'Latex', 'FontSize', 14);
subplot(212), plot((-512:511)/512, db(fftshift(abs(fft(b,1024)))), 'k');  
grid
ylim([-140 5])
xlabel('$$\omega/\pi \rightarrow$$','Interpreter', 'Latex', 'FontSize', 14);
title('\it{b) Amplitudengang}','Interpreter', 'Latex', 'FontSize', 14);
ylabel('$$|H_0(e^{j\omega})|/dB$$','Interpreter', 'Latex', 'FontSize', 14);
% create polyphase filter
N_samp=length(b);
MM=2*M;
E=zeros(64,512);
for k_poly=1:64
    E(k_poly,:)=b;
    for ind=1:64
        if ind==k_poly
            % keep coeff
        else
            E(k_poly,(ind:MM:N_samp))=0;
        end
    end
end
for k_poly=1:64
    for i=(64+k_poly):2*64:512
        E(k_poly,i)=-E(k_poly,i);  % invert coeff
    end
end
% create cos.-mod.-matrix
T_a=zeros(M,2*M);
T_s=zeros(M,2*M);
for k=1:M   % number of subbands
    for n=1:2*M   % number of coef.
        teta=(-1)^(k-1)*pi/4;
        coef=(N_FIR/(1));
        T_a(k,n)=2*cos(pi/M*((k-1)+0.5)*((n-1)-(coef)/2)+teta);  % Analyse
        T_s(k,n)=2*cos(pi/M*((k-1)+0.5)*((n-1)-(coef)/2)-teta);  % Synthese
    end
end
% create polyphase vector
e=[];       
e=E(1,:);
for k_poly=2:64
    e=[e; E(k_poly,:)];
end
% do the cos.-modulation
h=zeros(M,length(e(1,:)));
g=zeros(M,length(e(1,:)));
H=zeros(M,512);
f2=figure(2);
for k=1:M
    buf_a=zeros(length(e(1,:)),1);
    buf_s=zeros(length(e(1,:)),1);
    for ind=1:64
        buf_a(:)=buf_a(:)+(T_a(k,ind)*e(ind,:)');
        buf_s(:)=buf_s(:)+(T_s(k,ind)*e(ind,:)');
    end
    h(k,:)= buf_a(:);
    g(k,:)= buf_s(:);
    H(k,:)= fft(h(k,:),512);
    subplot(311),plot((-256:255)/256, db(fftshift(abs(H(k,:)))),'k');
    hold on
 end
grid;
hold off;
title('\it{a) Teilb\"ander}','Interpreter', 'Latex', 'FontSize', 14);
ylabel('$$|H_{i}(e^{j \omega})|/dB, i=0...31$$','Interpreter', 'Latex', 'FontSize', 14);
xlabel('$$\omega/\pi \rightarrow$$','Interpreter', 'Latex', 'FontSize', 14);
xlim([0 1])
ylim([-150 5])
% distortion function
F_dist=zeros(512,1);
for n=1:512
    F_dist(n)=F_dist(n)+sum(H(:,n));
end
subplot(312), plot((-256:255)/256,db(fftshift(abs(F_dist))),'k');%,'LineWidth',.75); 
grid
xlabel('$$\omega/\pi \rightarrow$$','Interpreter', 'Latex', 'FontSize', 14);
title('\it{b) Verzerrungsfunktion}','Interpreter', 'Latex', 'FontSize', 14);
ylabel('$$F_{dist}(e^{j \omega})/dB$$','Interpreter', 'Latex', 'FontSize', 14);
xlim([0 1])
% create alias components
W=zeros(64,512);      % twiddle operator
H_=zeros(32,31,512);  % analysefilter-alias components
G=zeros(32,512);      % synthesis filter
for alias=1:31        % M alias components
    for i=1:length(h(1,:))  
        W(alias,i) =exp(-2i*pi/M*(i-1)*(alias)); 
    end
    for filter=1:M
        H_(filter,alias,:)=(((fft(h(filter,:).*W(alias,:),512))));
        G(filter,:)=fft(g(filter,:),512);
    end
    
end
% alias function
F_alias=zeros(31,512);
buf=zeros(31,512);
for alias=1:31    
    for ind=1:512
        for filter=1:32
            buf(alias,ind)=buf(alias,ind)+(G(filter,ind).*H_(filter,alias,ind));
        end
    end
    F_alias(alias,:)=F_alias(alias,:)+abs(buf(alias,:));
end
F_alias_final=zeros(512,1);
for ind=1:512
    F_alias_final(ind)=sqrt(1/M*sum(F_alias(:,ind)));
end
subplot(313), plot((-256:255)/256,db(fftshift(F_alias_final)),'k');
grid
xlabel('$$\omega/\pi \rightarrow$$','Interpreter', 'Latex', 'FontSize', 14);
title('\it{c) Aliasfunktion}','Interpreter', 'Latex', 'FontSize', 14);
ylabel('$$F_{alias}(e^{j \omega})/dB$$','Interpreter', 'Latex', 'FontSize', 14);
xlim([0 1])

disp('done')