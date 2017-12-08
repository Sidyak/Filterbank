%----------------------------------------------------------------------------------------- 
%  S. K. Mitra.  Digital Signal Processing: A Computer-Based Approach. Second Edition, 
%                McGraw-Hill,2001; 
% This function creates an "optimal" lp-prototype filter  for pseudo-QMF
% bank with nbands bands of length N
%------------------------------------------------------------------------------------------ 
function [hopt,H]=opt_filter(N,nbands) 
stopedge=1/nbands; 
delta=0.001; 
passedge=1/(4*nbands); 
toll=0.000001; 
step=0.1*passedge; 
tcost=0; 
mpc=1; 
way=-1; 
pcost=10; 
flag=0; 
s_time=cputime;
flops=0;        % diese Init. fehlt im Buch
s_flops=flops; 
while flag==0 
    hopt=remez(N,[0,passedge,stopedge,1],[1,1,0,0],[5,1]); 
%    length(hopt) 
    figure(10) 
    stem(hopt); grid
    H=fft(hopt,4096); 
    HH=ovlp_ripple(H,nbands); % calc. cost fct. being minimized
    [tcost]=max(abs(HH-ones(max(size(HH)),1))); 
    if tcost > pcost 
        step=step/2; 
        way=-way; 
    end 
    if abs(pcost-tcost)<toll 
        flag=1; 
    end 
    pcost=tcost; 
    passedge=passedge+way*step; 
end 
final_time=cputime-s_time; 
total_flops=flops-s_flops;
 
%save hopt.mat hopt -ascii;
%save pass passedge -ascii;

disp('Optimal Filterdesign done!')
 