function [ b ] = MinimalPhaseFIR( A, delta, reziprok)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K.Radmacher, 03.09.14
%
% Funktion zur Berechnung eines minimalphasigen-, gültigen Halbbandfilters
% Da die hälfte der Koeff. wegfällt und keine Symmetrie mehr vorhanden ist,
% ist das resultierende Filter auch NICHT mehr LINEARPHASIG!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% zunächst gültiges Halbbandfilter machen (nicht negativen Freq.-Gang)
N=length(A);
x=0.1;
A((N+1)/2)=A((N+1)/2)+delta*(1+x); % Dämpfungsfaktor 5 % größer gewählt -> Sperrdämpf. sinkt!
% So skalieren, dass größter Koeff. =0.5, also max(T)=0.5
% bzw. H(e^jw)=0.5, für f=Fs/4
T=A*0.5/(0.5+delta*(1+x));
% dieses Filter wird immer geplottet
[h w]=freqz(T);
f345=figure(234); 
subplot(311), stem(T,'k','LineWidth',2); title('Impulse-response Prototyp, valid halfband Filter'); grid
subplot(312), hz1=zplane(T); title('zplane Prototyp'); grid
subplot(313), plot(w/(2*pi), abs(h),'k','LineWidth',2); 
title_string=sprintf('Freq-response Prototyp %d-ter Ordn. (nicht negativ \forall f)',length(T)-1);
title(title_string); grid
% change color and line für zplane
set(findobj(hz1, 'Type', 'line'), 'Color', 'k','LineWidth',2); 

% dann minimalphasig (alle NS im EK)
R(1,:)=roots(T); % Nullstellen berechnen
fprintf('\nBetrags-Nullstellen (Wenn =1 -> auf dem EK, KRITISCH für spektrale Faktorisierung!\n')
abs(R)
if(reziprok==1) % andere Methode -> Nullstellen außerhalb d. EK reziproken
    z=R(:)';
    for i=1:length(R)
       if abs(z(i))>=1 
           z(i)=1/z(i);
       end
    end
    k=sqrt(T((N+1)/2-1))/(prod(R));
    % Koeffizienten des minimalphasigen Prototypenfilters ermitteln
    [b,a] = zp2tf(z',zeros(1,length(z)),k);
    b=b/sum(b); % Frequenzgang auf 1 normieren
else
    % Nur Nullstellen im EK herausnehmen (nicht mehr lin. phasig)
    j=1;
    for i=1:length(R)
       if abs(R(i))<=1 
           z(j,1)=R(i);
           j=j+1;
       end
    end
    prodR=prod(R);
    % Skalierungs- Gainfaktor ermitteln
    gain=sqrt(T((N+1)/2-1))/(prodR);
    [b,a] = zp2tf(z,zeros(1,length(z)),real(gain));
    % Da Skalierungsfaktor k die H(z) nicht auf 1 Skaliert, 
    % nochmal anhand der Koeffizienten neu skalieren
    b=b/sum(b); % Frequenzgang auf 1 normieren
end

end

