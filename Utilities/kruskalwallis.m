for ii = 7:81
    cell(:,ii) = ch_28(:,ii);
end

for ii = 0:83;
    Y = fft(allflash(ii+1,:));
    L = 161
    Fs = 100
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);

    f = Fs*(0:(L/2))/L;
end
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%%%%%%%%%%%%%%%%%%
Fs = 32;
nfft = 2490;
X = fft(cell,nfft);
X = X(1:nfft/2);
mx = abs(X);
f = (0:nfft/2-1)*Fs/nfft;
plot(f,mx);


%%%%%%%%%%%%%%%%%%%%%%%%%%
cf = cfit(cell);
