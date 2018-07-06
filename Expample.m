y = randi(100,1,100);
y = complex(y); % turn to complex so as to produce s with nfft rows instead of just 1+nfft/2 rows.

% parameters
windowLen = 10;
hop = 3;
nfft = length(y);
noverLap = windowLen-hop;

% get STFT magnitude
s = spectrogram(y,windowLen,noverLap,nfft);
S = abs(s);

% STliFT
y_ = STliFt(S, length(y), windowLen, hop);
plot(1:100,y);hold on;plot(1:100,y_)
