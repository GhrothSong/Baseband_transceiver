
clf;

RowN=1000;
BitN=64;
SamplingRate=16;
ColN=BitN*SamplingRate;


x=randi([0 1],RowN,BitN); % Unpolar Sequence
x=x*2-1;      % Bipolar Sequence

% oversampling for the first bit of each row) 
y=x(:,1);
for n=1:log2(SamplingRate)
y=cat(2,y,y);
end                            

% oversampling for all bit sequences
for m=2:BitN
x1=x(:,m);
for n=1:log2(SamplingRate)
x1=cat(2,x1,x1);
end 
y=cat(2,y,x1);
end
% 

% Normalise the magnitude of each sample
y=y*sqrt(1/SamplingRate);  % normalise the bit enery to be one regardless of sampling rate 


% PSD 
Psd=zeros(1,ColN);
for m=1:RowN
    Psd1=abs(fft(y(m,:).*exp(1j*pi*(0:ColN-1)))/sqrt(ColN)).^2;  % shift zero freqency to the center of the spectrum otherwise you could use the commented -out line below )
% Psd1=abs(fft(y(m,:))/sqrt(ColN)).^2;
Psd=Psd+Psd1;
end
Psd=Psd/RowN;

% plot
plot((-ColN/2+1:ColN/2),Psd,'r')
axis([-ColN/2-1 ColN/2 0 1])
