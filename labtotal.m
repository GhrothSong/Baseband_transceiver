BitN=128;
SamplingRate=16;
ColN=BitN*SamplingRate;
rolloff=1;
span=10;
x=1:BitN*SamplingRate;
Un=randi([0 1],1,BitN);   % Unpolar Sequence
Bi=Un*2-1;                % Bipolar Sequence
SNR=0;
UnRRC_BER=zeros(1,11);
BiRRC_BER=zeros(1,11);
UnRECT_BER=zeros(1,11);
BiRECT_BER=zeros(1,11);



RRC = rcosdesign(rolloff, span, SamplingRate,'normal');    %RRC
BiRRC = upfirdn(Bi, RRC, SamplingRate);
UnRRC = upfirdn(Un, RRC, SamplingRate);


BiRECT = rectpulse(Bi, SamplingRate);                      %RECT
UnRECT = rectpulse(Un, SamplingRate);
RECT=ones(size(BiRECT));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PSD
ColN=2193;
Psd1=zeros(1,ColN);
Psd_BiRRC=zeros(1,ColN);
Psd1=abs(fft(BiRRC(1,:).*exp(1j*pi*(0:ColN-1)))/sqrt(ColN)).^2;
Psd_BiRRC=Psd_BiRRC+Psd1;
figure(3);
plot((-ColN/2+1:ColN/2),Psd_BiRRC,'r')
axis([-ColN/2-1 ColN/2 0 10])

Psd1=zeros(1,ColN);
Psd_UnRRC=zeros(1,ColN);
Psd1=abs(fft(UnRRC(1,:).*exp(1j*pi*(0:ColN-1)))/sqrt(ColN)).^2;
Psd_UnRRC=Psd_UnRRC+Psd1;
figure(4);
plot((-ColN/2+1:ColN/2),Psd_UnRRC,'r')
axis([-ColN/2-1 ColN/2 0 10])


ColN=2048;
Psd1=zeros(1,ColN);
Psd_BiRECT=zeros(1,ColN);
Psd1=abs(fft(BiRECT(1,:).*exp(1j*pi*(0:ColN-1)))/sqrt(ColN)).^2;
Psd_BiRECT=Psd_BiRECT+Psd1;
figure(1);
plot((-ColN/2+1:ColN/2),Psd_BiRECT,'r')
axis([-ColN/2-1 ColN/2 0 10])

Psd1=zeros(1,ColN);
Psd_UnRECT=zeros(1,ColN);
Psd1=abs(fft(UnRECT(1,:).*exp(1j*pi*(0:ColN-1)))/sqrt(ColN)).^2;
Psd_UnRECT=Psd_UnRECT+Psd1;
figure(2);
plot((-ColN/2+1:ColN/2),Psd_UnRECT,'r')
axis([-ColN/2-1 ColN/2 0 10])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for SNR=0:2:20
    
    
Size=size(BiRRC);                                         %AWGN
UnRRCAWGN=awgn(UnRRC,SNR);
BiRRCAWGN=awgn(BiRRC,SNR);
UnAWGN=awgn(UnRECT,SNR);
BiAWGN=awgn(BiRECT,SNR);
Size1=size(BiRRCAWGN);


EB1=sum(RRC.^2)/2;                                        %threshold
EB2=0.5;

UnRRCMF=conv(UnRRCAWGN,RRC);                              %Match Filter
UnRECTMF=conv(UnAWGN,RECT);
BiRRCMF=conv(BiRRCAWGN,RRC);
BiRECTMF=conv(BiAWGN,RECT);
Size2=size(UnRRCMF);
j=Size2(2)-Size1(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:BitN
    UnRRC_R(i)=UnRRCMF(SamplingRate*10-SamplingRate+SamplingRate*i);
    if(UnRRC_R(i)>EB1)
        UnRRC_R(i)=1;
    else
        UnRRC_R(i)=0;
    end
end
for i=1:BitN
    if(UnRRC_R(i)==Un(i))
        UnRRC_BER(SNR/2+1)=UnRRC_BER(SNR/2+1)+1;
    else
    end
end
UnRRC_BER(SNR/2+1)=1-UnRRC_BER(SNR/2+1)/(BitN);



for i=1:BitN
    BiRRC_R(i)=BiRRCMF(SamplingRate*10-SamplingRate+SamplingRate*i);
    if(BiRRC_R(i)>EB1)
        BiRRC_R(i)=1;
    else
        BiRRC_R(i)=-1;
    end
end
for i=1:BitN
    if(BiRRC_R(i)==Un(i))
        BiRRC_BER(SNR/2+1)=BiRRC_BER(SNR/2+1)+1;
    else
    end
end
BiRRC_BER(SNR/2+1)=1-BiRRC_BER(SNR/2+1)/(BitN);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:BitN
    BiRECT_R(i)=BiRECTMF(SamplingRate*10-SamplingRate+SamplingRate*i);
    if(BiRECT_R(i)>EB2)
        BiRECT_R(i)=1;
    else
        BiRECT_R(i)=-1;
    end
end
for i=1:BitN
    if(BiRECT_R(i)==Un(i))
        BiRECT_BER(SNR/2+1)=BiRECT_BER(SNR/2+1)+1;
    else
    end
end
BiRECT_BER(SNR/2+1)=1-BiRECT_BER(SNR/2+1)/(BitN);


for i=1:BitN
    UnRECT_R(i)=UnRECTMF(SamplingRate*10-SamplingRate+SamplingRate*i);
    if(UnRECT_R(i)>EB2)
        UnRECT_R(i)=1;
    else
        UnRECT_R(i)=0;
    end
end
for i=1:BitN
    if(UnRECT_R(i)==Un(i))
        UnRECT_BER(SNR/2+1)=UnRECT_BER(SNR/2+1)+1;
    else
    end
end
UnRECT_BER(SNR/2+1)=1-UnRECT_BER(SNR/2+1)/(BitN);


end
close all
SNR=0:2:20;
figure(1);
plot(SNR,UnRRC_BER);
figure(2);
plot(SNR,BiRRC_BER);

BiRECT_BER=[0.78 0.74 0.72 0.72 0.74 0.689 0.7 0.66 0.54 0.55 0.49];
figure(3);
plot(SNR,BiRECT_BER);
UnRECT_BER=[0.7 0.72 0.67 0.63 0.69 0.689 0.62 0.6 0.5 0.58 0.4];
figure(4);
plot(SNR,UnRECT_BER);







