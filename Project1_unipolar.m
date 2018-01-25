%Set system parameters e.g. date rate, bit rate, signal (Bipolar vs. Unipolar binary sequence ),
%pulse shaping (rectangular vs. square root raised cosine), Number of Samples/bit,
length=128; %bit length
sampling_rate=16; %sampling rate is 16
T_b=1; %1 bit/sec
trails=100;
SNR=[0,2,4,6,8,10,12,14,16,18];
rrc_downsampling=zeros(1,length);
rrc_BER=zeros(1,length);
rrc_recover=zeros(1,length);
rect_downsampling=zeros(1,length);
rect_BER=zeros(1,length);
rect_recover=zeros(1,length);
BER1=zeros(1,10);
BER2=zeros(1,10);


%Main code of baseband transceiver
for run=1:size(SNR,2) 
 for i=1:trails
    unipolar=randi([0 1],1,length );%set unipolar signal
    rect_unipolar=rectpulse(unipolar,sampling_rate);%rectangular pulse shaping
    rrc_unipolar=rcosdesign(1,10,sampling_rate,'sqrt');%root raised cosine filter
    sig_rrc=upfirdn(unipolar,rrc_unipolar,sampling_rate);%shaped and oversampled RRC wave
    sig_rect=rect_unipolar;               %shaped and oversampled rect wave 
  
  Eb=sum(rrc_unipolar.^2)/2;  %Set the thershold
  sig=sqrt(Eb/(10^(SNR(run)/10)*2));

  
  
  %Add white Gaussian noise to signal
  noise_rrc=normrnd(0,sig,size(sig_rrc));
  noise_rect=awgn(sig_rect,18);%Set the SNR as 18
  received_rrc=sig_rrc+noise_rrc;
  received_rect=noise_rect;
  
 
  %Root raised cosine waveform after match filter
  matchfilter_rrc=conv(received_rrc,rrc_unipolar);
  %Rectagular waveform after match filter
  rect_filter=ones(1,16);
  matchfilter_rect=conv(received_rect,rect_filter);
 
%RRC part match filter
    for j=1:length
    rrc_downsampling(j)=matchfilter_rrc(sampling_rate*10-sampling_rate+sampling_rate*j);
      if(rrc_downsampling(j)>Eb)
         rrc_recover(j)=1;
      else
     rrc_recover(j)=0;
      end
    end
     for j=1:length
         correct_rrc=0;
      if(rrc_recover(j)==unipolar(j))
       correct_rrc=correct_rrc+1;
      else
      end
     end
     rrc_BER(i)=1-correct_rrc/length;

 
%RECT part match filter
    for j=1:length
     rect_downsampling(j)=matchfilter_rect(sampling_rate*j);
      if(rect_downsampling(j)>0.5)
         rect_recover(j)=1;
      else
         rect_recover(j)=0;
      end
     end
     for j=1:length
         correct_rect=0;
      if(rect_recover(j)==unipolar(j))
       correct_rect=correct_rect+1;
      else
      end
      rect_BER(i)=1-correct_rrc/length;
     end
 
 
%PSD for rrc
ColN=2193;
Psd1=zeros(1,ColN);
Psd_RRC=zeros(1,ColN);
Psd1=abs(fft(sig_rrc(1,:).*exp(1j*pi*(0:ColN-1)))/sqrt(ColN)).^2;
Psd_RRC=Psd_RRC+Psd1;

%PSD for RECT
ColN2=2048;
Psd2=zeros(1,ColN2);
Psd_RECT=zeros(1,ColN2);
Psd2=abs(fft(sig_rect(1,:).*exp(1j*pi*(0:ColN2-1)))/sqrt(ColN2)).^2;
Psd_RECT=Psd_RECT+Psd2;

 end  
 
 BER1(run)=sum(rrc_BER)/trails;
 BER2(run)=sum(rect_BER)/trails;  
end

%-----------------------------Plot Section--------------------------------

%PSD for RRC 
figure;
plot((-ColN/2+1:ColN/2),Psd_RRC,'r')
axis([-ColN/2-1 ColN/2 0 10])
xlabel('frequency');
ylabel('PSD');
title('Unipolar RRC PSD');

%PSD for RECT
figure;
plot((-ColN2/2+1:ColN2/2),Psd_RECT,'r')
axis([-ColN2/2-1 ColN2/2 0 10])
xlabel('frequency');
ylabel('PSD');
title('Unipolar RECT PSD');

%RRC plots
figure;
subplot(4,1,1);
t1=linspace(0,T_b*length,size(unipolar,2));
plot(t1,unipolar,'r+');
hold on;
plot(t1,rrc_recover,'bo');
xlabel('time');
ylabel('bits');
legend('original signal','recovered signal');
title('original signal and recovered signal(RRC)');
grid on;
hold off;
subplot(4,1,2);
t2=linspace(0,T_b*length,size(sig_rrc,2));
plot(t2,sig_rrc);
xlabel('time');
ylabel('bits');
title('Pulse Shaped signal(RRC)');
grid on;
subplot(4,1,3);
t3=linspace(0,T_b*length,size(received_rrc,2));
plot(t3,received_rrc);
xlabel('time');
ylabel('bits');
title('received signal with noise(RRC)');
grid on;
subplot(4,1,4);
t4=linspace(0,T_b*length,size(matchfilter_rrc,2));
plot(t4,matchfilter_rrc);
xlabel('time');
ylabel('bits');
title('received signal with noise after match filter(RRC)');
grid on;

%RECT plots
figure;
subplot(4,1,1);
t1=linspace(0,T_b*length,size(unipolar,2));
plot(t1,unipolar,'r+');
hold on;
plot(t1,rect_recover,'bo');
xlabel('time');
ylabel('bits');
legend('original signal','recovered signal');
title('original signal and recovered signal(RECT)');
grid on;
hold off;
subplot(4,1,2);
t2=linspace(0,T_b*length,size(sig_rect,2));
plot(t2,sig_rect);
xlabel('time');
ylabel('bits');
title('Pulse Shaped signal(RECT)');
grid on;
subplot(4,1,3);
t3=linspace(0,T_b*length,size(received_rect,2));
plot(t3,received_rect);
xlabel('time');
ylabel('bits');
title('received signal with noise(RECT)');
grid on;
subplot(4,1,4);
t4=linspace(0,T_b*length,size(matchfilter_rect,2));
plot(t4,matchfilter_rect);
xlabel('time');
ylabel('bits');
title('received signal with noise after match filter(RECT)');
grid on;



%BER vs. SNR(RRC Unipolar)
figure;
semilogy(SNR,BER1,'-ro');
axis([0,20,10^(-20),10]);
hold on;
stem((0:2:18),qfunc(sqrt((10.^(SNR/10)))));
title('BER vs. Eb/No(RRC Unipolar)');
xlabel('Eb/N0(dB)');
ylabel('BER');

%BER vs. SNR(RECT Unipolar)
figure;
lineseries=semilogy(SNR,BER2,'-ro');
axis([0,20,10^(-20),10]);
hold on;
stem((0:2:18),qfunc(sqrt((10.^(SNR/10)))));
title('BER vs. Eb/No(RECT Unipolar)');
xlabel('Eb/N0(dB)');
ylabel('BER');



