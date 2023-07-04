clc
clear all
close all
figure(1)
data = table2array(readtable('Signal_R_F_1.txt')); %conerting table to array
y1=data(1:4500) %getting only upto 4500 data points
fs=4500%sampling frequency
t=0:4500-1
T1= y1; 
plot (t,T1);
xlabel('time')
ylabel('amplitude')
title('EEG signal plot (in time)')%eeg signal plot in tiime domain
%fft of the signal
figure(2)
xfft=fft(T1)
xfft=xfft(1:length(T1)/2+1)/length(T1)
xfft(2:end-1) = 2*xfft(2:end-1);
plot(abs(xfft));
xlabel("frequency")
ylabel("Amplitude")
title("frequency spectrum of eeg signal")
%adding noise to signal
figure(3)
w=200*cos(2*pi*12.5*t)
n=y1+w
subplot(3,1,1)
plot(t,n)
title("signal added with noise")
xlabel("time")
ylabel("amplitude")
N=7
f = 0.05
[b1, a1] = butter(N,0.32,'low');
[h, w] = freqz(b1, a1); 
subplot(3,1,2)
plot(abs(h),'b'); 
title("frequency response of low pass filter using butterworth");%here we use low pass filter to remove noise of freq 12.4hz 
ylabel('|H(w)|');
xlabel('Frequency');
grid on
y0 = filter(b1,a1,n);
subplot(3,1,3)
plot(y0)
title("Signal after filter noise")
xlabel("time")
ylabel("amplitude")
%DECOMPOSITION OF SIGNAL INTO 5 LEVELS
%DISCRETE WAVELET TRANSFORM, 5 level wavelet db2
figure(4)
waveletFunction = 'db8';
                [C,L] = wavedec(y0,8,waveletFunction);
       
                cD1 = detcoef(C,L,1);
                cD2 = detcoef(C,L,2);
                cD3 = detcoef(C,L,3);
                cD4 = detcoef(C,L,4);
                cD5 = detcoef(C,L,5); %GAMA
                cD6 = detcoef(C,L,6); %BETA
                cD7 = detcoef(C,L,7); %ALPHA
                cD8 = detcoef(C,L,8); %THETA
                cA8 = appcoef(C,L,waveletFunction,8); %DELTA
                D1 = wrcoef('d',C,L,waveletFunction,1);
                D2 = wrcoef('d',C,L,waveletFunction,2);
                D3 = wrcoef('d',C,L,waveletFunction,3);
                D4 = wrcoef('d',C,L,waveletFunction,4);
                D5 = wrcoef('d',C,L,waveletFunction,5); %GAMMA
                D6 = wrcoef('d',C,L,waveletFunction,6); %BETA
                D7 = wrcoef('d',C,L,waveletFunction,7); %ALPHA
                D8 = wrcoef('d',C,L,waveletFunction,8); %THETA
                A8 = wrcoef('a',C,L,waveletFunction,8); %DELTA
                
                Gamma = D5;
                 subplot(5,1,1); plot(1:1:length(Gamma),Gamma);
                 title('Decomposed signal GAMMA');
                 xlabel("time")
                 ylabel("amplitude")
               
                Beta = D6;
                subplot(5,1,2); plot(1:1:length(Beta), Beta); title('BETA');
                xlabel("time")
                ylabel("amplitude")
                
                Alpha = D7;
                subplot(5,1,3); plot(1:1:length(Alpha),Alpha); title('ALPHA'); 
                xlabel("time")
                ylabel("amplitude")
                
                Theta = D8;
                subplot(5,1,4); plot(1:1:length(Theta),Theta);title('THETA');
                D8 = detrend(D8,0);
                xlabel("time")
ylabel("amplitude")
                 Delta = A8;
                %figure, plot(0:1/fs:1,Delta);
                subplot(5,1,5);plot(1:1:length(Delta),Delta);title('DELTA');
                xlabel("time")
ylabel("amplitude")
D5 = detrend(D5,0);
xdft = fft(D5);
freq = 0:N/length(D5):N/2;
xdft = xdft(1:length(D5)/2+1);
figure(5);subplot(5,1,1);
plot(freq,abs(xdft));
title('GAMMA-FREQUENCY SPECTRUM');
xlabel("time")
ylabel("freq")
[~,I] = max(abs(xdft));
fprintf('Gamma:Maximum occurs at %3.2f Hz.\n',freq(I));
D6 = detrend(D6,0);
xdft2 = fft(D6);
freq2 = 0:N/length(D6):N/2;
xdft2 = xdft2(1:length(D6)/2+1);
% figure;
subplot(5,1,2);
plot(freq2,abs(xdft2));
title('BETA');
xlabel("time")
ylabel("freq")
[~,I] = max(abs(xdft2));
fprintf('Beta:Maximum occurs at %3.2f Hz.\n',freq2(I));
D7 = detrend(D7,0);
xdft3 = fft(D7);
freq3 = 0:N/length(D7):N/2;
xdft3 = xdft3(1:length(D7)/2+1);
% figure;
subplot(5,1,3);
plot(freq3,abs(xdft3));
title('ALPHA');
xlabel("time")
ylabel("freq")
[~,I] = max(abs(xdft3));
fprintf('Alpha:Maximum occurs at %f Hz.\n',freq3(I))
 xdft4 = fft(D8);
freq4 = 0:N/length(D8):N/2;
xdft4 = xdft4(1:length(D8)/2+1);
% figure;
subplot(5,1,4);
plot(freq4,abs(xdft4));
title('THETA');
xlabel("time")
ylabel("freq")
[~,I] = max(abs(xdft4));
fprintf('Theta:Maximum occurs at %f Hz.\n',freq4(I));
A8 = detrend(A8,0);
xdft5 = fft(A8);
freq5 = 0:N/length(A8):N/2;
xdft5 = xdft5(1:length(A8)/2+1);
% figure;
subplot(5,1,5);
plot(freq3,abs(xdft5));
title('DELTA');
xlabel("time")
ylabel("freq")
[~,I] = max(abs(xdft5));
fprintf('Delta:Maximum occurs at %f Hz.\n',freq5(I));
%feature extraction
%finding peaks of signal
pks_gamma=findpeaks(cD5)
pks_beta=findpeaks(cD6)
pks_alpha=findpeaks(cD7)
pks_theta=findpeaks(cD8)
pks_delta=findpeaks(cA8)
max_gamma=maxpeak(pks_gamma);
max_beta=maxpeak(pks_beta);
ma_alpha=maxpeak(pks_alpha);
max_theta=maxpeak(pks_theta);
max_delta=maxpeak(pks_delta);

%finding the standard deviation
sd1=std(pks_gamma);
sd2=std(pks_beta)
sd3=std(pks_alpha)
sd4=std(pks_theta)
sd5=std(pks_delta)
sd_v=[sd1,sd2,sd3,sd4,sd5]
p1=pks_delta(1:5)
p2=pks_gamma(1:5)
max_peakval=[max_gamma,max_beta,ma_alpha,max_theta,max_delta]
sd_v=[sd1,sd2,sd3,sd4,sd5]
x=[pks_gamma(1:5);pks_gamma(1:5)]

ans=findepilepsy(max_gamma,max_beta,ma_alpha,max_theta,max_delta);

if(ans==1)
    disp("epilepsy is present")
end
if(ans==0)
    disp('epilepsy is not present')
end
%function to get average value of peak of components signal
function y_peak = maxpeak(a)
 avg=0;
 count=0;
for i=1:length(a)
    y_peak=0;
    if(a(i)>0)
        count=count+1;
        avg=a(i)+avg;
    end
    y_peak=avg/count;
   
end
end
%function to find if epilepsy is present or not

function epilepsy = findepilepsy(n1, n2, n3, n4, n5)
epilepsy_present=[];
%This function calculates the maximum of the five values given as input
 max =  n5;
  if n1 > max
    epilepsy=0
    epilepsy_present=[epilepsy_present,epilepsy];
    return
  end
  if n2 > max
    epilepsy=0
        epilepsy_present=[epilepsy_present,epilepsy];
    return
  end

  if n3 > max
    epilepsy=0
        epilepsy_present=[epilepsy_present,epilepsy];
    return
  end
  if n4 > max
    epilepsy=0
        epilepsy_present=[epilepsy_present,epilepsy];
    return
  end
  epilepsy=1;
      epilepsy_present=[epilepsy_present,epilepsy];
end



           
          