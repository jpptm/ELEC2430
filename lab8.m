Lab 8: DFT and Digital Filters
Johanne Montano
Ismail Adebayo
Excercise 1: Signal Analysis using FFT

1.1 Reproducing the sampled signal
Option A: T = 0.05
Option B: T = 0.02
clear all;
% Constants
tau = 0.5;
w0 = 60;

% Sampling Periods
Ta = 0.05;
Tb = 0.02;

% Time duration of samples
ta = -2:Ta:2;
tb = -2:Tb:2;

% Defining each pulse Function
pulseA = rectangularPulse(-tau/2,tau/2,ta);
pulseB = rectangularPulse(-tau/2,tau/2,tb);

% Chuck values in for each option
xa = exp(-abs(ta)).*cos(w0.*ta).*pulseA;
xb = exp(-abs(tb)).*cos(w0.*tb).*pulseB;
% xb = xb(1:100)

% Plot the values
figure
plot(ta,xa);
title("Option A - Sampled function")
grid
xlabel("n")
ylabel("x[n]")

figure
plot(tb,xb);
grid
title("Option B - Sampled function")
xlabel("n")
ylabel("x[n]")
1.2 Computing the Discrete Fourier Transform
% Calculate discrete fourier transform for each option
Xka = fft(xa);
Xkb = fft(xb);
fa = linspace(0, 1/Ta, length(Xka));
fb = linspace(0, 1/Tb, length(Xkb));

% Normalise and chop graph by half to avoid repeating pulse    
%cut = (1/Tb)*(0:floor(length(Xkb)/2))/(length(Xkb));
%Xkb = Xkb(1:length(cut));

figure
plot(fa, abs(Xka));
xlim([0 1/Ta/2]);
grid
title("Option A - Magnitude plot")
xlabel("f")
ylabel("|X(w)|")

figure
plot(fb, abs(Xkb));
xlim([0 1/Tb/2]);
grid
title("Option B - Magnitude plot")
xlabel("f")
ylabel("|X(w)|")

Excercise 2: Filtering an Audio Signal
The audio files given for each option are as following:
Option A: bugsbunny1.wav
Option B: daffyduck1.wav
2.1 Loading the sample signal
[xa,Fsa] = audioread('bugsbunny1.wav');
[xb,Fsb] = audioread('daffyduck1.wav');
The sampling frequency for both options  is found to be 11.025Hz
% Confirming the sound at normal voice sampling frequency
%sound(xa,8192)
%sound(xb,8192)
Playing the sound at a higher sample rate than the rate in which it is recorded causes a noticable increase of pitch in the playback, which is related to the increase in frequency. Playing the sound a lower rate, the rate in which it was recorded means the normal voice is able to be heard.
2.2 Computing the Discrete Fourier Transform
% Using the fast fourier transform function
Xka = fft(xa);
Xkb = fft(xb);

% Normalising the sample frequency with an array half the length of the
% fourier transform
fa = Fsa*(0:floor(length(Xka)/2))/(length(Xka));
fb = Fsb*(0:floor(length(Xkb)/2))/(length(Xkb));

% Finds half the length of the fourier series (to the nearest whole number)
newXka = Xka(1:length(fa));
newXkb = Xkb(1:length(fb));  

% Plot of both Magnitude Responces
figure()
plot(fa,abs(newXka)), xlabel('Frequency (Hz)'), ylabel('|Xka|')
title('Option A: Magnitude Spectrum');
figure()
plot(fb,abs(newXkb)), xlabel('Frequency (Hz)'), ylabel('|Xkb|')
title('Option B: Magnitude Spectrum');
2.3 Butterworth Filter Design
% n is order of the butterworth filter, wn is the normalised bandwidth
    % Incresed order, increases the sharpness of the frequency responce
    % (makes it more square)
   
% Cutoff frequency (wc) - same as the bandwith for a low pass filter
wc = 2500; % 0 to 2.5kHz, has a bandwith of 2.5kHz
n = 10; % 10th Order butterworth filter

% Normalised cutoff frequency (wn), so that: 0 < wn < 1  
wna = wc/(Fsa/2); % Half the sample rate
wnb = wc/(Fsb/2);

% num is vector Y coefficient
% den is vector X coefficient
[numa,dena] = butter(n,wna); %[B,A] = butter(filterorder,normalisedcutoff)
[numb,denb] = butter(n,wnb);
2.4 Recursion
Option A
x0a = zeros(1,length(numa)-1); % Initial conditions
y0a = zeros(1,length(dena)-1);
dena = dena(2:1:n+1); % Ensures dimentions suit recur function
ta = 0:1:(length(xa)-1);
ya = recur(dena, numa, ta, xa', x0a, y0a); % Calculate approximation recursively
Option B
x0b = zeros(1,length(numb)-1); % Initial conditions
y0b = zeros(1,length(denb)-1);
denb = denb(2:1:n+1); % Ensures dimentions suit recur function
tb = 0:1:(length(xb)-1);
yb = recur(denb, numb, tb, xb', x0b, y0b); % Calculate approximation recursively
2.5 Comparison of Before/After filtering
Option A
% Play sound
%sound(ya,Fsa)

% Repeat 2.2
    % Using the fast fourier transform function
    Ya = fft(ya);
    % Normalising the sample frequency with an array half the length of the
    % fourier transform
    Yfa = Fsa*(0:(length(Ya)/2))/(length(Ya));
    % Finds half the length of the fourier series (to the nearest whole number)
    newYa = Ya(1:length(Yfa));  

% Display plots
figure()
subplot(2,1,1)
plot(fa,abs(newXka)),title('Option A: Before Butterworth filter is applied to the signal')
xlabel('f [Hz]'), ylabel('|Xka|'), grid on, grid minor

subplot(2,1,2)
plot(Yfa,abs(newYa)), title('Option A: After Butterworth filter applied')
xlabel('f [Hz]'), ylabel('|Ya|'), xline(wc,'--r','cutoff frequency')
grid on, grid minor
% Waveform Plots
figure()
subplot(2,1,1) % Before filtering
aDuration = audioinfo('bugsbunny1.wav').Duration; % Get song duration, for plotting against time
plot(linspace(0,aDuration,length(xa)),xa),title('Option A: Waveform before 2.5kHz filter')
xlabel('t (sec)'), ylabel('xa(t)')

subplot(2,1,2) % After filtering
plot(linspace(0,aDuration,length(ya)),ya), title('Option A: Waveform after 2.5kHz filter')
xlabel('t (sec)'), ylabel('ya(t)')
Option B
% Play sound
%sound(yb,Fsb)

% Repeat 2.2
    % Using the fast fourier transform function
    Yb = fft(yb);
    % Normalising the sample frequency with an array half the length of the
    % fourier transform
    Yfb = Fsa*(0:(length(Yb)/2))/(length(Yb));
    % Finds half the length of the fourier series (to the nearest whole number)
    newYb = Yb(1:length(Yfb));  

% Display plots
figure()
subplot(2,1,1)
plot(fb,abs(newXkb)),title('Option B: Before Butterworth filter is applied to the signal')
xlabel('f [Hz]'), ylabel('|Xkb|'), grid on, grid minor

subplot(2,1,2)
plot(Yfb,abs(newYb)), title('Option B: After Butterworth filter applied')
xlabel('f [Hz]'), ylabel('|Yb|'), xline(wc,'--r','cutoff frequency')
grid on, grid minor
% Waveform Plots
figure()
subplot(2,1,1) % Before filtering
bDuration = audioinfo('daffyduck1.wav').Duration; % Get song duration, for plotting against time
plot(linspace(0,bDuration,length(xb)),xb),title('Option B: Waveform before 2.5kHz filter')
xlabel('t (sec)'), ylabel('xb(t)')

subplot(2,1,2) % After filtering
plot(linspace(0,bDuration,length(yb)),yb), title('Option B: Waveform after 2.5kHz filter')
xlabel('t (sec)'), ylabel('yb(t)')
Result of 2.5kHz low-pass Butterworth filter
For both sound files the magnitude decays noticably at the 2.5kHz cutoff frequency. There was little noticeable difference in sound quality with only a slight decrease in depth and the audio could still be heard clearly. The waveform for each looked quite similar also, the option B file (daffyduck.wav) had a noticable reduction in the amplitude of peaks.
2.6 Repeating with a 1200Hz Bandwith
1200Hz - 2.3 Butterworth Filter Design
% n is order of the butterworth filter, wn is the normalised bandwidth
    % Incresed order, increases the sharpness of the frequency responce
    % (makes it more square)
   
% Cutoff frequency (wc) - same as the bandwith for a low pass filter
wc2 = 1200; % Bandwith of 1.2kHz
n = 10; % 10th Order butterworth filter

% Normalised cutoff frequency (wn), so that: 0 < wn < 1  
wna2 = wc2/(Fsa/2); % Half the sample rate
wnb2 = wc2/(Fsb/2);

% num is vector Y coefficient
% den is vector X coefficient
[numa2,dena2] = butter(n,wna2); %[B,A] = butter(filterorder,normalisedcutoff)
[numb2,denb2] = butter(n,wnb2);
1200Hz - 2.4 Recursion
Option A
x0a2 = zeros(1,length(numa2)-1); % Initial conditions
y0a2 = zeros(1,length(dena2)-1);
dena2 = dena2(2:1:n+1); % Ensures dimentions suit recur function
ta2 = 0:1:(length(xa)-1);
ya2 = recur(dena2, numa2, ta2, xa', x0a2, y0a2); % Calculate approximation recursively
Option B
x0b2 = zeros(1,length(numb2)-1); % Initial conditions
y0b2 = zeros(1,length(denb2)-1);
denb2 = denb2(2:1:n+1); % Ensures dimentions suit recur function
tb2 = 0:1:(length(xb)-1);
yb2 = recur(denb2, numb2, tb2, xb', x0b2, y0b2); % Calculate approximation recursively
1200Hz - 2.5 Comparison of Before/After filtering
Option A
% Play sound
%sound(ya2,Fsa)

% Repeat 2.2
    % Using the fast fourier transform function
    Ya2 = fft(ya2);
    % Normalising the sample frequency with an array half the length of the
    % fourier transform
    Yfa2 = Fsa*(0:(length(Ya2)/2))/(length(Ya2));
    % Finds half the length of the fourier series (to the nearest whole number)
    newYa2 = Ya2(1:length(Yfa));  

% Display plots
figure()
subplot(2,1,1)
plot(fa,abs(newXka)),title('Option A: Before Butterworth filter is applied to the signal')
xlabel('f [Hz]'), ylabel('|Xka|'), grid on, grid minor

subplot(2,1,2)
plot(Yfa2,abs(newYa2)), title('Option A: After Butterworth filter applied')
xlabel('f [Hz]'), ylabel('|Ya|'), xline(wc2,'--r','cutoff frequency')
grid on, grid minor
% Waveform Plots
figure()
subplot(2,1,1) % Before filtering
plot(linspace(0,aDuration,length(xa)),xa),title('Option A: Waveform before 1.2kHz filter')
xlabel('t (sec)'), ylabel('xa(t)')

subplot(2,1,2) % After filtering
plot(linspace(0,aDuration,length(ya2)),ya2), title('Option A: Waveform after 1.2kHz filter')
xlabel('t (sec)'), ylabel('ya2(t)')
Option B
% Play sound
%sound(yb2,Fsb)

% Repeat 2.2
    % Using the fast fourier transform function
    Yb2 = fft(yb2);
    % Normalising the sample frequency with an array half the length of the
    % fourier transform
    Yfb2 = Fsa*(0:(length(Yb2)/2))/(length(Yb2));
    % Finds half the length of the fourier series (to the nearest whole number)
    newYb2 = Yb2(1:length(Yfb2));  

% Display plots
figure()
subplot(2,1,1)
plot(fb,abs(newXkb)),title('Option B: Before 1.2kHz Butterworth filter is applied to the signal')
xlabel('f [Hz]'), ylabel('|Xkb|'), grid on, grid minor

subplot(2,1,2)
plot(Yfb2,abs(newYb2)), title('Option B: After 1.2kHz Butterworth filter applied')
xlabel('f [Hz]'), ylabel('|Yb|'), xline(wc2,'--r','cutoff frequency')
grid on, grid minor
% Waveform Plots
figure()
subplot(2,1,1) % Before filtering
bDuration = audioinfo('daffyduck1.wav').Duration; % Get song duration, for plotting against time
plot(linspace(0,bDuration,length(xb)),xb),title('Option B: Waveform before 1.2kHz filter')
xlabel('t (sec)'), ylabel('xb(t)')

subplot(2,1,2) % After filtering
plot(linspace(0,bDuration,length(yb)),yb), title('Option B: Waveform after 1.2kHz filter')
xlabel('t (sec)'), ylabel('yb2(t)')
Result of 1.2kHz low-pass Butterworth filter
There is a large decay when the frequency goes over the cutoff of 1200Hz. The lower bandwith causes the sound to have a shallow frequency range with little high notes, making the sound seem less deep and low quality.The sound seems muffled as if it is coming from another room with the door closed.

2.7 Adding Noise to the Signal
c = 0.2; % scaling parameter
wa = c*(rand(size(xa))-0.5); % Random noise for option A
wb = c*(rand(size(xb))-0.5); % Random noise for option B
A vibration needs to have large changes in amplitude over a short period of time (0.05ms to 20ms) to cause a sound resonance that can be heard. The amplitude of the noise is shifted down by 0.5 to ensure it is centred on the equalibrium and fluctuates between -c/2 and c/2. Making the length of the noise the same as the source allows for the noise to be added to the sound clip without matrix dimention errors and also makes sure the sound is preset throughout the entire clip.
2.7.1 Playback of Generated Noise and Noisy Sound Clip
%Playback
%sound(wa,Fsa)
%sound(wb,Fsb)

% Adding the noise to the sound clip
noisea = xa+wa;
noiseb = xb+wb;

%sound(noisea,Fsa)
%sound(noiseb,Fsb)
2.7.2 Filtering the Noisy Signal
% Cutoff frequency same as 2.3, so variables can be reused

% num is vector Y coefficient
% den is vector X coefficient
[numa3,dena3] = butter(n,wna); %[B,A] = butter(filterorder,normalisedcutoff)
[numb3,denb3] = butter(n,wnb);

% Recursion - Option A
x0a3 = zeros(1,length(numa3)-1); % Initial conditions
y0a3 = zeros(1,length(dena3)-1);
dena3 = dena3(2:1:n+1); % Ensures dimentions suit recur function
ya3 = recur(dena3, numa3, ta, xa', x0a3, y0a3); % Calculate approximation recursively

% Recursion - Option B
x0b3 = zeros(1,length(numb3)-1); % Initial conditions
y0b3 = zeros(1,length(denb3)-1);
denb3 = denb3(2:1:n+1); % Ensures dimentions suit recur function
yb3 = recur(denb3, numb3, tb, xb', x0b3, y0b3); % Calculate approximation recursively

% Playback Filtered Sound
%sound(ya3,Fsa)
%sound(yb3,Fsb)
2.7.3 Plotting the Spectra
% Finding the transfer functions

% Transform filtered signals
Yka3 = fft(ya3);
Yfa3 = Fsa*(0:(length(Yka3)/2))/(length(Yka3));  % Normalising the sample frequency
newYka3 = Yka3(1:length(Yfa3)); % Adjust length to match normalised sample freqency

Ykb3 = fft(yb3);
Yfb3 = Fsb*(0:(length(Ykb3)/2))/(length(Ykb3));  % Normalising the sample frequency
newYkb3 = Ykb3(1:length(Yfb3)); % Adjust length to match normalised sample freqency

% Transform noise
Wka3 = fft(wa);
newWka3 = Wka3(1:length(Yfa3));   % Adjust length to match normalised sample freqency
Wkb3 = fft(wb);
newWkb3 = Wkb3(1:length(Yfb3));   % Adjust length to match normalised sample freqency

% Transform noisy signal
Nka = fft(noisea);
newNka = Nka(1:length(Yfa3));   % Adjust length to match normalised sample freqency
Nkb = fft(noiseb);
newNkb = Nkb(1:length(Yfb3));   % Adjust length to match normalised sample freqency

% Plots
figure()
% Noise - A
subplot(3,2,1)
plot(Yfa3,abs(newWka3)),title('Option A: Noise')
xlabel('f [Hz]'), ylabel('|Wka|')

% Noise - B
subplot(3,2,2)
plot(Yfb3,abs(newWkb3)),title('Option B: Noise')
xlabel('f [Hz]'), ylabel('|Wkb|')

% Noisy signal - A
subplot(3,2,3)
plot(Yfa3,abs(newNka)), title('Option A: Noisy Signal')
xlabel('f [Hz]'), ylabel('|Nka|')

% Noisy signal - B
subplot(3,2,4)
plot(Yfb3,abs(newNkb)), title('Option B: Noisy Signal')
xlabel('f [Hz]'), ylabel('|Nkb|')

% Filtered signal - A
subplot(3,2,5)
plot(Yfa3,abs(newYka3)),title('Option A: Filtered Signal')
xlabel('f [Hz]'), ylabel('|Yka|')

% Filtered signal - B
subplot(3,2,6)
plot(Yfb3,abs(newYkb3)),title('Option B: Filtered Signal')
xlabel('f [Hz]'), ylabel('|Ykb|')
The high frequency noise is filtered away leaving only the low frequency noise which isn't as noticeable when mixed with the sound signal that is primarily low frequency.
