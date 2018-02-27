filtord = 63;                                       % Filter Order
Fs = 2400000;                                           % Sampling Frequency
Fn = Fs/2;                                          % Nyquist Frequency
Fc = 70000;                                            % Cutoff Frequency
blp = fir1(filtord, Fc/Fn);                         % Calculate Filter Coefficients

freqz(blp,1,[],Fs)                                  % Plot the filter

yB1 = filter(blp,1,x)    % pass the signal through the low pass filter
yN1 = downsample(yB1,10) % a downsampling factor of 10
simpleSA(yN1,2^14,2400)
Zdis = discrim(yN1)

filtord = 63;                                       % Filter Order
Fs = 200000;                                           % Sampling Frequency
Fn = Fs/2;                                          % Nyquist Frequency
Fc = 15000;                                            % Cutoff Frequency
blp = fir1(filtord, Fc/Fn);                         % Calculate Filter Coefficients

yB2 = filter(blp,1,Zdis)
yN2 = downsample(yB2,5)   %pass the signal through the second decimator

sound(yN2,48000)

f3=1/(2*pi*0.000075)
a1= exp(-2*pi*f3/48000)
a =[1, -a1]
b=[1-a1]
y = filter(b,a,yN2)

sound(y,48000)
