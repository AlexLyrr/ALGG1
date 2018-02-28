function [Px,F] = simpleSA(x,N,fs,min_dB,max_dB,color)
% [Px,F] = simpleSA(x,N,fs,min_dB,max_dB,color)
% Plot the estimated power spectrum of real and complex baseband signals
% using the toolbox function psd(). This replaces the use of psd().
% using the toolbox function psd(). This replaces the use of psd().
%
%      x = complex baseband data record, at least N samples in length
%      N = FFT length
%     fs = Sampling frequency in Hz
% min_dB = floor of spectral plot in dB
% max_dB = ceiling of spectral plot in dB
%  color = line color and type as a string, e.g., 'r' or 'r--
%
%////////// Optional output variables ///////////////////////////////////
% Px = power spectrum estimate values on two-sided frequency axis
%  F = the corresponding frequency axis values
%
% Mark Wickert, October 2008

% while length(x) < N
%     N = N/2;
% end 

warning off
if isreal(x)
    [Px,F] = pwelch(x,N,[],[],fs) %psd(x,N,fs);
else
    %[Px,F] = psd(x,[],[],N,fs);
    [Px,F] = pwelch(x,N,[],[],fs) %psd(x,N,fs);
    N = fix(length(F)/2);
    F = [F(end-N+1:end)-fs; F(1:N)];
    Px = [Px(end-N+1:end); Px(1:N)];
end
warning on

if nargout == 2
    return
end

if nargin == 3 || nargin == 5,
    plot(F,10*log10(Px))
else
    plot(F,10*log10(Px),color)
end

if nargin >= 5
    axis([F(1) F(end) min_dB max_dB]);
end

xlabel('Frequency (kHz)')
ylabel('Power Spectrum in dB')
grid on