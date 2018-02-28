function y = loadFile(filename)
%  y = loadFile(filename)
%
% reads  complex samples from the rtlsdr file 
%
% SDR>rtl_sdr -s f_in_sps -f f_in_Hz -g gain_dB capture.bin
% SDR>rtl_sdr -s 2400000 -f 70000000 -g 25 capture70R1k.bin

fid = fopen(filename,'rb');
y = fread(fid,'uint8=>double');

y = y-127;
y = y(1:2:end) + i*y(2:2:end);
