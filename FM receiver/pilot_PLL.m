function [theta, phi_error] = pilot_PLL(xr,fq,fs,loop_type,Bn,zeta)
% [theta, phi_error] = pilot_PLL(xr,fq,fs,loop_type,Bn,zeta)
%
%
% Mark Wickert, April 2014

T = 1/fs;
% Set the VCO gain in Hz/V  
Kv = 1.0;
% Design a lowpass filter to remove the double freq term
Norder = 5;
[b_lp,a_lp] = butter(Norder,2*(fq/2)/fs);
fstate = zeros(1,Norder); % LPF state vector

Kv = 2*pi*Kv; % convert Kv in Hz/v to rad/s/v

if loop_type == 1
    % First-order loop parameters
    fn = Bn;
    Kt = 2*pi*fn; % loop natural frequency in rad/s
elseif loop_type == 2
    % Second-order loop parameters
    fn = 1/(2*pi) * 2*Bn/(zeta + 1/(4*zeta)); % given Bn in Hz
    Kt = 4*pi*zeta*fn; % loop natural frequency in rad/s
    a = pi*fn/zeta;
else
    error('Loop type must be 1 or 2');
end

% Initialize integration approximation filters
filt_in_last = 0; filt_out_last = 0;
vco_in_last = 0; vco_out = 0; vco_out_last = 0;

% Initialize working and final output vectors
n = 0:length(xr)-1;
theta = zeros(size(xr));
ev = zeros(size(xr));
phi_error = zeros(size(xr));
% Normalize total power in an attemp to make the 19kHz sinusoid
% component have amplitude ~1.
%xr = xr/(2/3*std(xr));
% Begin the simulation loop
for k = 1:length(n)
    % Sinusoidal phase detector (simple multiplier)
    phi_error(k) = 2*xr(k)*sin(vco_out);
    % LPF to remove double frequency term
    [phi_error(k),fstate] = filter(b_lp,a_lp,phi_error(k),fstate);
    pd_out = phi_error(k);
    %pd_out = 0;
    % Loop gain
    gain_out = Kt/Kv*pd_out; % apply VCO gain at VCO
    % Loop filter
    if loop_type == 2
        filt_in = a*gain_out;
        filt_out = filt_out_last + T/2*(filt_in + filt_in_last);
        filt_in_last = filt_in;
        filt_out_last = filt_out;
        filt_out = filt_out + gain_out;
    else
        filt_out = gain_out;
    end
    % VCO
    vco_in = filt_out + fq/(Kv/(2*pi)); % bias to quiescent freq.
    vco_out = vco_out_last + T/2*(vco_in + vco_in_last);
    vco_in_last = vco_in;
    vco_out_last = vco_out;
    vco_out = Kv*vco_out; % apply Kv;
    % Measured loop signals
    ev(k) = filt_out;
    theta(k) = mod(vco_out,2*pi); % The vco phase mod 2pi
end