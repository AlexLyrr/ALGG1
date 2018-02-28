function [rx_symb_d,clk,track] = sccs_bit_sync(y,Ns)
%
%//////////////////////////////////////////////////////
% Symbol synchronization algorithm using SCCS
%//////////////////////////////////////////////////////
%     y = baseband NRZ data waveform
%    Ns = nominal number of samples per symbol
%
% Reworked from ECE 5675 Project
% Mark Wickert April 2014

rx_symb_d = zeros(fix(length(y)/Ns),1); % decimated symbol sequence for SEP
track = zeros(fix(length(y)/Ns),1);
bit_count = 0;
y_abs = zeros(size(y));
clk = zeros(size(y));
k = Ns+1; %initial 1-of-Ns symbol synch clock phase;
% Sample-by-sample processing required
for i=1:length(y), 
    %y_abs(i) = abs(round(real(y(i))));
    if i >= Ns % do not process first Ns samples
        % Collect timing decision unit (TDU) samples:
        y_abs(i) = abs(sum(y(i-Ns+1:i)));
        % Update sampling instant and take a sample
        % For causality reason the early sample is 'i',
        % the on-time or prompt sample is 'i-1', and  
        % the late sample is 'i-2'.
        if (k == 0)
            % Load the samples into the 3x1 TDU register w_hat.
            % w_hat[1] = late, w_hat[2] = on-time; w_hat[3] = early.
            w_hat = y_abs(i-2:i);
            bit_count = bit_count + 1;
            if w_hat(2) ~= 0,
                if w_hat(1) < w_hat(3),
                    k = Ns-1;
                    clk(i-2) = 1;
                    rx_symb_d(bit_count) = y(i-2-round(Ns/2));
                elseif w_hat(1) > w_hat(3)
                    k = Ns+1;
                    clk(i) = 1;
                    rx_symb_d(bit_count) = y(i-round(Ns/2));
                else
                    k = Ns;
                    clk(i-1) = 1;
                    rx_symb_d(bit_count) = y(i-1-round(Ns/2));
                end
            else
                k = Ns;
                clk(i-1) = 1;
                rx_symb_d(bit_count) = y(i-1-round(Ns/2));
            end
            track(bit_count) = mod(i,Ns);
        end
    end
    k = k - 1;
end
% Trim the final output to bit_count
rx_symb_d = rx_symb_d(1:bit_count);
end