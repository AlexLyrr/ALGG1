function disdata = discrim(x)

% function disdata = discrimf(x)
% x is the received signal in complex baseband form
% Mark Wickert

X=real(x);           % X is the real part of the received signal
Y=imag(x);           % Y is the imaginary part of the received signal
N=length(x);         % N is the length of X and Y
b=[1 -1];            % filter coefficients for discrete derivative
a=[1 0];             %    "
derY=filter(b,a,Y);  % derivative of Y, 
derX=filter(b,a,X);  %    "          X,

disdata=(X.*derY-Y.*derX)./(X.^2+Y.^2);
        

