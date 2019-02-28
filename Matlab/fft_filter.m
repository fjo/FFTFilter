%Setup
% clear all
close all

% Make some data
fs = 2048;
T = 10;
t = linspace(0,T-1/fs,T*fs);
x = sin(2*pi*100*t) + sin(2*pi*500*t);
xt = sin(2*pi*100*t);

% Plot the data
figure('position',[125   545   550   420])
title('Input signal')
pidx = 1:100;
plot(t(pidx),x(pidx))

% Design a filter to take out the higher tone
%   Order 24 will produce 25 coeff, group delay = N-1 /2

N     = 24;
Fpass = 200;  % Passband Frequency
Fstop = 400;  % Stopband Frequency
Wpass = 1;    % Passband Weight
Wstop = 1;    % Stopband Weight
dens  = 20;   % Density Factor

% Calculate the coefficients using the FIRPM function.
b  = firpm(N, [0 Fpass Fstop fs/2]/(fs/2), [1 1 0 0], [Wpass Wstop], ...
           {dens});
% Filter data in time domain
y = filter(b,1,x);

gdelay = N/2;
figure('position',[680   545   550   420])
title('Time Domain filtered output')
plot(x(pidx));
hold on
plot(y(pidx),'r');
plot(y(pidx+gdelay),'r--')

%params for fft filter
fftsize = 256;
overlap = N ;      

% Zero pad Overlap and apply filter
B = fft(b,fftsize);
xbuf = buffer(x,fftsize-overlap,0);
XBUF = bsxfun(@times,fft(xbuf,fftsize),transpose(B));
ybuf = ifft(XBUF);

figure('position',[1235   545   550   420])
title('First frame of Freq Domain filter')
plot(y(1:512))
hold on
plot(ybuf(:,1),'r')

% now do the add part of "overlap add"
yoa = zeros(1,prod(size(ybuf)));
yoa(1:fftsize) = ybuf(:,1);
for k = 2:size(ybuf,2)
   start = (k-1)*(fftsize-overlap)+1;
   stop = start + fftsize-1;
   yoa(start:stop) = yoa(start:stop) + transpose(ybuf(:,k));
end

%% This is a bunch of code I used to test and verify results
figure('position',[680   30   550   420])
plot(yoa(1:768))
hold on
plot(y(1:768),'r')

b1 = ifft(fft(x(1:232),256).*B);
b2 = ifft(fft(x(233:464),256).*B);

xbuf = buffer(x,232,0);
bb1 = ifft(fft(xbuf(:,1)',256).*B);
bb2 = ifft(fft(xbuf(:,2),256).*transpose(B));

%add together
bbc = zeros(1,512);
bbc(1:256) = bb1;
bbc(233:488) = bbc(233:488)+bb2.';

yc = zeros(1,512);
yc(1:256) = ybuf(:,1).';
cidx = fftsize-overlap+1:2*fftsize-overlap;
yc(cidx) = yc(cidx)+ybuf(:,2).';

plot(yc,'k--')

%% Downsample by 4
dsr = fftsize/4;
dfftsize = fftsize/4;
doverlap = overlap/4;
idx1 = 1:dsr;
idx2 = dsr+1:2*dsr;
idx3 = 2*dsr+1:3*dsr;
idx4 = 3*dsr+1:4*dsr;

% alais the signal!
XDSB = XBUF(idx1,:) + XBUF(idx2,:) + XBUF(idx3,:) +XBUF(idx4,:);
ybuf = real(ifft(XDSB))./4;

yoad = zeros(1,prod(size(ybuf)));
yoad(1:dfftsize) = ybuf(:,1);
for k = 2:size(ybuf,2)
   start = (k-1)*(dfftsize-doverlap)+1;
   stop = start + dfftsize-1;
   yoad(start:stop) = yoad(start:stop) + ybuf(:,k).';
end












