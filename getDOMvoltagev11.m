%% Disclaimer
%Please read%
%this script is messy but commented. :( also it might be broken 
%also it has some plot color issues?
%also good luck

%% README
%{

This file will calculate
1. Voltages
2. Velocities
3. Velocity Autocorrelation
4. Density of Modes

For a given number of devices at a given frequency (rate) and a given time
(timespan) entered below. THE DEVICES MUST MEASURE VOLTAGE.

The file input must be columns of voltage, one for each measurement (no
time)

IMPORTANT: You need the SIGNAL PROCESSING TOOLBOX (so go get it now)
%}

%% Update Me
%CALLING THE TXT%
clear all
fileID = fopen('720whitenoise4.txt');
txtvoltageunconcatenated = importdata('720whitenoise4.txt');
txtvoltage = sum(txtvoltageunconcatenated, 1); %summing up all the rows

%Initialize the sauce%

%noise = wgn(200000,1,0); %use if you want to generate a white noise signal
rate = 500000; %UPDATE THIS---
timespan = 2; %UPDATE THIS---
nosamples = rate*timespan;
nodevices = 1; %UPDATE THIS---

closestpow2 = nextpow2(2*nosamples-1);%not used?

%% Initialize 
voltagevTotal = zeros([nosamples 1]); %initialize the voltage vector (this is used)

velocityvTotal =zeros([nosamples 1]); %initialize the velocity vector

VACvTotal = zeros([2*nosamples 1]); %initialize the VAC vector

DoMvTotal = zeros([2*nosamples 1]); %initialize the DoM vector

%% Total Voltage, Velocity, VAC and DoM
%mix the secret geophone ingredients%
%we can take each column, do a calculation, then throw it add it to the
%total and throw it out to make things faster. 
for i = 1:nodevices
    voltvectorcurrent = txtvoltageunconcatenated(:,i); %take the ith column
    %voltvectorcurrent = voltvectorcurrent(1:1000); %cut out the bullshit
    %at the start and end- used for splicing the data
    voltvectorcurrent = voltvectorcurrent - mean(voltvectorcurrent); %subtract the mean manually 
    %no clue why this works but the other way doesnt theyre the same. 
    %Dont question it... We subtract mean and THEN put it into getDoM
    [voltagevcurr, velocityvcurr, VACvcurr, DOMvcurr]= getDoM(voltvectorcurrent, rate, timespan);
    
    
    voltagevTotal= voltagevTotal+voltagevcurr; %summing up all the voltages? 
    %why am I doing this? Why is this useful? Dont ask questions Adi just
    %do what Ted says
    
    velocityvTotal = velocityvTotal+velocityvcurr; %i dont think this is used
    
    VACvTotal = VACvTotal+VACvcurr;
    
    DoMvTotal = DoMvTotal+DOMvcurr;
end

%The main getDOM function below ensures <V> and <v> are centered. 

%% Voltages vs Time (input)

%basically, i dont know how to create new variables using a forloop so i
%might bite the bullet and do it manually (pain is beginning)
tv4 = timekeeper(nosamples, 1/rate);
for i = 1:nodevices
    voltvectorcurrent = txtvoltageunconcatenated(:,i);
    %voltvectorcurrent = voltvectorcurrent(1:1000); %can be used to cut the
    %bad data at the start and end of the voltage file. (will change length
    %so you need to change tv4 then of course)
    voltvectorcurrent = voltvectorcurrent - mean(voltvectorcurrent); %SUBTRACT MEAN
    [voltagevcurr, velocityvcurr, VACvcurr, DOMvcurr]= getDoM(voltvectorcurrent, rate, timespan);
    
    plot(tv4, voltagevcurr)
    hold on
end
hold off
%axis([0 0.2 -0.00001 0.00001])
title('Voltage vs Time') 
ylabel('Voltage (V)')
xlabel('Time (s)')
%legend('0% strain', '2% strain', '4% strain', '6% strain' , '8% strain')
%this legend command sucks and I need a better way to do it but I dont care
%enough
figure()

%% Velocities vs Time

tv5 = timekeeper(nosamples, 1/rate);
tv5transpose = transpose(tv5);
for i = 1:nodevices
    voltvectorcurrent = txtvoltageunconcatenated(:,i);
    %voltvectorcurrent = voltvectorcurrent(1:1000); %use for splicing
    voltvectorcurrent = voltvectorcurrent - mean(voltvectorcurrent); %SUBTRACT MEAN
    [voltagevcurr, velocityvcurr, VACvcurr, DOMvcurr]= getDoM(voltvectorcurrent, rate, timespan);
    
    plot(tv5, velocityvcurr)
    hold on
end
hold off
%axis([0 3 -0.005 0.005])
title('Velocity vs Time')
ylabel('Velocity (au)')
xlabel('Time (s)')
figure()

%% Veloctiy Autocorr vs Lag Time(half of duration) (Avg. in black)

tv6 = timekeeper(2*nosamples, 1/rate);
%tv6 = tv6 - timespan/2; %wait why is this happening? 
tv6transpose = transpose(tv6);
for i = 1:nodevices
    voltvectorcurrent = txtvoltageunconcatenated(:,i);
    %voltvectorcurrent = voltvectorcurrent(1:1000);
    voltvectorcurrent = voltvectorcurrent - mean(voltvectorcurrent);
    [voltagevcurr, velocityvcurr, VACvcurr, DOMvcurr]= getDoM(voltvectorcurrent, rate, timespan);
    
    plot(tv6, real(VACvcurr))
    hold on
end

%plot(tv6, VACvTotal/nodevices, "black", "LineWidth", 1) %plotting the average in black

hold off
title('V.A.C. vs Delay Time')
ylabel('Velocity Autocorrelation (average in black)')
xlabel('Time (s)')
xlim([0, timespan/2])
ylim([-1, 1])
figure()

%% Density of Modes vs Freq (avg. in black)

fv1 = 0.5*rate*(0:(nosamples/2))/nosamples; %boy i hope this works!
fv1transpose = transpose(fv1);
lenfv1 = length(fv1);
slope1 = zeros(lenfv1,1);
slope2 = zeros(lenfv1,1);


for i = 1:nodevices
    voltvectorcurrent = txtvoltageunconcatenated(:,i);
    %voltvectorcurrent = voltvectorcurrent(1:1000);
    voltvectorcurrent = voltvectorcurrent - mean(voltvectorcurrent);
    [voltagevcurr, velocityvcurr, VACvcurr, DOMvcurr]= getDoM(voltvectorcurrent, rate, timespan);
    
    P2 = abs(DOMvcurr/nosamples); %i stole this from matlab documentation on ffts
    P1 = P2(1:nosamples/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
    %okay so basically I dont really know what the three lines above do,
    %but I have reason to believe that it makes sense. Also it makes plots
    %that make Karen and Ted fairly happy and that is my job.
    
    %Update(v8): I now understand it. 
   
    
    loglog(fv1, P1)
    hold on
end


P2 = abs(DoMvTotal/nosamples);
P1 = P2(1:nosamples/2+1);
P1(2:end-1) = 2*P1(2:end-1);
    

for k = 1:lenfv1
    slope1(k,1) = 1/fv1(1,k) * exp(-12.5);
    slope2(k,1) = 1/power(fv1(1,k),2) * exp(11);
end

%loglog(fv1, P1/nodevices, "black")
%loglog(fv1(1:25),slope1(1:25), "black", 'Linewidth', 1.25)
loglog(fv1(250:end),slope2(250:end), "black", 'Linewidth', 1.25)


hold off

title('Density of Modes')
ylabel('Density of Modes')
xlabel('Frequency (Hz)')
xlim([0, rate/2])
%ylim([1/1000000000000000000000000, 1/1000000000000]) %just dont ask
figure()

%% Normalize the DoM by debye (avg. in black)
%my attempt to remove Debye scaling
%this section SHOULD NOT EXIST, so if there are too many bugs, just get rid
%of the whole damn thing and it will not be missed.

fv2 = 0.5*rate*(0:(nosamples/2))/nosamples;
lenfv2 = length(fv2);
debye = zeros(lenfv2, 1); %normalize by this

for i = 1:nodevices
    voltvectorcurrent = txtvoltageunconcatenated(:,i);
    %voltvectorcurrent = voltvectorcurrent(1:1000);
    voltvectorcurrent = voltvectorcurrent - mean(voltvectorcurrent);
    [voltagevcurr, velocityvcurr, VACvcurr, DOMvcurr]= getDoM(voltvectorcurrent, rate, timespan);
    
    P2 = abs(DOMvcurr/nosamples);
    P1 = P2(1:nosamples/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
    for k = 1:lenfv2
        debye(k,1) = power(fv2(1,k),0.5);
        P3 = P1./debye; %the normalization
    end
    
    loglog(fv2, P3)
    hold on
end


P2 = abs(DoMvTotal/nosamples);
P1 = P2(1:nosamples/2+1);
P1(2:end-1) = 2*P1(2:end-1);


for k = 1:lenfv2
    debye(k,1) = power(fv2(1,k),0.5);
    P3 = P1./debye; %the normalization
end
 
    
loglog(fv2, P3, "black")


hold off

title('DoM scaled for Debye Scaling')
ylabel('Density of Modes (normalized)')
xlabel('Frequency (Hz)')
xlim([0, rate/2])
%ylim([1/1000000000000000000000000, 1/1000000000000])
figure()

%% Voltage Semilogy for checking noise 
%take the fft of the voltage vector total, then conjugate, multiply and
%plot with timekeeper

%PSD is so useful for diagnosing issues, boy I wish someone wrote a script
%to plot all the PSDs next to each other. But I have no time. 

xdft = fft(voltagevTotal);
xdft = xdft(1:nosamples/2+1); %lets take a good length
psdx = (1/(rate*nosamples)) * abs(xdft).^2; %now
psdx(2:end-1) = 2*psdx(2:end-1);
fv3 = 0:rate/length(voltagevTotal):rate/2;


semilogy(fv3 , psdx)

%set(gca,'yscale','log')

title('PSD vs f') %do I understand this?
ylabel('Power Spectral Density')
xlabel('Frequency (Hz)')
xlim([0, rate/2])
ylim([1/1000000000000000 1])
figure()

%% Helper functions
%Misc Helper functions-----------------------------------------------------

%the timekeeper vector creator- the timekeeper keeps the time
%do not bother the timekeeper

function tvector = timekeeper(n, timedifference)
tvector = zeros([n-1 1]);

    for i = 1:n-1 
        tvector(i+1)= tvector(i)+timedifference;
    end
tvector = tvector';
end


%the main function
%4 outputs to play with
function[voltout, velocityout, VACout, DoMout] = getDoM(voltvector, ratescalar, durationscalar)

maxlensc = length(voltvector);

voltmeansc = mean(voltvector);
%voltout = voltvector - voltmeansc; %why does this not work? do not use i
%guess and i will subtract mean before putting into the func
voltout = voltvector;

%this is the first filteirng attempt (at voltage state)
%voltoutfft = fft(voltout);
%for i = 120:180
%    voltoutfft(i) = voltoutfft(i)/200;
%end
%voltoutfftinv = fft(voltoutfft);
%I am deciding for no reason at all to remove this filter because it
%definitely works as intended


velocityout = cumtrapz(voltout); %the integral is delivered
velocitymeansc = mean(velocityout);
velocityout = velocityout - velocitymeansc; %this works but the other one doesnt?

%we need to subtract mean before doing autocorr

%Filter number 2: you need SIGNAL PROCESSING TOOLBOX for this guy
velocityout = bandstop(velocityout, [47 68], ratescalar);
%works like a charm
%?

%now we pad with zeros of length equal to the vector
velocitylength = length(velocityout);
velocitypad = zeros([velocitylength 1]);
velocitytovac = cat(1,velocityout, velocitypad); %padded

%now to implement VAC in fourier space with Karen's method
velocityfft = fft(velocitytovac);
for i = 1:length(velocityfft)
    velocityfft(i) = velocityfft(i) * conj(velocityfft(i));
end

%Filter 3- works but is violent?
%now we can clean out the 60Hz noise
%for i = 225:260
%    velocityfft(i) = velocityfft(i)/1000;
%end
%I remove this one because she is too crazy

velocityfftinverse = ifft(velocityfft);
VACout = velocityfftinverse;

%normalize the VAC
VACout = VACout/VACout(1);

%VACout = xcorr(velocityout, 0.5*ratescalar*durationscalar,
%'normalized');%old VAC method for the newbies (I never used this) (real)

pow2 = nextpow2(length(VACout));
FFTVAC = fft(VACout); %removed padding

DoMout = real(FFTVAC);

%fuck- did i forget a factor of 2?
end