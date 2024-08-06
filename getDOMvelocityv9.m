%% disclaimer 
%Please read%
%The program is a mess- sorry I am as dissapointed as you are. 

%% README
%{
oioioi

This file will calculate
2. Velocities
3. Velocity Autocorrelation
4. Density of Modes

For a given number of devices at a given frequency (rate) and a given time
(timespan) entered below. THE DEVICES MUST MEASURE VELOCITY. (Hint: like a
vibrometer)

You need to input just a time vector line of all the velocities you need.
Thats all! (Multiple columns also works)

IMPORTANT: You need SIGNAL PROCESSING TOOLBOX for this. (go get that now
and come back)

The getDoMvoltage script is better commented, so it might help to take a
look at it. 
%}

%% Update Me
%CALLING THE TXT%
fileID = fopen('ScheckForest_noise2.txt');
txtvoltageunconcatenated = importdata('ScheckForest_noise2.txt'); %enter file here!
txtvoltage = sum(txtvoltageunconcatenated, 1);


%Initialize the sauce%

rate = 16000; %UPDATE THIS---
timespan = 4.24175; %UPDATE THIS---
nosamples = rate*timespan; %must be even because of an off-by-one thing
nodevices = 24; %UPDATE THIS---

closestpow2 = nextpow2(2*nosamples-1); %might not use this

%% Initialize
velocityvTotal =zeros([nosamples 1]);

VACvTotal = zeros([2*nosamples 1]);

DoMvTotal = zeros([2*nosamples 1]); %removed padding

%% Total Velocity, VAC and DoM
%mix the secret geophone ingredients%
%we can take each column, do a calculation, then throw it add it to the
%total and throw it out to make things faster. 
for i = 1:nodevices
    velvectorcurrent = txtvoltageunconcatenated(:,i); 
    %velvectorcurrent = velvectorcurrent(1:1000); %cut out the bullshit at
    %the start and end- use for splicing data only if needed
    velvectorcurrent = velvectorcurrent - mean(velvectorcurrent); %subtract the mean before getDoM 
    [velocityvcurr, VACvcurr, DOMvcurr]= getDoM(velvectorcurrent, rate, timespan);

    
    velocityvTotal = velocityvTotal+velocityvcurr;
    %%summing up all the velocities? 
    %why am I doing this? Why is this useful? Dont ask questions Adi just
    %do what Ted says
    
    VACvTotal = VACvTotal+VACvcurr;

    
    DoMvTotal = DoMvTotal+DOMvcurr;
end

%% Velocity vs Time (same as input)
tv5 = timekeeper(nosamples, 1/rate);
for i = 1:nodevices
    velvectorcurrent = txtvoltageunconcatenated(:,i);
    %voltvectorcurrent = voltvectorcurrent(1:1000);
    velvectorcurrent = velvectorcurrent - mean(velvectorcurrent);
    [velocityvcurr, VACvcurr, DOMvcurr]= getDoM(velvectorcurrent, rate, timespan);
    
    plot(tv5, velocityvcurr)
    hold on
end
hold off
%axis([0 3 -0.005 0.005])
title('Velocities vs Time')
ylabel('Velocity (au)')
xlabel('Time (s)')
figure()

%% Velocity Autocorr vs Lag Time (Equals half of timespan) (Avg. in black)
tv6 = timekeeper(nosamples+1, 1/rate);
tv6 = tv6 - timespan/2 ;
for i = 1:nodevices
    velvectorcurrent = txtvoltageunconcatenated(:,i);
    %voltvectorcurrent = voltvectorcurrent(1:1000);
    velvectorcurrent = velvectorcurrent - mean(velvectorcurrent);
    [velocityvcurr, VACvcurr, DOMvcurr]= getDoM(velvectorcurrent, rate, timespan);
    
    %plot(tv6, VACvcurr)
    hold on
end

plot(tv6, VACvTotal/nodevices, "black", "LineWidth", 1)

hold off
title('V.A.C. vs Delay Time')
ylabel('Velocity Autocorrelation')
xlabel('Time (s)')
xlim([-timespan/2, timespan/2])
figure()

%% plot the DoMs with total in black
fv1 = 0.5*rate*(0:(nosamples/2))/nosamples;

for i = 1:nodevices
    velvectorcurrent = txtvoltageunconcatenated(:,i);
    %voltvectorcurrent = voltvectorcurrent(1:1000);
    velvectorcurrent = velvectorcurrent - mean(velvectorcurrent);
    [velocityvcurr, VACvcurr, DOMvcurr]= getDoM(velvectorcurrent, rate, timespan);
    
    P2 = abs(DOMvcurr/nosamples);
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
    
loglog(fv1, P1/nodevices, "black")


hold off

title('Density of Modes v/s Frequency')
ylabel('Density of Modes')
xlabel('Frequency (Hz)')
xlim([0, 4000])
%ylim([1/1000000000000000000000000, 1/1000000000000])
figure()

%% Normalize the DoM by debye (avg. in black)
%my attempt to remove Debye scaling
%this section SHOULD NOT EXIST, so if there are too many bugs, just get rid
%of the whole damn thing and it will not be missed.

fv2 = 0.5*rate*(0:(nosamples/2))/nosamples; %check this pls
lenfv2 = length(fv2);
debye = zeros(lenfv2, 1); %normalize by this

for i = 1:nodevices
    velvectorcurrent = txtvoltageunconcatenated(:,i);
    %voltvectorcurrent = voltvectorcurrent(1:1000);
    velvectorcurrent = velvectorcurrent - mean(velvectorcurrent);
    [velocityvcurr, VACvcurr, DOMvcurr]= getDoM(velvectorcurrent, rate, timespan);
    
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
xlim([0, 1300])
%ylim([1/1000000000000000000000000, 1/1000000000000])
figure()

%% Helper functions
%Misc Helper functions-----------------------------------------------------

%the timekeeper vector creator- do not mess with the timekeeper

function tvector = timekeeper(n, timedifference)
tvector = zeros([n-1 1]);

    for i = 1:n-1 
        tvector(i+1)= tvector(i)+timedifference;
    end
tvector = tvector';
end


%the main function
%4 outputs to play with
function[velocityout, VACout, DoMout] = getDoM(velocityvector, ratescalar, durationscalar)

maxlensc = length(velocityvector);
voltmeansc = mean(velocityvector);

velocityout = velocityvector;
velocitymeansc = mean(velocityout);
velocityout = velocityout - velocitymeansc; %I need to input with mean
%already subtracted in this version

%Filter number 1: you need SIGNAL PROCESSING TOOLBOX for this guy
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

%now we can clean out the 60Hz noise? holy shit

%this is agony
%for i = 225:260
%    velocityfft(i) = velocityfft(i)/1000;
%end
%okay so this filter works but its pretty aggressive and blunt. its like a
%hammer that will bonk the relevant part of the signal down
%she is crazy and must go

velocityfftinverse = ifft(velocityfft);
VACout = velocityfftinverse;

%normalize the VAC
VACout = VACout/VACout(1);

%VACout = xcorr(velocityout, 0.5*ratescalar*durationscalar , 'normalized');
%old VAC method for the newbies (I never used this) (real)

pow2 = nextpow2(length(VACout)); %removed padding
FFTVAC = fft(VACout);

DoMout = real(FFTVAC);

%fuck- did I forget a factor of 2?
end