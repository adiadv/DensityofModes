%% disclaimer 
%Please read%
%The program is a mess- sorry I am as dissapointed as you are. 

%% README
%{

This file will calculate
2. Velocities
3. Velocity Autocorrelation
4. Density of Modes

For a given number of devices at a given frequency (rate) and a given time
(timespan) entered below. THE DEVICES MUST MEASURE VELOCITY.
%}

%% Update Me
%CALLING THE TXT%
fileID = fopen('voltageDoM0610');
txtvoltageunconcatenated = importdata('voltageDoM0610.txt'); %enter file here!
txtvoltage = sum(txtvoltageunconcatenated, 1);


%Initialize the sauce%

rate = 4000; %UPDATE THIS---
timespan = 8.00275; %UPDATE THIS---
nosamples = rate*timespan; %must be even because of an off-by-one thing
nodevices = 1; %UPDATE THIS---

closestpow2 = nextpow2(2*nosamples-1); %might not use this

%% Initialize
%forgive the horrible notation in the next 20 lines, it will never be used
%again.
initialiazervoltvector = txtvoltage; 
[velocityvinitialize, VACvinitialize, DOMvinitialize]= getDoM(initialiazervoltvector, rate, timespan);

geophonevelocitylength = length(velocityvinitialize);
velocityvTotal =zeros([nosamples 1]);

geophoneVAClength = length(VACvinitialize);
VACvTotal = zeros([nosamples+1 1]);

geophoneDOMlength = length(DOMvinitialize);
DoMvTotal = zeros([nosamples+1 1]); %removed padding

%% Total Velocity, VAC and DoM
%mix the secret geophone ingredients%
%we can take each column, do a calculation, then throw it add it to the
%total and throw it out to make things faster. 
for i = 1:nodevices
    velvectorcurrent = txtvoltageunconcatenated(:,i); 
    %velvectorcurrent = velvectorcurrent(1:1000); %cut out the bullshit at the start and end
    velvectorcurrent = velvectorcurrent - mean(velvectorcurrent); %subtract the mean manually 
    [velocityvcurr, VACvcurr, DOMvcurr]= getDoM(velvectorcurrent, rate, timespan);

    
    velocityvTotal = velocityvTotal+velocityvcurr;

    
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
tv6 = timekeeper(nosamples, 1/rate);
tv6 = tv6 - timespan/2;
for i = 1:nodevices
    velvectorcurrent = txtvoltageunconcatenated(:,i);
    %voltvectorcurrent = voltvectorcurrent(1:1000);
    velvectorcurrent = velvectorcurrent - mean(velvectorcurrent);
    [velocityvcurr, VACvcurr, DOMvcurr]= getDoM(velvectorcurrent, rate, timespan);
    
    plot(tv6, VACvcurr)
    hold on
end

plot(tv6, VACvTotal/nodevices, "black", "LineWidth", 1)

hold off
title('V.A.C. vs Delay Time')
ylabel('Velocity Autocorrelation')
xlabel('Time (s)')
%axis([0 5 -0.006 0.021])
figure()

%% plot the DoMs with total in black
fv1 = rate*(0:(nosamples/2))/nosamples;

for i = 1:nodevices
    velvectorcurrent = txtvoltageunconcatenated(:,i);
    %voltvectorcurrent = voltvectorcurrent(1:1000);
    velvectorcurrent = velvectorcurrent - mean(velvectorcurrent);
    [velocityvcurr, VACvcurr, DOMvcurr]= getDoM(velvectorcurrent, rate, timespan);
    
    P2 = abs(DOMvcurr/nosamples);
    P1 = P2(1:nosamples/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
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
xlim([0, 1300])
%ylim([1/1000000000000000000000000, 1/1000000000000])
figure()

%% Normalize the DoM by debye (avg. in black)
fv2 = rate*(0:(nosamples/2))/nosamples; %check this pls
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
%velocitymeansc = mean(velocityout);
%velocityout = velocityout - velocitymeansc; %I need to input with mean
%already subtracted in this version

VACout = xcorr(velocityout, 0.5*ratescalar*durationscalar, 'normalized');

pow2 = nextpow2(length(VACout)); %removed padding
FFTVAC = fft(VACout);

DoMout = real(FFTVAC);

end