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
fileID = fopen('testnakul');
txtvoltageunconcatenated = importdata('testnakul.txt'); %enter file here!
txtvoltage = sum(txtvoltageunconcatenated, 2);


%Initialize the sauce%

rate = 2500; %UPDATE THIS---
timespan = 1.5996; %UPDATE THIS---
nosamples = rate*timespan;
nodevices = 1; %UPDATE THIS---

closestpow2 = nextpow2(2*nosamples-1); %might not use this

%% Initializer mess
%forgive the horrible notation in the next 20 lines, it will never be used
%again.
initialiazervoltvector = txtvoltage; 
[velocityvinitialize, VACvinitialize, DOMvinitialize]= getDoM(initialiazervoltvector);

geophonevelocitylength = length(velocityvinitialize);
velocityvTotal =zeros([nosamples 1]);

geophoneVAClength = length(VACvinitialize);
VACvTotal = zeros([2*nosamples-1 1]);

geophoneDOMlength = length(DOMvinitialize);
DoMvTotal = zeros([2*nosamples-1 1]); %removed padding

%% Total velocity, VAC and DOM
%mix the secret geophone ingredients%
%we can take each column, do a calculation, then throw it add it to the
%total and throw it out to make things faster. 
for i = 1:nodevices
    velvectorcurrent = txtvoltageunconcatenated(:,i); 
    %velvectorcurrent = velvectorcurrent(1:1000); %cut out the bullshit at the start and end
    velvectorcurrent = velvectorcurrent - mean(velvectorcurrent); %subtract the mean manually 
    %no clue why this works but the other way doesnt theyre the same. 
    %Dont question it..
    [velocityvcurr, VACvcurr, DOMvcurr]= getDoM(velvectorcurrent);

    
    velocityvTotal = velocityvTotal+velocityvcurr;

    
    VACvTotal = VACvTotal+VACvcurr;

    
    DoMvTotal = DoMvTotal+DOMvcurr;
end

%% plot the velocities
%Ted's work- starting here-------------------------------------------------

%Ted2- v vs t but the same as above (i'll just copy paste the stuff above)
tv5 = timekeeper(nosamples, 1/rate);
for i = 1:nodevices
    velvectorcurrent = txtvoltageunconcatenated(:,i);
    %voltvectorcurrent = voltvectorcurrent(1:1000);
    velvectorcurrent = velvectorcurrent - mean(velvectorcurrent);
    [velocityvcurr, VACvcurr, DOMvcurr]= getDoM(velvectorcurrent);
    
    plot(tv5, velocityvcurr)
    hold on
end
hold off
%axis([0 3 -0.005 0.005])
title('Velocities vs Time')
ylabel('Velocity (m/s)')
xlabel('Time (s)')
figure()

%% plot the VACs with total in black
%Ted3- VAC vs t, but the same again... How original...
tv6 = timekeeper(2*nosamples-1, 1/rate);
tv6 = tv6 - timespan;
for i = 1:nodevices
    velvectorcurrent = txtvoltageunconcatenated(:,i);
    %voltvectorcurrent = voltvectorcurrent(1:1000);
    velvectorcurrent = velvectorcurrent - mean(velvectorcurrent);
    [velocityvcurr, VACvcurr, DOMvcurr]= getDoM(velvectorcurrent);
    
    plot(tv6, VACvcurr)
    hold on
end

plot(tv6, VACvTotal/nodevices, "black", "LineWidth", 1)

hold off
title('V.A.C. vs Delay Time')
ylabel('Velocity Autocorrelation (average in black)')
xlabel('Time (s)')
%axis([0 5 -0.006 0.021])
figure()

%% plot the DoMs with total in black
%Ted4- Plotting the DoM with total in black
fv1 = 0.5*rate*(0:(nosamples/2))/nosamples;

for i = 1:nodevices
    velvectorcurrent = txtvoltageunconcatenated(:,i);
    %voltvectorcurrent = voltvectorcurrent(1:1000);
    velvectorcurrent = velvectorcurrent - mean(velvectorcurrent);
    [velocityvcurr, VACvcurr, DOMvcurr]= getDoM(velvectorcurrent);
    
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

title('DOM (fft)')
ylabel('Density of Modes')
xlabel('Frequency (Hz)')
xlim([0, 1300])
ylim([1/1000000000000000000000000, 1/1000000000000])
figure()

%% Checker plotters
%checker plotting%---------------------------------------------------------
%uncomment the plot you want to check
%be sure to check timekeeper lengths


%velocityvTotal
%{
%velocityvTotal----------
plot(tv, velocityvTotal, ".")
axis([0 3 -0.005 0.005]) %axis scales
title('Velocity vs Time')
ylabel('Velocity (m/s)')
figure()
%}

%VACvTotal
%{
%VACvTotal---------------
tv2 = timekeeper(15999, 1/rate);
plot(tv2, VACvTotal, ".")
title('V.A.C. vs Delay Time')
ylabel('Velocity Autocorrelation')
axis([0 5 -0.006 0.021])
figure()
%}

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
function[velocityout, VACout, DoMout] = getDoM(velocityvector)

maxlensc = length(velocityvector);
voltmeansc = mean(velocityvector);

velocityout = velocityvector;
%velocitymeansc = mean(velocityout);
%velocityout = velocityout - velocitymeansc; %this works but the other one doesnt?

VACout = xcorr(velocityout);

pow2 = nextpow2(length(VACout)); %removed padding
FFTVAC = fft(VACout);

DoMout = real(FFTVAC);

end