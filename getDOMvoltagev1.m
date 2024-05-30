%% disclaimer
%Please read%
%this script is messy but commented. :( also it might be broken 

%% README
%{

This file will calculate
1. Voltages
2. Velocities
3. Velocity Autocorrelation
4. Density of Modes

For a given number of devices at a given frequency (rate) and a given time
(timespan) entered below. THE DEVICES MUST MEASURE VOLTAGE.
%}

%% Update Me
%CALLING THE TXT%
fileID = fopen('SB_21_new.txt');
txtvoltageunconcatenated = importdata('SB_21_new.txt');
txtvoltage = sum(txtvoltageunconcatenated, 2); %summing up all the rows

%Initialize the sauce%


rate = 4000; %UPDATE THIS---
timespan = 1; %UPDATE THIS---
nosamples = rate*timespan;
nodevices = 24; %UPDATE THIS---

closestpow2 = nextpow2(2*nosamples-1);%not used?

%% initializer section 
%forgive the horrible notation in the next 20 lines, it will never be used
%again.
initialiazervoltvector = txtvoltage; 
[voltagevinitialize, velocityvinitialize, VACvinitialize, DOMvinitialize]= getDoM(initialiazervoltvector);

geophonevoltagelength = length(voltagevinitialize); %how many voltages (these are never used?)
voltagevTotal = zeros([nosamples 1]); %initialize the voltage vector (this is used)

geophonevelocitylength = length(velocityvinitialize); %how many velocities
velocityvTotal =zeros([nosamples 1]); %initialize the velocity vector

geophoneVAClength = length(VACvinitialize); %how many VACs
VACvTotal = zeros([2*nosamples-1 1]); %initialize the VAC vector

geophoneDOMlength = length(DOMvinitialize);%how many DoMs
DoMvTotal = zeros([2*nosamples-1 1]); %initialize the DoM vector

%% finding total voltage, velocity, VAC and DOM
%mix the secret geophone ingredients%
%we can take each column, do a calculation, then throw it add it to the
%total and throw it out to make things faster. 
for i = 1:nodevices
    voltvectorcurrent = txtvoltageunconcatenated(:,i); %take the ith column
    %voltvectorcurrent = voltvectorcurrent(1:1000); %cut out the bullshit at the start and end
    voltvectorcurrent = voltvectorcurrent - mean(voltvectorcurrent); %subtract the mean manually 
    %no clue why this works but the other way doesnt theyre the same. 
    %Dont question it..
    [voltagevcurr, velocityvcurr, VACvcurr, DOMvcurr]= getDoM(voltvectorcurrent);
    
    
    voltagevTotal= voltagevTotal+voltagevcurr; %summing up all the voltages? 
    %why am I doing this? Why is this useful? 
    
    velocityvTotal = velocityvTotal+velocityvcurr; %i dont think this is used

    
    VACvTotal = VACvTotal+VACvcurr;

    
    DoMvTotal = DoMvTotal+DOMvcurr;
end

%The main getDOM function below ensures <V> and <v> are centered. 

%% plot all voltages 
%Ted's work- starting here-------------------------------------------------
%Ted1- V vs t, 24 colors on the same graph (i cant even name 24 colors)

%basically, i dont know how to create new variables using a forloop so i
%might bite the bullet and do it manually (pain is beginning)
tv4 = timekeeper(nosamples, 1/rate);
for i = 1:nodevices
    voltvectorcurrent = txtvoltageunconcatenated(:,i);
    %voltvectorcurrent = voltvectorcurrent(1:1000); %can be used to cut the
    %bad data at the start and end of the voltage file. (will change length
    %so you need to change tv4 then of course)
    voltvectorcurrent = voltvectorcurrent - mean(voltvectorcurrent);
    [voltagevcurr, velocityvcurr, VACvcurr, DOMvcurr]= getDoM(voltvectorcurrent);
    
    plot(tv4, voltagevcurr)
    hold on
end
hold off
%axis([0 0.2 -0.00001 0.00001])
title('24 Voltages vs Time') 
ylabel('Voltage (V)')
xlabel('Time (s)')
figure()

%% plot all velocities
%Ted2- v vs t but the same as above (i'll just copy paste the stuff above)
tv5 = timekeeper(nosamples, 1/rate);
for i = 1:nodevices
    voltvectorcurrent = txtvoltageunconcatenated(:,i);
    %voltvectorcurrent = voltvectorcurrent(1:1000);
    voltvectorcurrent = voltvectorcurrent - mean(voltvectorcurrent);
    [voltagevcurr, velocityvcurr, VACvcurr, DOMvcurr]= getDoM(voltvectorcurrent);
    
    plot(tv5, velocityvcurr)
    hold on
end
hold off
%axis([0 3 -0.005 0.005])
title('24 Velocities vs Time')
ylabel('Velocity (m/s)')
xlabel('Time (s)')
figure()

%% plot all VACs with total in black
%Ted3- VAC vs t, but the same again...
tv6 = timekeeper(2*nosamples-1, 1/rate);
tv6 = tv6 - timespan; %wait why is this happening? 
for i = 1:nodevices
    voltvectorcurrent = txtvoltageunconcatenated(:,i);
    %voltvectorcurrent = voltvectorcurrent(1:1000);
    voltvectorcurrent = voltvectorcurrent - mean(voltvectorcurrent);
    [voltagevcurr, velocityvcurr, VACvcurr, DOMvcurr]= getDoM(voltvectorcurrent);
    
    plot(tv6, VACvcurr)
    hold on
end

plot(tv6, VACvTotal/nodevices, "black", "LineWidth", 1) %plotting the average in black

hold off
title('V.A.C. vs Delay Time')
ylabel('Velocity Autocorrelation (average in black)')
xlabel('Time (s)')
%axis([0 5 -0.006 0.021])
figure()

%% plot all DoMs with total in black
%Ted4- You guessed it. getDOMx24
fv1 = rate*(0:(nosamples/2))/nosamples; %boy i hope this works!

for i = 1:nodevices
    voltvectorcurrent = txtvoltageunconcatenated(:,i);
    %voltvectorcurrent = voltvectorcurrent(1:1000);
    voltvectorcurrent = voltvectorcurrent - mean(voltvectorcurrent);
    [voltagevcurr, velocityvcurr, VACvcurr, DOMvcurr]= getDoM(voltvectorcurrent);
    
    P2 = abs(DOMvcurr/nosamples); %i stole this from matlab documentation on ffts
    P1 = P2(1:nosamples/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
    semilogx(fv1, P1)
    hold on
end


P2 = abs(DoMvTotal/nosamples);
P1 = P2(1:nosamples/2+1);
P1(2:end-1) = 2*P1(2:end-1);
    
semilogx(fv1, P1/nodevices, "black")


hold off

title('DOM plotted against frequency (average in freq)')
ylabel('Density of Modes')
xlabel('Frequency (Hz)')
%axis([0 60 -0.0005 0.0009])
figure()

%% Voltage Semilogy for checking noise 
%Ted 5- psd using the voltages
%take the fft of the voltage vector total, then conjugate, multiply and
%plot with timekeeper


xdft = fft(voltagevTotal);
xdft = xdft(1:nosamples/2+1); %lets take a good length
psdx = (1/(rate*nosamples)) * abs(xdft).^2; %now
psdx(2:end-1) = 2*psdx(2:end-1);
fv2 = 0:rate/length(voltagevTotal):rate/2;


semilogy(fv2,psdx, ".")

%set(gca,'yscale','log')

title('PSD vs f') %do I understand this?
ylabel('Power Spectral Density')
axis([0 2000 1/10000000000000000 1/10000000000])
figure()

%% Checker plots
%checker plotting%---------------------------------------------------------
%uncomment the plot you want to check
%be sure to check timekeeper lengths
%I think this section is broken or outdated, but that doesn't really matter

%voltagevTotal
%{
%voltagevTotal--------------
tv = timekeeper(nosamples, 1/rate);
plot(tv, voltagevTotal, ".")
axis([0 3 -0.00001 0.00001])
title('Voltage vs Time')
ylabel('Voltage (V)')
figure()
%}

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
tv2 = timekeeper(2*nosamples-1, 1/rate);
plot(tv2, VACvTotal, ".")
title('V.A.C. vs Delay Time')
ylabel('Velocity Autocorrelation')
axis([0 5 -0.006 0.021])
figure()
%}

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
function[voltout, velocityout, VACout, DoMout] = getDoM(voltvector)

maxlensc = length(voltvector);

voltmeansc = mean(voltvector);
%voltout = voltvector - voltmeansc; %why does this not work? do not use i
%guess
voltout = voltvector;

velocityout = cumsum(voltout); %the integral is delivered
%velocitymeansc = mean(velocityout);
%velocityout = velocityout - velocitymeansc; %this works but the other one doesnt?

%as of now, we are subtracting the mean manually

VACout = xcorr2(velocityout);

pow2 = nextpow2(length(VACout));
FFTVAC = fft(VACout); %removed padding

DoMout = real(FFTVAC);

end