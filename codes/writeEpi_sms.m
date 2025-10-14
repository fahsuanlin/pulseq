% this is a demo low-performance EPI sequence;
% it doesn't use ramp-samping and is only good for educational purposes.
%
close all; clear all;

seq=mr.Sequence();              % Create a new sequence object
fov=220e-3; Nx=64; Ny=64;       % Define FOV and resolution
thickness=3e-3;                 % slice thinckness
Nslices=10;

% Set system limits
lims = mr.opts('MaxGrad',32,'GradUnit','mT/m',...
               'MaxSlew',130,'SlewUnit','T/m/s', ...
               'rfRingdownTime', 30e-6, 'rfDeadTime', 100e-6);


% Create 90 degree slice selection pulse and gradient
[rf, gz] = mr.makeSincPulse(pi/2,'system',lims,'Duration',3e-3,...
    'SliceThickness',thickness,'apodization',0.5,'timeBwProduct',4,...
    'use', 'excitation');


BW = 4/(3e-3); % RF bandwidth, which is equal to RF timebandwidthproduct/RFdur.
RFdur = 3e-3;
SliceLeakRatio = 1; % meaning adjacent 2 slices have no gap and no overlapping.
MBfac = 2; % 1: single slice; 2: multi(2) slices.
Nslices = Nslices/MBfac;
% >>>>>>>>>>> Phase modulation for multi-band. >>>>>>>>>>>>>>>>
if MBfac == 2
%    if MBhalfFOVshift == 1
%         mbFreqSept = SliceLeakRatio * gz.amplitude * thickness * Nslices; % mbFreqSept = SliceLeakRatio * BW * Nslices;
%         % N = round(RFdur/ sys.rfRasterTime); t = (1 : N) * sys.rfRasterTime; tt = t - (mr.calcRfCenter(rf) * 0.5); % Something wrong here...!!!
%         N = round(RFdur/ sys.rfRasterTime); t = (1 : N) * sys.rfRasterTime; tt = t - (RFdur * 0.5); % Something wrong here...!!!
%         MB1 = exp(1i * 2 * pi * (-0.5) * mbFreqSept * tt);
%         MB2 = exp(1i * 2 * pi * (+0.5) * mbFreqSept * tt);
%         rfsignal_neg = rf.signal;
%         rfsignal_neg(1 : N) = squeeze(rf.signal(1 : N)) .* (-MB1(:) + MB2(:));
%         rfsignal_pos = rf.signal;
%         rfsignal_pos(1 : N) = squeeze(rf.signal(1 : N)) .* ( MB1(:) + MB2(:));
%    else
        mbFreqSept = SliceLeakRatio * BW * Nslices; % mbFreqSept = SliceLeakRatio * BW * Nslices;
        N = round(RFdur / lims.rfRasterTime); t = (1 : N) * lims.rfRasterTime; tt = t - (RFdur * 0.5);
        MB1 = exp(1i * 2 * pi * (0.0) * mbFreqSept * tt);
        MB2 = exp(1i * 2 * pi * (1.0) * mbFreqSept * tt);
        rf.signal(1 : N) = squeeze(rf.signal(1 : N)) .* (MB1(:) + MB2(:));
%    end
end
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


% Define other gradients and ADC events
deltak=1/fov;
kWidth = Nx*deltak;
dwellTime = 4e-6; % I want it to be divisible by 2
readoutTime = Nx*dwellTime;
flatTime=ceil(readoutTime*1e5)*1e-5; % round-up to the gradient raster
gx = mr.makeTrapezoid('x',lims,'Amplitude',kWidth/readoutTime,'FlatTime',flatTime);
adc = mr.makeAdc(Nx,'Duration',readoutTime,'Delay',gx.riseTime+flatTime/2-(readoutTime-dwellTime)/2);

% Pre-phasing gradients
preTime=8e-4;
gxPre = mr.makeTrapezoid('x',lims,'Area',-gx.area/2,'Duration',preTime); % removed -deltak/2 to aligh the echo between the samples
gzReph = mr.makeTrapezoid('z',lims,'Area',-gz.area/2,'Duration',preTime);
gyPre = mr.makeTrapezoid('y',lims,'Area',-Ny/2*deltak,'Duration',preTime);

% Phase blip in shortest possible time
dur = ceil(2*sqrt(deltak/lims.maxSlew)/10e-6)*10e-6;
gy = mr.makeTrapezoid('y',lims,'Area',deltak,'Duration',dur);

MBhalfFOVshift=1;
if(MBhalfFOVshift)
    %CAIPI
    caipi_fov_shift=2; %shift in FOV, 2 means 1/2 FOV shift, 3 means 1/3 FOV shift, ...
    deltak_sms=1/(Nslices*thickness)/caipi_fov_shift;
    deltak_sms_max=1/(Nslices*thickness)/caipi_fov_shift*(caipi_fov_shift-1);
    dur_sms = ceil(2*sqrt(deltak_sms_max/lims.maxSlew)/10e-6)*10e-6;
    if(dur_sms>dur)
        gy = mr.makeTrapezoid('y',lims,'Area',deltak,'Duration',dur_sms);
        gz_sms = mr.makeTrapezoid('z',lims,'Area',deltak_sms,'Duration',dur_sms);
        gz_sms_rep = mr.makeTrapezoid('z',lims,'Area',-deltak_sms_max,'Duration',dur_sms);
    else
        gz_sms = mr.makeTrapezoid('z',lims,'Area',deltak_sms,'Duration',dur);
        gz_sms_rep = mr.makeTrapezoid('z',lims,'Area',-deltak_sms_max,'Duration',dur);
    end;
end;

% Define sequence blocks
% seq.addBlock(mr.makeDelay(1)); % older scanners like Trio may need this
                                 % dummy delay to keep up with timing
for s=1:Nslices
    rf.freqOffset=gz.amplitude*thickness*(s-1-(Nslices-1)/2);
    seq.addBlock(rf,gz);
    seq.addBlock(gxPre,gyPre,gzReph);

    caipi_counter=1;
    for i=1:Ny
        seq.addBlock(gx,adc);           % Read one line of k-space
        if MBhalfFOVshift == 0
            seq.addBlock(gy);               % Phase blip
        else
            if(mod(caipi_counter,caipi_fov_shift)==0)
                seq.addBlock(gy,gz_sms_rep);           % Phase blip
            else
                seq.addBlock(gy,gz_sms);               % Phase blip
            end;
            caipi_counter=caipi_counter+1;
        end;
        gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
    end
    if s==1
        TR_1slice=seq.duration; % note the actual TR per slice
    end
end

%% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% Plot sequence waveforms
%seq.plot();

seq.plot('stacked',1,'timeRange',[0 TR_1slice], 'timeDisp','ms'); % niceer plot for the 1st sclice

%% trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

% plot k-spaces
figure; plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
hold; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold; plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.');

seq.write('epi_sms_caipi3.seq');   % Output sequence for scanner
% seq.sound(); % simulate the seq's tone
