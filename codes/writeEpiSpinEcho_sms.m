%change the original example with multi-slice capbility
close all; clear all;

% Define FOV and resolution
fov=256e-3; Nx=64; Ny=64; 
thickness=3e-3;
sliceGap=1.5e-3;             % slice gap im mm
Nslices=4;

% Set system limits
lims = mr.opts('MaxGrad', 32, 'GradUnit', 'mT/m',...
    'MaxSlew', 130, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 20e-6, ...
    'rfDeadTime', 100e-6, 'adcDeadTime', 20e-6);  
% Create a new sequence object
seq=mr.Sequence(lims);

% Create 90 degree slice selection pulse and gradient
timeBwProduct=4;
RFdur=3e-3;
[rf, gz] = mr.makeSincPulse(pi/2,lims,'Duration',RFdur,...
    'SliceThickness',thickness,'apodization',0.5,'timeBwProduct',timeBwProduct,...
    'use','excitation');

BW = timeBwProduct/(RFdur); % RF bandwidth, which is equal to RF timebandwidthproduct/RFdur.
SliceLeakRatio = 1; % meaning adjacent 2 slices have no gap and no overlapping.
MBfac = 2; % 1: single slice; 2: multi(2) slices.
MBFOVshift=2; %shift in FOV, 0 means no shift, 2 means 1/2 FOV shift, 3 means 1/3 FOV shift, ...

% >>>>>>>>>>> Phase modulation for multi-band. >>>>>>>>>>>>>>>>
if MBfac >= 2
        Nslices = Nslices/MBfac;
        mbFreqSept = SliceLeakRatio * BW * Nslices; % mbFreqSept = SliceLeakRatio * BW * Nslices;
        N = round(RFdur / lims.rfRasterTime); t = (1 : N) * lims.rfRasterTime; tt = t - (RFdur * 0.5);
%         MB1 = exp(1i * 2 * pi * (0.0) * mbFreqSept * tt);
%         MB2 = exp(1i * 2 * pi * (1.0) * mbFreqSept * tt);
%         rf.signal(1 : N) = squeeze(rf.signal(1 : N)) .* (MB1(:) + MB2(:));
        rf_signal0=rf.signal(1 : N);
        tmp=zeros(N,1);
        for mb_idx=1:MBfac
            MB = exp(1i * 2 * pi * (mb_idx-1) * mbFreqSept * tt);
            tmp = tmp+ rf_signal0(:) .* MB(:);
        end;
        rf.signal(1 : N)=tmp;
end



% Define other gradients and ADC events
deltak=1/fov;
kWidth = Nx*deltak;
%readoutTime = 3.2e-4;
ro_os=2.0;
readoutTime = round(lims.adcRasterTime*50*Nx*ro_os/1e-6)*1e-6;
gx = mr.makeTrapezoid('x',lims,'FlatArea',kWidth,'FlatTime',readoutTime);
adcSamples=Nx;
adc = mr.makeAdc(Nx,lims,'Duration',gx.flatTime,'Delay',gx.riseTime);

% Pre-phasing gradients
preTime=8e-4;
%gxPre = mr.makeTrapezoid('x',lims,'Area',-gx.area/2-deltak/2,'Duration',preTime);
gzReph = mr.makeTrapezoid('z',lims,'Area',-gz.area/2,'Duration',preTime);
%gyPre = mr.makeTrapezoid('y',lims,'Area',-Ny/2*deltak,'Duration',preTime);
% we need no minus for in-plane prephasers because of the spin-echo (position reflection in k-space)
gxPre = mr.makeTrapezoid('x',lims,'Area',gx.area/2-deltak/2,'Duration',preTime);
gyPre = mr.makeTrapezoid('y',lims,'Area',Ny/2*deltak,'Duration',preTime);

% Phase blip in shortest possible time
dur = ceil(2*sqrt(deltak/lims.maxSlew)/10e-6)*10e-6;
gy = mr.makeTrapezoid('y',lims,'Area',deltak,'Duration',dur);

% Refocusing pulse with spoiling gradients
% rf180 = mr.makeBlockPulse(pi,lims,'Duration',500e-6,'use','refocusing');
gzSpoil = mr.makeTrapezoid('z',lims,'Area',gz.area*2,'Duration',3*preTime);
[rf180, gz180] = mr.makeSincPulse(pi,lims,'Duration',6e-3,...
    'SliceThickness',thickness,'apodization',0.5,'timeBwProduct',4,'PhaseOffset',pi/2,'use','refocusing');

% Calculate delay time %% MZ: I thisk this is very wrong!
TE=80e-3;
durationToCenter = (Nx/2+0.5)*mr.calcDuration(gx) + Ny/2*mr.calcDuration(gy);
rfCenterInclDelay=rf.delay + mr.calcRfCenter(rf);
rf180centerInclDelay=rf180.delay + mr.calcRfCenter(rf180);
delayTE1=TE/2 - mr.calcDuration(gz) + rfCenterInclDelay - preTime -preTime - mr.calcDuration(gzSpoil) - rf180centerInclDelay;
delayTE2=TE/2 - mr.calcDuration(rf180) + rf180centerInclDelay - mr.calcDuration(gzSpoil) - durationToCenter;


if(MBFOVshift)
    %CAIPI
    deltak_sms=1/(Nslices*thickness)/MBFOVshift;
    deltak_sms_max=1/(Nslices*thickness)/MBFOVshift*(MBFOVshift-1);
    dur_sms = ceil(2*sqrt(deltak_sms_max/lims.maxSlew)/10e-6)*10e-6;
    if(dur_sms>dur)
        gy = mr.makeTrapezoid('y',lims,'Area',deltak,'Duration',dur_sms);
        gy_parts = mr.splitGradientAt(gy, dur_sms/2, lims);
        [gy_blipup, gy_blipdown,~]=mr.align('right',gy_parts(1),'left',gy_parts(2),gx);
        gy_blipdownup=mr.addGradients({gy_blipdown, gy_blipup}, lims);

        gz_sms = mr.makeTrapezoid('z',lims,'Area',deltak_sms,'Duration',dur_sms);
        gz_sms_rep = mr.makeTrapezoid('z',lims,'Area',-deltak_sms_max,'Duration',dur_sms);

%         gz_sms_parts = mr.splitGradientAt(gz_sms, dur_sms/2, lims);
%         [gz_sms_blipup, gz_sms_blipdown,~]=mr.align('right',gz_sms_parts(1),'left',gz_sms_parts(2),gx);
% 
%         gz_sms_rep_parts = mr.splitGradientAt(gz_sms_rep, dur_sms/2, lims);
%         [gz_sms_rep_blipup, gz_sms_rep_blipdown,~]=mr.align('right',gz_sms_rep_parts(1),'left',gz_sms_rep_parts(2),gx);
% 
%         gz_sms_blipdownup=mr.addGradients({gz_sms_blipdown, gz_sms_blipup}, lims);    
%         gz_sms_blipdownuprep=mr.addGradients({gz_sms_blipdown, gz_sms_rep_blipup}, lims);    
%         gz_sms_blipdownrepup=mr.addGradients({gz_sms_rep_blipdown, gz_sms_blipup}, lims); 
    else
        gz_sms = mr.makeTrapezoid('z',lims,'Area',deltak_sms,'Duration',dur);
        gz_sms_rep = mr.makeTrapezoid('z',lims,'Area',-deltak_sms_max,'Duration',dur);

%         gz_sms_parts = mr.splitGradientAt(gz_sms, dur/2, lims);
%         [gz_sms_blipup, gz_sms_blipdown,~]=mr.align('right',gz_sms_parts(1),'left',gz_sms_parts(2),gx);
% 
%         gz_sms_rep_parts = mr.splitGradientAt(gz_sms_rep, dur/2, lims);
%         [gz_sms_rep_blipup, gz_sms_rep_blipdown,~]=mr.align('right',gz_sms_rep_parts(1),'left',gz_sms_rep_parts(2),gx);
% 
%         gz_sms_blipdownup=mr.addGradients({gz_sms_blipdown, gz_sms_blipup}, lims);    
%         gz_sms_blipdownuprep=mr.addGradients({gz_sms_blipdown, gz_sms_rep_blipup}, lims); 
%         gz_sms_blipdownrepup=mr.addGradients({gz_sms_rep_blipdown, gz_sms_blipup}, lims); 
    end;
end;


slicePositions=(thickness+sliceGap)*((0:(Nslices-1)) - (Nslices-1)/2);

for s=1:Nslices
    % Define sequence blocks
    rf.freqOffset=gz.amplitude*slicePositions(s);
    rf.phaseOffset=-2*pi*rf.freqOffset*mr.calcRfCenter(rf); % compensate for the slice-offset induced phase
    seq.addBlock(rf,gz);
    seq.addBlock(gzReph);

    seq.addBlock(gxPre,gyPre);
    
    seq.addBlock(mr.makeDelay(delayTE1));
    seq.addBlock(gzSpoil);
    
    %seq.addBlock(rf180);
    rf180.freqOffset=gz.amplitude*slicePositions(s);
    rf180.phaseOffset=-2*pi*rf180.freqOffset*mr.calcRfCenter(rf180); % compensate for the slice-offset induced phase    
    seq.addBlock(rf180,gz180);
    seq.addBlock(gzSpoil);

    seq.addBlock(mr.makeDelay(delayTE2));
    for i=1:Ny
        seq.addBlock(gx,adc);           % Read one line of k-space
        if MBFOVshift == 0
            seq.addBlock(gy);               % Phase blip
        else
            if(mod(i,MBFOVshift)~=0)
               seq.addBlock(gy,gz_sms);               % Phase blip
            else
               seq.addBlock(gy,gz_sms_rep);               % Phase blip
            end;
        end;
        gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
    end

    
    seq.addBlock(mr.makeDelay(1e-4));
end;
%% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% export and visualization
seq.setDefinition('FOV', [fov fov thickness]);
seq.setDefinition('Name', 'epise');
seq.write('epi_se_sms.seq');   % Output sequence for scanner
seq.plot();             % Plot sequence waveforms

%% calculate trajectory 
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();
%[ktraj_adc, ktraj, t_excitation, t_refocusing, t_adc] = seq.calculateKspace();

%% plot k-spaces
figure; plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
hold; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis

figure; plot(ktraj(1,:),ktraj(2,:),'b',...
             ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display

%% sanity checks
TE_check=(t_refocusing(1)-t_excitation(1))*2;
fprintf('intended TE=%.03f ms, actual spin echo TE=%.03fms\n', TE*1e3, TE_check*1e3); 

seq.paperPlot('blockRange',[1 81]);

