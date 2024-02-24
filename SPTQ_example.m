% Small example for the use of the SPTQ method to determine the 23Na TQ
% signal from the FID only.
% Here, we used a simulated FID with relaxation times of a 6% agar + 154Mm
% NaCl sample. 
% see [Paper doi] for detailed desctiption of the SP method


% % % Can be replaced by experimental FID
% 23Na T2 relaxation times for 6% agarose + 154mM NaCl solution at 9.4T
T2s = 32.0e-3; % ms
T2f =  3.7e-3; % ms

ts = 0:1e-4:400e-3; % sampling times
FID = 0.4*exp(-ts./T2s) + 0.6*exp(-ts./T2f); % simulated FID at sampling times
% % % 

% calculating the SPTQ signal
[intFID, fid, TQ, tevos] = getSPTQ(ts, FID);


% Plot intFID, fid and TQ signals
figure(); hold on;
plot(tevos, intFID*100, 'linewidth', 2);
plot(tevos, fid*100, 'linewidth', 2);
plot(tevos, TQ*100, 'linewidth', 2);
pbaspect([1 1 1]);
ylabel('Normalized signal [%]')
xlabel("Evolution time [ms]");
legend("SQ(\tau_{evo})+TQ(\tau_{evo}) = \int FID(t,\tau_{evo})dt", "SQ(\tau_{evo})=FID(0,\tau_{evo})", "TQ(\tau_{evo}) = \int FID(t,\tau_{evo})dt-FID(0,\tau_{evo})", 'Location','northeast')
title("Signal of SPTQ method");

set(gca,'FontSize',12);
set(gca,'linewidth',1.5);
grid on;
box on;
% set(gca,'Fontweight','bold');



% Determine the 23Na TQ signal using only the FID
% Input:    ts: sampling times
%          FID: FID sampled at ts
% Output: intFID: integral over FID for all evolution times
%            fid: FID(0:NumTevoPoints-1)
%             TQ: TQ signal
function [intFID, fid, TQ, tevos] = getSPTQ(ts, FID)
    FID = FID/real(FID(1)); % normalize FID to avoid rounding errors
    NumTevoPoints = floor(length(FID)/2); 
    intFID  = zeros([1 NumTevoPoints]);
    fid = zeros([1 NumTevoPoints]);
    tevos = ts(1:NumTevoPoints); 

    for i = 1:1:NumTevoPoints
        fid_i  = (FID(i:end-(NumTevoPoints-i) )); 
        
        intFID_i  = sum(fid_i); % integral of fid
        fid0_i = fid_i(1); % first point of fid ~ int(spec)
        fid(i) = fid0_i; 
        intFID(i)  = intFID_i;     
    end

    intFID  =  intFID/intFID(1) ; % normalize intFID
    fid = fid/fid(1);   % normalize fid

    TQ = intFID-fid; % TQ signal is given by normalized difference
end

