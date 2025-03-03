function [IE,ph_alt] = hw_fun_sim_offset(v,ORFreq,Gmax,Gmean,FA,RF_dur,RF_sep,isUnipolar)
G_amp = Gmax; meanGz = Gmean; FlipAngleDeg = FA;

g = GetGamma;  % rad/s/mT            
zmax = 3.0 / 100; % cm to m conversion (maximum spin offset from the labelling plane to simulate)

Pa = [-0.02 0]'; Pb = [0 0]'; % Position of vessels "A" and "B" within the labelling plane in m
Ps = -0.05*2:0.001:0.05*2; Ps(2,:) = 0; 
cycle = 3; % In this cycle we label at the position corresponding to vessel A and control at vessel B
z_offset = 0; % Position the labelling at isocentre for simplicity

% Define T1 and T2
T1 = 1.932; T2 = 0.275;  % @ 3T, ref: Stanisz, MRM 2005

RF_amp = 0.04 * 1e-4 * 1000; % Gauss to mT conversion
                             % Wong = 0.04 G
RF_shape = 'hanning';  
RF_shape_params = []; % No extra parameters needed for Hanning pulses
dt = 10e-6; pulse_samples = round(RF_dur/dt); 
RF_pulse = Hann_window(pulse_samples)*RF_amp; RF_pulse(2,:) = 0;
Act_FlipAngle_degs = FlipAngle(RF_pulse,dt);
RF_amp = RF_amp * FA / Act_FlipAngle_degs;

dt = RF_dur/200; % /20
FinalMz = zeros(size(ORFreq,2),size(Ps,2)); IE = FinalMz; jj = 1;

for ii = 1:size(ORFreq,2)
    parfor kk = 1:size(Ps,2)
        disp(['Velocity step ' num2str(jj) ', Off-reson step ' num2str(ii)]);
        % Run the simulation
        [M, t, P, G, RF] = test_VEPCASL_seq_meanGz(v, zmax, z_offset, meanGz, G_amp, RF_shape, ...
                                                    RF_shape_params, RF_amp, RF_dur, ...
                                                    RF_sep, Pa, Pb, cycle, ...
                                                    Ps(:,kk), dt, T1,T2, ORFreq(ii),isUnipolar);

        % Record the final z magnetisation
        FinalMz(ii,kk) = M(3,end);
        IE(ii,kk) = InvEff(FinalMz(ii,kk),2*zmax/v(jj), zmax/v(jj), T1);
    end
end
ph_alt = (Ps(1,:) - Pa(1)) / (Pb(1) - Pa(1)) * pi;


end