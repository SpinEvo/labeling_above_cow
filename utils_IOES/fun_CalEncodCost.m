function [cost,condition,min_half_lambda,ang_locs_tagmode,hadamard_full,e] = fun_CalEncodCost(vesloc,modmat_phase,modmat_uni_lookup,offset)

hadamard_full = fun_HW_ChooseHadamard(vesloc); % random
hadamard_full(end+1,:) = -hadamard_full(1,:); % tag all
hadamard_full(end+1,:) = hadamard_full(1,:); % control all
hadamard_full(1,:) = [];

motion = 4; 
pad_size = 1024;
matrix_dim = pad_size;
units = (2*pi/matrix_dim);                                      %units for converting coordinates in to spatial frequencies
k_centre = [(matrix_dim/2)+1 (matrix_dim/2)+1];                 %coordinates of the centre of kspace
delta_k_padded = 1/(matrix_dim);                                %delta k in 1/mm of the zero padded matrix
lamda_limit = motion*4;                                         %minimum possible wavelength [mm]
max_mask_radius = 1/(delta_k_padded*lamda_limit);               %maximum radius of the mask, depends on the minimum possible wavelength [#pixels = 1/((mm^-1)*mm)]
half_FOV_padded = (matrix_dim)/2;
xy = -half_FOV_padded:(2*half_FOV_padded)/matrix_dim:half_FOV_padded; %make a vector (in mm) to fill the axes of later figures
half_kFOV_padded = (matrix_dim)/2;
kxy = -(half_kFOV_padded*delta_k_padded):delta_k_padded:(half_kFOV_padded*delta_k_padded);

position_weighting = zeros(matrix_dim,matrix_dim);
mask = zeros(matrix_dim,matrix_dim);
for g = 1:matrix_dim
    for f = 1:matrix_dim
        position_weighting(g,f) = sqrt((g-k_centre(1,1))^2 + (f-k_centre(1,2))^2);
        mask(g,f) = position_weighting(g,f);
        if mask(g,f) > max_mask_radius
            mask(g,f) = 0;
        elseif mask(g,f) > 100
            mask(g,f) = 0;
        end
    end
end

for g = 1:matrix_dim
    for f = 1:matrix_dim
        if mask(g,f) > 0
            mask(g,f) = 1;  
        end
    end
end
p = (-position_weighting/max(max(position_weighting)) + ones(matrix_dim,matrix_dim))/2;
[pixves] = convert_vessels_mm2pix(vesloc,pad_size);
h_phi = (hadamard_full+1)*(pi/2);
h_phi_offres = zeros(size(h_phi));
%remove the off resonance phase that will accrue between pulses due to the field inhomogeneity
for n = 1:size(h_phi,1)
    h_phi_offres(n,:) = h_phi(n,:) - offset;  % no offset
end
%Turn the encodings in to complex numbers (e^itheta)
h_offres = exp(1i*h_phi_offres);
hadamard_wo_ctrl_and_tissue = h_offres;
%Pre allocate variable names to speed up processing time in matlab
r = zeros(1,size(hadamard_wo_ctrl_and_tissue,1));
c = zeros(1,size(hadamard_wo_ctrl_and_tissue,1));
frequency = zeros(size(hadamard_wo_ctrl_and_tissue,1),2);
phase = zeros(1,size(hadamard_wo_ctrl_and_tissue,1));
encoding = zeros(size(hadamard_wo_ctrl_and_tissue,1),size(vesloc,2));
pp = zeros(size(hadamard_wo_ctrl_and_tissue,1),size(vesloc,2));
for m = 1:size(hadamard_wo_ctrl_and_tissue,1)
    Mv = zeros(matrix_dim,matrix_dim);
    for n = 1:size(vesloc,2)
        Mv(pixves(2,n),pixves(1,n)) = hadamard_wo_ctrl_and_tissue(m,n);
    end
    k_vessel_centre=fftshift(fft2(fftshift(Mv)));
    k = abs(k_vessel_centre)/max(max(abs(k_vessel_centre)));
    weighted_kspace_complex = p + k; weighted_kspace_complex = weighted_kspace_complex.*mask;
    [r(m),c(m)] = find(weighted_kspace_complex == max(max(weighted_kspace_complex)),1);
    frequency(m,1:2) = ([c(m) r(m)] - k_centre)*units; 
    phase(m) = angle(k_vessel_centre(r(m),c(m)));
    
    phase(m) = phase(m) + pi;
    
    for d=1:size(vesloc,2)
        pp(m,d) = frequency(m,1:2)*vesloc(1:2,d) + phase(m);
        if or(pp(m,d)>3.1416, pp(m,d)<-3.1415)
            pp(m,d) = mod(pp(m,d),2*3.1416);
            if and(pp(m,d)>3.1416,pp(m,d)<2*3.1416)
                pt = mod(pp(m,d),3.1416);
                pp(m,d) = pt-3.1415;
            end
        end
        encoding(m,d) = interp1(modmat_phase, modmat_uni_lookup,pp(m,d),'linear','extrap');
    end
end
static_tissue = ones(size(hadamard_wo_ctrl_and_tissue,1),1);
e = [encoding,static_tissue]; % calculate condition number
k_mm = frequency/(2*pi);
tagang = zeros(1,size(frequency,1));
ab = zeros(1,size(frequency,1));
for n = 1:size(frequency,1)
    tagang(n) = atan2(k_mm(n,1),k_mm(n,2));
    %convert to degrees
    tagang(n) = torad2deg(tagang(n));
    ab(n) = -phase(n)/(2*pi*sqrt(k_mm(n,1).^2 + k_mm(n,2).^2));
end

%Create pairs of the ab tag/control locations if standard or PCASL encoding
%is used
encs = zeros(2,size(frequency,1));
half_lambda = zeros(1,size(frequency,1));
for n = 1:size(frequency,1)
    half_lambda(:,n) = 1/(2*sqrt(k_mm(n,1).^2 + k_mm(n,2).^2));
    encs(:,n) = [ab(n)-half_lambda(:,n), ab(n)];
end

tagmode = 2*ones(1,size(frequency,1));
ang_locs_tagmode = [tagang;encs;tagmode];
ang_locs_tagmode = ang_locs_tagmode';

condition = cond(e); min_half_lambda = min(half_lambda(:));
cost_1 = 0.5*(1-1/condition)^2; cost_2 = 0.5*(1/(min(half_lambda(:))/motion-1))^2;
cost = cost_1 + cost_2;

end