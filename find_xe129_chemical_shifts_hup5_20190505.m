% find_xe129_chemical_shifts_hup5_20190505

%% modified code from the 13C dataset
% define chemical shifts
cs_gas = 0;
cs_blood = 210;
cs_tissue = 196;

% find the index maximum peak (pyruvate)
[~,ref_index] = max(squeeze(sum(sum(abs(complex_img),1),2)));

% find ppm vector
n = 1:Np;
sw_ppm = Bandwidth/Freq*1e6;
delta_cs = sw_ppm/length(n);
ref_ppm = 0;
cs = (ref_index-n)*delta_cs+ref_ppm;

% find corresponding indices

if show_real_mrsi
    cs_gas = cs(intersect(find(cs<18),find(cs>-18)));
    cs_blood = cs(intersect(find(cs<223),find(cs>205)));
    cs_tissue = cs(intersect(find(cs<199),find(cs>134)));
    
    gas_idx = intersect(find(cs<18),find(cs>-18)); % gas
    blood_idx = intersect(find(cs<223),find(cs>205)); % blood
    tissue_idx = intersect(find(cs<199),find(cs>134)); % tissue
else
    cs_gas = cs(intersect(find(cs<40),find(cs>-40)));
    cs_blood = cs(intersect(find(cs<230),find(cs>212)));
    cs_tissue = cs(intersect(find(cs<196),find(cs>165)));
    
    
    gas_idx = intersect(find(cs<40),find(cs>-40)); % gas
    blood_idx = intersect(find(cs<230),find(cs>210)); % blood
    tissue_idx = intersect(find(cs<194),find(cs>165)); % tissue
end

peaks = [gas_idx blood_idx tissue_idx];
peaks_dis = [blood_idx tissue_idx];
cs_peaks = [cs_blood cs_tissue];
baseline_idx = setdiff(1:Np,peaks); 
