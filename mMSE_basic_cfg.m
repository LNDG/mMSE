cfg = [];
cfg.timwin        = 0.5; % sliding window size
cfg.m             = 2; % pattern length
cfg.r             = 0.5; % similarity criterion 0.5
cfg.timescales    = 1:42; %1:40; % scale list
cfg.recompute_r = 'perscale_toi_sp';  
cfg.coarsegrainmethod = 'filtskip';  % pointavg filtskip
cfg.filtmethod = 'lp'; % low pass filter for pointskip 
cfg.mem_available = 16e9; % in bytes, 8e9 default
cfg.allowgpu = true;
cfg.toi  = -0.5:0.05:1; % set this according to your trial length

ft_entropyanalysis(cfg, data);  % 256 Hz data from the ft_preprocessing function
