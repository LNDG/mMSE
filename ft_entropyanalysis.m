function mse = ft_entropyanalysis(cfg, data)

% 181019 JQK | line 131, 249, 259, 261 287 commented
%            | line 163 verbosity changed
%            | save r parameter
%            | add option for HPF & BPF
%            | also use coarsegrainmethod for scale 1 (e.g., HPF)
% 180226 JQK | 364 ff. change encoding from sc to s to allow for nonlinear scale encoding
% 190130 JQK | filter entire time series first, then temporally segment
%            | do not require data.trialinfo, get trialsize from .trial
% 190405 JQK | integrated changes from previous scripts: adpative padding
%               length; no LPF at scale 1; changed indenting; removed
%               outdated code

% FT_ENTROPYANALYSIS performs entropy and time-entropy analysis
% on time series data over multiple trials
%
% Use as
%   mse = ft_entropyanalysis(cfg, data)
%
% The input data should be organised in a structure as obtained from
% the FT_PREPROCESSING function. The configuration
% depends on the type of computation that you want to perform.
%
% cfg is a configuration structure that should contain
%
%  cfg.toi        = vector 1 x numtoi, the times on which the analysis
%                    windows should be centered (in seconds), or TODO a string
%                    such as '50%' or 'all' (default).  Both string options
%                    use all timepoints available in the data, but 'all'
%                    centers an entropy estimate on each sample, whereas
%                    the percentage specifies the degree of overlap between
%                    the shortest time windows from cfg.timwin.
%  cfg.timwin     = vector 1 x numfoi, length of time window (in seconds)
%  cfg.timescales = vector 1 x numtimescales, the time scales to compute MSE for.
%                   Scale 1 is the fastest scale, i.e. sample entropy at the native
%                   sampling rate of the signal. Slower scales are achieved
%                   by coarse graining the data. The highest scales achievable
%                   is determined by pattern length m and the time window timwin: at
%                   least m+1 samples need to be present in the time window
%                   for MSE computation at this scale.
%  cfg.coarsegrainmethod = string, method used for coarse%  graining:'filt_skip'
%                    (default) (filter, then skip points) or 'pointavg'
%                    (average groups of timepoints)
%  cfg.filtmethod = string, method used for filtering: {lp, hp, bp, no}
%  cfg.m          = pattern length for MSE computation, default is 2
%  cfg.r          = similarity criterion, set as a fraction of the time
%                   series SD. Default is 0.5.
%  cfg.recompute_r = recompute r parameter. 'perscale' or 'perscale_toi_sp'
%  cfg.mem_available = Memory available to perform computations (default
%                     8e9 bytes).
%  cfg.allowgpu     = 1 to use gpu if available, 0 to force
%                     cpu computation (default 1). if a gpu is found,
%                     available memory on that gpu is used.
%
%
% The configuration can optionally contain
%   cfg.option3   = value, explain it here (default is automatic)
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also <<give a list of function names, all in capitals>>

% Copyright (C) 2018, MPIB Berlin, Niels Kloosterman
%
% Here comes the Revision tag, which is auto-updated by the version control system
% $Id$

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the initial part deals with parsing the input options and data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function

% the ft_preamble function works by calling a number of scripts from
% fieldtrip/utility/private that are able to modify the local workspace

ft_defaults                   % this ensures that the path is correct and that the ft_defaults global variable is available
ft_preamble init              % this will reset ft_warning and show the function help if nargin==0 and return an error
ft_preamble debug             % this allows for displaying or saving the function name and input arguments upon an error
ft_preamble loadvar    data % this reads the input data in case the user specified the cfg.inputfile option
ft_preamble provenance data % this records the time and memory usage at the beginning of the function
ft_preamble trackconfig       % this converts the cfg structure in a config object, which tracks the cfg options that are being used

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
    % do not continue function execution in case the outputfile is present and the user indicated to keep it
    return
end

% ensure that the input data is valid for this function, this will also do
% backward-compatibility conversions of old data that for example was
% read from an old *.mat file
data = ft_checkdata(data, 'datatype', {'raw+comp', 'raw'}, 'feedback', 'yes', 'hassampleinfo', 'yes');

% % TODO check if the input cfg is valid for this function
% cfg = ft_checkconfig(cfg, 'renamed',     {'blc', 'demean'});
% cfg = ft_checkconfig(cfg, 'renamed',     {'blcwindow', 'baselinewindow'});

% ensure that the required options are present
cfg = ft_checkconfig(cfg, 'required', {'toi', 'timescales'});

% ensure that the options are valid
% cfg = ft_checkopt(cfg, 'normalized_r', 'double', {0, 1});
cfg = ft_checkopt(cfg, 'recompute_r', 'char', {'perscale_toi_sp', 'per_scale'});
cfg = ft_checkopt(cfg, 'coarsegrainmethod', 'char', {'filtskip', 'pointavg', 'bp', 'hp'});

% get the options
cfg.trials        = ft_getopt(cfg, 'trials',     'all', 1);
cfg.channel       = ft_getopt(cfg, 'channel',    'all');
toi               = ft_getopt(cfg, 'toi'); % time points for mse, e.g. cfg.toi  = -0.75:0.05:1.5;
timescales        = ft_getopt(cfg, 'timescales'); % time scales, depends on sample rate and winsize
timwin            = ft_getopt(cfg, 'timwin', 0.5); % e.g. 0.5 s
m                 = ft_getopt(cfg, 'm', 2); % pattern length, e.g. 2
r                 = ft_getopt(cfg, 'r', 0.5); % similarity criterion, 0.5
recompute_r       = ft_getopt(cfg, 'recompute_r', 'perscale_toi_sp'); % recompute r for each scale (1)
coarsegrainmethod = ft_getopt(cfg, 'coarsegrainmethod', 'filtskip'); % coarsening_filt_skip or coarsening_avg
filtmethod        = ft_getopt(cfg, 'filtmethod', 'lp'); % coarsening_filt_skip or coarsening_avg
mem_available     = ft_getopt(cfg, 'mem_available', 8e9); % 8 GB
allowgpu          = ft_getopt(cfg, 'allowgpu', 1); % 8 GB

gpuavailable = gpuDeviceCount;
if allowgpu && gpuavailable
    fprintf('GPU device found. Running things there\n')
    gpu = gpuDevice;
    mem_available = gpu.AvailableMemory * 0.6; % only use % of available mem, other vars also required there
end

% select channels and trials of interest, by default this will select all channels and trials
tmpcfg = keepfields(cfg, {'trials', 'channel', 'showcallinfo'});
data = ft_selectdata(tmpcfg, data);
% restore the provenance information
%[cfg, data] = rollback_provenance(cfg, data);

% preallocate matrices
nchan = length(data.label);
nscales = length(timescales);
ntoi = size(toi,2);
sampen = nan(nchan, nscales, ntoi);

r_estimate = nan(nchan, nscales, ntoi, nscales); % dimord chan nsc ntoi nstartingpts

for s = 1:numel(timescales) %  loop through timescales
    sc = timescales(s);
    
    % apply filtering here
    switch filtmethod
        case 'lp'
            if sc == 1
                data_filt = data;
            else
                fs = data.fsample;
                nyquist = (fs/2);
                fcLowPass = (1/sc)*nyquist;
                if fcLowPass == nyquist
                    fcLowPass = fcLowPass-1;
                end
                [B,A] = butter(6,fcLowPass/nyquist);
                cfg.freq(1,s) = fcLowPass;
                
                padlength = ceil(size(data.trial{1},2)./2); % use half the length of trial 1 as padding (JQK)
                x_pad = cellfun(@(a) ft_preproc_padding(a, 'mean', padlength), data.trial, 'UniformOutput', false );  % add padding
                x_pad = cellfun(@transpose, x_pad, 'UniformOutput', false); % transpose for filtfilt: time x chan
                resamp_x_pad = cellfun(@(x_pad) filtfilt(B,A,x_pad), x_pad, 'UniformOutput', false );  % filter data
                resamp_x_pad = cellfun(@transpose, resamp_x_pad, 'UniformOutput', false); % transpose back : chan x time again
                resamp_x = cellfun(@(resamp_x_pad) ft_preproc_padding(resamp_x_pad, 'remove', padlength), resamp_x_pad, 'UniformOutput', false );  % remove padding
                % create data_filt structure
                data_filt = data;
                data_filt.trial = resamp_x;
                clear resamp_* x_pad;
            end
        case 'bp'
            fs = data.fsample;
            nyquist = fs/2;
            fcLowPass = (1/sc)*nyquist;
            if fcLowPass == nyquist
                fcLowPass = fcLowPass-1;
            end
            if s == numel(timescales)
                fcHighPass = 0.5;
            else
                fcHighPass = (1/(timescales(s+1)))*nyquist;
            end
            [B,A] = butter(6,fcLowPass/nyquist);            % define low-pass filter: https://de.mathworks.com/help/signal/ref/butter.html
            [D,C] = butter(6,fcHighPass/nyquist, 'high');   % define high-pass filter
            cfg.freq(1,s) = fcLowPass;
            cfg.freq(2,s) = fcHighPass;
            
            padlength = ceil(size(data.trial{1},2)./2); % use half the length of trial 1 as padding (JQK)
            x_pad = cellfun(@(a) ft_preproc_padding(a, 'mean', padlength), data.trial, 'UniformOutput', false );    % add padding
            x_pad = cellfun(@transpose, x_pad, 'UniformOutput', false);                                                 % transpose for filtfilt: time x chan
            resamp_x_pad = cellfun(@(x_pad) filtfilt(B,A,x_pad), x_pad, 'UniformOutput', false );                       % low-pass filter data
            if sc == 1
            else
                resamp_x_pad = cellfun(@(resamp_x_pad) filtfilt(D,C,resamp_x_pad), resamp_x_pad, 'UniformOutput', false );  % high-pass filter data
            end
            resamp_x_pad = cellfun(@transpose, resamp_x_pad, 'UniformOutput', false);                                   % transpose back : chan x time again
            resamp_x = cellfun(@(resamp_x_pad) ft_preproc_padding(resamp_x_pad, 'remove', padlength), ...                % remove padding
                resamp_x_pad, 'UniformOutput', false );
            %figure; hold on; plot(resamp_x{1}(1,:)); plot(data.trial{1}(1,:))
            % create data_filt structure
            data_filt = data;
            data_filt.trial = resamp_x;
            clear resamp_* x_pad;
        case 'hp'
            fs = data.fsample;
            nyquist = (fs/2);
            fcHighPass = (1/(sc+1))*nyquist;
            [D,C] = butter(6,fcHighPass/nyquist, 'high');   % define high-pass filter
            cfg.freq(1,s) = fcHighPass;
            
            padlength = ceil(size(data.trial{1},2)./2); % use half the length of trial 1 as padding (JQK)
            x_pad = cellfun(@(a) ft_preproc_padding(a, 'mean', padlength), data.trial, 'UniformOutput', false );    % add padding
            x_pad = cellfun(@transpose, x_pad, 'UniformOutput', false);                                                 % transpose for filtfilt: time x chan
            resamp_x_pad = cellfun(@(x_pad) filtfilt(D,C,x_pad), x_pad, 'UniformOutput', false );                       % low-pass filter data
            resamp_x_pad = cellfun(@transpose, resamp_x_pad, 'UniformOutput', false);                                   % transpose back : chan x time again
            resamp_x = cellfun(@(resamp_x_pad) ft_preproc_padding(resamp_x_pad, 'remove', padlength), ...                % remove padding
                resamp_x_pad, 'UniformOutput', false );
            %figure; hold on; plot(resamp_x{1}(1,:)); plot(data_sel.trial{1}(1,:))
            % create data_filt structure
            data_filt = data;
            data_filt.trial = resamp_x;
            clear resamp_* x_pad;
        case 'no'
            data_filt = data;
    end
    
    if strcmp(recompute_r, 'per_scale') % recompute r for each scale or: sc toi sp
        % per_scale
        % per_toi
        % pertoi_sp (fixed per scale)
        % perscale_toi_sp (run til now)
        % perscale_toi
        
        r_new = r * std(cell2mat(data_filt.trial),1,2);
    end
    
    for itoi = 1:ntoi
        
        fprintf('Scale %d of %d; Time %d of %d\n', s, length(timescales),itoi, ntoi)
        
        % select time window of interest from each trial
        tmpcfg=[];
        tmpcfg.toilim = [toi(itoi)-timwin*0.5 toi(itoi)+timwin*0.5];
        tmpcfg.showcallinfo = 'no';
        data_sel = ft_redefinetrial(tmpcfg, data_filt);
        
        % only take trials that have the whole interval
        tmpcfg = [];
        tmpcfg.minlength = timwin;
        tmpcfg.showcallinfo = 'no';
        data_sel = ft_redefinetrial(tmpcfg, data_sel);
        
        % need 40 samples for mse calc, 3 smp per trial for scale 42: 40/3 = 13.3 trials, make 15
        ntrials = size(data_sel.trial,2);
        if ntrials < 1
            warning('Time point %g: Not enough trials remain', toi(itoi))
            break % subsequent time points will also not work
        end
        
        if strcmp(recompute_r, 'per_toi') % not per scale
            
            % select time window of interest from each trial
            tmpcfg=[];
            tmpcfg.toilim = [toi(itoi)-timwin*0.5 toi(itoi)+timwin*0.5];
            data_sel_unfilt = ft_redefinetrial(tmpcfg, data);
            
            % only take trials that have the whole interval
            tmpcfg = [];
            tmpcfg.minlength = timwin;
            data_sel_unfilt = ft_redefinetrial(tmpcfg, data_sel_unfilt);
            
            % need 40 samples for mse calc, 3 smp per trial for scale 42: 40/3 = 13.3 trials, make 15
            ntrials = size(data_sel_unfilt.trial,2);
            if ntrials < 1
                warning('Time point %g: Not enough trials remain', toi(itoi))
                break % subsequent time points will also not work
            end
            
            % calculate similarity criterion
            r_new = r * std(cell2mat(data_sel_unfilt.trial),1,2);
            nchan = size(data_sel_unfilt.trial{1},1);
        elseif strcmp(recompute_r, 'perscale_toi')
            % calculate similarity criterion
            r_new = r * std(cell2mat(data_sel.trial),1,2);
            nchan = size(data_sel.trial{1},1);
        end
        
        % do point skipping
        cg_data = {};
        switch coarsegrainmethod
            case 'filtskip'
                nloops = sc;
                cg_data = cell(nloops,1); % make cell: cg_data{istart}{trials}(chan-by-time)
                resamp_x = data_sel.trial;
                for is = 1:nloops % loop over starting points here!
                    cg_data{is} = cellfun(@(resamp_x) resamp_x(:, is:(sc-1+1):end), resamp_x, 'UniformOutput', false );  % add padding% Filter
                end
                clear resamp_x;
            case 'pointavg' % original point averaging coarse graining, no loop over starting points
                if sc == 1 % no coarse graining for native sampling rate
                    cg_data{1} = data_sel.trial; %only keep trial data
                    nloops = 1; % no loop across starting points
                else % coarse-grain time series at this time scale
                    nloops = 1; % no loop across starting points
                    nchan = size(data_sel.trial{1},1);
                    for itrial = 1:length(data_sel.trial)
                        num_cg_tpts = floor(length(data_sel.trial{itrial})/sc); % number of coarse-grained time points
                        cg_data{1}{itrial} = zeros(nchan, num_cg_tpts); % preallocate cg_data matrix
                        for t = 1:num_cg_tpts % calculate coarse_grained time series
                            cg_data{1}{itrial}(:,t) = mean( data_sel.trial{itrial}(:, (t-1)*sc + [1:sc]) ,2);
                        end
                    end
                end
        end
        
        % after coarsegraining, loop mse computation across starting points
        allcont = zeros(sc, nchan, m+1); % start_chan_m
        for istart = 1:nloops
            
            if max(cellfun(@(x) size(x,2), cg_data{istart})) == m % TODO check this at start
                fprintf('Coarse grained trials below %d + 1 samples, skipping remaining starting points\n', m)
                break
            end
            
            %  concatenate trials and convert to single
            y = single(cell2mat(cg_data{istart}));
            
            % collect trial bounds and create mask with valid time points for pats
            trl_bounds = cumsum(cellfun(@(x) size(x,2), cg_data{istart}))';
            trl_mask = true(size(y,2),1);
            if allowgpu && gpuavailable
                trl_mask = gpuArray(trl_mask);
            end
            trl_mask([trl_bounds-1; trl_bounds]) = false;
            
            %  break if n data points < 100 (See Grandy et al., 2016)
            ndatapoints = length(trl_mask); % TODO check this at start
            if ndatapoints < 100
                fprintf('N data points < 100, breaking\n')
                break
            end
            
            %  Calculate sample entropy of coarse grained time series
            if strcmp(recompute_r, 'perscale_toi_sp')
                r_new = r * std(y,1,2);
            end
            
            % keep the estimated r parameter
            r_estimate(:, s, itoi, istart) = r_new; % dimord chan nsc ntoi nstartingpts
            
            % chunk y to keep memory usage in check
            max_numel = mem_available/4; % single = 4 bytes
            chunk_size = round(max_numel/numel(y));
            n_chunks = ceil(size(y,2)/chunk_size);
            temp = 1;
            chunk_borders = zeros(n_chunks,2);
            for ic = 1:n_chunks
                chunk_borders(ic,:) = [temp temp+chunk_size];
                temp = temp+chunk_size-1; % chunks need to overlap to avoid missing pats at chunk borders
            end
            chunk_borders(end) = size(y,2);
            clear temp
            
            %fprintf('starting point %d\n', istart)
            cont = zeros(m+1, n_chunks, nchan);
            y_chunk1 = shiftdim(y', -1 ); % insert singleton dim
            r_new2 = shiftdim(r_new, -2);
            if allowgpu && gpuavailable
                cont = gpuArray(cont);
                y_chunk1 = gpuArray(y_chunk1);
                r_new2 = gpuArray(r_new2);
            end
            
            fprintf('%d chunks: ', n_chunks)
            for ic = 1:n_chunks
                fprintf('%d ', ic)
                
                y_inds = transpose(chunk_borders(ic,1):chunk_borders(ic,2));
                
                y_chunk2 = permute(y_chunk1(1,y_inds,:), [2 1 3]); % insert singleton dim
                if allowgpu && gpuavailable
                    y_chunk2 = gpuArray(y_chunk2);
                end
                
                ymat = bsxfun(@le, abs(bsxfun(@minus, y_chunk1, y_chunk2 )), r_new2 );
                for ichan=1:nchan % loop since triu only supports 2D
                    ymat(:,:,ichan) = triu(ymat(:,:,ichan), chunk_borders(ic,1));
                end
                
                for k = 1:m+1
                    if k >= m % TODO try for m > 2
                        cont(k,ic,:) = sum(reshape(ymat(trl_mask(y_inds(1:end-2)), trl_mask, :), [], nchan));
                    end
                    if k < m+1
                        ymat = ymat & circshift(ymat, [-1 -1 0]);
                    end
                end
                clear ymat y_inds y_chunk2
            end
            
            allcont(istart, :, :) = gather(squeeze(sum(cont,2))'); % sum over chunks. dimord: start_chan_m
            %fprintf('\n')
        end % cg starting points
        
        allcont = sum(allcont,1); % sum counts over starting points
        
        if ndatapoints < 100
            fprintf('N data points < 100, breaking\n')
            break
        end
        
        %  calculate sample entropy
        for ichan=1:nchan
            if allcont(1,ichan,m+1) == 0 || allcont(1,ichan,m) == 0
                fprintf('zero patterns found!\n')
                %         nlin_sc = size(pnts,1); % ori THG code
                %         mse(s) = -log(1/((nlin_sc)*(nlin_sc -1)));
                npossiblepats = length(find(trl_mask));
                sampen(ichan,s,itoi) = -log(1/(npossiblepats*(npossiblepats-1)));
            else
                sampen(ichan,s,itoi) = -log(allcont(1,ichan,m+1)./allcont(1,ichan,m)); % same as log(cont(m)/cont(m+1))
            end
        end
        
    end % for toi
end % for timescales

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with the output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mse = [];
mse.label = data.label;
mse.fsample = bsxfun(@rdivide, data.fsample, timescales); % sample rates after coarse graining
mse.timescales = 1000 ./ mse.fsample; % by convention
mse.time = toi;
mse.dimord = 'chan_timescales_time'; 
mse.sampen = sampen;
mse.r = squeeze(nanmean(r_estimate,4)); % average across starting points
if isfield(data, 'trialinfo')
  mse.trialinfo = data.trialinfo;
end
mse.config = cfg;

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug               % this clears the onCleanup function used for debugging in case of an error
ft_postamble trackconfig         % this converts the config object back into a struct and can report on the unused fields
ft_postamble previous   data   % this copies the data.cfg structure into the cfg.previous field. You can also use it for multiple inputs, or for "varargin"
ft_postamble provenance mse  % this records the time and memory at the end of the function, prints them on screen and adds this information together with the function name and MATLAB version etc. to the output cfg
ft_postamble history    mse  % this adds the local cfg structure to the output data structure, i.e. dataout.cfg = cfg
ft_postamble savevar    mse  % this saves the output data structure to disk in case the user specified the cfg.outputfile option
