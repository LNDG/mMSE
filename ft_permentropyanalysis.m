function mse = ft_permentropyanalysis(cfg, data)

% 181019 JQK | line 131, 249, 259, 261 287 commented
%            | line 163 verbosity changed
%            | save r parameter
%            | add option for HPF & BPF
%            | also use coarsegrainmethod for scale 1 (e.g., HPF)
% 181026 JQK | 364 ff. change encoding from sc to s to allow for nonlinear scale encoding
% 191030 JQK | filter entire time series first, then temporally segment
%            | do not require data.trialinfo, get trialsize from .trial
% 190401 JQK | removed floor() for frequencies, at scale 1 use no LPF
% 190402 JQK | replaced sample entropy with permutation entropy
% 190405 JQK | corrected BP setting: at scale 1 only use HP
%               vs. LP

% FT_PERMENTROPYANALYSIS performs temporally-resolved permutation entropy
% on time series data over multiple trials
%
% Use as
%   mse = ft_permentropyanalysis(cfg, data)
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
%  cfg.t          = time delay, default is 1
%  cfg.nrm          = normalize permutation entropy to 1? (0/1; default: yes(1))
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
% Copyright (C) 2019, MPIB Berlin, Julian Q. Kosciessa
% uses code from Ouyang et al., (2013)
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
cfg = ft_checkopt(cfg, 'nrm', 'double', {0, 1});
cfg = ft_checkopt(cfg, 'coarsegrainmethod', 'char', {'filtskip', 'pointavg', 'bp', 'hp'});

% get the options
cfg.trials        = ft_getopt(cfg, 'trials',     'all', 1);
cfg.channel       = ft_getopt(cfg, 'channel',    'all');
toi               = ft_getopt(cfg, 'toi'); % time points for mse, e.g. cfg.toi  = -0.75:0.05:1.5;
timescales        = ft_getopt(cfg, 'timescales'); % time scales, depends on sample rate and winsize
timwin            = ft_getopt(cfg, 'timwin', 0.5); % e.g. 0.5 s
m                 = ft_getopt(cfg, 'm', 2); % pattern length, e.g. 2
t_lag             = ft_getopt(cfg, 't_lag', 1); % temporal lag, e.g. 1
nrm               = ft_getopt(cfg, 'nrm', 1); % normalize permutation entropy?
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
            if sc == 1 % only HPF
               resamp_x_pad = cellfun(@(x_pad) filtfilt(D,C,x_pad), x_pad, 'UniformOutput', false );  % high-pass filter data
            else
                resamp_x_pad = cellfun(@(x_pad) filtfilt(B,A,x_pad), x_pad, 'UniformOutput', false );                       % low-pass filter data
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

    for itoi = 1:ntoi
        
        fprintf('Scale %d of %d; Time %d of %d\n', s, length(timescales),itoi, ntoi)
        
        % select time window of interest from each trial
        tmpcfg=[];
        tmpcfg.toilim = [toi(itoi)-timwin*0.5 toi(itoi)+timwin*0.5];
        data_sel = ft_redefinetrial(tmpcfg, data_filt);

        % only take trials that have the whole interval
        tmpcfg = [];
        tmpcfg.minlength = timwin; 
        data_sel = ft_redefinetrial(tmpcfg, data_sel);

        % need 40 samples for mse calc, 3 smp per trial for scale 42: 40/3 = 13.3 trials, make 15
        ntrials = size(data_sel.trial,2);
        if ntrials < 1
            warning('Time point %g: Not enough trials remain', toi(itoi))
            break % subsequent time points will also not work
        end
            
        cg_data = {};
        switch coarsegrainmethod
            case 'filtskip'
                nloops = sc;
                cg_data = cell(nloops,1); % make cell: cg_data{istart}{trials}(chan-by-time)
                for is = 1:nloops % loop over starting points here!
                    resamp_x = data_sel.trial;
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
        tmp_permen = []; % start_chan_m
        for istart = 1:nloops
            
            if max(cellfun(@(x) size(x,2), cg_data{istart})) == m % TODO check this at start
                fprintf('Coarse grained trials below %d + 1 samples, skipping remaining starting points\n', m)
                break
            end
            
            %  concatenate trials and convert to single
            y = single(cell2mat(cg_data{istart}));
                        
            %  Calculate permutation entropy of coarse grained time series
            
            for ichan = 1:size(y,1)
            
                X = y(ichan,:);                         % select current coarse grained signal
                lX = length(X);                         % amount of samples in single-channel time series
                permlist = perms(1:m);                  % get all possible permutations of order m
                c(1:length(permlist))=0;                % preallocate output matrix
                 for j=1:lX-t_lag*(m-1)                     % step through all allowed patterns
                     [~,iv]=sort(X(j:t_lag:j+t_lag*(m-1)));     % sort in ascending order; step with step size t
                     for jj=1:length(permlist)          % check for all possible patterns
                         if (abs(permlist(jj,:)-iv))==0 % check whether there is a match for current pattern
                             c(jj) = c(jj) + 1 ;
                         end
                     end
                 end
                %hist = c;                   % save counter of how many matches were observed for each pattern
                c=c(c~=0);                  % keep only pattern that have occurred at all to avoid NaN situation
                p = c/sum(c);               % take ratio of matches to overall matches (normalizes to one across all patterns)
                                            % should be identical to : p = c./(lX-((m-1)*t)) (e.g. see Li et al., 2013)
                pe = -sum(p .* log(p));     % generate a single score of entropy across patterns
                %If specified, normalize PE for total number of patterns (see Ouyang, G., Li, J., Liu, X., & Li, X. (2013). Dynamic characteristics of absence EEG recordings with multiscale permutation entropy analysis.Epilepsy Research, 104(3), 246?252.
                if nrm==1
                    pattern_count=numel(permlist(:,1));%gives row count, reflecting discrete patterns from "perm_patterns."
                    pe = pe/log(pattern_count); % normalize such that matched amount of patterns leads to PE = 1 (Li et al., 2013)
                end
                tmp_permen(istart,ichan) = pe;
            end % channel
        end % cg starting points
        permen(:,s,itoi) = squeeze(nanmean(tmp_permen,1)); % average over chunks.
        clear tmp_permen;
        
    end % for timescales
end % for toi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with the output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mse = [];
mse.label = data.label;
mse.fsample = bsxfun(@rdivide, data.fsample, timescales); % sample rates after coarse graining
mse.timescales = 1000 ./ mse.fsample; % by convention
mse.time = toi;
mse.dimord = 'chan_timescales_time';
mse.permen = permen;
mse.config = cfg;

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug               % this clears the onCleanup function used for debugging in case of an error
ft_postamble trackconfig         % this converts the config object back into a struct and can report on the unused fields
ft_postamble previous   data   % this copies the data.cfg structure into the cfg.previous field. You can also use it for multiple inputs, or for "varargin"
ft_postamble provenance mse  % this records the time and memory at the end of the function, prints them on screen and adds this information together with the function name and MATLAB version etc. to the output cfg
ft_postamble history    mse  % this adds the local cfg structure to the output data structure, i.e. dataout.cfg = cfg
ft_postamble savevar    mse  % this saves the output data structure to disk in case the user specified the cfg.outputfile option
