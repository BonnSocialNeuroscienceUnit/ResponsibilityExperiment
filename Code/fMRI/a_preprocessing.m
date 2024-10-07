%% Batch script to preprocess fMRI data using SPM12
% Does DICOM conversion, realignment and unwarping (no fieldmap use
% though), normalisation, smoothing, of the specified runs and T1 image 
% folders, in all participant folders found in the studyPath data folder.

% Paths
spm( 'defaults', 'fmri' );
spm_jobman( 'initcfg' );
spmDir = spm('dir'); % this works if SPM12 is in the path.

% Works through all participants found in studyPath:
studyPath = '/Volumes/LaCie_RAID_A/Documents/exp/SoDec/';

% Settings
preprocLogFile = fopen( 'Responsibility_preprocessing_SPM12_add_date_here.txt', 'w' );
smoothingFWHM = [6 6 6]; smoothingPrefix = 's6';
inputDataType = 'IMA'; % or 'f'
runNames = {'rutledge1','rutledge2'};
T1dirName = '5_T1'; % Not a cell array, and I assume there's only 1 T1 per participant

%% Now go to work
participantList = dir( studyPath );
failedParticipants = zeros(1,length(participantList));
timeItTook = zeros(1,length(participantList));

for participant = 1 : length( participantList ) % for all participants:
    
    isdir1 = ~strcmp( participantList( participant ).name, '.' ) ;
    isdir2 = ~strcmp( participantList( participant ).name, '..' ) ;
    if ( participantList( participant ).isdir == 1 ) && isdir1 && isdir2 % if suitable folder
        
        participantDir = participantList( participant ).name;
        fprintf( preprocLogFile, 'Start with : %s/r/n', participantDir );
        
        name = participantDir;
        
        %%%%_______________________________________________________________%%%%%
        idx = 1; % start index for matlabbatch
        clear matlabbatch

        %% Convert DICOM to NIFTI if we need to
        if strcmpi(inputDataType,'IMA')
          for r = 1:length(runNames)
            thisRunDir = [studyPath participantDir filesep runNames{r}];
            dataFiles{r} = spm_select('FPList', thisRunDir, '^.*\.IMA$'); % get the files
            if ~isempty(dataFiles{r}) % there are datafiles
              matlabbatch{idx}.spm.util.import.dicom.data = cellstr(dataFiles{r});
              matlabbatch{idx}.spm.util.import.dicom.root = 'flat';
              matlabbatch{idx}.spm.util.import.dicom.outdir = cellstr(thisRunDir);
              matlabbatch{idx}.spm.util.import.dicom.protfilter = '.*';
              matlabbatch{idx}.spm.util.import.dicom.convopts.format = 'nii';
              matlabbatch{idx}.spm.util.import.dicom.convopts.icedims = 0;
              idx = idx+1;
            end % there are datafiles
          end % each run
          
          % Now realign and unwarp the to-be-converted files ! -> need to use cfg_dep
          idx2 = 0; % run counter for the realignment
          for r = 1:length(runNames)
            dataFiles{r} = spm_select('FPList', thisRunDir, '^.*\.IMA$'); % get the files
            if ~isempty(dataFiles{r}) % there are datafiles
              idx2 = idx2+1;
              % we specify the files in the following way so as to have realignment parameters in each run folder:
              matlabbatch{idx}.spm.spatial.realignunwarp.data(idx2).scans = cfg_dep('DICOM Import: Converted Images', substruct(...
                '.','val', '{}',{idx2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
              matlabbatch{idx}.spm.spatial.realignunwarp.data(idx2).pmscan = '';
            end % we have data
          end % runs
          
        else % No need to convert from IMA, so get the f*.nii directly
          % Realign and unwarp the existing files
          idx2 = 0; % just for the realignment!
          for r = 1:length(runNames)
            thisRunDir = [studyPath participantDir filesep runNames{r}];
            dataFiles{r} = spm_select('FPList', thisRunDir, '^f.*\.nii$'); % get the nifti files
            if ~isempty(dataFiles{r}) % there are datafiles
              idx2 = idx2+1;
              % we specify the files in the following way so as to have realignment parameters in each run folder:
              matlabbatch{idx}.spm.spatial.realignunwarp.data(idx2).scans = cellstr(dataFiles{r});
              matlabbatch{idx}.spm.spatial.realignunwarp.data(idx2).pmscan = '';
            end % we have data
          end % runs
        end % DICOM (IMA) to Nifti conversion or not

        % Continue specifying the realignment options
        matlabbatch{idx}.spm.spatial.realignunwarp.eoptions.quality = 0.9;
        matlabbatch{idx}.spm.spatial.realignunwarp.eoptions.sep = 4;
        matlabbatch{idx}.spm.spatial.realignunwarp.eoptions.fwhm = 5;
        matlabbatch{idx}.spm.spatial.realignunwarp.eoptions.rtm = 0;
        matlabbatch{idx}.spm.spatial.realignunwarp.eoptions.einterp = 2;
        matlabbatch{idx}.spm.spatial.realignunwarp.eoptions.ewrap = [0 0 0];
        matlabbatch{idx}.spm.spatial.realignunwarp.eoptions.weight = '';
        matlabbatch{idx}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
        matlabbatch{idx}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
        matlabbatch{idx}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
        matlabbatch{idx}.spm.spatial.realignunwarp.uweoptions.jm = 0;
        matlabbatch{idx}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
        matlabbatch{idx}.spm.spatial.realignunwarp.uweoptions.sot = [];
        matlabbatch{idx}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
        matlabbatch{idx}.spm.spatial.realignunwarp.uweoptions.rem = 1;
        matlabbatch{idx}.spm.spatial.realignunwarp.uweoptions.noi = 5;
        matlabbatch{idx}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
        matlabbatch{idx}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
        matlabbatch{idx}.spm.spatial.realignunwarp.uwroptions.rinterp = 4;
        matlabbatch{idx}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0];
        matlabbatch{idx}.spm.spatial.realignunwarp.uwroptions.mask = 1;
        matlabbatch{idx}.spm.spatial.realignunwarp.uwroptions.prefix = 'u';
        
        % Next, normalise, using the files created during realignment
        idx = idx+1;
        matlabbatch{idx}.spm.spatial.normalise.estwrite.subj.vol(1) = cfg_dep('Realign & Unwarp: Unwarped Mean Image', substruct('.','val', '{}',{idx-1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','meanuwr'));
        for ii = 1:idx2 % for all runs that have data
            str = sprintf('Realign & Unwarp: Unwarped Images (Run %d)',ii);
            matlabbatch{idx}.spm.spatial.normalise.estwrite.subj.resample(ii) = cfg_dep(str, substruct(...
              '.','val', '{}',{idx-1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{ii}, '.','uwrfiles'));
        end        
        matlabbatch{idx}.spm.spatial.normalise.estwrite.subj.resample(ii+1) = cfg_dep('Realign & Unwarp: Unwarped Mean Image', substruct('.','val', '{}',{idx-1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','meanuwr'));
        matlabbatch{idx}.spm.spatial.normalise.estwrite.subj.vol(1) = cfg_dep('Realign & Unwarp: Unwarped Mean Image', substruct('.','val', '{}',{idx-1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','meanuwr'));
        for ii = 1:idx2 % for all sessions that have data
            str = sprintf('Realign & Unwarp: Unwarped Images (Run %d)',ii);
            matlabbatch{idx}.spm.spatial.normalise.estwrite.subj.resample(ii) = cfg_dep(str, substruct(...
              '.','val', '{}',{idx-1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{ii}, '.','uwrfiles'));
        end        
        matlabbatch{idx}.spm.spatial.normalise.estwrite.subj.resample(ii+1) = cfg_dep('Realign & Unwarp: Unwarped Mean Image', substruct('.','val', '{}',{idx-1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','meanuwr'));
        matlabbatch{idx}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
        matlabbatch{idx}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
        matlabbatch{idx}.spm.spatial.normalise.estwrite.eoptions.tpm = {fullfile(spmDir, 'tpm', 'TPM.nii')};
        matlabbatch{idx}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
        matlabbatch{idx}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
        matlabbatch{idx}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
        matlabbatch{idx}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
        matlabbatch{idx}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70
            78 76 85];
        matlabbatch{idx}.spm.spatial.normalise.estwrite.woptions.vox = [2 2 2];
        matlabbatch{idx}.spm.spatial.normalise.estwrite.woptions.interp = 4;
        matlabbatch{idx}.spm.spatial.normalise.estwrite.woptions.prefix = 'w';

        % Smoothe the normalised images
        idx = idx+1;
        matlabbatch{idx}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Estimate & Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{idx-1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
        matlabbatch{idx}.spm.spatial.smooth.fwhm = smoothingFWHM;
        matlabbatch{idx}.spm.spatial.smooth.dtype = 0;
        matlabbatch{idx}.spm.spatial.smooth.im = 0;
        matlabbatch{idx}.spm.spatial.smooth.prefix = smoothingPrefix;
        
        % Done with the functional images. Now process the T1 image.

        % T1 processing: get datafiles, convert to NII if needed, then normalise
        dataFolderT1 = [studyPath participantDir filesep T1dirName];

        if strcmpi(inputDataType,'IMA') % then as for the functional data, get the T1 image DICOMs
          idx = idx+1;
          dataT1 = spm_select('FPList',dataFolderT1,'^.*\.IMA$');  % char array, many DICOMs for 1 T1 scan
          if isempty(dataT1); disp(['No T1 DICOM images found in ' dataFolderT1]); end
          matlabbatch{idx}.spm.util.import.dicom.data = cellstr(dataT1);
          matlabbatch{idx}.spm.util.import.dicom.root = 'flat';
          matlabbatch{idx}.spm.util.import.dicom.outdir = cellstr(dataFolderT1);
          matlabbatch{idx}.spm.util.import.dicom.protfilter = '.*';
          matlabbatch{idx}.spm.util.import.dicom.convopts.format = 'nii';
          matlabbatch{idx}.spm.util.import.dicom.convopts.icedims = 0;

          idx = idx+1;
          matlabbatch{idx}.spm.spatial.normalise.estwrite.subj.vol(1) = cfg_dep('DICOM Import: Converted Images', substruct('.','val', '{}',{idx-1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
          matlabbatch{idx}.spm.spatial.normalise.estwrite.subj.resample(1) = cfg_dep('DICOM Import: Converted Images', substruct('.','val', '{}',{idx-1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));

        else % get the s*.nii
          idx = idx+1;
          dataT1 = spm_select('FPList',dataFolderT1,'^s.*\.nii$');  % char array, only 1 T1 scan
          if isempty(dataT1); disp(['No T1 Nifti image found in ' dataFolderT1]); end
          matlabbatch{idx}.spm.spatial.normalise.estwrite.subj.vol = cellstr(dataT1);
          matlabbatch{idx}.spm.spatial.normalise.estwrite.subj.resample = cellstr(dataT1);
%           matlabbatch{idx}.spm.spatial.normalise.estwrite.subj.vol(1) = cellstr(dataT1);
%           matlabbatch{idx}.spm.spatial.normalise.estwrite.subj.resample(1) = cellstr(dataT1);
        end
        
        % Continue specifying the normalisation options
        matlabbatch{idx}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
        matlabbatch{idx}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
        matlabbatch{idx}.spm.spatial.normalise.estwrite.eoptions.tpm = {fullfile(spmDir, 'tpm', 'TPM.nii')};
        matlabbatch{idx}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
        matlabbatch{idx}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
        matlabbatch{idx}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
        matlabbatch{idx}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
        matlabbatch{idx}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70
            78 76 85];
        matlabbatch{idx}.spm.spatial.normalise.estwrite.woptions.vox = [1 1 1];
        matlabbatch{idx}.spm.spatial.normalise.estwrite.woptions.interp = 4;
        matlabbatch{idx}.spm.spatial.normalise.estwrite.woptions.prefix = 'w';
        
        % Run, keep track of time
        tic
        try % if crashes for one person, does not stop the whole batch
          % spm_jobman( 'interactive' , matlabbatch ); % For debugging
          spm_jobman( 'run' , matlabbatch );
        catch
          failedParticipants(participant) = 1;
        end
        timeItTook(participant) = toc;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end % for each suitable participant folder
    
end % for all participants

%% Finished, report
if sum(failedParticipants)>0 % some did not work
  disp('Failed:')
  disp(char(participantList(find(failedParticipants)).name))
else
  disp('All good!')
  disp(['This took, for each participant, in seconds: ' num2str(round(timeItTook))])
end
fprintf(preprocLogFile, 'Done!' );
fclose( preprocLogFile );
