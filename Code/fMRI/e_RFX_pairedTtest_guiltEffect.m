% Script to run SPM RFX paired t tests
%
% a direct comparison between negative partner outcomes of lotteries
% chosen by the participant vs the same outcomes resulting from lotteries
% chosen by the partner, using paired t-test.

%% Base settings
baseDir = '/Volumes/LaCie_RAID_A/Documents/exp/SoDec/Responsibility/';
subj = [ 301 302 304 305 306 307 308 309 310 311 312 313 314 315 316 317 318 319 320 321 322 323 324 325 326 327 328 329 330 331 332 333 347 334 335 336 337 338 339 340 ];

%% FFX contrasts
%  1 Decision safe social                             
%  2 Decision safe partner                            
%  3 Decision safe solo                               
%  4 Decision risky social                            
%  5 Decision risky partner                           
%  6 Decision risky solo                              
%  7 Outcome safe social                              
%  8 Outcome safe partner                             
%  9 Outcome safe solo                                
% 10 Outcome other pos social                         
% 11 Outcome other pos partner                        
% 12 Outcome self pos solo                            
% 13 Outcome other neg social                         
% 14 Outcome other neg partner                        
% 15 Outcome self neg solo                            
% 16 Guilt effect - Outcome other neg social > partner

%% Contrasts to test
% make RFX paired t-tests for the critical guilt contrast
% Guilt: risky > safe    Reg No: 13 vs 14
i=1; cons{i,1} = 'Guilt - neg social > neg partner';           cons{i,2} = [13 14];  cons{i,3} = [baseDir 'FFX/main/'];

% Explicit mask for model
explicitMask = '/Volumes/LaCie_RAID_A/Documents/exp/SoDec/Responsibility/RFX/FullFact/main/outcome/risky>safe_0p05FWE_clust.nii';

% Thresholding and SVC if desired
thresh = 0.001;
minClustSize = 1;
SVCimage = '/Users/johannesschultz/sciebo/NIFTIbrainMasks/Anatomy/Hammers_anterior_insula_left_86+88+92.nii';

% Go through the specified models
for c = 1 %1:size(cons,1)
  % set up
  clear matlabbatch
  matlabbatch{1}.spm.stats.factorial_design.dir = {[baseDir 'RFX/ttests/' cons{c,1} filesep]};
  for s = 1:length(subj)
    matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(s).scans = {
      fullfile(cons{c,3},sprintf('SODEC_FMRI_%d',subj(s)),sprintf('con_%04d.nii,1',cons{c,2}(1)));...
      fullfile(cons{c,3},sprintf('SODEC_FMRI_%d',subj(s)),sprintf('con_%04d.nii,1',cons{c,2}(2)));...
      };
  end
  matlabbatch{1}.spm.stats.factorial_design.des.pt.gmsca = 0;
  matlabbatch{1}.spm.stats.factorial_design.des.pt.ancova = 0;
  matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
  matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
  matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
  matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
  if isempty(explicitMask)
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
  else
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {explicitMask};
  end
  matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
  matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
  matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
  
  % estimate
  matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
  matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
  matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
  
  % create contrast
  matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
  matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = cons{c,1};
  matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
  matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
  matlabbatch{3}.spm.stats.con.delete = 0;
  
  % pick this contrast and threshold
  matlabbatch{4}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
  matlabbatch{4}.spm.stats.results.conspec.titlestr = '';
  matlabbatch{4}.spm.stats.results.conspec.contrasts = 1;
  matlabbatch{4}.spm.stats.results.conspec.threshdesc = 'none';
  matlabbatch{4}.spm.stats.results.conspec.thresh = thresh;
  matlabbatch{4}.spm.stats.results.conspec.extent = minClustSize;
  matlabbatch{4}.spm.stats.results.conspec.conjunction = 1;
  if isempty(SVCimage)
    matlabbatch{4}.spm.stats.results.conspec.mask.none = 1;
  else
    matlabbatch{4}.spm.stats.results.conspec.mask.image.name = {SVCimage};
  end
  matlabbatch{4}.spm.stats.results.conspec.mask.image.mtype = 0;
  matlabbatch{4}.spm.stats.results.units = 1;
  matlabbatch{4}.spm.stats.results.export{1}.ps = false;

  % run all this; to start it by hand, run this: spm_jobman('interactive',matlabbatch);
  % spm_jobman('interactive',matlabbatch);
  spm_jobman('run',matlabbatch);
  
  % then if desired, do small volume correction (using GUI) for left insula using this:
  if ~isempty(SVCimage)
    xY.def = 'mask';
    xY.spec = SVCimage;
    load([matlabbatch{1}.spm.stats.factorial_design.dir{1} filesep 'SPM.mat'])
    hReg = [];
    [TabDat,xSVC] = spm_VOI(SPM,xSPM,hReg,xY); % Applies SVC and updates results table in figure
  end
end
