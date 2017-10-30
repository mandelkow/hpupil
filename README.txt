%HPUPIL [2B5] Bounding box for pupil from monkey video.
%
% Pupil(t,[x,y,dx,dy,t]) = hpupil(VideoFile,'param',value,...)
%
% EXAMPLES:
%   P = hpupil2('EyeVideo.avi');
%   % Process EyeVideo.avi, H = P(:,4) = pupil hight in each frame.
%
%   P = hpupil2('EyeVideo.avi','SKIP',10,'DISP',1);
%   % Read EyeVideo.avi and process every 10th frame with graphic display.
%
% For eyeliner correction let pupil cover less than 50% of the window
% width.
%
% AUTHOR: Hendrik.Mandelkow@nih.gov. 2014-07-21
%
% REFERENCE:
% Chang C, Leopold DA, Schölvinck ML, Mandelkow H, Picchioni D, Liu X, Ye FQ,
% Turchi JN, Duyn JH. Tracking brain arousal fluctuations with fMRI. Proc Natl Acad
% Sci U S A. 2016 Apr 19;113(16):4518-23. doi: 10.1073/pnas.1520613113. Epub 2016
% Apr 5. PubMed PMID: 27051064; PubMed Central PMCID: PMC4843437.
