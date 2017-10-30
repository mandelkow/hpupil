function Box = hpupil(VidFile,varargin)
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

% PREC: hSusani_ImgFN
% AUTH: HM, 2014-07-21, v1a
% AUTH: HM, 2014-07-21, v2a: Use parfor.
% AUTH: HM, 2015-09-14, v2a1: Add eyeliner correction.
% AUTH: HM, 2016-10-04, v2b2: Upgrade to Matlab 2016 VideoReader.
% AUTH: HM, 2017-08-16, v2b3: Add hpar and EYELINER parameter.
% AUTH: HM, 2017-08-16, v2b4: Add IMFUN parameter.
% AUTH: HM, 2017-08-16, v2b5: Replace SKIP and DISP with hpar(varargin)

% Pass varargin = {'parameter',value,...} to change any of these parameters:
SKIP = 1; % analyse video frames 1:SKIP:end
DISP = false; % graphic display off / on (slow!)
Thr = 25; % intensity threshold segregating dark pupil, querried later
EYELINER = false; % try to eliminate dark eyeliner
IMFUN = @(x) min(x,[],2); % e.g.
IMFUN = []; % pre-process image?
hpar(varargin{:});

% try NP = gcp('nocreate'); NP = NP.NumWorkers;
% catch, NP = 0;
% end
% img = harrcrop(img,Do.Crop);
% img = eval(sprintf('img(%s);',Do.Crop));
% [Do.CropThr,img] = harrcropthr(img,Do.CropThr,2);

%% Get AVI file
if nargout<1
	ExID = strtok(VidFile,'.');
	DataFile = [ExID,'_pup.mat']
	if exist(DataFile,'file'),
		% warning(['Overwrite file (or Ctrl-C)? ',DataFile]); pause;
		warning(['File exists ',DataFile,' - SKIP!']); Box = []; return;
	end
end

%% Read AVI file
fprintf('\nOpen video file: %s\n\n',VidFile);
Vid = VideoReader(VidFile) % +disp info
Fs = Vid.FrameRate;
NF = Vid.Duration*Vid.FrameRate; % replaces Vid.NumberOfFrames
img = Vid.readFrame;
img = img(:,:,1); % 1st color channel only
% NOTE: All img coordinates are ij = +y+x
H = figure('name',['hPupil - ',VidFile]); movegui(H,'northwest');
H = get(imshow(img),'parent');
% title('Aaaw, cute monkey...!','fontsize',16,'fontangle','italic','fontweight','bold','color','k');
xlabel('Draw crop box.','fontsize',14,'color','b','fontangle','italic');
Crop = round(getrect(H)); % [x0 y0 dx dy] = [m0 n0 dm dn] for Inm
img = img(Crop(2)+[0:Crop(4)-1],Crop(1)+[0:Crop(3)-1],1);
IMSPAR = {'Border','tight','Colormap',jet(64),'InitialMagn',200};
H = get(imshow(img,[0 127],IMSPAR{:}),'parent');
colorbar; hold on;
try
	Box = getfield(regionprops(img<Thr,'BoundingBox'),'BoundingBox');
	h = rectangle('position',Box,'edgecolor','m','linewidth',1);
	h(2) = rectangle('position',Box,'curvature',[1 1],'edgecolor','y','linewidth',1);
catch
end
tmp = input('Input upper threshold (default Thr=25): ');
if tmp, Thr = tmp; end; Thr
Box(floor(NF/SKIP),1:5) = 0;
fprintf('Read %u video frames (SKIP=%u): %6u',NF,SKIP,1);
% IMSPAR = {'Border','tight','Colormap',gray(64),'InitialMagn',200};
cmap = gray(64); cmap(1,:) = [0 0 1]; cmap(end,:) = [1 0 0];
IMSPAR = {'Border','tight','Colormap',cmap,'InitialMagn',200};
tic;
for nb = 1:size(Box,1),
	t = (nb-1)*SKIP/Fs;
	Vid.CurrentTime = t;
	tmp = Vid.readFrame;
	img = tmp(Crop(2)+[0:Crop(4)-1],Crop(1)+[0:Crop(3)-1],1);
	if EYELINER % Eyeliner correction
		img(img>120) = 120; % 127 -> 126 % aetetic
		[b,bi] = max(sum(img<Thr,2)); img(1:bi,:) = 127;
		img(mean(img<Thr,2)>0.5,:) = 127; % ***
		img(mean(img(:,[1:5,end-4:end])<Thr,2)>0.5,:) = 127; % ***
  end
  Dimg = img;
  if ~isempty(IMFUN)
    img = IMFUN(img);
  end
	try
		B = regionprops(img<Thr,'BoundingBox');
    B = cat(1,B.BoundingBox);
		B = sortrows(B,-4); B = B(1,:); % Find the box largest in hight.
		Box(nb,:) = [B,t];
	catch
		continue;
	end
	if DISP,
		delete(H); % clf;
		% H = get(imshow(img,[0 127],IMSPAR{:}),'parent'); hold on;
		H = get(imshow(Dimg,[Thr 127],IMSPAR{:}),'parent'); hold on;
		h = rectangle('position',Box(nb,1:4),'edgecolor','m','linewidth',1);
		h(2) = rectangle('position',Box(nb,1:4),'curvature',[1 1],'edgecolor','y','linewidth',1);
		drawnow;
	end
	fprintf('\b\b\b\b\b\b%6u',(nb-1)*SKIP+1);
end
fprintf('\tDONE.\n');
toc
clear Vid
%%
if nargout<1
	disp(['Save file: ',DataFile]);
	save(DataFile,'Box');
	Box = DataFile;
end

function hpar(varargin)
% Eval input in parameter-value pairs.
while ~isempty(varargin) && ~isempty(varargin{1})
  try
    assignin('caller',varargin{1},varargin{2});
  catch
    error('varargin must be pairs of {''parameter'',value,...}');
  end
  varargin(1:2) = [];
end
