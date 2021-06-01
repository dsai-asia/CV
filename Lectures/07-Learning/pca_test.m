function [Locs] = pca_test( filename, MeanImg, PCs, Vars )

k = size(PCs,2);
if ( size(Vars,1) < size(Vars,2) )
  Vars = Vars';
end;
Stds = sqrt( Vars(1:k) );

nr = size( MeanImg, 1 );
nc = size( MeanImg, 2 );
M = nr*nc;
MeanImg = MeanImg';
MeanImg = MeanImg(:);

Img = imread( filename );
if size(Img,3) == 3
  Img = rgb2gray( Img );
end

rowstep = round(nr/10);
colstep = round(nc/10);
scalefactor = 0.85;
nlocsinit = size( Img, 1 ) / rowstep * size( Img, 2 ) / colstep * 5;

% Loop over image scales
nlocs = 0;
Locs = zeros( nlocsinit, 5 );
iscale = 1;
while( size(Img, 1) >= nr & size( Img, 2 ) >= nc )
  % Skip tiny patches unlikely to be a real positive
  if ( iscale > 5 )
    Imgthis = Img;
    imgrows = size(Img,1);
    imgcols = size(Img,2);
    fprintf(1,'Scale %dx%d row', imgrows, imgcols );
    fflush(1);
    % Loop over image rows
    for row = 1:rowstep:(imgrows-nr+1)
      if ( mod(row-1,rowstep*10)==0 )
          fprintf(1,' %d...', row );
          fflush(1);
      end;
      % Loop over image columns
      for col = 1:colstep:(imgcols-nc+1)
        % Extract and normalize the image patch
        testimg = mean( Img(row:row+nr-1,col:col+nc-1,:), 3 );
        testimg = testimg';
        testimg = testimg(:);
        testimg = testimg - mean(testimg);
        std_dev = std(testimg);
        % Ignore patches without much contrast (as measured by std deviation)
        if ( std_dev < 25 )
          continue;
        end;
        testimg = testimg / std_dev;
        testimg = testimg - MeanImg;
        % Project the patch into the eigenspace and reconstruct
        %dist = norm(PCs*(PCs'*testimg)-testimg);
        dist = norm((PCs'*testimg)./Stds);
        % Save this location/score
        nlocs = nlocs + 1;
        Locs(nlocs,1) = row;
        Locs(nlocs,2) = col;
        Locs(nlocs,3) = imgrows;
        Locs(nlocs,4) = imgcols;
        Locs(nlocs,5) = dist;
      end;
    end;
    fprintf(1,'\n');
    fflush(1);
  end;
  Img = imscale( Img, scalefactor );
  iscale = iscale + 1;
end;
Locs = Locs(1:nlocs,:);

%-------------------------------------------------------
%                     Functions
%-------------------------------------------------------

% Scale an image by the given factor with bilinear interpolation
function Imgout = imscale( Img, factor )
  % Get x/y locations of the desired pixels in source image
  nr = size( Img, 1 );
  nc = size( Img, 2 );
  np = size( Img, 3 );
  nrnew = round( nr*factor );
  ncnew = round( nc*factor );
  Imgout = zeros( nrnew, ncnew, np );
  for p = 1:np
      Imgout(:,:,p) = _imscale( Img(:,:,p), factor );
  end
      
function Img = _imscale( Img, factor )
  nr = size( Img, 1 );
  nc = size( Img, 2 );
  nrnew = round( nr*factor );
  ncnew = round( nc*factor );
  xs = nc/ncnew;
  ys = nr/nrnew;
  xlocs = xs * (1:ncnew) - (xs-1)/2;
  xlocs(find(xlocs<1)) = 1;
  xlocs(find(xlocs>nc)) = nc;
  xlocs = repmat(xlocs,nrnew,1);
  ylocs = ys * (1:nrnew)' - (ys-1)/2;
  ylocs(find(ylocs<1)) = 1;
  ylocs(find(ylocs>nr)) = nr;
  ylocs = repmat(ylocs,1,ncnew);
  % Get indices into the vectorized image
  xlocs = xlocs(:);
  ylocs = ylocs(:);
  inds = floor(xlocs-1)*nr+floor(ylocs);
  % Get left-, up-, and left-up-shifted versions of the image
  Il = zeros(nr,nc);
  Il(:,1:nc-1) = Img(:,2:nc);
  Il(:,nc) = Img(:,nc);
  Iu = zeros(nr,nc);
  Iu(1:nr-1,:) = Img(2:nr,:);
  Iu(nr,:) = Img(nr,:);
  Ilu = zeros(nr,nc);
  Ilu(1:nr-1,1:nc-1) = Img(2:nr,2:nc);
  Ilu(nr,:) = Img(nr,:);
  Ilu(:,nc) = Img(:,nc);
  % Get weights
  w1 = (floor(1+xlocs)-xlocs) .* (floor(1+ylocs)-ylocs);
  w2 = (1-(floor(1+xlocs)-xlocs)) .* (floor(1+ylocs)-ylocs);
  w3 = (floor(1+xlocs)-xlocs) .* (1-(floor(1+ylocs)-ylocs));
  w4 = (1-(floor(1+xlocs)-xlocs)) .* (1-(floor(1+ylocs)-ylocs));
  % Get weighted image
  Img = round( double(Img(inds)).*w1 + Il(inds).*w2 + ...
               Iu(inds).*w3 + Iu(inds).*w4 );
  Img = reshape(Img,nrnew,ncnew);

