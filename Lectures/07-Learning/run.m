
% PCA training and testing for face detection

test_filename = 'Jigme.png';

% Train on 300 example face images in Examples

fprintf(1,'Training...\n');
fflush(1);
[MeanImg,PCs_all,Vars] = pca_train();
nrpatch = size(MeanImg,1);
ncpatch = size(MeanImg,2);

% Pick a number k of eigenfeatures to use.  You have to experiment with
% this.  If too big, you're just paying attention to the details of the
% training set.  If too small, you won't be fully modeling the set of
% possible face appearances.

k = 10;

% Test on the test image

fprintf(1,'Testing...\n');
fflush(1);
[Locs] = pca_test( test_filename, MeanImg, PCs_all(:,1:k), Vars );

% Sort locations by distance to the eigenspace

[junk, rindices] = sort( Locs(:,5) );
Locs = Locs(rindices,:);

% Convert best location in the scaled image to a location in the original image

Img = imread( test_filename );
if size(Img,3) == 3
  Img = rgb2gray( Img );
end
nr = size( Img, 1 );
nc = size( Img, 2 );

nloc=20;
if size( Img, 3 ) == 3
    R=Img(:,:,1);
    B=Img(:,:,2);
    G=Img(:,:,3);
else
    R=Img;
    B=Img;
    G=Img;
end
for loc = 1:nloc
    besty = Locs(loc,1);
    bestx = Locs(loc,2);
    bestys = Locs(loc,3);
    bestxs = Locs(loc,4);
    bestdist = Locs(loc,5);
    bestx1 = floor( (bestx-1) / bestxs * nc );
    besty1 = floor( (besty-1) / bestys * nr );
    bestx2 = ceil( bestx1 + ncpatch / bestxs * nc );
    besty2 = ceil( besty1 + nrpatch / bestys * nr );
    if ( loc==1 )
        fprintf( 1, 'Best location is (%d,%d) to (%d,%d) (dist=%f)\n', ...
                 bestx1, besty1, bestx2, besty2, bestdist );
    end;

    if bestx1 < 1, bestx1 = 1; end;
    R(besty1,bestx1:bestx2)=255;
    R(besty2,bestx1:bestx2)=255;
    R(besty1:besty2,bestx1)=255;
    R(besty1:besty2,bestx2)=255;
    G(besty1,bestx1:bestx2)=0;
    G(besty2,bestx1:bestx2)=0;
    G(besty1:besty2,bestx1)=0;
    G(besty2:besty2,bestx2)=0;
    B(besty1,bestx1:bestx2)=0;
    B(besty2,bestx1:bestx2)=0;
    B(besty1:besty2,bestx1)=0;
    B(besty2:besty2,bestx2)=0;

end;

% Display the image and the patch closest to the eigenspace

RGBImg = uint8( zeros( nr, nc, 3 ));
RGBImg(:,:,1) = R;
RGBImg(:,:,2) = G;
RGBImg(:,:,3) = B;
imshow(RGBImg);

