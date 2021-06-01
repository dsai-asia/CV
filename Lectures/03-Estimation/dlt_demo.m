function [Htrue, H, Hnorm] = dlt_demo

% This function is another tutorial on how to do computer vision
% type stuff with Octave/Matlab.  It shows how to implement
% Hartley and Zisserman's "DLT" algorithm for estimating a homography
% from point correspondences.  It also shows how to test an algorithm
% using randomly-generated synthetic data.  When implementing algorithms
% you should ALWAYS test the algorithm on synthetic data you know the
% right answer for.  If you don't, you might never catch bugs in your
% implementation.

% This function is written for Octave.  You will have to tweak some of
% the function calls to get it running on Matlab.  Mostly just the plot
% functions and the random number generator seed functions.

% Matt Dailey, June 2007

% Initalize the random number generator with a fixed seed to get
% repeatable results.  This is for Octave; see "help rand" for
% corresponding Matlab call

rand( 'seed', 0 );
randn( 'seed', 10 );

% Let's generate 12 points uniformly distributed over image 1

n = 12;
w = 640;
h = 480;
x = [ w 0 ; 0 h ] * rand( 2, n );

% Generate a random homography

Htrue = rand_homography( x, w, h );

xp = Htrue * [ x; ones( 1, n ) ];
xp = xp ./ repmat( xp(3,:), 3, 1 );

% Add Gaussian noise and round off to generate the measurements

xmeas = [ round( x + randn( 2, n )) ; ones( 1, n ) ];
xpmeas = [ round( xp(1:2,:) + randn( 2, n )); ones( 1, n ) ];

% Plot the measured data.  axis() arguments might have to change for Matlab

figure(1); plot( xmeas(1,:), xmeas(2,:), 'ro' );
title( 'Points in image 1' ); axis([1 640 1 480], 'ij' );
figure(2); plot( xpmeas(1,:), xpmeas(2,:), 'ro' );
title( 'Points in image 2' ); axis([1 640 1 480], 'ij' );
hold on;
for i = 1:12
    plot( [xmeas(1,i),xpmeas(1,i)], [xmeas(2,i),xpmeas(2,i)], 'g-' );
end;
hold off;

% Run the DLT without normalization

H = dlt( xmeas, xpmeas );

% Use H to predict xp from x and plot on image 2

xpp = H * xmeas; xpp = xpp ./ repmat( xpp(3,:), 3, 1 );
figure(2); hold on;
plot( xpp(1,:), xpp(2,:), 'bx' );
hold off;

% Run the DLT with normalization

Hnorm = dlt_norm( xmeas, xpmeas );

% Use Hnorm to predict xp from x and plot on image 2

xppn = Hnorm * xmeas; xppn = xppn ./ repmat( xppn(3,:), 3, 1 );
figure(2); hold on;
plot( xppn(1,:), xppn(2,:), 'kx' );
hold off;

% Calculate one-way transfer error for both estimates of H

err = sum( sum( ( xpp(1:2,:)-xpmeas(1:2,:) ).^2 ));
err_norm = sum( sum( ( xppn(1:2,:)-xpmeas(1:2,:) ).^2 ));

fprintf( 1, 'Error %f for basic DLT and %f for normalized DLT.\n', ...
         err, err_norm );

% Usually if there's not much perspective distortion, the unnormalized
% algorithm works reasonably well.  To see how well it does under more
% extreme perspective distortions, try commenting out line 172 below
% and see how the error changes.

%------------------------------------------------------------------------

function H = dlt( x, xp )

  n = size( x, 2 );
  if n < 4
    error( 'DLT requires at least 4 points' );
  end;
  if ( size( x, 1 ) ~= 3 | size( xp, 1 ) ~= 3 )
    error( 'DLT requres homogeneous coordinates' );
  end;

  A = [];

  for i = 1:n

    xip = xp( 1, i );
    yip = xp( 2, i );
    wip = xp( 3, i );

    xi = x( :, i );

    Ai = [ 0, 0, 0,    -wip * xi',   yip * xi' ;
           wip * xi',     0, 0, 0,  -xip * xi' ];

    A = [ A ; Ai ];
  end;

  [U,D,V] = svd( A );

  % In Octave, the SVD is sorted with decreasing singular values
  % so we want the last column of V

  H = reshape( V(:,9), 3, 3 )';
  H = H / H(3,3);

%------------------------------------------------------------------------

function H = dlt_norm( x, xp )

  T = normalize_transform( x );
  Tp = normalize_transform( xp );

  Htilde = dlt( T * x, Tp * xp );
  H = inv( Tp ) * Htilde * T;
  H = H / H(3,3);

%------------------------------------------------------------------------

function T = normalize_transform( x )

  % Transform taking x's centroid to the origin

  Ttrans = [ 1 0 -mean( x(1,:) ) ; 0 1 -mean( x(2,:) ) ; 0 0 1 ];

  % Calculate appropriate scaling factor

  x = Ttrans * x;
  lengths = sqrt( sum( x(1:2,:).^2 ));
  s = sqrt(2) / mean(lengths);

  % Transform scaling x to an average length of sqrt(2)

  Tscale = [ s 0 0 ; 0 s 0 ; 0 0 1 ];

  % Compose the transforms

  T = Tscale * Ttrans;

%------------------------------------------------------------------------

function H = rand_homography( x, w, h )

  % Start with a completely random homography

  H = rand( 3, 3 );

  % Make sure it's not too much perspective warping

  H = H / H(3,3);
  if norm( H(3,1:2) ) > 0.01
    %H(3,1:2) = H(3,1:2) / norm( H(3,1:2) ) * 0.01;  % Comment this out for fun
  end; 

  % Figure out where that takes our points x

  n = size( x, 2 );
  xp = H * [ x; ones( 1, n ) ];
  xp = xp ./ repmat( xp(3,:), 3, 1 );

  % Scale so that the points are spread nicely

  xspread = max( xp(1,:) ) - min( xp(1,:) );
  yspread = max( xp(2,:) ) - min( xp(2,:) );
  Hscale = [ 0.9 * 1/xspread * w, 0, 0 ; 0, 0.9 * 1/yspread * h, 0 ; 0 0 1 ];
  xp = Hscale * xp;
  xp = xp ./ repmat( xp(3,:), 3, 1 );

  % Translate so that the points are inside the image

  xmin = min( xp(1,:) );
  ymin = min( xp(2,:) );
  Htrans = [ 1, 0, -xmin + 0.05 * w; 0 1 -ymin + 0.05 * h; 0 0 1 ];

  % Compose the transforms

  H = Htrans * Hscale * H;
  H = H / H(3,3);

