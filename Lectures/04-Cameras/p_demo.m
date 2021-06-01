function p_demo

  % Demonstration of computing a camera matrix P and decomposing
  % it as K [R t] in Octave/Matlab

  pkg load optim;

  rand('seed',12);
  randn('seed',13);

  % Get some random synthetic data

  [Xw,P,K,R,t,x,avg_dist_to_cam] = get_synthetic_data();
  plot_img_data( x, 1 );

  % Get the camera matrix

  Pest = compute_p( x, Xw );
  [Cest,Test,Rest,Kest] = decompose_p( Pest );

  % Plot the world points, actual camera, and estimated camera from above

  plot_camera_pts( R, t, Rest, Test, Xw, 2 );

  % Calculate camera center error

  err_c = norm( -Rest'*Test - (-R'*t) );
  fprintf( 1, 'Camera center error %g (%g relative to cam dist)\n', ...
           err_c, mean( err_c ) / avg_dist_to_cam );
end;

%-----------------------------------------------------------------------------

function R = rot_y( theta )

  ct = cos(theta); st = sin(theta);
  R = [ ct 0 -st ; 0 1 0 ; st 0 ct ];

end;

%-----------------------------------------------------------------------------

function [Xw,P,K,R,t,x,avg_dist_to_cam] = get_synthetic_data

  % Random data around (0,0,5)' in 3D

  zmean = 5;
  muX = [0;0;zmean];
  n = 30;
  Xw = [ randn(3,n) + repmat(muX,1,n); ones(1,n) ];

  % Actual intrinsic parameters

  K = [ 320 0 320 ; 0 320 240 ; 0 0 1 ];

  % Actual extrinsic parameters

  R = eye(3);
  C = [0;0;0];
  t = -R * C;
  P = K * [ R t ];

  % Calculate ideal projections into the image

  x = P * Xw; x = x ./ repmat( x(3,:), 3, 1 );

  % Add Gaussian noise and discretize

  x(1:2,:) = round( x(1:2,:) + randn(2,n) );

  avg_dist_to_cam = mean( sqrt( sum(( Xw(1:3,:) - repmat( C, 1, n )).^2)));

end;

%-----------------------------------------------------------------------------

function err = get_reproj_error( x, P, Xw );

  xt = P * Xw;
  xt = xt ./ repmat( xt(3,:), 3, 1 );
  err = sqrt( mean( sum( ( x(1:2,:) - xt(1:2,:) ).^2 )));

end;

%-----------------------------------------------------------------------------

function plot_img_data( x, fignum )
  plot(x(1,:),x(2,:),'r+');
  axis([0 640 0 480],'ij');
  title(sprintf('Image data'));

end;

%-----------------------------------------------------------------------------
