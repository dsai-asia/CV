function plot_camera_pts( R1, t1, R2, t2, Xw, fignum )

  % Plot a camera and set of points

  % Matt Dailey, June 2007

  % Uses the help of EGT toolbox's f_3Dcamera camera drawing function

  % The f_3Dcamera function takes a 3D homography, expected
  % to be a rigid transform, with the translation vector
  % being the camera center in world coordinates and the
  % rotation being the camera-to-world rotation.

  % Convert R, t to the form required by f_3Dcamera

  Xmax = max( Xw(1:3,:)' );
  Xmin = min( Xw(1:3,:)' );
  R1camtow = R1';
  T1caminw = -R1'*t1;
  R2camtow = R2';
  T2caminw = -R2'*t2;
  Xmax = max( [ Xmax ; T2caminw' ; T1caminw' ] );
  Xmin = min( [ Xmin ; T2caminw' ; T1caminw' ] );

  % Try to adjust the axes

  Xwidth = Xmax - Xmin;
  maxw = max( Xwidth );
  camscale = maxw/20;
  Xmid = ( Xmin + Xmax ) / 2;
  Xmin = floor( Xmid - repmat( 1.2*(maxw/2), 1, 3 ));
  Xmax = ceil( Xmid + repmat( 1.2*(maxw/2), 1, 3 ));

  % Plot the data and cameras

  figure( fignum );
  clf;
  axis([Xmin(1) Xmax(1) Xmin(2) Xmax(2) Xmin(3) Xmax(3)]);
  hold on;
  title( 'Reconstructed camera' );
  f_3Dcamera( [ R1camtow, T1caminw ; 0 0 0 1 ], 'r', camscale );
  f_3Dcamera( [ R2camtow, T2caminw ; 0 0 0 1 ], 'g', camscale );
  plot3(Xw(1,:),Xw(2,:),Xw(3,:),'r+');
  view( 0, 0 );
  hold off;

