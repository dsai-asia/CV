function P = compute_p( x, Xw )

  % Estimate the camera projection P best fitting
  % 3D points X and observed image points W

  % Matt Dailey, June 2007

  % Check args

  if nargin ~= 2
    error( 'Usage: compute_p( x, Xw )' );
  end;
  n = size( Xw, 2 );
  if size( x, 2 ) ~= n
    error( 'x and Xw should contain the same number of points' );
  end;
  if size( Xw, 1 ) ~= 4
    error( 'Xw should be 4xn' );
  end;
  if size( x, 1 ) ~= 3
    error( 'x should be 3xn' );
  end;

  % Normalize homogeneous coords

  Xw = Xw ./ repmat( Xw(4,:), 4, 1 );
  x = x ./ repmat( x(3,:), 3, 1 );

  % Compute isotropic scaling transforms to normalize x and Xw

  T = normalize_transform( x );
  U = normalize_transform( Xw );

  xtilde = T * x;
  Xwtilde = U * Xw;

  % Get the linear estimate of P

  Plin = dlt( xtilde, Xwtilde );
  fprintf( 1, 'Reprojection error from linear estimate of P: %g\n', ...
           reprojection_error( x, Xw, inv(T)*Plin*U ));

  % Get the nonlinear estimate of P

  Xobs = reshape( repmat(Xwtilde,2,1), 4, 2*n )'; % Independent variables
  Ydes = reshape( xtilde(1:2,:), 2*n, 1 );        % Desired outputs
  P0 = Plin(:);                                   % Parameters
  Yact = p_errf(Xobs,P0);
  %fprintf( 1, 'Error from linear estimate %g\n', sum( (Ydes-Yact).^2 ));
  [f,P1] = leasqr( Xobs, Ydes, P0, 'p_errf' );
  %fprintf( 1, 'Error from nonlinear estimate %g\n', sum( (Ydes-f).^2 ));
  Pnonlin = reshape( P1, 3, 4 );

  % Denormalize

  P = inv( T ) * Pnonlin * U;
  fprintf( 1, 'Reprojection error from nonlinear estimate of P: %g\n', ...
           reprojection_error( x, Xw, P ));

%------------------------------------------------------------------------

function P = dlt( x, Xw )

  n = size( Xw, 2 );
  if n < 6
    error( 'DLT for P requires at least 6 points' );
  end;

  A = [];

  for i = 1:n

    xi = x( 1, i );
    yi = x( 2, i );
    wi = x( 3, i );
    Xwi = Xw( :, i );

    Ai = [ 0, 0, 0, 0,   -wi * Xwi',    yi * Xwi' ;
            wi * Xwi',   0, 0, 0, 0,   -xi * Xwi' ];

    A = [ A ; Ai ];
  end;

  [U,D,V] = svd( A );

  % The solution P minimizing |Ap| subject to |p|=1 is the last column of V

  P = reshape( V(:,12), 4, 3 )';
  P = P / norm(P(3,1:3));

%------------------------------------------------------------------------

function err = reprojection_error( x, Xw, P )

  xtilde = P * Xw;
  xtilde = xtilde ./ repmat( xtilde(3,:), 3, 1 );
  err = sqrt( mean( sum( ( x(1:2,:) - xtilde(1:2,:) ).^2 )));

%------------------------------------------------------------------------

function T = normalize_transform( X )

  % Are we working on 2D or 3D points?

  k = size(X,1)-1;

  % Transform taking X's centroid to the origin

  Ttrans = [ eye( k ), -mean( X(1:k,:)' )' ; zeros( 1, k ), 1 ];

  % Calculate appropriate scaling factor

  X = Ttrans * X;
  lengths = sqrt( sum( X(1:k,:).^2 ));
  s = sqrt(k) / mean(lengths);

  % Transform scaling x to an average length of sqrt(2)

  Tscale = [ s * eye(k), zeros(k,1) ; zeros(1,k), 1 ];

  % Compose the transforms

  T = Tscale * Ttrans;

%------------------------------------------------------------------------

