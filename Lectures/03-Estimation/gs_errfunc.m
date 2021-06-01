function Y = gs_errfunc( X, P )

  % Check argument consistency, extract H and xhat from P,
  % and extract x1 and x2 from X.

  H = reshape( P(1:9), 3, 3 );
  xhat = P(10:end);
  n = size( xhat, 1 ) / 2;
  if size( X, 1 ) ~= 4*n
    error( 'Unexpected dimensionality for X' );
  end;
  xhat = [ reshape( xhat, 2, n ) ; ones( 1, n ) ];
  x1 = X(1:2*n);
  x2 = X(2*n+1:end);
  x1 = reshape( x1, 2, n );
  x2 = reshape( x2, 2, n );

  % Calculate H*xhat and the error between estimates and observations

  xhatp = H * xhat; xhatp = xhatp ./ repmat( xhatp(3,:), 3, 1 );
  error1 = x1 - xhat(1:2,:);
  error2 = x2 - xhatp(1:2,:);

  % Unroll them into the residual vector Y

  Y = [ error1(:) ; error2(:) ];

  %err = sum( sum( Y.^2 ))

