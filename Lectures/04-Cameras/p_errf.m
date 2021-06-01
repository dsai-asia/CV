function Y = p_errf( X, P )

  % Error function for nonlinear least squares estimation of camera matrix P

  n = size(X,1) / 2;
  P = reshape( P, 3, 4 );
  X = X(1:2:end,:)';
  Y = P * X;
  Y = Y(1:2,:) ./ repmat(Y(3,:),2,1);
  Y = Y(:);

