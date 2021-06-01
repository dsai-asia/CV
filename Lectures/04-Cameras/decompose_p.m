function [C,T,R,K] = decompose_p( P )

  % Decompose a camera matrix into K, R, and t

  % Rassarin Chinnachodteeranun, May 2007
  % Modified by Matt Dailey, June 2007

  % Camera center is the right null space of P

  [U,S,V]=svd(P);
  C = V(1:3,4)/V(4,4);

  % K and R are obtained using the RQ factorization

  [K,R] = rq(P(:,1:3));
  K = K/K(3,3);

  % If focal lengths come out negative in K, fix them

  fix_t = eye(3);
  if K(1,1) < 0
    fix_t(1,1) = -1;
  end
  if K(2,2) < 0
    fix_t(2,2) = -1;
  end
  K = K * fix_t;
  R = fix_t * R;

  % If R is oriented backwards, fix it

  if det(R) < 0
    R = -R;
  end;

  % Translation vector is -R*C

  T = -R*C;

end

