
% Simple example of computing a homography from 4 exact point correspondences

% Matt Dailey

% The points in the left image

X = [ 200 250 210 255 ;
      150 145 170 162 ;
        1   1   1   1 ];

% The points in the right image

Xprime = [ 200 250 200 250 ;
           150 150 170 170 ;
             1   1   1   1 ];

% We have Xprime = H * X

% Let's loop over the points and extract the linear equations in H

A = [];

for iPoint = 1:size(X,2)
    x = X(:,iPoint);
    xprime = Xprime(:,iPoint);

    % The equation is
    % xprime(1) = ( h11 x(1) + h12 x(2) + h13 ) / ( h31 x(1) + h32 x(2) + h33 )

    % Put this in the form ( coeffs * (h11 h12 h13 h21 h22 h23 h31 h32 h33)' ):

    coeffs = ...
    [ -x(1) -x(2) -1 0 0 0 xprime(1) * x(1) xprime(1) * x(2) xprime(1) ];

    A = [ A ; coeffs ];

    % Do the same for the equation
    % xprime(2) = ( h21 x(1) + h22 x(2) + h23 ) / ( h31 x(1) + h32 x(2) + h33 )

    coeffs = ...
    [ 0 0 0 -x(1) -x(2) -1 xprime(2) * x(1) xprime(2) * x(2) xprime(2) ];

    A = [ A ; coeffs ];

end

% The solution is the right null space of A

h = null(A);

% Form the 9 parameters of H into a matrix

H = reshape( h, 3, 3 )';

% Generate 100 random points

Xtest = [ ( rand(2,100)' * diag([640,480]) )' ; ones(1,100) ];

% Plot them

figure(1);
plot( Xtest(1,:), Xtest(2,:), 'rx' );
hold on;
plot(X(1,:),X(2,:),'bo');
axis([0,640,0,480],'ij');
title('Left image points X');
hold off;

% Transfer the test points

Xtestprime = H * Xtest;
Xtestprime = Xtestprime ./ repmat(Xtestprime(3,:),3,1);

Xprimeest = H * X;
Xprimeest = Xprimeest ./ repmat(Xprimeest(3,:),3,1);

% Plot them

figure(2);
plot( Xtestprime(1,:), Xtestprime(2,:), 'rx' );
hold on;
plot(Xprimeest(1,:),Xprimeest(2,:),'bo');
axis([0,640,0,480],'ij');
title('Right image points X');
hold off;

