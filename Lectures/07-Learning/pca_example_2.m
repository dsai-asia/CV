
n = 500;
Data = [ 10, 0 ; 0 , 3 ] * randn( 2, n );

p1 = [ 63/12*30 ; 107/2.2 ];
p2 = [ 72/12*30 ; 170/2.2 ];
mu = mean( [ p1' ; p2' ] )';
theta = atan2( p2(2)-mu(2), p2(1)-mu(1) );
R = [ sin( theta ) , -cos( theta ) ; cos( theta ) , sin( theta ) ];
Data = R * Data;
[U,V] = eig( Data*Data' );
Data = Data + repmat( mu, 1, n );
X1 = mu;
X2 = 20*U(:,2);
if X2(1) < 0, X2 = -X2; end;
X2 = X2 + mu;
apts = [ 18 , 18 ; 1, -1 ];
apts = R * apts + repmat( mu, 1, 2 );

hold off;
plot( Data(1,:), Data(2,:), 'rx' );
xlabel('Height (cm)');
ylabel('Weight (kg)');
title('Height vs. weight in a sample of 500 people');
axis('equal');
hold on;
plot( [X1(1),X2(1)],[X1(2),X2(2)],'b-' );
plot( [X2(1),apts(1)],[X2(2),apts(2)],'b-' );
plot( [X2(1),apts(3)],[X2(2),apts(4)],'b-' );

print( 'height-weight-2.eps', '-deps', '-color' );

