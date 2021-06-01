
n = 500;
Data = [ 10, 0 ; 0 , 3 ] * randn( 2, n );

p1 = [ 63/12*30 ; 107/2.2 ];
p2 = [ 72/12*30 ; 170/2.2 ];
mu = mean( [ p1' ; p2' ] )';
theta = atan2( p2(2)-mu(2), p2(1)-mu(1) );
Data = [ sin( theta ) , -cos( theta ) ; cos( theta ) , sin( theta ) ] * Data;
Data = Data + repmat( mu, 1, n );

plot( Data(1,:), Data(2,:), 'rx' );
xlabel('Height (cm)');
ylabel('Weight (kg)');
title('Height vs. weight in a sample of 500 people');

print( 'height-weight.eps', '-deps', '-color' );

