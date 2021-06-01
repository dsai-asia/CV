function runsvm

% Generate data

mu = cell(2,2);
covar = cell(2,2);

mu{1,1} = [ 1; 1 ];
covar{1,1} = [ .1 0 ; 0 .1 ];
mu{1,2} = [ 2; 2 ];
covar{1,2} = [ .1 0 ; 0 .1 ];

mu{2,1} = [ 1; 2 ];
covar{2,1} = [ .1 0 ; 0 .1 ];
mu{2,2} = [ 2; 1 ];
covar{2,2} = [ .1 0 ; 0 .1 ];

N = 100;
Xtrain = cell(1,2);
Xtest = cell(1,2);
Xtrain{1} = [ myrandgauss(mu{1,1},covar{1,1},N),...
              myrandgauss(mu{1,2},covar{1,2},N) ];
Xtrain{2} = [ myrandgauss(mu{2,1},covar{2,1},N),...
              myrandgauss(mu{2,2},covar{2,2},N) ];
Xtest{1} = [ myrandgauss(mu{1,1},covar{1,1},N),...
             myrandgauss(mu{1,2},covar{1,2},N) ];
Xtest{2} = [ myrandgauss(mu{2,1},covar{2,1},N),...
             myrandgauss(mu{2,2},covar{2,2},N) ];

% Plot training data

figure(1);
plot( Xtrain{1}(1,:), Xtrain{1}(2,:), 'rx' );
hold on;
plot( Xtrain{2}(1,:), Xtrain{2}(2,:), 'gx' );
title('Training data');
hold off;


% Create training data file

fp = fopen( 'svmtrain.txt', 'w' );
for i = 1:size(Xtrain{1},2)
    fprintf( fp, '1 1:%g 2:%g\n', Xtrain{1}(:,i) );
end;
for i = 1:size(Xtrain{2},2)
    fprintf( fp, '2 1:%g 2:%g\n', Xtrain{2}(:,i) );
end;
fclose( fp );

% Create test data file

fp = fopen( 'svmtest.txt', 'w' );
for i = 1:size(Xtest{1},2)
    fprintf( fp, '1 1:%g 2:%g\n', Xtest{1}(:,i) );
end;
for i = 1:size(Xtest{2},2)
    fprintf( fp, '2 1:%g 2:%g\n', Xtest{2}(:,i) );
end;
fclose( fp );

% Run svmtrain with nu-SVC, radial basis function kernel, default gamma, c

fprintf( 1, 'Training SVM...\n' ); fflush( 1 );
system( 'svm-train -s 1 -t 2 svmtrain.txt svmmodel.txt' );
fprintf( 1, 'Testing SVM...\n' ); fflush( 1 );
system( 'svm-predict svmtest.txt svmmodel.txt svmtest.out.txt' );

% Read classification results

fp = fopen( 'svmtest.out.txt', 'r' );
Xtest_pred = fscanf( fp, '%d' );
fclose( fp );

% Plot classification results

figure(2);
plot( Xtest{1}(1,:), Xtest{1}(2,:), 'rx' );
hold on;
plot( Xtest{2}(1,:), Xtest{2}(2,:), 'gx' );
title('Test data');
pred1 = Xtest_pred(1:2*N);
pred2 = Xtest_pred(2*N+1:end);
err1 = find(pred1==2);
err2 = find(pred2==1);
plot(Xtest{1}(1,err1), Xtest{1}(2,err1), 'ro' );
plot(Xtest{2}(1,err2), Xtest{2}(2,err2), 'go' );
hold off;


% ---------------------------------------------------------------------------

function X = myrandgauss( mu, covar, N )

    X = chol( covar )' * randn( 2, N ) + repmat( mu, 1, N );

