function [X,mu,Sigma,err] = ekf_demo()

% Simulate the EKF on a problem in stereo vision
% This version measures X directly by triangulation

% Set up and display simulation parameters

X = gen_sim_data();
[Tr,R,f] = cameras();
Z = get_obs( X, Tr, R, f );
figure(1);
plot_sim_data( X, Z, [], Tr, R );

% Initial given distribution

[mu{1}, Sigma{1}, Sigma_x] = initial_distribution();
err = mu{1}-X{1};
err = norm(err(1:3));
for T = 2:size(Z,2)
    muprev = mu{T-1};
    [pred_x, Jf] = sysf( muprev );
    [pred_z, H, Sigma_z] = sysh( pred_x, R, Tr, f );
    pred_sig = Jf * Sigma{T-1} * Jf' + Sigma_x;
    residual = Z{T} - pred_z;
    e(T)=norm(residual);
    K = pred_sig * H' * inv( H * pred_sig * H' + Sigma_z );
    Sigma{T} = ( eye( size(K,1) ) - K * H ) * pred_sig;
    munew = pred_x + K * residual;
    mu{T} = munew;
    xnew = X{T};
    err(T) = norm(munew(1:3)-xnew(1:3));
end;
figure(2); 
plot_sim_data( X, Z, mu, Tr, R );

% ------------------------------- FUNCTIONS -----------------------------

% System equation

function [xnew,J] = sysf( xold )
    x = xold(1);
    y = xold(2);
    z = xold(3);
    p = xold(4);
    w = xold(5);
    v = xold(6);
    delta_t = 1/30;
    xnew = [ x + delta_t * v * cos(p) * cos(w)
             y + delta_t * v * sin(p)
             z + delta_t * v * cos(p) * sin(w)
             p
             w
             v ];

    % Jacobian

    grad_xnew = [ 1, 0, 0, delta_t * v * cos(w) * -sin(p), ...
                  delta_t * v * cos(p) * -sin(w), ...
                  delta_t * cos(p) * cos(w) ];
    grad_ynew = [ 0, 1, 0, delta_t * v * cos(p), 0, delta_t * sin(p) ];
    grad_znew = [ 0, 0, 1, delta_t * v * sin(w) * -sin(p), ...
                  delta_t * v * cos(p) * cos(w), ...
                  delta_t * cos(p) * sin(w) ];

    J = [ grad_xnew
          grad_ynew
          grad_znew
          0 0 0 1 0 0
          0 0 0 0 1 0
          0 0 0 0 0 1 ];


% Measurement equation

function [Z,H,Sigma_z] = sysh( X, R, T, f )
    H = [ eye(3), zeros(3) ];
    Z = H * X;
    Sigma_img = diag( [ 2 2 2 2 ].^2 );
    [imgobs, Jac_img] = old_sysh( X, R, T, f );
    Jac_img = Jac_img(:,1:3);
    % See "Backward Transport of Covariance", Hartley and Zisserman p. 127
    Sigma_z = inv( Jac_img' * inv( Sigma_img ) * Jac_img );


% Test correctness of jacobian of f using finite differences

function test_jacf( X, Jf )
    FX = sysf( X );
    diff = 1e-6;
    for i = 1:size(X,1)
        X(i) = X(i) - diff;
        FX1 = sysf( X );
        X(i) = X(i) + 2*diff;
        FX2 = sysf( X );
        X(i) = X(i) - diff;
        J(1:size(X,1),i) = ( FX2 - FX1 ) / ( diff * 2 );
    end;
    Jf
    J


% Simulate an object traveling in the X-Z plane in a circle with
% constant acceleration

function X = gen_sim_data()
    X{1} = [ 0 0 0 0 0 0 ]';
    % The object takes about 10 sec to travel around a circle in XZ plane
    t = 1;
    while 1
        t = t+1;
        prevX = X{t-1};         % Last state
        v = prevX(6);           % Last velocity
        v = v + 0.01;           % Accelerate a bit
        dist = v/30;            % Sampling at 30 Hz
        travel = atan2( prevX(1), 1-prevX(3) );
        travel = travel + dist;
        newX = [ sin(travel) 0 1-cos(travel) 0 travel v ]';
        X{t} = newX;
        if newX(1) >= 0 && prevX(1) < 0
            break;
        end;
    end;


% Plot data X and camera positions T

function plot_sim_data( X, Z, Xest, T, R )
    for i = 1:size(X,2)
        thisx = X{i};
        plotx(i) = thisx(1);
        ploty(i) = thisx(2);
        plotz(i) = thisx(3);
    end;
    hold off;
    plot( plotz, plotx, '@15' );
    hold on;
    for i = 1:size(X,2)
        thisx = Z{i};
        plotx(i) = thisx(1);
        ploty(i) = thisx(2);
        plotz(i) = thisx(3);
    end;
    plot( plotz, plotx, '@44' );
    if ~isempty( Xest )
        for i = 1:size(X,2)
            thisx = Xest{i};
            plotx(i) = thisx(1);
            ploty(i) = thisx(2);
            plotz(i) = thisx(3);
        end;
        plot( plotz, plotx, '21' );
    end;
    for i = 1:size(T,2)
        this = R{i}'*-T{i};
        plot( this(3), this(1), '@21' );
    end;
    axis([-0.75 4.2 -1.2 1.2],"equal");
    hold off;


% Generate rotation, translation, and focal lengths for two cameras

function [T,R,f] = cameras
    R{1} = [ 1 0 0 ; 0 1 0 ; 0 0 1 ];
    T{1} = [ 0.2 0 -4 ]';
    R{2} = [ 1 0 0 ; 0 1 0 ; 0 0 1 ];
    T{2} = [ -0.2 0 -4 ]';
    f = 320/tan(30/180*pi);         % 60 degree horizontal field of view


% Return the prior state distribution

function [mu, Sigma, Sigma_x] = initial_distribution
    mu = [ 0 0 0 0 0 0 ]';
    Sigma = diag( [ 0.1 0.1 0.1 10/180*pi 10/180*pi 0.10].^2 );
    deltat = 1/30;
    Sigma_x = diag( ( deltat * [ 0.1 0.1 0.1 pi pi 0.5 ] ).^2 );


% Return simulated observations

function Z = get_obs( X, T, R, f )
    for i = 1:size(X,2)
        xthis = X{i};
        ptdir = [];
        for j = 1:size(T,2)
            trans = T{j};
            rot = R{j};
            camobs = rot*xthis(1:3)+trans;
            camobs = [-f 0 0; 0 -f 0; 0 0 1]*camobs;
            camobs = camobs(1:2) / camobs(3);
            camobs = round( camobs + randn(2,1) * 2 );
            % Get a point and direction in world coords
            world1 = rot' * ( [ camobs/f ; -1 ] - trans );
            world2 = rot' * ( - trans );
            dir = world1 - world2;
            dir = dir / norm( dir );
            ptdir = [ ptdir, [ world2 ; dir ] ];
        end;
        % Triangulate
        p1 = ptdir(1:3,1);
        p2 = ptdir(1:3,2);
        d1 = ptdir(4:6,1);
        d2 = ptdir(4:6,2);
        tnear = -((d2'*d1 )*((p1-p2)'*d1)+(p2-p1)'*d2) / (1-(d2'*d1)^2);
        p2near = p2 + d2 * tnear;
        proj = (p2near-p1)'*d1 * d1 + p1;
        mid = mean([p2near,proj],2);
        Z{i} = mid;
    end;


% Measurement equation

function [z, Jh] = old_sysh( State, Rot, Trans, f )
    z = [];
    Jh = [];
    for j = 1:size(Rot,2)
        R = Rot{j};
        T = Trans{j};
        obs = [ ( R(1,:) * State(1:3) + T(1) ) * -f / ...
                  ( R(3,:) * State(1:3) + T(3) ) ;
                ( R(2,:) * State(1:3) + T(2) ) * -f / ...
                  ( R(3,:) * State(1:3) + T(3) ) ];
        z = [ z ; obs ];

        % Jacobian

        grad_xobs = [
           ( ( R(3,:) * State(1:3) + T(3) ) * R(1,1) * -f - ...
             ( R(1,:) * State(1:3) + T(1) ) * -f * R(3,1) ) / ...
           ( R(3,:) * State(1:3) + T(3) )^2, ...
           ( ( R(3,:) * State(1:3) + T(3) ) * R(1,2) * -f - ...
             ( R(1,:) * State(1:3) + T(1) ) * -f * R(3,2) ) / ...
           ( R(3,:) * State(1:3) + T(3) )^2, ...
           ( ( R(3,:) * State(1:3) + T(3) ) * R(1,3) * -f - ...
             ( R(1,:) * State(1:3) + T(1) ) * -f * R(3,3) ) / ...
           ( R(3,:) * State(1:3) + T(3) )^2, ...
           0 0 0 ];

        grad_yobs = [
           ( ( R(3,:) * State(1:3) + T(3) ) * R(2,1) * -f - ...
             ( R(2,:) * State(1:3) + T(2) ) * -f * R(3,1) ) / ...
           ( R(3,:) * State(1:3) + T(3) )^2, ...
           ( ( R(3,:) * State(1:3) + T(3) ) * R(2,2) * -f - ...
             ( R(2,:) * State(1:3) + T(2) ) * -f * R(3,2) ) / ...
           ( R(3,:) * State(1:3) + T(3) )^2, ...
           ( ( R(3,:) * State(1:3) + T(3) ) * R(2,3) * -f - ...
             ( R(2,:) * State(1:3) + T(2) ) * -f * R(3,3) ) / ...
           ( R(3,:) * State(1:3) + T(3) )^2, ...
           0 0 0 ];

       Jh = [ Jh ; grad_xobs ; grad_yobs ];

    end;


