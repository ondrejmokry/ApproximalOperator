%% How does the function phi look for g = |x| ?

% Date: 18/08/2020
% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@mensa.cz

clear
rng(2)

% global parameters
redundancy = 2; % redundancy of the transform
alpha      = 2; % constant of the tight frame
dim        = [ redundancy*2, 2 ]; % dimenions of the spaces
B          = rand(dim(1)); % a basis of the first space

% make B orthogonal again!
B(:,1) = B(:,1)/norm(B(:,1));
for j = 2:dim(1)
    ortogonalizator = 0;
    for k = 1:j-1
        ortogonalizator = ortogonalizator + (B(:,k)'*B(:,j))*B(:,k);
    end
    B(:,j) = B(:,j) - ortogonalizator;
    B(:,j) = B(:,j)/norm(B(:,j));
end

% projecting the basis onto R^dim(2) and incorporating alpha
SYN = sqrt(alpha)*B(1:dim(2),:); % synthesis matrix
ANA = SYN'; % analysis matrix

%% the functions to be compared

% function f
f = @(x) norm(ANA*x, 1);

% function phi
phi = @(x) func(x, dim(1), alpha, ANA);

%% computing the function values
values = -2:0.2:2; % parameter for the XY grid to plot the function

% compute the values
[ X, Y ] = meshgrid(values, values);
PHI = NaN(size(X));
F = NaN(size(X));
for i = 1:length(values)
    for j = 1:length(values)
        PHI(i,j) = phi([X(i,j); Y(i,j)]);
        F(i,j)   = f  ([X(i,j); Y(i,j)]);
    end
end

%% plot

% plot the frame atoms
figure
for i = 1:dim(1)
    line([0 SYN(1,i)], [0 SYN(2,i)],'color','b')
    direction = SYN(:,i)/norm(SYN(:,i));
    text(SYN(1,i) + 0.1*direction(1),...
         SYN(2,i) + 0.1*direction(2),...
         ['generator \#',num2str(i)])
end
xlim([-1.5 1.5])
ylim([-1.5 1.5])
xlabel('$$x[1]$$')
ylabel('$$x[2]$$')
grid on
box on
axis square

% plot the 3D figure
figure
hold on
surf(X,Y,PHI,'FaceColor','r','FaceAlpha',0.5)
surf(X,Y,F,'FaceColor','b','FaceAlpha',0.5)
xlabel('$$x[1]$$')
ylabel('$$x[2]$$')
legend('$$\phi(\mathbf{x})$$','$$f(\mathbf{x}) = g(A\mathbf{x})$$')
view(3)
grid on

%% infimal convolution
function value = func(x, dimension, alpha, T)
    Tx = T*x;
    cvx_begin quiet
        cvx_precision best
        variable y(dimension)
        minimize( alpha*norm( y, 1 ) + 1/2 *( Tx-y )'*( Tx-y ) )
        subject to
            T'*(Tx - y) == 0
    cvx_end
    value = cvx_optval/alpha;
end