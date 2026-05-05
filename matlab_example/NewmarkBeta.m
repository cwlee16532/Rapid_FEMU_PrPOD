%
%   Time-integration: Newmart-Beta Method
%
function [GU,GV,GA] = NewmarkBeta(GMM,GDM,GSM,GU,GV,GA,GF,Nmk)

t  = Nmk.t;
dt = Nmk.dt;

if isempty(GDM)
    GDM = 0;%sparse(zeros(size(GSM)));
end


GA(:,1) = GMM\(GF(:,1)-GSM*GU(:,1)-GDM*GV(:,1));
% GA(:,1) = zeros(size(GMM,1),1); % because GF(:,1) = GU(:,1) = GV(:,1) = 0

de = 1/2;
al = 1/4;

a(8) = 1/(al*dt^2);
a(1) = de/(al*dt);
a(2) = 1/(al*dt);
a(3) = 1/(2*al) - 1;
a(4) = de/al - 1;
a(5) = dt/2*(de/al - 2);
a(6) = dt*(1-de);
a(7) = de*dt;

GSMp = GSM + a(8)*GMM + a(1)*GDM;
h = waitbar(0,'seismic analysis');
for ii = 1:length(t)-1
    % tic
    GFp = GF(:,ii+1) + GMM*( a(8)*GU(:,ii) + a(2)*GV(:,ii) + a(3)*GA(:,ii) ) ...
                     + GDM*( a(1)*GU(:,ii) + a(4)*GV(:,ii) + a(5)*GA(:,ii) );
    GU(:,ii+1) = GSMp\GFp;
    GA(:,ii+1) = a(8)*( GU(:,ii+1) - GU(:,ii) ) - a(2)*GV(:,ii) - a(3)*GA(:,ii);
    GV(:,ii+1) = GV(:,ii) + a(6)*GA(:,ii) + a(7)*GA(:,ii+1);
    waitbar(ii/(length(t)-1))
    % toc
end




