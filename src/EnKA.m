function Xu=EnKA(X,y,Yp,alpha,R,dostoch)
%% Xu=EnKA(X,y,Yp,alpha,R,dostoch)
% Efficient implementation of the ensemble Kalman analysis step.
% Based on Evensen et al. (2022) and Emerick (2016).
% Assumes that R is diagonal and specified as a vector or scalar.
No=size(y,1);
Ne=size(X,2);


if isscalar(R)
    r=R.*ones(No,1);
    R=diag(r);
elseif size(R,1)==No&&numel(R)==No
    r=R;
    R=diag(R);
elseif size(R,1)==No&&size(R, 2)==No
    r=diag(R);
else
    error('Check R');
end

Xm=mean(X,2);
Xa=X-Xm;
Ypm=mean(Yp,2);
Ypa=Yp-Ypm;

if No>Ne


ris=1./sqrt(r);
Z=ris.*Ypa; % Normalize predicted anomalies via broadcasting.
[~,S,V]=svd(Z,'econ');
sigmas=diag(S);
tv=sum(sigmas);
svr=cumsum(sigmas)./sum(sigmas);
svthresh=0.99;
Nr=find(svr>=svthresh,1,'first');
if tv==0 % In case of complete ensemble collapse 
    % e.g.EnKF for FSCA after melt.
    Xu=X; 
    return;
end
Vr=V(:,1:Nr);
Sr2=sigmas(1:Nr).^2;
D=1./(Sr2+alpha*Ne); % Nr x 1
if dostoch % Classic stochastic EnKF
    Ew=sqrt(alpha)*randn(No,Ne);
    innoes=(ris.*(y-Yp))+Ew;
    pinnoes=Z'*innoes;
    X1=Vr'*pinnoes;
    X2=diag(D)*X1;
    W=Vr*X2;
    Xu=X+Xa*W;
else % DEnKF
    % Update mean
    innoms=ris.*(y-Ypm);
    pinnoms=Z'*innoms;
    X1=Vr'*pinnoms;
    X2=D.*X1;
    wm=Vr*X2;
    Xum=Xm+(Xa*wm);

    % Update anomaly
    Wa=Vr*(diag(D.*Sr2)*Vr');
    Xua=Xa-0.5*(Xa*Wa);

    % Combine for final update
    Xu=Xum+Xua;
end

else

    Ypat=Ypa';
    C_XY=(Xa*Ypat)./Ne;
    C_YY=(Ypa*Ypat)./Ne;
    K=C_XY/(C_YY+alpha.*R);

    if dostoch
        perts=sqrt(alpha).*sqrt(R)*randn(No,Ne);
        Y=y+perts;
        Xu=X+K*(Y-Yp);
    else
        Innm=y-mean(Yp,2);
        Xm=mean(X,2);
        Xum=Xm+K*Innm;
        Xua=Xa-0.5.*(K*Ypa);
        Xu=Xum+Xua;
    end

end

end
