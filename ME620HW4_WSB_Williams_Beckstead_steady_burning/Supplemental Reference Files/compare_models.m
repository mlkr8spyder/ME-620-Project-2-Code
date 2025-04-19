% Andrew Rettenmaier
% Compute WSB, WDB, Beckstead burning rate models

clear all
close all
clc


% Choose what to compare
type = 1;

% Pressure Response
if type == 1
    P  = [1,7,20,70,80,90];
    T0 = [293,303,313,663,673,683]
end
% ZN stability from Ward et al
if type == 2
    P   = [1 7 10 20 70 80];
    T0  = [200:50:450];
end
WSB_HMX;
beckstead;

%% Stability

% For WSB
for ii = 1:length(P)-1
    r(ii,:) = diff(Tsd(ii,:))./diff(T0);
    rte(ii,:)    = (Tsd(ii,3) - Tsd(ii,1))/(T0(3) - T0(1));
    k(ii,:) = (Tsd(ii,1:end-1) - T0(1:end-1)).*diff(log(rd(ii,:)))./diff(T0);
end
for jj = 1:length(T0)-1
    up(:,jj) = diff(log(rd(:,jj)))./diff(log(P'*101325*10));
    mu(:,jj) = (1./(Tsd(1:end-1,jj) - T0(jj))).*(diff(Tsd(:,jj))./diff(log(P'*101325*10)));
end
delta = up.*r - mu.*k;
% For WDB
for ii = 1:length(P)-1
    rdb(ii,:) = diff(Ts_DB(ii,:))./diff(T0);
    kdb(ii,:) = (Ts_DB(ii,1:end-1) - T0(1:end-1)).*diff(log(md_DB(ii,:)))./diff(T0);
end
for jj = 1:length(T0)-1
    updb(:,jj) = diff(log(rd_DB(:,jj)))./diff(log(P'*101325*10));
    mudb(:,jj) = (1./(Ts_DB(1:end-1,jj) - T0(jj))).*(diff(Ts_DB(:,jj))./diff(log(P'*101325*10)));
end
deltadb = updb.*rdb - mudb.*kdb;

% For Beckstead
for ii = 1:length(P)-1
    rbe(ii,:) = diff(Tsb(ii,:))./diff(T0);
    kbe(ii,:) = (Tsb(ii,1:end-1) - T0(1:end-1)).*diff(log(mb(ii,:)))./diff(T0);
end
for jj = 1:length(T0)-1
    upbe(:,jj) = diff(log(rb(:,jj)))./diff(log(P'*101325*10));
    mube(:,jj) = (1./(Tsb(1:end-1,jj) - T0(jj))).*(diff(Tsb(:,jj))./diff(log(P'*101325*10)));
end
deltabe = upbe.*rbe - mube.*kbe;

%%
if type == 1
    
    f = linspace(10,5000,1000);
    w = 2*pi*f;
    %% For 7 atm
    p = find(P==7,1,'first');
    t = find(T0==293,1,'first');
    %% WSB
    Omega = w*8e-4/rd(p,t).^2;
    lamb  = 0.5 + 0.5*(1 + 4*1i*Omega).^0.5;
    z     = -0.5 + 0.5*(1 + 4*1i*Omega).^0.5;
    
    up1    = up(p,t);
    mu1    = mu(p,t);
    k1     = k(p,t);
    r1     = r(p,t);
    delta1 = up1*r1 - mu1*k1;
    
    A      = k1/r1;
    B      = 1/k1;
    
    k2      = 1.07;
    r2      = 0.05;
    delta2 = -0.061;
    A2      = 21;
    B2      = .93;
    up2    = .8;
    
    for ii = 1:length(P) - 1
        for jj = 1:length(T0) - 1
            Rp(ii,jj,:)    = (up(ii,jj) + delta(ii,jj).*(lamb-1))./(1 + (r(ii,jj) - k(ii,jj)./lamb).*(lamb - 1));
        end
    end
    
    %% For Beckstead
    Omegabe = w*8e-4/rb(p,t)^2;
    lamb  = 0.5 + 0.5*(1 + 4*1i*Omegabe).^0.5;
    z     = -0.5 + 0.5*(1 + 4*1i*Omegabe).^0.5;
    
    upbe(p,t)
    deltabe(p,t)
    rbe(p,t)
    kbe(p,t)
    for ii = 1:length(P) - 1
        for jj = 1:length(T0) - 1
            Rpbe(ii,jj,:)    = (upbe(ii,jj) + deltabe(ii,jj).*(lamb-1))./(1 + (rbe(ii,jj) - kbe(ii,jj)./lamb).*(lamb - 1));
        end
    end
    %% WDB
    Omegadb = w*8e-4/rd_DB(p,t)^2;
    lamb  = 0.5 + 0.5*(1 + 4*1i*Omegadb).^0.5;
    z     = -0.5 + 0.5*(1 + 4*1i*Omegadb).^0.5;
    
    for ii = 1:length(P) - 1
        for jj = 1:length(T0) - 1
            Rpdb(ii,jj,:)    = (updb(ii,jj) + deltadb(ii,jj).*(lamb-1))./(1 + (rdb(ii,jj) - kdb(ii,jj)./lamb).*(lamb - 1));
        end
    end
    %% For 70 atm
    p68 = find(P==70,1,'first');
    t68 = find(T0==293,1,'first');
    %% WSB
    Omega68 = w*8e-4/rd(p68,t68).^2;
    lamb  = 0.5 + 0.5*(1 + 4*1i*Omega68).^0.5;
    z     = -0.5 + 0.5*(1 + 4*1i*Omega68).^0.5;
    
    up1    = up(p,t);
    mu1    = mu(p,t);
    k1     = k(p,t);
    r1     = r(p,t);
    delta1 = up1*r1 - mu1*k1;
    
    A      = k1/r1;
    B      = 1/k1;
    
    k2      = 1.07;
    r2      = 0.05;
    delta2 = -0.061;
    A2      = 21;
    B2      = .93;
    up2    = .8;
    
    for ii = 1:length(P) - 1
        for jj = 1:length(T0) - 1
            Rp68(ii,jj,:)    = (up(ii,jj) + delta(ii,jj).*(lamb-1))./(1 + (r(ii,jj) - k(ii,jj)./lamb).*(lamb - 1));
        end
    end
    
    %% For Beckstead
    Omegabe68 = w*8e-4/rb(p68,t68)^2;
    lamb  = 0.5 + 0.5*(1 + 4*1i*Omegabe68).^0.5;
    z     = -0.5 + 0.5*(1 + 4*1i*Omegabe68).^0.5;
    
    for ii = 1:length(P) - 1
        for jj = 1:length(T0) - 1
            Rpbe68(ii,jj,:)    = (upbe(ii,jj) + deltabe(ii,jj).*(lamb-1))./(1 + (rbe(ii,jj) - kbe(ii,jj)./lamb).*(lamb - 1));
        end
    end
    %% WDB
    Omegadb68 = w*8e-4/rd_DB(p68,t68)^2;
    lamb  = 0.5 + 0.5*(1 + 4*1i*Omegadb68).^0.5;
    z     = -0.5 + 0.5*(1 + 4*1i*Omegadb68).^0.5;
    
    for ii = 1:length(P) - 1
        for jj = 1:length(T0) - 1
            Rpdb68(ii,jj,:)    = (updb(ii,jj) + deltadb(ii,jj).*(lamb-1))./(1 + (rdb(ii,jj) - kdb(ii,jj)./lamb).*(lamb - 1));
        end
    end
elseif type == 2
    %% Stability
    kp  = 1:.001:3;
    rp   = (kp-1).^2./(kp+1);
    kp   = [0 1 kp];
    rp   = [0 0 rp];
    for ii = 1:length(P)-1
        r200(ii,:) = [r(ii,1) rdb(ii,1) rbe(ii,1)];
        r250(ii,:) = [r(ii,2) rdb(ii,2) rbe(ii,2)];
        r300(ii,:) = [r(ii,3) rdb(ii,3) rbe(ii,3)];
        r350(ii,:) = [r(ii,4) rdb(ii,4) rbe(ii,4)];
        r400(ii,:) = [r(ii,5) rdb(ii,5) rbe(ii,5)];
        %r400(ii,:) = [r(ii,6) rdb(ii,6) rbe(ii,6)];
        
        k200(ii,:) = [k(ii,1) kdb(ii,1) kbe(ii,1)];
        k250(ii,:) = [k(ii,2) kdb(ii,2) kbe(ii,2)];
        k300(ii,:) = [k(ii,3) kdb(ii,3) kbe(ii,3)];
        k350(ii,:) = [k(ii,4) kdb(ii,4) kbe(ii,4)];
        k400(ii,:) = [k(ii,5) kdb(ii,5) kbe(ii,5)];
        %k400(ii,:) = [k(ii,6) kdb(ii,6) kbe(ii,6)];
    end
    k1              = [k300(1,1) k300(1,2) k300(1,3)];
    k5              = [k300(2,1) k300(2,2) k300(2,3)];
    k10             = [k300(3,1) k300(3,2) k300(3,3)];
    k20             = [k300(4,1) k300(4,2) k300(4,3)];
    k70             = [k300(5,1) k300(5,2) k300(5,3)];
    
    r1              = [r300(1,1) r300(1,2) r300(1,3)];
    r5              = [r300(2,1) r300(2,2) r300(2,3)];
    r10             = [r300(3,1) r300(3,2) r300(3,3)];
    r20             = [r300(4,1) r300(4,2) r300(4,3)];
    r70             = [r300(5,1) r300(5,2) r300(5,3)];
end

%% Plotting

if type == 2
    figure(1);
    plot(k200(4,3),r200(4,3),'ok',k250(4,3),r250(4,3),'^k',k300(4,3),r300(4,3),'sk',k350(4,3),r350(4,3),'<k',k400(4,3),r400(4,3),'dk');hold on;

    plot(k1,r1,'ok',k5,r5,'^k',k10,r10,'sk',k20,r20,'<k',k70,r70,'dk','MarkerFaceColor','k');hold on;
    
    legend('P=20 atm, T_0=200 K','P=20 atm, T_0=250 K','P=20 atm, T_0=300 K','P=20 atm, T_0=350 K','P=20 atm, T_0=400 K','P=1 atm, T_0=300 K','P=5 atm, T_0=300 K','P=14 atm, T_0=300 K','P=20 atm, T_0=300 K','P=70 atm, T_0=300 K','Location','best');
    plot(k200(4,2),r200(4,3),'ok',k250(4,2),r250(4,2),'^k',k300(4,2),r300(4,2),'sk',k350(4,2),r350(4,2),'<k',k400(4,2),r400(4,2),'dk');hold on;
    plot(k200(4,1),r200(4,1),'ok',k250(4,1),r250(4,1),'^k',k300(4,1),r300(4,1),'sk',k350(4,1),r350(4,1),'<k',k400(4,1),r400(4,1),'dk');hold on;

    plot(kp,rp,'k');hold off;
    xlabel('k');ylabel('r');
    xlim([0 3]);
    figset(16,2,'b',10,1,800,800);
elseif type == 1
    fT14 = [314.106 314.106 460.699 818.333 841.03 841.03 1267.754 1267.754 1267.754 322.818];
    RpT14 =[1.488 1.349 .516 .595 .714 .833 .754 .575 .397 1.052];
    fT68 = [350.43 500.106 912.969 864.357 1910.988 1963.991 2018.465];
    RpT68 =[1.052 1.21 1.111 1.27 1.429 1.726 2.202];
    
    figure(2);
    subplot(2,1,1);
%     semilogx(fT14,RpT14,'ok',fT68,RpT68,'sk','MarkerFaceColor','k');hold on;
%     legend('T-burner (7 atm)','T-burner (70 atm)');
    semilogx(f,abs(reshape(Rp(p,t,:),1,length(Omega))),'k');
    hold on;
    semilogx(f,abs(reshape(Rp68(p68,t68,:),1,length(Omega68))),'-.k');
%     
%     semilogx(f,real(reshape(Rpbe(p,t,:),1,length(Omegabe))),'r');
%     semilogx(f,real(reshape(Rpbe68(p68,t68,:),1,length(Omegabe68))),'-.r');
%     
    semilogx(f,abs(reshape(Rpdb(p,t,:),1,length(Omegadb))),'b');
    semilogx(f,abs(reshape(Rpdb68(p68,t68,:),1,length(Omegadb68))),'-.b');hold off;
    legend('Rp WSB 7 atm','Rp WSB 70 atm', 'Rp WDB 7 atm','Rp WDB 70 atm');
    xlabel('Frequency, Hz');
    ylabel('Re[ R_p ]');
    %figset(16,2,'b',10,2,800,800);
    subplot(2,1,2);
%     semilogx(fT14, angle(RpT14), 'ok', fT68, angle(RpT68), 'sk', 'MarkerFaceColor','k');
    hold on;
    semilogx(f,angle(reshape(Rp(p,t,:),1,length(Omega))),'k');
    semilogx(f,angle(reshape(Rp68(p68,t68,:),1,length(Omega68))),'-.k');
    semilogx(f,angle(reshape(Rpdb(p,t,:),1,length(Omegadb))),'b');
    semilogx(f,angle(reshape(Rpdb68(p68,t68,:),1,length(Omegadb68))),'-.b');hold off;
%     legend('T-burner (7 atm)','T-burner (70 atm)')
    xlabel('Frequency, Hz');
    ylabel('Phase[ R_p ]');
end


