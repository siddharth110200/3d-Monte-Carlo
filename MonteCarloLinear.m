global hbar integfochi Tchig sample1 lambda;
rng('shuffle');
%FIXED PARAMETERS
hbar=137.06;
param=256; %Number of evaluation points
fine=0.001; %cutoff probability value input into smilei tables
integfochi=h5read('radiation_tables.h5','/integfochi');%actual values on sample1, as defined on Tables page
Tchig=h5read('multiphoton_Breit_Wheeler_tables.h5','/integration_dt_dchi');%actual values on sample1, as defined on Tables page
xirad=h5read('radiation_tables.h5','/xi');%actual values on X,Y, as defined on Tables page
xibreit=h5read('multiphoton_Breit_Wheeler_tables.h5','/xi');%actual values on X,Y, as defined on Tables page
xilowrad=h5read('radiation_tables.h5','/min_photon_chi_for_xi');%log10 of actualvalues stored on sample1
xilowbreit=h5read('multiphoton_Breit_Wheeler_tables.h5','/min_particle_chi_for_xi');%log10 of actual values stored
xihigh=linspace(-4,3,param);%log10 of actual upper limits stored
sample1=logspace(-4,3,param);%range for chi\pm
Y=repmat(sample1,param,1);
X=double.empty(param,0);
for iter=1:param
    append=logspace(xilowrad(iter),xihigh(iter),param);
    append=append(:);
    X=[X,append];
end
radinter=scatteredInterpolant(X(:),Y(:),xirad(:),'linear','nearest');
breitinter=scatteredInterpolant(X(:),Y(:),xibreit(:),'linear','nearest');
timeperiod=(w1(0)+w2(0))*pi/(w1(0)*w2(0));

%INPUT DATA
a=450;
lambda=1;%Global since code obtained through modification
tfinal=10*timeperiod;%Total time for integration
temperature=0.1;
automatic=true;%User input, or randomly selected particles
autoinitial=100;%Number of electrons to begin with when automatic

%INITIALISATION
if(~automatic)
    %elec=[xelec,yelec,zelec,pxelec,pyelec,pzelec,telec];
    elec=zeros(1,7);%Electron initial conditions array: 1: x,2:y,3:z,4:px,5:py,6:pz,7:Optical depth
    elec(4)=5;%Setting initial values, starting with 1 electron, no positron/photon
    elec(5)=5;
    elec(6)=10;
    posit=double.empty(0,7);%Positron initial conditions array
    phot=double.empty(0,7);%Photon initial conditions array
    Nelec=1;%Number of electrons
    Nposit=0;%Number of positrons
    Nphot=0;%Number of photons
    elecsol={};%Cell array storing electron solutions
    positsol={};%Cell array for positron solutions
    photsol={};%Cell array for photon solutions
    tinelec=[0];%Introduction/Production time of electrons. Firsst electron introduced at t=0
    tinposit=[];%Positron production times
    tinphot=[];%Photon production times
else
    Nelec=autoinitial;
    Nposit=0;
    Nphot=0;
    elecsol={};
    positsol={};
    photsol={};
    tinelec=zeros(Nelec,1);
    tinposit=[];
    tinphot=[];
    posit=double.empty(0,7);
    phot=double.empty(0,7);
    for iter=1:Nelec
        elec(iter,:)=maxjutt(temperature);
    end
end

%SIMULATION
iter=0;%Counts iterations completed
while(any(~(cellfun(@(solution) solution.x(end), [elecsol,positsol,photsol])==tfinal))||(iter==0))
    iter=iter+1;
%for iter=1:10000
    disp(string(iter)+" iterations "+string(Nelec)+" electrons "+string(min(cellfun(@(solution) solution.x(end), [elecsol,positsol,photsol]))/tfinal)+" time elapsed");
    for k=1:Nelec
        target=-log(rand);
        %opts=odeset('OutputFcn',@odetpbar);
        opts=odeset('Events',@(t,arr) eve(t,arr,target));
        if(k>size(elecsol,2) && tinelec(k)<tfinal)
            elecsol{k}=ode45(@(t,arr) delec(a,t,arr),[tinelec(k) tfinal],elec(k,:),opts);%Solves ODE for new electron
        elseif(k<=size(elecsol,2))
            if(elecsol{k}.x(end)<tfinal)
                elecsol{k}=odextend(elecsol{k},@(t,arr) delec(a,t,arr),tfinal,elec(k,:),opts);%Extends ODE solution for old electron
            end
        end
        if(tinelec(k)<tfinal)
            if(elecsol{k}.x(end)<tfinal)
                Nphot=Nphot+1;%Increases photon count during emission 
                temp=elecsol{k}.y(:,end);%Temporary variable storing electron solution just created/extended
                tempet=eta(a,elecsol{k}.x(end),temp(1),temp(2),temp(3),velx(temp(4),temp(5),temp(6)),vely(temp(4),temp(5),temp(6)),velz(temp(4),temp(5),temp(6)));
                temprand=(1-fine)*rand+fine;
                chif=fzero(@(lim)radinter(lim,tempet/2)-temprand-1E-6,[0,tempet/2]);
                modp=2*chif*sqrt(1+temp(4)^2+temp(5)^2+temp(6)^2)/tempet;
                transfer=modp*normalize(temp(4:6),'norm');%Assigning stuff using the Monte Carlo scheme
                phot=[phot;[temp(1),temp(2),temp(3),transfer(1),transfer(2),transfer(3),0]];%Transferring momentum: Entering initial conditions of new photon into array
                tinphot=[tinphot;elecsol{k}.x(end)];%Entering photon creation time into array 
                elec(k,:)=[temp(1),temp(2),temp(3),temp(4)-transfer(1),temp(5)-transfer(2),temp(6)-transfer(3),0];%Updating electron initial conditions array: Unnecessary, but anyway
            end
        end
    end
    for k=1:Nposit
        target=-log(rand);
        opts=odeset('Events',@(t,arr) eve(t,arr,target));
        if(k>size(positsol,2)&& tinposit(k)<tfinal)
            positsol{k}=ode45(@(t,arr) dposit(a,t,arr),[tinposit(k) tfinal],posit(k,:),opts);
        elseif(k<=size(positsol,2))
            if(positsol{k}.x(end)<tfinal)
                positsol{k}=odextend(positsol{k},@(t,arr) dposit(a,t,arr),tfinal,posit(k,:),opts);
            end
        end
        if(tinposit(k)<tfinal)
            if(positsol{k}.x(end)<tfinal)
                Nphot=Nphot+1;
                temp=positsol{k}.y(:,end);
                tempet=eta(a,positsol{k}.x(end),temp(1),temp(2),temp(3),velx(temp(4),temp(5),temp(6)),vely(temp(4),temp(5),temp(6)),velz(temp(4),temp(5),temp(6)));
                temprand=(1-fine)*rand+fine;
                chif=fzero(@(lim)radinter(lim,tempet/2)-temprand-1E-6,[0,tempet/2]);
                modp=2*chif*sqrt(1+temp(4)^2+temp(5)^2+temp(6)^2)/tempet;
                transfer=modp*normalize(temp(4:6),'norm');
                phot=[phot;[temp(1),temp(2),temp(3),transfer(1),transfer(2),transfer(3),0]];
                tinphot=[tinphot;positsol{k}.x(end)];
                posit(k,:)=[temp(1),temp(2),temp(3),temp(4)-transfer(1),temp(5)-transfer(2),temp(6)-transfer(3),0];
            end
        end
    end
    for k=1:Nphot
        photonenergy=sqrt(phot(k,4)^2+phot(k,5)^2+phot(k,6)^2);
        target=-log(rand);
        opts=odeset('Events',@(t,arr) eve(t,arr,target));
        if(k>size(photsol,2)&& tinphot(k)<tfinal)
            photsol{k}=ode45(@(t,arr) dphot(a,t,arr),[tinphot(k),tfinal],phot(k,:),opts);
        elseif(k<=size(photsol,2))
            if(photsol{k}.x(end)<tfinal)
                photsol{k}=odextend(photsol{k},@(t,arr) dphot(a,t,arr),tfinal,phot(k,:),opts);
            end
        end
        if(photonenergy>2)
            if(tinphot(k)<tfinal)
                if(photsol{k}.x(end)<tfinal)
                    Nelec=Nelec+1;
                    Nposit=Nposit+1;
                    temp=photsol{k}.y(:,end);
                    tempchig=chig(a,photsol{k}.x(end),temp(1),temp(2),temp(3),temp(4)/photonenergy,temp(5)/photonenergy,temp(6)/photonenergy,photonenergy/hbar);
                    temprand=(1-fine)*rand+fine;
                    chimf=fzero(@(lim)breitinter(lim,tempchig)-temprand/2-1E-6,[0,tempchig]);
                    chipf=tempchig-chimf;
                    energye=1+(photonenergy-2)*chimf/tempchig;
                    energyp=1+(photonenergy-2)*chipf/tempchig;
                    modpe=sqrt(energye^2-1);
                    modpp=sqrt(energyp^2-1);
                    transfere=modpe*normalize(temp(4:6),'norm');
                    transferp=modpp*normalize(temp(4:6),'norm');
                    elec=[elec;[temp(1),temp(2),temp(3),transfere(1),transfere(2),transfere(3),0]];
                    posit=[posit;[temp(1),temp(2),temp(3),transferp(1),transferp(2),transferp(3),0]];
                    tinelec=[tinelec;photsol{k}.x(end)];
                    tinposit=[tinposit;photsol{k}.x(end)];
                    phot(k,:)=[temp(1),temp(2),temp(3),temp(4)-(transfere(1)+transferp(1)),temp(5)-(transfere(2)+transferp(2)),temp(6)-(transfere(3)+transferp(3)),0];
                end
            end
        end
    end
    %keep=[keep,chif/(tempet/2)];
end

%PLOTTING SOLUTIONS.
% for k=1:2%First 2 electrons
%     plot3(elecsol{k}.y(1,:),elecsol{k}.y(2,:),elecsol{k}.y(3,:))%Postion space. Use 4,5,6 for momentum space.
%     hold on
% end
%tobe(:,1)=cellfun(@(solution) solution.y(4,end), elecsol)
%tobe(:,2)=cellfun(@(solution) sqrt(solution.y(5,end)^2+solution.y(6,end)^2), elecsol)
alsoplot1=cellfun(@(solution) atan(solution.y(4,end)/sqrt(solution.y(5,end)^2+solution.y(6,end)^2)), elecsol);
alsoplot2=cellfun(@(solution) atan(solution.y(4,end)/sqrt(solution.y(5,end)^2+solution.y(6,end)^2)), positsol);
alsoplot3=cellfun(@(solution) atan(solution.y(4,end)/sqrt(solution.y(5,end)^2+solution.y(6,end)^2)), photsol);
disp("Electron mean: "+string(nanmean(alsoplot1))+" Electron SD: "+string(nanstd(alsoplot1))+" Positron mean: "+string(nanmean(alsoplot2))+" Positron SD: "+string(nanstd(alsoplot2))+" Photon mean: "+string(nanmean(alsoplot3))+" Photon SD: "+string(nanstd(alsoplot3)));
alsoplot4=cellfun(@(solution) sqrt(1+solution.y(4,end)^2+solution.y(5,end)^2+solution.y(6,end)^2), elecsol);
%[histw,intervals]=histwc(alsoplot2,alsoplot5,50);bar(intervals, histw)
%[histw,intervals]=histwc(alsoplot2,alsoplot2*0+100,50);bar(intervals, histw)
%clear all;load('1.6_3.mat','elecsol');
% disp("1)"+ string(mean(cellfun(@(solution) atan(solution.y(4,end)/sqrt(solution.y(5,end)^2+solution.y(6,end)^2)), elecsol)))+","+"2)"+string(std(cellfun(@(solution) atan(solution.y(4,end)/sqrt(solution.y(5,end)^2+solution.y(6,end)^2)), elecsol)))+","+"3)"+string(sum(cellfun(@(solution) sqrt(1+solution.y(4,end)^2+solution.y(5,end)^2+solution.y(6,end)^2)*atan(solution.y(4,end)/sqrt(solution.y(5,end)^2+solution.y(6,end)^2)), elecsol))/sum(cellfun(@(solution) sqrt(1+solution.y(4,end)^2+solution.y(5,end)^2+solution.y(6,end)^2), elecsol)))+","+"4)"+string(sqrt(sum(cellfun(@(solution) sqrt(1+solution.y(4,end)^2+solution.y(5,end)^2+solution.y(6,end)^2)*atan(solution.y(4,end)/sqrt(solution.y(5,end)^2+solution.y(6,end)^2))^2, elecsol))/sum(cellfun(@(solution) sqrt(1+solution.y(4,end)^2+solution.y(5,end)^2+solution.y(6,end)^2), elecsol))-sum(cellfun(@(solution) sqrt(1+solution.y(4,end)^2+solution.y(5,end)^2+solution.y(6,end)^2)*atan(solution.y(4,end)/sqrt(solution.y(5,end)^2+solution.y(6,end)^2)), elecsol))^2/sum(cellfun(@(solution) sqrt(1+solution.y(4,end)^2+solution.y(5,end)^2+solution.y(6,end)^2), elecsol))^2)))
% 
% disp("1)"+ string(mean(cellfun(@(solution) atan(solution.y(4,end)/sqrt(solution.y(5,end)^2+solution.y(6,end)^2)), positsol)))+","+"2)"+string(std(cellfun(@(solution) atan(solution.y(4,end)/sqrt(solution.y(5,end)^2+solution.y(6,end)^2)), positsol)))+","+"3)"+string(sum(cellfun(@(solution) sqrt(1+solution.y(4,end)^2+solution.y(5,end)^2+solution.y(6,end)^2)*atan(solution.y(4,end)/sqrt(solution.y(5,end)^2+solution.y(6,end)^2)), positsol))/sum(cellfun(@(solution) sqrt(1+solution.y(4,end)^2+solution.y(5,end)^2+solution.y(6,end)^2), positsol)))+","+"4)"+string(sqrt(sum(cellfun(@(solution) sqrt(1+solution.y(4,end)^2+solution.y(5,end)^2+solution.y(6,end)^2)*atan(solution.y(4,end)/sqrt(solution.y(5,end)^2+solution.y(6,end)^2))^2, positsol))/sum(cellfun(@(solution) sqrt(1+solution.y(4,end)^2+solution.y(5,end)^2+solution.y(6,end)^2), positsol))-sum(cellfun(@(solution) sqrt(1+solution.y(4,end)^2+solution.y(5,end)^2+solution.y(6,end)^2)*atan(solution.y(4,end)/sqrt(solution.y(5,end)^2+solution.y(6,end)^2)), positsol))^2/sum(cellfun(@(solution) sqrt(1+solution.y(4,end)^2+solution.y(5,end)^2+solution.y(6,end)^2), positsol))^2)))
% 
% disp("1)"+ string(nanmean(cellfun(@(solution) atan(solution.y(4,end)/sqrt(solution.y(5,end)^2+solution.y(6,end)^2)), photsol)))+","+"2)"+string(nanstd(cellfun(@(solution) atan(solution.y(4,end)/sqrt(solution.y(5,end)^2+solution.y(6,end)^2)), photsol)))+","+"3)"+string(nansum(cellfun(@(solution) sqrt(solution.y(4,end)^2+solution.y(5,end)^2+solution.y(6,end)^2)*atan(solution.y(4,end)/sqrt(solution.y(5,end)^2+solution.y(6,end)^2)), photsol))/nansum(cellfun(@(solution) sqrt(solution.y(4,end)^2+solution.y(5,end)^2+solution.y(6,end)^2), photsol)))+","+"4)"+string(sqrt(nansum(cellfun(@(solution) sqrt(solution.y(4,end)^2+solution.y(5,end)^2+solution.y(6,end)^2)*atan(solution.y(4,end)/sqrt(solution.y(5,end)^2+solution.y(6,end)^2))^2, photsol))/nansum(cellfun(@(solution) sqrt(solution.y(4,end)^2+solution.y(5,end)^2+solution.y(6,end)^2), photsol))-nansum(cellfun(@(solution) sqrt(solution.y(4,end)^2+solution.y(5,end)^2+solution.y(6,end)^2)*atan(solution.y(4,end)/sqrt(solution.y(5,end)^2+solution.y(6,end)^2)), photsol))^2/nansum(cellfun(@(solution) sqrt(solution.y(4,end)^2+solution.y(5,end)^2+solution.y(6,end)^2), photsol))^2)))
%xt = get(gca, 'XTick');                                 % 'XTick' Values
%set(gca, 'XTick', xt, 'XTickLabel', xt/10)            
%histogram(alsoplot,'Normalization','pdf')
storeindex=[];
avgenergy=[];
energ={};
for time=0:tfinal/100:tfinal
    time/tfinal
    for i=1:Nelec
        if(tinelec(i)<time)
            storeindex=[storeindex,i];
        end
    end
    energ={};
    for i=1:size(storeindex,2)
        energ{i}=elecsol{storeindex(i)};
    end
    avgenergy=[avgenergy,std(cellfun(@(solution) atan(deval(solution,time,4)/sqrt(deval(solution,time,5)^2+deval(solution,time,6)^2)), energ))];
end
 %set(gco,'LineWidth',1)
%FUNCTIONS
function out=w1(t) %FREQUENCY vs TIME, 1st wave
out=1.2e-8;
end
function out=w2(t)%FREQUENCY vs TIME, 2nd wave
out=2.1e-8;
end
function out=maxjutt(theta)
temprand=rand;
maxrand=fzero(@(fin) (integral(@(gamma) (gamma.*sqrt(gamma.^2-1).*exp(-gamma./theta)./(theta.*besselk(2,1/theta))),1,fin)-temprand),[1,1000]);
randp=sqrt(maxrand^2-1);
wavelength=2*pi/sqrt(w1(0)*w2(0));
elevation=asin(2*rand-1);
azimuth=2*pi*rand;
randx=wavelength*(2*rand-1);
randy=wavelength*(2*rand-1);
randz=wavelength*(2*rand-1);
[randpx,randpy,randpz]=sph2cart(azimuth,elevation,1);
out=[randx,randy,randz,randp*randpx,randp*randpy,randp*randpz,0];
end
function out=xforce(a0,t,x,~,~,~,~,vz)
global lambda;
out=(1/2).*a0.*vz.*(sin (((-1).*t+x).*w1(t)).*w1(t)+lambda.* sin ((t+x).*w2(t)).*w2(t));
end
function out=yforce(~,~,~,~,~,~,~,~)
out=0;
end
function out=zforce(a0,t,x,~,~,vx,~,~)
global lambda;
out=(-1/2).*a0.*(((-1)+vx).*sin (((-1).*t+x).*w1(t)).*w1(t)+lambda.*(1+vx).*sin ((t+x).*w2(t)).*w2(t));
end
function out=chig(a0,t,x,~,~,vx,~,vz,freq) %Photon velocity and freq!
global hbar lambda;
out=(hbar^2)*freq*abs(sqrt((1/4).*a0.^2.*(((-1)+vx).^2.*sin(((-1).*t+x).*w1(t)).^2.*w1(t).^2+2.*lambda.*((-1)+vx.^2+2.*vz.^2).*sin(((-1).*t+x).*w1(t)).*sin((t+x).*w2(t)).*w1(t).*w2(t)+lambda.^2.*(1+vx).^2.*sin((t+x).*w2(t)).^2.*w2(t).^2)));
end

function out=eta(a0,t,x,~,~,vx,vy,vz)
global hbar lambda;
out=hbar*sqrt((1/4).*a0.^2.*(((-1)+vx).^2.*sin(((-1).*t+x).*w1(t)).^2.*w1(t).^2+2.*lambda.*((-1)+vx.^2+2.*vz.^2).*sin(((-1).*t+x).*w1(t)).*sin((t+x).*w2(t)).*w1(t).*w2(t)+lambda.^2.*(1+vx).^2.*sin((t+x).*w2(t)).^2.*w2(t).^2))/sqrt(1-vx^2-vy^2-vz^2);
end

function out=velx(px,py,pz)
out=px/sqrt(1+px^2+py^2+pz^2);
end
function out=vely(px,py,pz)
out=py/sqrt(1+px^2+py^2+pz^2);
end
function out=velz(px,py,pz)
out=pz/sqrt(1+px^2+py^2+pz^2);
end
function out=elecx(~,~,~,~,~,px,py,pz)
out=velx(px,py,pz);
end
function out=elecy(~,~,~,~,~,px,py,pz)
out=vely(px,py,pz);
end
function out=elecz(~,~,~,~,~,px,py,pz)
out=velz(px,py,pz);
end
function out=positx(~,~,~,~,~,px,py,pz)
out=velx(px,py,pz);
end
function out=posity(~,~,~,~,~,px,py,pz)
out=vely(px,py,pz);
end
function out=positz(~,~,~,~,~,px,py,pz)
out=velz(px,py,pz);
end
function out=photx(~,~,~,~,~,px,py,pz)
out=px/sqrt(px^2+py^2+pz^2);
end
function out=photy(~,~,~,~,~,px,py,pz)
out=py/sqrt(px^2+py^2+pz^2);
end
function out=photz(~,~,~,~,~,px,py,pz)
out=pz/sqrt(px^2+py^2+pz^2);
end
function out=elecpx(a0,t,x,~,~,px,py,pz)
out=xforce(a0,t,x,0,0,velx(px,py,pz),vely(px,py,pz),velz(px,py,pz));
end
function out=elecpy(a0,t,x,~,~,px,py,pz)
out=yforce(a0,t,x,0,0,velx(px,py,pz),vely(px,py,pz),velz(px,py,pz));
end
function out=elecpz(a0,t,x,~,~,px,py,pz)
out=zforce(a0,t,x,0,0,velx(px,py,pz),vely(px,py,pz),velz(px,py,pz));
end
function out=positpx(a0,t,x,~,~,px,py,pz)
out=-xforce(a0,t,x,0,0,velx(px,py,pz),vely(px,py,pz),velz(px,py,pz));
end
function out=positpy(a0,t,x,~,~,px,py,pz)
out=-yforce(a0,t,x,0,0,velx(px,py,pz),vely(px,py,pz),velz(px,py,pz));
end
function out=positpz(a0,t,x,~,~,px,py,pz)
out=-zforce(a0,t,x,0,0,velx(px,py,pz),vely(px,py,pz),velz(px,py,pz));
end
function out=photpx(~,~,~,~,~,~,~,~)
out=0;
end
function out=photpy(~,~,~,~,~,~,~,~)
out=0;
end
function out=photpz(~,~,~,~,~,~,~,~)
out=0;
end

function out=elect(a0,t,x,~,~,px,py,pz)
global hbar integfochi sample1;
intervx=velx(px,py,pz);
intervy=vely(px,py,pz);
intervz=velz(px,py,pz);
n=eta(a0,t,x,0,0,intervx,intervy,intervz);
out=sqrt(3)*n*interp1(sample1,integfochi,n,'linear','extrap')./(2*pi*hbar^2*sqrt(1+px^2+py^2+pz^2));
end

function out=positt(a0,t,x,~,~,px,py,pz)
global hbar integfochi sample1;
intervx=velx(px,py,pz);
intervy=vely(px,py,pz);
intervz=velz(px,py,pz);
n=eta(a0,t,x,0,0,intervx,intervy,intervz);
out=sqrt(3)*n*interp1(sample1,integfochi,n,'linear','extrap')./(2*pi*hbar^2*sqrt(1+px^2+py^2+pz^2));
end

function out=phott(a0,t,x,y,z,px,py,pz)%Photon momenta!
global hbar Tchig sample1;
intervx=px/sqrt(px^2+py^2+pz^2);
intervy=py/sqrt(px^2+py^2+pz^2);
intervz=pz/sqrt(px^2+py^2+pz^2);
freq=sqrt(px^2+py^2+pz^2)/hbar;
interchig=chig(a0,t,x,y,z,intervx,intervy,intervz,freq);
out=interp1(sample1,Tchig,interchig,'linear','extrap')/(hbar^3*freq*pi*sqrt(3)*interchig);
end
function out=delec(a0,t,arr)
    out=[elecx(a0,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6)),elecy(a0,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6)),elecz(a0,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6)),elecpx(a0,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6)),elecpy(a0,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6)),elecpz(a0,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6)),elect(a0,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6))];
    out=out(:);
end
function out=dposit(a0,t,arr)
    out=[positx(a0,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6)),posity(a0,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6)),positz(a0,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6)),positpx(a0,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6)),positpy(a0,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6)),positpz(a0,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6)),positt(a0,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6))];
    out=out(:);
end
function out=dphot(a0,t,arr)
    out=[photx(a0,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6)),photy(a0,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6)),photz(a0,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6)),photpx(a0,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6)),photpy(a0,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6)),photpz(a0,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6)),phott(a0,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6))];
    out=out(:);
end
function [value,isterminal,direction]=eve(~,arr,target)
    value=arr(7)-target;
    isterminal=1;
    direction=0;%-1: Where event function is decreasing
end