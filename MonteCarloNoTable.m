global hbar;
hbar=137.06;
%lambda/a0 is pair production length/time scale
a=10;
w1=1;
w2=1;

%modp is transfer momentum magnitude; p is producer momentum
% xelec=0;yelec=0;zelec=0;pxelec=0;pyelec=0;pzelec=0;telec=0;elec=[xelec,yelec,zelec,pxelec,pyelec,pzelec,telec];
elec=zeros(1,7);
posit=double.empty(0,7);
phot=double.empty(0,7);
Nelec=1;
Nposit=0;
Nphot=0;
elecsol={};
positsol={};
photsol={};
tinelec=[0];
tinposit=[];
tinphot=[];
tfinal=5;
%reset tau, momenta at end of k-loop
%  x=linspace(0,0.5);
%  plot(x,F(1,x))
% x=logspace(-1,5);
% loglog(x,arrayfun(@(g)integral(@(m) ((integral(@(s)airy(s),nthroot(g/(m*(g-m)),3)^2,999)+((2-g^2/(m*(g-m)))/(nthroot(g/(m*(g-m)),3)^2))*airy(1,nthroot(g/(m*(g-m)),3)^2))),0,g,'ArrayValued',1)/g^2,x))
fplot(@(chi)F(1300,chi)/chi,[0,650])
drawnow
for iter=1:1
    for k=1:Nelec
        target=-log(rand);
        opts=odeset('OutputFcn',@odetpbar);
        %opts=odeset('Events',@(t,arr) eve(t,arr,target));
        if(k>size(elecsol,2))
            elecsol{k}=ode45(@(t,arr) delec(a,w1,w2,t,arr),[tinelec(k) tfinal],elec(k,:),opts);
        else
            elecsol{k}=odextend(elecsol{k},@(t,arr) delec(a,w1,w2,t,arr),tfinal,elec(k,:),opts);
        end
        Nphot=Nphot+1;
        temp=elecsol{k}.y(:,end);
        tempet=eta(a,w1,w2,elecsol{k}.x(end),temp(1),temp(2),temp(3),velx(temp(4),temp(5),temp(6)),vely(temp(4),temp(5),temp(6)),velz(temp(4),temp(5),temp(6)));
        temprand=rand;
        chif=fzero(@(lim) integral(@(chi)spect(chi,a,w1,w2,elecsol{k}.x(end),temp(1),temp(2),temp(3),velx(temp(4),temp(5),temp(6)),vely(temp(4),temp(5),temp(6)),velz(temp(4),temp(5),temp(6))),0,lim)-temprand*elect(a,w1,w2,elecsol{k}.x(end),temp(1),temp(2),temp(3),temp(4),temp(5),temp(6)),[eps,tempet/2-eps]);
        modp=2*chif*sqrt(1+temp(4)^2+temp(5)^2+temp(6)^2)/tempet;
        transfer=modp*normalize(temp(4:6),'norm');
        phot=[phot;[temp(1),temp(2),temp(3),transfer(1),transfer(2),transfer(3),0]];
        tinphot=[tinphot;elecsol{k}.x(end)];
        elec(k,:)=[temp(1),temp(2),temp(3),temp(4)-transfer(1),temp(5)-transfer(2),temp(6)-transfer(3),0];
    end
    for k=1:Nposit
        target=-log(rand);
        opts=odeset('Events',@(t,arr) eve(t,arr,target));
        if(k>size(positsol,2))
            positsol{k}=ode45(@(t,arr) dposit(a,w1,w2,t,arr),[tinposit(k) tfinal],posit(k,:),opts);
        else
            positsol{k}=odextend(positsol{k},@(t,arr) dposit(a,w1,w2,t,arr),tfinal,posit(k,:),opts);
        end
        Nphot=Nphot+1;
        temp=positsol{k}.y(:,end);
        tempet=eta(a,w1,w2,positsol{k}.x(end),temp(1),temp(2),temp(3),velx(temp(4),temp(5),temp(6)),vely(temp(4),temp(5),temp(6)),velz(temp(4),temp(5),temp(6)));
        temprand=rand;
        chif=fzero(@(lim) integral(@(chi)spect(chi,a,w1,w2,positsol{k}.x(end),temp(1),temp(2),temp(3),velx(temp(4),temp(5),temp(6)),vely(temp(4),temp(5),temp(6)),velz(temp(4),temp(5),temp(6))),0,lim)-temprand*positt(a,w1,w2,positsol{k}.x(end),temp(1),temp(2),temp(3),temp(4),temp(5),temp(6)),[eps,tempet/2-eps]);
        modp=2*chif*sqrt(1+temp(4)^2+temp(5)^2+temp(6)^2)/tempet;
        transfer=modp*normalize(temp(4:6),'norm');
        phot=[phot;[temp(1),temp(2),temp(3),transfer(1),transfer(2),transfer(3),0]];
        tinphot=[tinphot;positsol{k}.x(end)];
        posit(k,:)=[temp(1),temp(2),temp(3),temp(4)-transfer(1),temp(5)-transfer(2),temp(6)-transfer(3),0];
    end
    for k=1:Nphot
        photonenergy=sqrt(phot(k,4)^2+phot(k,5)^2+phot(k,6)^2);
        target=-log(rand);
        opts=odeset('Events',@(t,arr) eve(t,arr,target));
        if(k>size(photsol,2))
            photsol{k}=ode45(@(t,arr) dphot(a,w1,w2,t,arr),[tinphot(k) tfinal],phot(k,:),opts);
        else
            photsol{k}=odextend(photsol{k},@(t,arr) dphot(a,w1,w2,t,arr),tfinal,phot(k,:),opts);
        end
        if(photonenergy>2)
            Nelec=Nelec+1;
            Nposit=Nposit+1;
            temp=photsol{k}.y(:,end);
            tempchig=chig(a,w1,w2,photsol{k}.x(end),temp(1),temp(2),temp(3),temp(4)/photonenergy,temp(5)/photonenergy,temp(6)/photonenergy,photonenergy/hbar);
            temprand=rand;
            chimf=fzero(@(lim) integral(@(chim) dNdtdchi(a,w1,w2,photsol{k}.x(end),temp(1),temp(2),temp(3),temp(4)/photonenergy,temp(5)/photonenergy,temp(6)/photonenergy,photonenergy/hbar,chim),0,lim,'ArrayValued',1)-temprand*phott(a,w1,w2,photsol{k}.x(end),temp(1),temp(2),temp(3),temp(4),temp(5),temp(6)),[eps,tempchig-eps]);
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
            phot(k,:)=[temp(1),temp(2),temp(3),temp(4)-(transfere(1)+transferp(1)),temp(5)-(transfere(2)+transferp(2)),temp(6)-(transfere(3)+trandferp(3)),0];
        end
    end
end

function out=xforce(a0,w1,w2,t,x,~,~,~,~,vz)
out=a0*vz*(w2*sin(w2*(t+x))-w1*sin(w1*(t-x)))/2;
end
function out=yforce(~,~,~,~,~,~,~,~,~,~)
out=0;
end
function out=zforce(a0,w1,w2,t,x,~,~,vx,~,~)
out=-a0*((1+vx)*w2*sin(w2*(t+x))+(1-vx)*w1*sin(w1*(t-x)))/2;
end
function out=chig(a0,w1,w2,t,x,~,~,vx,~,vz,freq) %Photon velocity and freq!
global hbar;
out=(hbar^2)*freq*a0*sqrt(((1+vx)*w2*sin(w2*(t+x)))^2+((1-vx)*w1*sin(w1*(t-x)))^2+2*w1*w2*(1-vx^2-2*vz^2)*sin(w1*(t-x))*sin(w2*(t+x)))/2;
end
function out=xcomb(a0,w1,w2,t,x,y,z,vx,vy,vz,freq,chim) %Photon velocity and freq!
inter=chig(a0,w1,w2,t,x,y,z,vx,vy,vz,freq);
out=((inter/(chim*(inter-chim)))^2)^(1/3);
end
function out=dNdtdchi(a0,w1,w2,t,x,y,z,vx,vy,vz,freq,chim) %Photon velocity and freq!
global hbar;
interchig=chig(a0,w1,w2,t,x,y,z,vx,vy,vz,freq);
interx=xcomb(a0,w1,w2,t,x,y,z,vx,vy,vz,freq,chim);
out=(integral(@(s)airy(s),interx,inf)+(2-interchig*(interx^(3/2)))*airy(1,interx)/interx)/(hbar^3*freq*interchig);
end
function out=eta(a0,w1,w2,t,x,~,~,vx,vy,vz)
global hbar;
out=hbar*a0*sqrt(((1+vx)*w2*sin(w2*(t+x)))^2+((1-vx)*w1*sin(w1*(t-x)))^2+2*w1*w2*(1-vx^2-2*vz^2)*sin(w1*(t-x))*sin(w2*(t+x)))/(2*sqrt(1-vx^2-vy^2-vz^2));
end
function out=y(et,chi)
out=2*chi./(3*et*(et-2*chi));
end
function out=F(et,chi)
% out=4*chi.^2.*y(et,chi).*besselk(2/3,y(et,chi))/et^2+(1-2*chi/et).*y(et,chi).*arrayfun(@(lim)integral(@(t) besselk(5/3,t),lim,inf),y(et,chi));
out=0;
if(chi<et/2)
    if(chi>0)
        out=4*chi.^2.*y(et,chi).*besselk(2/3,y(et,chi))/et^2+(1-2*chi/et).*y(et,chi).*arrayfun(@(lim)integral(@(t) besselk(5/3,t),lim,inf),y(et,chi));
    end
end
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
function out=elecx(~,~,~,~,~,~,~,px,py,pz)
out=velx(px,py,pz);
end
function out=elecy(~,~,~,~,~,~,~,px,py,pz)
out=vely(px,py,pz);
end
function out=elecz(~,~,~,~,~,~,~,px,py,pz)
out=velz(px,py,pz);
end
function out=positx(~,~,~,~,~,~,~,px,py,pz)
out=velx(px,py,pz);
end
function out=posity(~,~,~,~,~,~,~,px,py,pz)
out=vely(px,py,pz);
end
function out=positz(~,~,~,~,~,~,~,px,py,pz)
out=velz(px,py,pz);
end
function out=photx(~,~,~,~,~,~,~,px,py,pz)
out=px/sqrt(px^2+py^2+pz^2);
end
function out=photy(~,~,~,~,~,~,~,px,py,pz)
out=py/sqrt(px^2+py^2+pz^2);
end
function out=photz(~,~,~,~,~,~,~,px,py,pz)
out=pz/sqrt(px^2+py^2+pz^2);
end
function out=elecpx(a0,w1,w2,t,x,~,~,px,py,pz)
out=xforce(a0,w1,w2,t,x,0,0,velx(px,py,pz),vely(px,py,pz),velz(px,py,pz));
end
function out=elecpy(a0,w1,w2,t,x,~,~,px,py,pz)
out=yforce(a0,w1,w2,t,x,0,0,velx(px,py,pz),vely(px,py,pz),velz(px,py,pz));
end
function out=elecpz(a0,w1,w2,t,x,~,~,px,py,pz)
out=zforce(a0,w1,w2,t,x,0,0,velx(px,py,pz),vely(px,py,pz),velz(px,py,pz));
end
function out=positpx(a0,w1,w2,t,x,~,~,px,py,pz)
out=-xforce(a0,w1,w2,t,x,0,0,velx(px,py,pz),vely(px,py,pz),velz(px,py,pz));
end
function out=positpy(a0,w1,w2,t,x,~,~,px,py,pz)
out=-yforce(a0,w1,w2,t,x,0,0,velx(px,py,pz),vely(px,py,pz),velz(px,py,pz));
end
function out=positpz(a0,w1,w2,t,x,~,~,px,py,pz)
out=-zforce(a0,w1,w2,t,x,0,0,velx(px,py,pz),vely(px,py,pz),velz(px,py,pz));
end
function out=photpx(~,~,~,~,~,~,~,~,~,~)
out=0;
end
function out=photpy(~,~,~,~,~,~,~,~,~,~)
out=0;
end
function out=photpz(~,~,~,~,~,~,~,~,~,~)
out=0;
end
function out=spect(chi,a0,w1,w2,t,x,~,~,vx,vy,vz)
global hbar;
n=eta(a0,w1,w2,t,x,0,0,vx,vy,vz);
out=sqrt(3)*sqrt(1-vx^2-vy^2-vz^2)*n*F(n,chi)./(2*pi*chi*hbar^2);
if(chi==n/2)
    out=0;
end
end

function out=elect(a0,w1,w2,t,x,~,~,px,py,pz)
tempet=eta(a0,w1,w2,t,x,0,0,velx(px,py,pz),vely(px,py,pz),velz(px,py,pz));
out=integral(@(chi)spect(chi,a0,w1,w2,t,x,0,0,velx(px,py,pz),vely(px,py,pz),velz(px,py,pz)),0, tempet/2);
if(tempet==0)
    out=0;
end
end
function out=positt(a0,w1,w2,t,x,y,z,px,py,pz)
tempet=eta(a0,w1,w2,t,x,0,0,velx(px,py,pz),vely(px,py,pz),velz(px,py,pz));
out=integral(@(chi)spect(chi,a0,w1,w2,t,x,y,z,velx(px,py,pz),vely(px,py,pz),velz(px,py,pz)),0, tempet/2);
if(tempet==0)
    out=0;
end
end
function out=phott(a0,w1,w2,t,x,y,z,px,py,pz)%Photon momenta!
global hbar;
intervx=px/sqrt(px^2+py^2+pz^2);
intervy=py/sqrt(px^2+py^2+pz^2);
intervz=pz/sqrt(px^2+py^2+pz^2);
freq=sqrt(px^2+py^2+pz^2)/hbar;
interchig=chig(a0,w1,w2,t,x,y,z,intervx,intervy,intervz,freq);
out=integral(@(chim) dNdtdchi(a0,w1,w2,t,x,y,z,intervx,intervy,intervz,freq,chim),0,interchig,'ArrayValued',1);
end
function out=delec(a0,w1,w2,t,arr)
    out=[elecx(a0,w1,w2,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6)),elecy(a0,w1,w2,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6)),elecz(a0,w1,w2,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6)),elecpx(a0,w1,w2,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6)),elecpy(a0,w1,w2,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6)),elecpz(a0,w1,w2,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6)),elect(a0,w1,w2,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6))];
    out=out(:);
end
function out=dposit(a0,w1,w2,t,arr)
    out=[positx(a0,w1,w2,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6)),posity(a0,w1,w2,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6)),positz(a0,w1,w2,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6)),positpx(a0,w1,w2,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6)),positpy(a0,w1,w2,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6)),positpz(a0,w1,w2,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6)),positt(a0,w1,w2,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6))];
    out=out(:);
end
function out=dphot(a0,w1,w2,t,arr)
    out=[photx(a0,w1,w2,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6)),photy(a0,w1,w2,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6)),photz(a0,w1,w2,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6)),photpx(a0,w1,w2,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6)),photpy(a0,w1,w2,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6)),photpz(a0,w1,w2,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6)),phott(a0,w1,w2,t,arr(1),arr(2),arr(3),arr(4),arr(5),arr(6))];
    out=out(:);
end
function [value,isterminal,direction]=eve(~,arr,target)
    value=arr(7)-target;
    isterminal=1;
    direction=0;%-1: Where event function is decreasing
end