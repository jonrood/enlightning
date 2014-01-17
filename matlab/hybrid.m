function hybrid(timeSteps,X,Y,srcPosX,srcPosY,...
    srcAmp,srcAlpha,maxFreq,ppw,cfl,bc,rc)

    %timeSteps is number of time steps to run
    %X,Y are number of cells of domain
    %srcPosX,srcPosY are the center location of source
    %srcAmp is density amplitude of source
    %srcAlpha is spread of point source
    %maxFreq is maximum frequency to account for
    %ppw is points per width of maxFreq wavelength
    %cfl condtion
    %bc, 1 for reflecting boundaries, 2 for absorbing
    %rc is hybrid threshold for subscheme switching
    
    global w_n w p T H0 tag;

    alpha=360;
    c=343;
    mygamma=1.402;
    p0=101325;
    p_ref=p0;
    RH=0.2;
    T0=293.16;
    T_ref=T0;
    rho0=1.21;
    R=287.06;
    c_v=720.4;
    c_p=1010;
    T_star_n=3352;
    T_star_o=2239;
    mu=1.846e-5;
    muB=0.6*mu;
    kappa=0.02624;
    
    w=zeros(X,Y,6);
    H0=zeros(X,Y);
    p=zeros(X,Y);
    T=zeros(X,Y);
    tag=zeros(X,Y);

    lambda=c/maxFreq;
    dx=lambda/ppw;
    dy=lambda/ppw;
    dt=cfl*(dx/c);
    
    zeta=10.79586*(1-(273.16/T0))-5.02808*...
        log10(T0/273.16)+(1.50474e-4)*...
        (1-10^(-8.29692*((T0/273.16)-1)))-...
        (4.2873e-4)*(1-10^(-4.76955*...
        ((273.16/T0)-1)))-2.2195983;
    p_vp=p_ref*10^zeta;
    h=((10^(-2))*RH*p_vp*T0)/p0;
    eta=(4.17)*(((T_ref/T0)^(1/3))-1);
    f_n=(p0/p_ref)*(24+(4.04e6)*h*...
        ((0.02+100*h)/(0.391+100*h)));
    f_o=(p0/p_ref)*(((T_ref/T0)^(1/2))*...
        (9+2.8e4*h*exp(-eta)));
    tau_n=1/(2*pi*f_n);
    tau_o=1/(2*pi*f_o);
    
    if dt>((tau_n)/4)
        dt=(tau_n)/4;
    end
  
    p(:,:)=p0;
    T(:,:)=T0; 
    w(:,:,1)=rho0;
    w(:,:,2)=0;
    w(:,:,3)=0;
    w(:,:,4)=0;
    w(:,:,5)=rho0*T0;
    w(:,:,6)=rho0*T0;
    w_n=w;
    
    updateState(p0,T0,rho0,R,c_v,c_p,mygamma,c);
    
    addSource(X,Y,srcPosX,srcPosY,...
        srcAmp,srcAlpha);
    
    rungeKutta(X,Y,dx,dy,muB,mu,...
            tau_n,tau_o,T_star_n,T_star_o,R,kappa,...
            alpha,dt,bc,rc,p0,T0,rho0,mygamma,c,c_v,c_p);
    H0(:,:)=0;      
        
    for step=2:timeSteps
        rungeKutta(X,Y,dx,dy,muB,mu,...
            tau_n,tau_o,T_star_n,T_star_o,R,kappa,...
            alpha,dt,bc,rc,p0,T0,rho0,mygamma,c,c_v,c_p);

        if mod(step,50)==0
            %subplot(1,2,1);
            surf(p-p0);
            axis([1 Y 1 X -srcAmp/10 srcAmp/10]);
            xlabel('y');
            ylabel('x');
            title('p');

            %subplot(1,2,2);
            %image(tag.*256);
            %xlabel('y');
            %ylabel('x');
            %title('tags');
            drawnow;
        end
    end
end

function rungeKutta(X,Y,dx,dy,muB,mu,...
    tau_n,tau_o,T_star_n,T_star_o,R,kappa,alpha,...
    dt,bc,rc,p0,T0,rho0,mygamma,c,c_v,c_p)
    global w w_n;
    
    K=rhs(X,Y,dx,dy,muB,mu,...
    tau_n,tau_o,T_star_n,T_star_o,R,kappa,alpha,rc);
    w_n=w+dt*K;
    updateState(p0,T0,rho0,R,c_v,c_p,mygamma,c);
    
    K=rhs(X,Y,dx,dy,muB,mu,...
    tau_n,tau_o,T_star_n,T_star_o,R,kappa,alpha,rc);
    w_n=(3/4)*w+(1/4)*w_n+(1/4)*dt*K;
    updateState(p0,T0,rho0,R,c_v,c_p,mygamma,c);
    
    K=rhs(X,Y,dx,dy,muB,mu,...
    tau_n,tau_o,T_star_n,T_star_o,R,kappa,alpha,rc);
    w_n=(1/3)*w+(2/3)*w_n+(2/3)*dt*K;
    updateState(p0,T0,rho0,R,c_v,c_p,mygamma,c);
    
    if bc==1
        bcReflect(X,Y);
    elseif bc==2
        bcAbsorb(X,Y,p0,T0,rho0);
    end
    
    w=w_n;
end

function K=rhs(X,Y,dx,dy,muB,mu,...
    tau_n,tau_o,T_star_n,T_star_o,R,kappa,alpha,rc)
    global w_n H0 p T tag;
    
    K=zeros(X,Y,6);
    F=zeros(7,6);
    G=zeros(7,6);
    H=zeros(1,6);
    tag(:,:)=0;
    xi=1;
    ep=0.9*rc*xi*xi/(1-0.9*rc*rc);
    
    for i=4:X-3
        for j=4:Y-3
            ix=4;
            for q=-3:3
               F(ix+q,1)=w_n(i+q,j,1)*...
                   (w_n(i+q,j,2)/w_n(i+q,j,1));
               G(ix+q,1)=w_n(i,j+q,1)*...
                   (w_n(i,j+q,3)/w_n(i,j+q,1));

               F(ix+q,2)=w_n(i+q,j,1)*...
                   ((w_n(i+q,j,2)/w_n(i+q,j,1))^2)+...
                   p(i+q,j);
               G(ix+q,2)=w_n(i,j+q,1)*...
                   (w_n(i,j+q,3)/w_n(i,j+q,1))*...
                   (w_n(i,j+q,2)/w_n(i,j+q,1));

               F(ix+q,3)=w_n(i+q,j,1)*...
                   (w_n(i+q,j,2)/w_n(i+q,j,1))*...
                   (w_n(i+q,j,3)/w_n(i+q,j,1));
               G(ix+q,3)=w_n(i,j+q,1)*...
                   ((w_n(i,j+q,3)/w_n(i,j+q,1))^2)+...
                   p(i,j+q);

               F(ix+q,4)=w_n(i+q,j,1)*...
                   (w_n(i+q,j,2)/w_n(i+q,j,1))*...
                   (w_n(i+q,j,4)/w_n(i+q,j,1));
               G(ix+q,4)=w_n(i,j+q,1)*...
                   (w_n(i,j+q,3)/w_n(i,j+q,1))*...
                   (w_n(i,j+q,4)/w_n(i,j+q,1));

               F(ix+q,5)=w_n(i+q,j,1)*...
                   (w_n(i+q,j,2)/w_n(i+q,j,1))*...
                   (w_n(i+q,j,5)/w_n(i+q,j,1));
               G(ix+q,5)=w_n(i,j+q,1)*...
                   (w_n(i,j+q,3)/w_n(i,j+q,1))*...
                   (w_n(i,j+q,5)/w_n(i,j+q,1));

               F(ix+q,6)=w_n(i+q,j,1)*...
                   (w_n(i+q,j,2)/w_n(i+q,j,1))*...
                   (w_n(i+q,j,6)/w_n(i+q,j,1));
               G(ix+q,6)=w_n(i,j+q,1)*...
                   (w_n(i,j+q,3)/w_n(i,j+q,1))*...
                   (w_n(i,j+q,6)/w_n(i,j+q,1));
            end
    
            H(1)=H0(i,j);

            u_dprime_x=((w_n(i-1,j,2)/w_n(i-1,j,1))-...
                2*(w_n(i,j,2)/w_n(i,j,1))+...
                (w_n(i+1,j,2)/w_n(i+1,j,1)))/(dx^2);
            u_dprime_y=((w_n(i,j-1,2)/w_n(i,j-1,1))-...
                2*(w_n(i,j,2)/w_n(i,j,1))+...
                (w_n(i,j+1,2)/w_n(i,j+1,1)))/(dy^2);
            v_dprime_xy=...
                ((w_n(i+1,j+1,3)/w_n(i+1,j+1,1))-...
                (w_n(i+1,j-1,3)/w_n(i+1,j-1,1))-...
                (w_n(i-1,j+1,3)/w_n(i-1,j+1,1))+...
                (w_n(i-1,j-1,3)/w_n(i-1,j-1,1)))/...
                (4*dx*dy);
            phi_xx_prime_x=(4.0/3.0)*...
                u_dprime_x-(2.0/3.0)*v_dprime_xy;
            phi_xy_prime_y=u_dprime_y+v_dprime_xy;

            H(2)=muB*(u_dprime_x+v_dprime_xy)+...
                mu*(phi_xx_prime_x+phi_xy_prime_y);

            v_dprime_y=((w_n(i,j-1,3)/w_n(i,j-1,1))-...
                2*(w_n(i,j,3)/w_n(i,j,1))+...
                (w_n(i,j+1,3)/w_n(i,j+1,1)))/(dy^2);
            v_dprime_x=((w_n(i-1,j,3)/w_n(i-1,j,1))-...
                2*(w_n(i,j,3)/w_n(i,j,1))+...
                (w_n(i+1,j,3)/w_n(i+1,j,1)))/(dx^2);
            u_dprime_yx=...
                ((w_n(i+1,j+1,2)/w_n(i+1,j+1,1))-...
                (w_n(i+1,j-1,2)/w_n(i+1,j-1,1))-...
                (w_n(i-1,j+1,2)/w_n(i-1,j+1,1))+...
                (w_n(i-1,j-1,2)/w_n(i-1,j-1,1)))/...
                (4*dx*dy);
            phi_yx_prime_x=u_dprime_yx+v_dprime_x;
            phi_yy_prime_y=(4.0/3.0)*...
                v_dprime_y-(2.0/3.0)*u_dprime_yx;

            H(3)=muB*(v_dprime_y+u_dprime_yx)+...
                mu*(phi_yx_prime_x+phi_yy_prime_y);

            u_prime_x=((w_n(i+1,j,2)/w_n(i+1,j,1))-...
                (w_n(i-1,j,2)/w_n(i-1,j,1)))/(2*dx);
            u_prime_y=((w_n(i,j+1,2)/w_n(i,j+1,1))-...
                (w_n(i,j-1,2)/w_n(i,j-1,1)))/(2*dy);
            v_prime_y=((w_n(i,j+1,3)/w_n(i,j+1,1))-...
                (w_n(i,j-1,3)/w_n(i,j-1,1)))/(2*dy);
            v_prime_x=((w_n(i+1,j,3)/w_n(i+1,j,1))-...
                (w_n(i-1,j,3)/w_n(i-1,j,1)))/(2*dx);
            t_prime_x=(T(i+1,j)-T(i-1,j))/(2*dx);
            t_prime_y=(T(i,j+1)-T(i,j-1))/(2*dy);
            t_dprime_x=(T(i-1,j)-2*T(i,j)+...
                T(i+1,j))/(dx^2);
            t_dprime_y=(T(i,j-1)-2*T(i,j)+...
                T(i,j+1))/(dy^2);
            T_n=w_n(i,j,5)/w_n(i,j,1);
            T_o=w_n(i,j,6)/w_n(i,j,1);
            Mn=(1.0/tau_n)*(T(i,j)-T_n);
            Mo=(1.0/tau_o)*(T(i,j)-T_o);
            tmp1=T_star_n/(T_n);
            tmp2=T_star_o/(T_o);
            c_vn=0.78*R*(tmp1^2)*exp(-tmp1);
            c_vo=0.21*R*(tmp2^2)*exp(-tmp2);
            A_n=((T(i,j)/T_n)-1)*c_vn;
            A_o=((T(i,j)/T_o)-1)*c_vo;
            phi_xx=(4.0/3.0)*u_prime_x-(2.0/3.0)*v_prime_y;
            phi_xy_yx=u_prime_y+v_prime_x;
            phi_yy=(4.0/3.0)*v_prime_y-(2.0/3.0)*u_prime_x;
            sigma=(muB/T(i,j))*((u_prime_x+v_prime_y)^2)+...
                (mu/(2*T(i,j)))*...
                (phi_xx^2+2*(phi_xy_yx^2)+phi_yy^2)+...
                (kappa/(T(i,j)^2))*...
                (t_prime_x^2+t_prime_y^2)+...
                (w_n(i,j,1)/T(i,j))*...
                (A_n*Mn+A_o*Mo);
            H(4)=sigma-...
                ((w_n(i,j,1)/T_n)*c_vn*Mn+...
                (w_n(i,j,1)/T_o)*c_vo*Mo)+...
                (-((kappa*(t_prime_x^2))/(T(i,j)^2))+...
                ((kappa*t_dprime_x)/T(i,j))-...
                ((kappa*(t_prime_y^2))/(T(i,j)^2))+...
                ((kappa*t_dprime_y)/T(i,j)));

            H(5)=(w_n(i,j,1)/tau_n)*(T(i,j)-T_n);

            H(6)=(w_n(i,j,1)/tau_o)*(T(i,j)-T_o);

            tmp1=w_n(i,j,1)-w_n(i-1,j,1);
            tmp2=w_n(i+1,j,1)-w_n(i,j,1);
            rx=(2*abs(tmp1*tmp2)+ep)/(tmp1^2+tmp2^2+ep);
            if (rx<0.9999)
                tag(i,j)=1;
            end

            tmp1=w_n(i,j,1)-w_n(i,j-1,1);
            tmp2=w_n(i,j+1,1)-w_n(i,j,1);
            ry=(2*abs(tmp1*tmp2)+ep)/(tmp1^2+tmp2^2+ep);
            if (ry<0.9999)
                tag(i,j)=1;
            end

            if tag(i,j)==0
                ix=4;
                a_1=0.79926643;
                a_2=-0.18941314;
                a_3=0.02651995;
                for k=1:6
                    K(i,j,k)=...
                        (-a_3*F(ix-3,k)-a_2*F(ix-2,k)-...
                        a_1*F(ix-1,k)+a_1*F(ix+1,k)+...
                        a_2*F(ix+2,k)+a_3*F(ix+3,k))/-dx...
                        +(-a_3*G(ix-3,k)-a_2*G(ix-2,k)-...
                        a_1*G(ix-1,k)+a_1*G(ix+1,k)+...
                        a_2*G(ix+2,k)+a_3*G(ix+3,k))/-dy...
                        +H(k);
                end
            else
                for k=1:6
                    ii=i-1;
                    ix=4-1;
                    fgm2=0.5*(F(ix-2,k)+alpha*w_n(ii-2,j,k));
                    fgm1=0.5*(F(ix-1,k)+alpha*w_n(ii-1,j,k));
                    fg0 =0.5*(F(ix  ,k)+alpha*w_n(ii  ,j,k));
                    fgp1=0.5*(F(ix+1,k)+alpha*w_n(ii+1,j,k));
                    fgp2=0.5*(F(ix+2,k)+alpha*w_n(ii+2,j,k));
                    mf1=wenoM(fgm2,fgm1,fg0,fgp1,fgp2);

                    ii=ii+1;
                    ix=ix+1;
                    fgm2=0.5*(F(ix-2,k)-alpha*w_n(ii-2,j,k));
                    fgm1=0.5*(F(ix-1,k)-alpha*w_n(ii-1,j,k));
                    fg0 =0.5*(F(ix  ,k)-alpha*w_n(ii  ,j,k));
                    fgp1=0.5*(F(ix+1,k)-alpha*w_n(ii+1,j,k));
                    fgp2=0.5*(F(ix+2,k)-alpha*w_n(ii+2,j,k));
                    pf1=wenoP(fgm2,fgm1,fg0,fgp1,fgp2);

                    fgm2=0.5*(F(ix-2,k)+alpha*w_n(ii-2,j,k));
                    fgm1=0.5*(F(ix-1,k)+alpha*w_n(ii-1,j,k));
                    fg0 =0.5*(F(ix  ,k)+alpha*w_n(ii  ,j,k));
                    fgp1=0.5*(F(ix+1,k)+alpha*w_n(ii+1,j,k));
                    fgp2=0.5*(F(ix+2,k)+alpha*w_n(ii+2,j,k));
                    mf2=wenoM(fgm2,fgm1,fg0,fgp1,fgp2);

                    ii=ii+1;
                    ix=ix+1;
                    fgm2=0.5*(F(ix-2,k)-alpha*w_n(ii-2,j,k));
                    fgm1=0.5*(F(ix-1,k)-alpha*w_n(ii-1,j,k));
                    fg0 =0.5*(F(ix  ,k)-alpha*w_n(ii  ,j,k));
                    fgp1=0.5*(F(ix+1,k)-alpha*w_n(ii+1,j,k));
                    fgp2=0.5*(F(ix+2,k)-alpha*w_n(ii+2,j,k));
                    pf2=wenoP(fgm2,fgm1,fg0,fgp1,fgp2);

                    jj=j-1;
                    ix=4-1;
                    fgm2=0.5*(G(ix-2,k)+alpha*w_n(i,jj-2,k));
                    fgm1=0.5*(G(ix-1,k)+alpha*w_n(i,jj-1,k));
                    fg0 =0.5*(G(ix  ,k)+alpha*w_n(i,jj  ,k));
                    fgp1=0.5*(G(ix+1,k)+alpha*w_n(i,jj+1,k));
                    fgp2=0.5*(G(ix+2,k)+alpha*w_n(i,jj+2,k));
                    mg1=wenoM(fgm2,fgm1,fg0,fgp1,fgp2);

                    jj=jj+1;
                    ix=ix+1;
                    fgm2=0.5*(G(ix-2,k)-alpha*w_n(i,jj-2,k));
                    fgm1=0.5*(G(ix-1,k)-alpha*w_n(i,jj-1,k));
                    fg0 =0.5*(G(ix  ,k)-alpha*w_n(i,jj  ,k));
                    fgp1=0.5*(G(ix+1,k)-alpha*w_n(i,jj+1,k));
                    fgp2=0.5*(G(ix+2,k)-alpha*w_n(i,jj+2,k));
                    pg1=wenoP(fgm2,fgm1,fg0,fgp1,fgp2);

                    fgm2=0.5*(G(ix-2,k)+alpha*w_n(i,jj-2,k));
                    fgm1=0.5*(G(ix-1,k)+alpha*w_n(i,jj-1,k));
                    fg0 =0.5*(G(ix  ,k)+alpha*w_n(i,jj  ,k));
                    fgp1=0.5*(G(ix+1,k)+alpha*w_n(i,jj+1,k));
                    fgp2=0.5*(G(ix+2,k)+alpha*w_n(i,jj+2,k));
                    mg2=wenoM(fgm2,fgm1,fg0,fgp1,fgp2);

                    jj=jj+1;
                    ix=ix+1;
                    fgm2=0.5*(G(ix-2,k)-alpha*w_n(i,jj-2,k));
                    fgm1=0.5*(G(ix-1,k)-alpha*w_n(i,jj-1,k));
                    fg0 =0.5*(G(ix  ,k)-alpha*w_n(i,jj  ,k));
                    fgp1=0.5*(G(ix+1,k)-alpha*w_n(i,jj+1,k));
                    fgp2=0.5*(G(ix+2,k)-alpha*w_n(i,jj+2,k));
                    pg2=wenoP(fgm2,fgm1,fg0,fgp1,fgp2);

                    K(i,j,k)=(mf1+pf1-mf2-pf2)/dx+...
                        (mg1+pg1-mg2-pg2)/dy+H(k);
                end
            end
        end
    end 
end

function vp=wenoP(fgm2,fgm1,fg0,fgp1,fgp2)
	ep=10e-6;
	
	beta0=(13/12)*(fg0-2*fgp1+fgp2)^2+...
        (1/4)*(3*fg0-4*fgp1+fgp2)^2;
	beta1=(13/12)*(fgm1-2*fg0+fgp1)^2+...
        (1/4)*(fgm1-fgp1)^2;
	beta2=(13/12)*(fgm2-2*fgm1+fg0)^2+...
        (1/4)*(fgm2-4*fgm1+3*fg0)^2;
	
	zeta0=0.1/((ep+beta0)^2);
	zeta1=0.6/((ep+beta1)^2);
	zeta2=0.3/((ep+beta2)^2);
	
	sigma=zeta0+zeta1+zeta2;
	
	omega0=zeta0/sigma;
	omega1=zeta1/sigma;
	omega2=zeta2/sigma;
    
    f0=(1/6)*(11*fg0-7*fgp1+2*fgp2);
	f1=(1/6)*(2*fgm1+5*fg0-fgp1);
	f2=(1/6)*(-fgm2+5*fgm1+2*fg0);
	
	vp=omega0*f0+omega1*f1+omega2*f2;
end

function vm=wenoM(fgm2,fgm1,fg0,fgp1,fgp2)
	ep=10e-6;
	
    beta0=(13/12)*(fg0-2*fgp1+fgp2)^2+...
        (1/4)*(3*fg0-4*fgp1+fgp2)^2;
    beta1=(13/12)*(fgm1-2*fg0+fgp1)^2+...
        (1/4)*(fgm1-fgp1)^2;
    beta2=(13/12)*(fgm2-2*fgm1+fg0)^2+...
        (1/4)*(fgm2-4*fgm1+3*fg0)^2;
	
	zeta0=0.3/((ep+beta0)^2);
	zeta1=0.6/((ep+beta1)^2);
	zeta2=0.1/((ep+beta2)^2);
	
	sigma=zeta0+zeta1+zeta2;
	
	omega0=zeta0/sigma;
	omega1=zeta1/sigma;
	omega2=zeta2/sigma;
    
    f0=(1/6)*(2*fg0+5*fgp1-fgp2);
    f1=(1/6)*(-fgm1+5*fg0+2*fgp1);
    f2=(1/6)*(2*fgm2-7*fgm1+11*fg0);
	
	vm=omega0*f0+omega1*f1+omega2*f2;
end

function updateState(p0,T0,rho0,R,c_v,c_p,mygamma,c)
    global w_n p T;
    
    p=(c^2).*((w_n(:,:,1)-rho0)+((mygamma-1)...
        ./(2.*rho0)).*((w_n(:,:,1)-rho0).^2)...
       +(w_n(:,:,1).*(T./T0)./c_p)...
        .*(w_n(:,:,4)./w_n(:,:,1)))+p0;
    T=T0.*exp((w_n(:,:,4)./w_n(:,:,1)...
       -R.*log(rho0./w_n(:,:,1)))./c_v);
end

function addSource(X,Y,srcPosX,srcPosY,...
    srcAmp,srcAlpha)
    global H0;

    for i=1:X
        for j=1:Y
            H0(i,j)=srcAmp*exp(-log(2)/((srcAlpha)^2)...
               *(((i-srcPosX)^2)+((j-srcPosY)^2)));
        end
    end
end

function bcReflect(X,Y)
    global w_n;
    
    w_n(1,:,:)=w_n(4,:,:);
    w_n(2,:,:)=w_n(4,:,:);
    w_n(3,:,:)=w_n(4,:,:);
    w_n(3,:,2)=-w_n(4,:,2);
    
    w_n(:,1,:)=w_n(:,4,:);
    w_n(:,2,:)=w_n(:,4,:);
    w_n(:,3,:)=w_n(:,4,:);
    w_n(:,3,3)=-w_n(:,4,3);
    
    w_n(X,:,:)=w_n(X-3,:,:);
    w_n(X-1,:,:)=w_n(X-3,:,:);
    w_n(X-2,:,:)=w_n(X-3,:,:);
    w_n(X-2,:,2)=-w_n(X-3,:,2);
    
    w_n(:,Y,:)=w_n(:,Y-3,:);
    w_n(:,Y-1,:)=w_n(:,Y-3,:);
    w_n(:,Y-2,:)=w_n(:,Y-3,:);
    w_n(:,Y-2,3)=-w_n(:,Y-3,3);
end

function bcAbsorb(X,Y,p0,T0,rho0)
    global w_n p T;
    
    numpoints=8;
    x_s=1; y_s=1;
    envAlpha=4;
    
    for i=1:numpoints
        tmp=(-exp(-log(2)/(envAlpha^2)*(i-1-x_s)^2)+1);
        w_n(i,:,1)=(w_n(i,:,1)-rho0)*tmp+rho0;
        w_n(i,:,2)=w_n(i,:,2)*tmp;
        w_n(i,:,3)=w_n(i,:,3)*tmp;
        w_n(i,:,4)=w_n(i,:,4)*tmp;
        w_n(i,:,5)=(w_n(i,:,5)-T0*rho0)*tmp+T0*rho0;
        w_n(i,:,6)=(w_n(i,:,6)-T0*rho0)*tmp+T0*rho0;
        p(i,:)=(p(i,:)-p0)*tmp+p0;
        T(i,:)=(T(i,:)-T0)*tmp+T0;
    end
    
    for j=1:numpoints
        tmp=(-exp(-log(2)/(envAlpha^2)*(j-1-y_s)^2)+1);
        w_n(:,j,1)=(w_n(:,j,1)-rho0)*tmp+rho0;
        w_n(:,j,2)=w_n(:,j,2)*tmp;
        w_n(:,j,3)=w_n(:,j,3)*tmp;
        w_n(:,j,4)=w_n(:,j,4)*tmp;
        w_n(:,j,5)=(w_n(:,j,5)-T0*rho0)*tmp+T0*rho0;
        w_n(:,j,6)=(w_n(:,j,6)-T0*rho0)*tmp+T0*rho0;
        p(:,j)=(p(:,j)-p0)*tmp+p0;
        T(:,j)=(T(:,j)-T0)*tmp+T0;
    end

    for i=X:-1:X-numpoints+1
        tmp=(-exp(-log(2)/(envAlpha^2)*(X-i-x_s)^2)+1);
        w_n(i,:,1)=(w_n(i,:,1)-rho0)*tmp+rho0;
        w_n(i,:,2)=w_n(i,:,2)*tmp;
        w_n(i,:,3)=w_n(i,:,3)*tmp;
        w_n(i,:,4)=w_n(i,:,4)*tmp;
        w_n(i,:,5)=(w_n(i,:,5)-T0*rho0)*tmp+T0*rho0;
        w_n(i,:,6)=(w_n(i,:,6)-T0*rho0)*tmp+T0*rho0;
        p(i,:)=(p(i,:)-p0)*tmp+p0;
        T(i,:)=(T(i,:)-T0)*tmp+T0;
    end
    
    for j=Y:-1:Y-numpoints+1
        tmp=(-exp(-log(2)/(envAlpha^2)*(Y-j-y_s)^2)+1);
        w_n(:,j,1)=(w_n(:,j,1)-rho0)*tmp+rho0;
        w_n(:,j,2)=w_n(:,j,2)*tmp;
        w_n(:,j,3)=w_n(:,j,3)*tmp;
        w_n(:,j,4)=w_n(:,j,4)*tmp;
        w_n(:,j,5)=(w_n(:,j,5)-T0*rho0)*tmp+T0*rho0;
        w_n(:,j,6)=(w_n(:,j,6)-T0*rho0)*tmp+T0*rho0;
        p(:,j)=(p(:,j)-p0)*tmp+p0;
        T(:,j)=(T(:,j)-T0)*tmp+T0;
    end
end