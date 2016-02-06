function monoClarinet(M,Fs,At,Av,Dt,Dv,St,Rt)
    y=0;
    for i=1:length(M(:,1))
        t=M(i,5):1/Fs:M(i,6);
        f= (2^((M(i,3)-69)/12))*(440);
        Ae=linspace(0,Av,At*length(t));
        De=linspace(Av,Dv,Dt*length(t));
        Se=linspace(Dv,Dv,St*length(t));
        Re=linspace(Dv,0,Rt*length(t));
        ADSR=[Ae De Se Re];
        ADSR=padarray(ADSR,[0 (length(t)-length(ADSR))],0,'post');
        w=2*pi*f*t;
        a=M(i,4)/127;
        y = [y (a*(sin(w)+0.75*sin(3*w)+0.5*sin(5*w)+0.14*sin(7*w)+0.5*sin(9*w)+0.14*sin(11*w)+0.17*sin(13*w))).*ADSR(1,:)];
    end
    figure(10);
    plot(ADSR);
    soundsc(y,Fs);
end