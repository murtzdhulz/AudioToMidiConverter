function polyClarinet(M,Fs,At,Av,Dt,Dv,St,Rt)
    tt = M(1,5):1/Fs:M(length(M(:,6)),6);
    y = zeros(1,length(tt));
    for i=1:length(M(:,1))
        t=M(i,5):1/Fs:M(i,6);
        ind=round(M(i,5)*Fs);
        f= (2^((M(i,3)-69)/12))*(440);
        Ae=linspace(0,Av,At*length(t));
        De=linspace(Av,Dv,Dt*length(t));
        Se=linspace(Dv,Dv,St*length(t));
        Re=linspace(Dv,0,Rt*length(t));
        ADSR=[Ae De Se Re];
        ADSR=padarray(ADSR,[0 (length(t)-length(ADSR))],0,'post');
        w=2*pi*f*t;
        a=M(i,4)/127;
        z=(a*(sin(w)+0.75*sin(3*w)+0.5*sin(5*w)+0.14*sin(7*w)+0.5*sin(9*w)+0.14*sin(11*w)+0.17*sin(13*w))).*ADSR(1,:);
        for j=1:length(t)
            y(j+ind) = y(j+ind)+z(j);
        end
    end
    figure(10);
    plot(ADSR);
    soundsc(y,Fs);
end