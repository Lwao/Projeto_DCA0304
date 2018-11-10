function F=nlsistema(s)
    x=s(1);
    y=s(2);
    z=s(3);
    
    
    %F(1) = 1/2*sin(x*y)-(y/(4*pi))-(x/2)
    %F(2) = (1-1/(4*pi))*((exp(2*x))-exp(1))-((exp(1)*y)/pi)-2*exp(1)*x
    
    F(1) = y+z-exp(-x)
    F(2) = x+z-exp(-z)
    F(3) = x+y-exp(-z)
    
end
