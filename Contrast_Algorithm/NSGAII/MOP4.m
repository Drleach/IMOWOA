%%  闲鱼：深度学习与智能算法
%%  唯一官方店铺：https://mbd.pub/o/author-aWWbm3BtZw==
%%  微信公众号：强盛机器学习，关注公众号获得更多免费代码！
function z=MOP4(x)

    a=0.8;
    
    b=3;
    
    z1=sum(-10*exp(-0.2*sqrt(x(1:end-1).^2+x(2:end).^2)));
    
    z2=sum(abs(x).^a+5*(sin(x)).^b);
    
    z=[z1 z2]';

end