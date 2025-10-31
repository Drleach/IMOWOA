%% 清理之前的数据
% 清除所有数据
clear all;
% 清除窗口输出
clc;
cec20 = CEC20_MM_Functions();
% 测试12个测试函数能否正常运行
for fun=1:24
    [func_name,obj_num,dim,xl,xu,repoint,N_ops,fobj] = cec20.get_function(fun);
    input = unifrnd(xl,xu);
    disp([fun,input]);
    a = fobj(input);
    disp(size(a));

end