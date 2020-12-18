function [f, g, h] = LSC(x, I_fno)

% final checked

[ps, D]=size(x);
global initial_flag
persistent o p s w m c
persistent M M_10 M_30 M_50 M_100 M2_10 M2_30 M2_50 M2_100

%% ------------------  Fully-separable Functions --------------------
%------------------------------------------------------------------------------
% f1:  D = 1000
%------------------------------------------------------------------------------
if(I_fno == 1)
    if initial_flag == 0
        load Function01
        initial_flag = 1;
    end
    z = x - repmat(o,size(x,1),1);
    f = sum(z.*sin(z),2);
    h1 = sum(z-100*cos(0.5*z)+100,2);
    h2 = sum(-z+100*cos(0.5*z)-100,2);
    h = [h1,h2];
    g = zeros(ps,1);
end

%% ------------  Partially Additive Separable Functions I -----------
%------------------------------------------------------------------------------
% f2: Shifted and Rotated Rastrigin's Function
% D = 1000
%------------------------------------------------------------------------------
if (I_fno == 2)
    if initial_flag == 0
        load Function02
        initial_flag = 1;
    end
    y = x - repmat(o,size(x,1),1);
    f = 0;
    g1 = 0;
    g2 = 0;
    ldim = 1;
    for i=1:length(s)
        if s(i) == 10
            z = (M_10 * y(:, p(ldim:ldim+s(i)-1))')';
            fit = rastrigin(z);
            ldim = ldim + s(i);
        elseif s(i) == 30
            z = (M_30 * y(:, p(ldim:ldim+s(i)-1))')';
            fit = rastrigin(z);
            ldim = ldim + s(i);
        elseif s(i) == 50
            z = (M_50 * y(:, p(ldim:ldim+s(i)-1))')';
            fit = rastrigin(z);
            ldim = ldim + s(i);
        elseif s(i) == 100
            z = (M_100 * y(:, p(ldim:ldim+s(i)-1))')';
            fit = rastrigin(z);
            ldim = ldim + s(i);
        end
        f = f + w(i)*fit;
        g1 = g1 + (-sum(abs(z),2)) + 4; 
        g2 = g2 + sum(z.^2,2) - 4;
    end
    z = y(:, p(ldim:end));
    f = f + rastrigin(z);
    g1 = g1 + (-sum(abs(z), 2)) + 4; 
    g2 = g2 + sum(z.^2, 2) - 4;
    g = [g1, g2];
    h =zeros(ps,1);
end

%------------------------------------------------------------------------------
% f3: Shifted and Rotated Ackley's Function
% D = 1000
%------------------------------------------------------------------------------
if (I_fno == 3)
    if initial_flag == 0
        load Function03
        initial_flag = 1;
    end
    y = x - repmat(o,size(x,1),1);
    f = 0;
    g = 0;
    h = 0;
    ldim = 1;
    for i=1:length(s)
        if s(i) == 10
            z = (M_10 * y(:, p(ldim:ldim+s(i)-1))')';
            ldim = ldim + s(i);
        elseif s(i) == 30
            z = (M_30 * y(:, p(ldim:ldim+s(i)-1))')';
            ldim = ldim + s(i);
        elseif s(i) == 50
            z = (M_50 * y(:, p(ldim:ldim+s(i)-1))')';
            ldim = ldim + s(i);
        elseif s(i) == 100
            z = (M_100 * y(:, p(ldim:ldim+s(i)-1))')';
            ldim = ldim + s(i);
        end
        if i <= 7
            fit = ackley(z);
            f = f + w(i)*fit;
        end
        h = h + sum(z.^2,2) - 4;
        g = g - abs(z(:,1)) + sum(z(:,2:end).^2,2) + 1;
    end
    z = y(:, p(sum(s(1:7))+1:end));
    f = f + ackley(z);
end

%% -----------  Partially Additive Separable Functions II -----------
%------------------------------------------------------------------------------
% f4: Maximum Function
% D = 1000
%------------------------------------------------------------------------------
if (I_fno == 4)
    if initial_flag == 0
        load Function04
        initial_flag = 1;
    end
    y = x - repmat(o,size(x,1),1);
    f = 0;
    h = 0;
    ldim = 1;
    for i=1:length(s)
        z = y(:, p(ldim:ldim+s(i)-1));
        fit = max(z, [], 2);
        ldim = ldim + s(i);
        f = f + w(i)*fit;
        hh = schwefel12(z);
        h = h + hh;
    end
    g = zeros(ps,1);
end

%------------------------------------------------------------------------------
% f5:  Schwefel¡¯s Problem 1.2
% D = 1000
%------------------------------------------------------------------------------
if (I_fno == 5)
    if initial_flag == 0
        load Function05
        initial_flag = 1;
    end
    y = x - repmat(o,size(x,1),1);
    f = 0;
    g = 0;
    ldim = 1;
    for i=1:length(s)
        if s(i) == 10
            z = (M_10 * y(:, p(ldim:ldim+s(i)-1))')';
            fit = schwefel12(z);
            ldim = ldim + s(i);
        elseif s(i) == 30
            z = (M_30 * y(:, p(ldim:ldim+s(i)-1))')';
            fit = schwefel12(z);
            ldim = ldim + s(i);
        elseif s(i) == 50
            z = (M_50 * y(:, p(ldim:ldim+s(i)-1))')';
            fit = schwefel12(z);
            ldim = ldim + s(i);
        elseif s(i) == 100
            z = (M_100 * y(:, p(ldim:ldim+s(i)-1))')';
            fit = schwefel12(z);
            ldim = ldim + s(i);
        end
        f = f + w(i)*fit;
        g = g + sum(z.^2-5000.*cos(0.1.*pi.*z)-4000,2);
    end
    h = zeros(ps,1);
end

%------------------------------------------------------------------------------
% f6: Rosenbrock's Function
% D = 1000
%------------------------------------------------------------------------------
if (I_fno == 6)
    if initial_flag == 0
        load Function06
        initial_flag = 1;
    end
    y = x - repmat(o,size(x,1),1);
    f = 0;
    g1 = 0;
    g2 = 0;
    ldim = 1;
    for i=1:length(s)
        if s(i) == 10
            z = y(:, p(ldim:ldim+s(i)-1));
            fit = rosenbrock(z);
            u = (M_10*z')';
            v = (M2_10*z')';
            ldim = ldim + s(i);
        elseif s(i) == 30
            z = y(:, p(ldim:ldim+s(i)-1));
            fit = rosenbrock(z);
            u = (M_30*z')';
            v = (M2_30*z')';
            ldim = ldim + s(i);
        elseif s(i) == 50
            z = y(:, p(ldim:ldim+s(i)-1));
            fit = rosenbrock(z);
            u = (M_50*z')';
            v = (M2_50*z')';
            ldim = ldim + s(i);
        elseif s(i) == 100
            z = y(:, p(ldim:ldim+s(i)-1));
            fit = rosenbrock(z);
            u = (M_100*z')';
            v = (M2_100*z')';
            ldim = ldim + s(i);
        end
        f = f + w(i)*fit;
        g1 = g1 + sum(u.^2-50.*cos(2.*pi.*u)-40,2);
        g2 = g2 + sum(v.^2-50.*cos(2.*pi.*v)-40,2);
    end
    g = [g1, g2];
    h = zeros(ps,1);
end

%------------------------------------------------------------------------------
% f7: Absolute Function
% D = 1000
%------------------------------------------------------------------------------
if(I_fno == 7)
    if initial_flag == 0
        load Function07
        initial_flag = 1;
    end
    y = x - repmat(o,size(x,1),1);
    f = 0;
    g = 0;
    h = 0;
    ldim = 1;
    for i=1:length(s)
        if s(i) == 10
            z = (M_10 * y(:, p(ldim:ldim+s(i)-1))')';
            fit = absolute(z);
            ldim = ldim + s(i);
        elseif s(i) == 30
            z = (M_30 * y(:, p(ldim:ldim+s(i)-1))')';
            fit = absolute(z);
            ldim = ldim + s(i);
        elseif s(i) == 50
            z = (M_50 * y(:, p(ldim:ldim+s(i)-1))')';
            fit = absolute(z);
            ldim = ldim + s(i);
        elseif s(i) == 100
            z = (M_100 * y(:, p(ldim:ldim+s(i)-1))')';
            fit = absolute(z);
            ldim = ldim + s(i);
        end
        f = f + w(i)*fit;
        g = g + sum(z.^2,2);
        h = h + (cos(fit) + sin(fit)).^2-exp(cos(fit) + sin(fit))-1+exp(1);
    end
    g = g - 100 * D;
end

%------------------------------------------------------------------------------
% f8: Schwefel Function
% D = 1000
%------------------------------------------------------------------------------
if(I_fno == 8)
    if initial_flag == 0
        load Function08
        initial_flag = 1;
    end
    y = x - repmat(o,size(x,1),1);
    f = 0;
    g1 = 0;
    g2 = 0;
    g3 = 0;
    ldim = 1;
    for i=1:length(s)
        if s(i) == 10
            z = (M_10 * y(:, p(ldim:ldim+s(i)-1))')';
            fit = schwefel(z);
            ldim = ldim + s(i);
        elseif s(i) == 30
            z = (M_30 * y(:, p(ldim:ldim+s(i)-1))')';
            fit = schwefel(z);
            ldim = ldim + s(i);
        elseif s(i) == 50
            z = (M_50 * y(:, p(ldim:ldim+s(i)-1))')';
            fit = schwefel(z);
            ldim = ldim + s(i);
        elseif s(i) == 100
            z = (M_100 * y(:, p(ldim:ldim+s(i)-1))')';
            fit = schwefel(z);
            ldim = ldim + s(i);
        end
        f = f + w(i)*fit;
        g1 = g1 + sum(z.^2, 2);
        g2 = g2 + sum(sin((1/50).*pi.*z), 2);
        gg = 1;
        for j=1:s(i)
            gg = gg.*(cos(z(:,j)/sqrt(j)));
        end
        g3 = g3 + sum(z.^2,2)/4000 - gg + 1;
    end
    g1 = 1/(100*D) * g1 - 50;
    g2 = 50/D * g2;
    g3 = -g3 + 30;
    g = [g1, g2, g3];
    h = zeros(ps,1);
end

%% --------------------  Overlapping Functions  ---------------------

%------------------------------------------------------------------------------
% f9: D = 1000
%------------------------------------------------------------------------------
if (I_fno == 9)
    if initial_flag == 0
        load Function09
        initial_flag = 1;
    end
    z = x - repmat(o,size(x,1),1);
    f = rosenbrock(z);
    g1 = rastrigin(z)-100;
    g2 = sum(z,2)-2*D;
    g3 = -sum(z,2) + 5;
    g = [g1,g2,g3];
    h =zeros(ps,1);
end

%------------------------------------------------------------------------------
% f10: D = 900
%------------------------------------------------------------------------------
if(I_fno == 10)
    if initial_flag == 0
        load Function10
        c = cumsum(s);
        initial_flag = 1;
    end
    y = x - repmat(o,size(x,1),1);
    f = 0;
    h1 = 0;
    h2 = 0;
    for i=1:length(s)
        if i == 1
            ldim = 1;
        else
            ldim = c(i-1) - ((i-1)*m) + 1;
        end
        udim = c(i) - ((i-1)*m);
        z = y(:, p(ldim:udim));
        f = f + max(z,[],2);
        h1 = h1 + schwefel12(z);
        h2 = h2 + sum((z(:,1:s(i)-1)-z(:,2:s(i))).^2, 2);
    end
    h = [h1, h2];
    g = zeros(ps,1);
end

%% ----------------  Fully Non-Separable Functions ------------------

%------------------------------------------------------------------------------
% f11: D = 1000
% Objective function is nonseperable
% Constraint function is separable
%------------------------------------------------------------------------------
if(I_fno == 11)
    if initial_flag == 0
        load Function11
        initial_flag = 1;
    end
    y = x - repmat(o,size(x,1),1);
    g1 = 0;
    g2 = 0;
    g3 = 0;
    ldim = 1;
    for i=1:length(s)
        if s(i) == 10
            u = (M_10 * y(:, p(ldim:ldim+s(i)-1))')';
            ldim = ldim + s(i);
        elseif s(i) == 30
            u = (M_30 * y(:, p(ldim:ldim+s(i)-1))')';
            ldim = ldim + s(i);
        elseif s(i) == 50
            u = (M_50 * y(:, p(ldim:ldim+s(i)-1))')';
            ldim = ldim + s(i);
        elseif s(i) == 100
            u = (M_100 * y(:, p(ldim:ldim+s(i)-1))')';
            ldim = ldim + s(i);
        end
        g1 = g1 + rastrigin(u);
        g2 = g2 + sum(u,2);
    end
    g1 = 1/100 * g1 - 21;
    g3 = 105 - g2;
    g2 = g2 - 2*D;
    z = (M*y')';
    f = rosenbrock(z);
    g = [g1, g2, g3];
    h = zeros(ps,1);
end

%------------------------------------------------------------------------------
% f12: D = 1000
%------------------------------------------------------------------------------
if(I_fno == 12)
    if initial_flag == 0
        load Function12
        initial_flag = 1;
    end
    z = x - repmat(o,size(x,1),1);
    z = (M*z')';
    f  = rastrigin(z);
    g1 = sum(-z.*sin(2*z),2);
    g2 = sum(z.*sin(z),2);
    g = [g1,g2];
    h =zeros(ps,1);
end

end

%------------------------------------------------------------------------------
% Rastrigin's Function
%------------------------------------------------------------------------------
function fit = rastrigin(x)
    fit = sum(x.^2-10.*cos(2.*pi.*x)+10,2);
end

%------------------------------------------------------------------------------
% Ackley's Function
%------------------------------------------------------------------------------
function fit = ackley(x)
    [~, D]=size(x);
    f = sum(x.^2,2);
    fit = 20-20.*exp(-0.2.*sqrt(f./D))-exp(sum(cos(2.*pi.*x),2)./D)+exp(1);
end

%------------------------------------------------------------------------------
% Schwefel's Problem 1.2
%------------------------------------------------------------------------------
function fit = schwefel12(x)
    [~, D]=size(x);
    fit = 0;
    for i=1:D
        fit = fit + sum(x(:,1:i),2).^2;
    end
end

%------------------------------------------------------------------------------
% Rosenbrock's Function
%------------------------------------------------------------------------------
function fit = rosenbrock(x)
    [~, D]=size(x);
    fit = sum(100*(x(:,1:D-1).^2-x(:,2:D)).^2+(x(:,1:D-1)-1).^2,2);
end

%------------------------------------------------------------------------------
% Absolute Function
%------------------------------------------------------------------------------
function fit = absolute(x)
    fit = sum(abs(x), 2);
end

%------------------------------------------------------------------------------
% Schwefel Function
%------------------------------------------------------------------------------
function fit = schwefel(x)
    [~, D]=size(x);
    f = x.*sin(sqrt(abs(x)));
    fit = 418.9829*D - sum(f, 2);
end