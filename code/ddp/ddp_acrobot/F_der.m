function [Fx, Fu, Fxx, Fux, Fxu] = F_der(x, u, param)
    x1 = x(1);
    x2 = x(2);
    x3 = x(3);
    x4 = x(4);

    T =  param.T;
    a1 = param.a1;
    a2 = param.a2;
    a3 = param.a3;
    a4 = param.a4;
    a5 = param.a5;
    b1 = param.b1;
    b2 = param.b2;

    Fx = NaN(4, 4);
    Fu = NaN(4, 1);
    Fxx = NaN(4,4,4);
    Fux = NaN(4,1,4);

    %Fx
    Fx(1,1) = 1;
    Fx(1,2) = 0;
    Fx(1,3) = T;
    Fx(1,4) = 0;
    Fx(2,1) = 0;
    Fx(2,2) = 1;
    Fx(2,3) = 0;
    Fx(2,4) = T;
    Fx(3,1) = (T*(a3*a4*cos(x1) - a2*a5*cos(x1 + x2)*cos(x2)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2);
    Fx(3,2) = (2*T*a2^2*cos(x2)*sin(x2)*(a3*u + a3*a4*sin(x1) + a2*u*cos(x2) - (a2^2*x3^2*sin(2*x2))/2 + a3*b1*x3 - a3*b2*x4 - a2*b2*x4*cos(x2) + 2*a2*a3*x3*sin(x2) + a2*a3*x4*sin(x2) - a2*a5*sin(x1 + x2)*cos(x2) - a2*a3*x3^2*sin(x2)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2)^2 - (2*T*a2*(u*sin(x2) + a5*cos(x1 + 2*x2) - 2*a3*x3*cos(x2) - a3*x4*cos(x2) - b2*x4*sin(x2) + a3*x3^2*cos(x2) + a2*x3^2*cos(2*x2)))/(a2^2*cos(2*x2) - 2*a1*a3 + a2^2 + 2*a3^2);
    Fx(3,3) = (T*(a3*b1 + 2*a2*a3*sin(x2) - 2*a2^2*x3*cos(x2)*sin(x2) - 2*a2*a3*x3*sin(x2)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2) + 1;
    Fx(3,4) = -(T*(a3*b2 - a2*a3*sin(x2) + a2*b2*cos(x2)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2);
    Fx(4,1) = -(T*(a3*a5*cos(x1 + x2) - a1*a5*cos(x1 + x2) + a3*a4*cos(x1) - a2*a5*cos(x1 + x2)*cos(x2) + a2*a4*cos(x1)*cos(x2)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2);
    Fx(4,2) = (2*T*(2*a2*u*sin(x2) + 2*a2^2*x3^2*cos(2*x2) + (a2*a4*cos(x1 - x2))/2 + a2*a5*cos(x1 + 2*x2) - 2*a2^2*x3*cos(2*x2) - a2^2*x4*cos(2*x2) + a1*a5*cos(x1 + x2) - (a2*a4*cos(x1 + x2))/2 - a3*a5*cos(x1 + x2) - 2*a2*a3*x3*cos(x2) - a2*a3*x4*cos(x2) + a2*b1*x3*sin(x2) - 2*a2*b2*x4*sin(x2) + a1*a2*x3^2*cos(x2)))/(a2^2*cos(2*x2) - 2*a1*a3 + a2^2 + 2*a3^2) - (2*T*a2^2*cos(x2)*sin(x2)*(a1*u + a3*a4*sin(x1) + 2*a2*u*cos(x2) - a2^2*x3^2*sin(2*x2) - a1*b2*x4 + a3*b1*x3 + a2^2*x3*sin(2*x2) + (a2^2*x4*sin(2*x2))/2 - a1*a5*sin(x1 + x2) + a3*a5*sin(x1 + x2) + a2*b1*x3*cos(x2) - 2*a2*b2*x4*cos(x2) + 2*a2*a3*x3*sin(x2) + a2*a3*x4*sin(x2) - a2*a5*sin(x1 + x2)*cos(x2) - a1*a2*x3^2*sin(x2) + a2*a4*cos(x2)*sin(x1)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2)^2;
    Fx(4,3) = -(T*(a3*b1 + a2^2*sin(2*x2) + 2*a2*a3*sin(x2) - 2*a2^2*x3*sin(2*x2) + a2*b1*cos(x2) - 2*a1*a2*x3*sin(x2)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2);
    Fx(4,4) = (T*(a1*b2 - a2*a3*sin(x2) - a2^2*cos(x2)*sin(x2) + 2*a2*b2*cos(x2)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2) + 1;

    %Fu
    Fu(1,1) = 0;
    Fu(2,1) = 0;
    Fu(3,1) = (T*(a3 + a2*cos(x2)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2);
    Fu(4,1) = -(T*(a1 + 2*a2*cos(x2)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2);

    %Fxx
    Fxx(1,1,1) = 0;
    Fxx(1,1,2) = 0;
    Fxx(1,1,3) = 0;
    Fxx(1,1,4) = 0;
    Fxx(1,2,1) = 0;
    Fxx(1,2,2) = 0;
    Fxx(1,2,3) = 0;
    Fxx(1,2,4) = 0;
    Fxx(1,3,1) = 0;
    Fxx(1,3,2) = 0;
    Fxx(1,3,3) = 0;
    Fxx(1,3,4) = 0;
    Fxx(1,4,1) = 0;
    Fxx(1,4,2) = 0;
    Fxx(1,4,3) = 0;
    Fxx(1,4,4) = 0;
    Fxx(2,1,1) = 0;
    Fxx(2,1,2) = 0;
    Fxx(2,1,3) = 0;
    Fxx(2,1,4) = 0;
    Fxx(2,2,1) = 0;
    Fxx(2,2,2) = 0;
    Fxx(2,2,3) = 0;
    Fxx(2,2,4) = 0;
    Fxx(2,3,1) = 0;
    Fxx(2,3,2) = 0;
    Fxx(2,3,3) = 0;
    Fxx(2,3,4) = 0;
    Fxx(2,4,1) = 0;
    Fxx(2,4,2) = 0;
    Fxx(2,4,3) = 0;
    Fxx(2,4,4) = 0;
    Fxx(3,1,1) = -(T*(a3*a4*sin(x1) - a2*a5*sin(x1 + x2)*cos(x2)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2);
    Fxx(3,1,2) = (2*T*a2*a5*sin(x1 + 2*x2))/(a2^2*cos(2*x2) - 2*a1*a3 + a2^2 + 2*a3^2) + (2*T*a2^2*cos(x2)*sin(x2)*(a3*a4*cos(x1) - a2*a5*cos(x1 + x2)*cos(x2)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2)^2;
    Fxx(3,1,3) = 0;
    Fxx(3,1,4) = 0;
    Fxx(3,2,1) = (2*T*a2*a5*sin(x1 + 2*x2))/(a2^2*cos(2*x2) - 2*a1*a3 + a2^2 + 2*a3^2) + (2*T*a2^2*cos(x2)*sin(x2)*(a3*a4*cos(x1) - a2*a5*cos(x1 + x2)*cos(x2)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2)^2;
    Fxx(3,2,2) = (2*T*a2*(2*a5*sin(x1 + 2*x2) - u*cos(x2) + b2*x4*cos(x2) - 2*a3*x3*sin(x2) - a3*x4*sin(x2) + a3*x3^2*sin(x2) + 2*a2*x3^2*sin(2*x2)))/(a2^2*cos(2*x2) - 2*a1*a3 + a2^2 + 2*a3^2) - (2*T*a2^2*sin(x2)^2*(a3*u + a3*a4*sin(x1) + a2*u*cos(x2) - (a2^2*x3^2*sin(2*x2))/2 + a3*b1*x3 - a3*b2*x4 - a2*b2*x4*cos(x2) + 2*a2*a3*x3*sin(x2) + a2*a3*x4*sin(x2) - a2*a5*sin(x1 + x2)*cos(x2) - a2*a3*x3^2*sin(x2)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2)^2 - (4*T*a2^3*sin(2*x2)*(u*sin(x2) + a5*cos(x1 + 2*x2) - 2*a3*x3*cos(x2) - a3*x4*cos(x2) - b2*x4*sin(x2) + a3*x3^2*cos(x2) + a2*x3^2*cos(2*x2)))/(a2^2*cos(2*x2) - 2*a1*a3 + a2^2 + 2*a3^2)^2 + (2*T*a2^2*cos(x2)^2*(a3*u + a3*a4*sin(x1) + a2*u*cos(x2) - (a2^2*x3^2*sin(2*x2))/2 + a3*b1*x3 - a3*b2*x4 - a2*b2*x4*cos(x2) + 2*a2*a3*x3*sin(x2) + a2*a3*x4*sin(x2) - a2*a5*sin(x1 + x2)*cos(x2) - a2*a3*x3^2*sin(x2)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2)^2 + (8*T*a2^4*cos(x2)^2*sin(x2)^2*(a3*u + a3*a4*sin(x1) + a2*u*cos(x2) - (a2^2*x3^2*sin(2*x2))/2 + a3*b1*x3 - a3*b2*x4 - a2*b2*x4*cos(x2) + 2*a2*a3*x3*sin(x2) + a2*a3*x4*sin(x2) - a2*a5*sin(x1 + x2)*cos(x2) - a2*a3*x3^2*sin(x2)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2)^3 - (2*T*a2^3*cos(x2)*sin(x2)*(u*sin(x2) + a5*cos(x1 + 2*x2) - 2*a3*x3*cos(x2) - a3*x4*cos(x2) - b2*x4*sin(x2) + a3*x3^2*cos(x2) + a2*x3^2*cos(2*x2)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2)^2;
    Fxx(3,2,3) = (2*T*a2^2*cos(x2)*sin(x2)*(a3*b1 + 2*a2*a3*sin(x2) - 2*a2^2*x3*cos(x2)*sin(x2) - 2*a2*a3*x3*sin(x2)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2)^2 - (2*T*a2*(a3*x3*cos(x2) - a3*cos(x2) + a2*x3*(2*cos(x2)^2 - 1)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2);
    Fxx(3,2,4) = (T*a2*(a3*cos(x2) + b2*sin(x2)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2) - (2*T*a2^2*cos(x2)*sin(x2)*(a3*b2 - a2*a3*sin(x2) + a2*b2*cos(x2)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2)^2;
    Fxx(3,3,1) = 0;
    Fxx(3,3,2) = (2*T*a2^2*cos(x2)*sin(x2)*(a3*b1 + 2*a2*a3*sin(x2) - 2*a2^2*x3*cos(x2)*sin(x2) - 2*a2*a3*x3*sin(x2)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2)^2 - (2*T*a2*(a3*x3*cos(x2) - a3*cos(x2) + a2*x3*(2*cos(x2)^2 - 1)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2);
    Fxx(3,3,3) = -(T*a2*(2*a3*sin(x2) + 2*a2*cos(x2)*sin(x2)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2);
    Fxx(3,3,4) = 0;
    Fxx(3,4,1) = 0;
    Fxx(3,4,2) = (T*a2*(a3*cos(x2) + b2*sin(x2)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2) - (2*T*a2^2*cos(x2)*sin(x2)*(a3*b2 - a2*a3*sin(x2) + a2*b2*cos(x2)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2)^2;
    Fxx(3,4,3) = 0;
    Fxx(3,4,4) = 0;
    Fxx(4,1,1) = (T*(a3*a4*sin(x1) - a1*a5*sin(x1 + x2) + a3*a5*sin(x1 + x2) - a2*a5*sin(x1 + x2)*cos(x2) + a2*a4*cos(x2)*sin(x1)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2);
    Fxx(4,1,2) = - (T*(a2*a4*sin(x1 - x2) + 2*a2*a5*sin(x1 + 2*x2) + 2*a1*a5*sin(x1 + x2) - a2*a4*sin(x1 + x2) - 2*a3*a5*sin(x1 + x2)))/(a2^2*cos(2*x2) - 2*a1*a3 + a2^2 + 2*a3^2) - (2*T*a2^2*cos(x2)*sin(x2)*(a3*a5*cos(x1 + x2) - a1*a5*cos(x1 + x2) + a3*a4*cos(x1) - a2*a5*cos(x1 + x2)*cos(x2) + a2*a4*cos(x1)*cos(x2)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2)^2;
    Fxx(4,1,3) = 0;
    Fxx(4,1,4) = 0;
    Fxx(4,2,1) = - (2*T*((a2*a4*sin(x1 - x2))/2 + a2*a5*sin(x1 + 2*x2) + a1*a5*sin(x1 + x2) - (a2*a4*sin(x1 + x2))/2 - a3*a5*sin(x1 + x2)))/(a2^2*cos(2*x2) - 2*a1*a3 + a2^2 + 2*a3^2) - (2*T*a2^2*cos(x2)*sin(x2)*(a3*a5*cos(x1 + x2) - a1*a5*cos(x1 + x2) + a3*a4*cos(x1) - a2*a5*cos(x1 + x2)*cos(x2) + a2*a4*cos(x1)*cos(x2)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2)^2;
    Fxx(4,2,2) = (2*T*(2*a2*u*cos(x2) - 4*a2^2*x3^2*sin(2*x2) + (a2*a4*sin(x1 - x2))/2 - 2*a2*a5*sin(x1 + 2*x2) + 4*a2^2*x3*sin(2*x2) + 2*a2^2*x4*sin(2*x2) - a1*a5*sin(x1 + x2) + (a2*a4*sin(x1 + x2))/2 + a3*a5*sin(x1 + x2) + a2*b1*x3*cos(x2) - 2*a2*b2*x4*cos(x2) + 2*a2*a3*x3*sin(x2) + a2*a3*x4*sin(x2) - a1*a2*x3^2*sin(x2)))/(a2^2*cos(2*x2) - 2*a1*a3 + a2^2 + 2*a3^2) + (4*T*a2^2*sin(2*x2)*(2*a2*u*sin(x2) + 2*a2^2*x3^2*cos(2*x2) + (a2*a4*cos(x1 - x2))/2 + a2*a5*cos(x1 + 2*x2) - 2*a2^2*x3*cos(2*x2) - a2^2*x4*cos(2*x2) + a1*a5*cos(x1 + x2) - (a2*a4*cos(x1 + x2))/2 - a3*a5*cos(x1 + x2) - 2*a2*a3*x3*cos(x2) - a2*a3*x4*cos(x2) + a2*b1*x3*sin(x2) - 2*a2*b2*x4*sin(x2) + a1*a2*x3^2*cos(x2)))/(a2^2*cos(2*x2) - 2*a1*a3 + a2^2 + 2*a3^2)^2 - (2*T*a2^2*cos(x2)^2*(a1*u + a3*a4*sin(x1) + 2*a2*u*cos(x2) - a2^2*x3^2*sin(2*x2) - a1*b2*x4 + a3*b1*x3 + a2^2*x3*sin(2*x2) + (a2^2*x4*sin(2*x2))/2 - a1*a5*sin(x1 + x2) + a3*a5*sin(x1 + x2) + a2*b1*x3*cos(x2) - 2*a2*b2*x4*cos(x2) + 2*a2*a3*x3*sin(x2) + a2*a3*x4*sin(x2) - a2*a5*sin(x1 + x2)*cos(x2) - a1*a2*x3^2*sin(x2) + a2*a4*cos(x2)*sin(x1)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2)^2 + (2*T*a2^2*sin(x2)^2*(a1*u + a3*a4*sin(x1) + 2*a2*u*cos(x2) - a2^2*x3^2*sin(2*x2) - a1*b2*x4 + a3*b1*x3 + a2^2*x3*sin(2*x2) + (a2^2*x4*sin(2*x2))/2 - a1*a5*sin(x1 + x2) + a3*a5*sin(x1 + x2) + a2*b1*x3*cos(x2) - 2*a2*b2*x4*cos(x2) + 2*a2*a3*x3*sin(x2) + a2*a3*x4*sin(x2) - a2*a5*sin(x1 + x2)*cos(x2) - a1*a2*x3^2*sin(x2) + a2*a4*cos(x2)*sin(x1)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2)^2 - (8*T*a2^4*cos(x2)^2*sin(x2)^2*(a1*u + a3*a4*sin(x1) + 2*a2*u*cos(x2) - a2^2*x3^2*sin(2*x2) - a1*b2*x4 + a3*b1*x3 + a2^2*x3*sin(2*x2) + (a2^2*x4*sin(2*x2))/2 - a1*a5*sin(x1 + x2) + a3*a5*sin(x1 + x2) + a2*b1*x3*cos(x2) - 2*a2*b2*x4*cos(x2) + 2*a2*a3*x3*sin(x2) + a2*a3*x4*sin(x2) - a2*a5*sin(x1 + x2)*cos(x2) - a1*a2*x3^2*sin(x2) + a2*a4*cos(x2)*sin(x1)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2)^3 + (2*T*a2^2*cos(x2)*sin(x2)*(2*a2*u*sin(x2) + 2*a2^2*x3^2*cos(2*x2) + (a2*a4*cos(x1 - x2))/2 + a2*a5*cos(x1 + 2*x2) - 2*a2^2*x3*cos(2*x2) - a2^2*x4*cos(2*x2) + a1*a5*cos(x1 + x2) - (a2*a4*cos(x1 + x2))/2 - a3*a5*cos(x1 + x2) - 2*a2*a3*x3*cos(x2) - a2*a3*x4*cos(x2) + a2*b1*x3*sin(x2) - 2*a2*b2*x4*sin(x2) + a1*a2*x3^2*cos(x2)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2)^2;
    Fxx(4,2,3) = (2*T*a2*(b1*sin(x2) - 2*a3*cos(x2) - 2*a2*(2*cos(x2)^2 - 1) + 2*a1*x3*cos(x2) + 4*a2*x3*(2*cos(x2)^2 - 1)))/(2*a2^2*cos(x2)^2 - 2*a1*a3 + 2*a3^2) - (2*T*a2^2*cos(x2)*sin(x2)*(a3*b1 + a2^2*sin(2*x2) + 2*a2*a3*sin(x2) - 2*a2^2*x3*sin(2*x2) + a2*b1*cos(x2) - 2*a1*a2*x3*sin(x2)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2)^2;
    Fxx(4,2,4) = (2*T*a2^2*cos(x2)*sin(x2)*(a1*b2 - a2*a3*sin(x2) - a2^2*cos(x2)*sin(x2) + 2*a2*b2*cos(x2)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2)^2 - (T*a2*(a3*cos(x2) - a2 + 2*b2*sin(x2) + 2*a2*cos(x2)^2))/(a2^2*cos(x2)^2 - a1*a3 + a3^2);
    Fxx(4,3,1) = 0;
    Fxx(4,3,2) = (T*a2*(b1*sin(x2) - 2*a3*cos(x2) - 2*a2*(2*cos(x2)^2 - 1) + 2*a1*x3*cos(x2) + 4*a2*x3*(2*cos(x2)^2 - 1)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2) - (2*T*a2^2*cos(x2)*sin(x2)*(a3*b1 + a2^2*sin(2*x2) + 2*a2*a3*sin(x2) - 2*a2^2*x3*sin(2*x2) + a2*b1*cos(x2) - 2*a1*a2*x3*sin(x2)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2)^2;
    Fxx(4,3,3) = (2*T*a2*(a1*sin(x2) + 2*a2*cos(x2)*sin(x2)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2);
    Fxx(4,3,4) = 0;
    Fxx(4,4,1) = 0;
    Fxx(4,4,2) = (2*T*a2^2*cos(x2)*sin(x2)*(a1*b2 - a2*a3*sin(x2) - a2^2*cos(x2)*sin(x2) + 2*a2*b2*cos(x2)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2)^2 - (T*a2*(a3*cos(x2) - a2 + 2*b2*sin(x2) + 2*a2*cos(x2)^2))/(a2^2*cos(x2)^2 - a1*a3 + a3^2);
    Fxx(4,4,3) = 0;
    Fxx(4,4,4) = 0;

    %Fux
    Fux(1,1,1) = 0;
    Fux(1,1,2) = 0;
    Fux(1,1,3) = 0;
    Fux(1,1,4) = 0;
    Fux(2,1,1) = 0;
    Fux(2,1,2) = 0;
    Fux(2,1,3) = 0;
    Fux(2,1,4) = 0;
    Fux(3,1,1) = 0;
    Fux(3,1,2) = (T*a2*sin(x2)*(a1*a3 + a2^2*cos(x2)^2 - a3^2 + 2*a2*a3*cos(x2)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2)^2;
    Fux(3,1,3) = 0;
    Fux(3,1,4) = 0;
    Fux(4,1,1) = 0;
    Fux(4,1,2) = -(2*T*a2*sin(x2)*(a3 + a2*cos(x2))*(a1 - a3 + a2*cos(x2)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2)^2;
    Fux(4,1,3) = 0;
    Fux(4,1,4) = 0;


    %Fxu
    Fxu(1,1,1) = 0;
    Fxu(1,2,1) = 0;
    Fxu(1,3,1) = 0;
    Fxu(1,4,1) = 0;
    Fxu(2,1,1) = 0;
    Fxu(2,2,1) = 0;
    Fxu(2,3,1) = 0;
    Fxu(2,4,1) = 0;
    Fxu(3,1,1) = 0;
    Fxu(3,2,1) = (T*a2*sin(x2)*(a1*a3 + a2^2*cos(x2)^2 - a3^2 + 2*a2*a3*cos(x2)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2)^2;
    Fxu(3,3,1) = 0;
    Fxu(3,4,1) = 0;
    Fxu(4,1,1) = 0;
    Fxu(4,2,1) = -(2*T*a2*sin(x2)*(a3 + a2*cos(x2))*(a1 - a3 + a2*cos(x2)))/(a2^2*cos(x2)^2 - a1*a3 + a3^2)^2;
    Fxu(4,3,1) = 0;
    Fxu(4,4,1) = 0;
end