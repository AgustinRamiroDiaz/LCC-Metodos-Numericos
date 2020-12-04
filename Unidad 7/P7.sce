function y = Lk(x,k)
    [m, n] = size(x)
    r = x
    r(k) = []
    p = poly(r,"x","roots")
    pk = prod(x(k) * ones(r) - r)
    y = p / pk
endfunction

// TODO: SE PUEDE OPTIMIZAR CALCULANDO LK ADENTRO 
function p = interpolacionLagrange(x, y)
    [m,n] = size(x)
    p = 0
    for k = 1:n
        p = p + Lk(x, k) * y(k)
    end
endfunction


function dd = diferenciasDivididas(x, y)
    [m, n] = size(x)
    if n == 1
        dd = y(1)
    else
        dd = (diferenciasDivididas(x(2:n), y(2:n)) - diferenciasDivididas(x(1:n-1), y(1:n-1))) / (x(n) - x(1))
    end
endfunction


function p = interpolacionNewton(x, y)
    [m,n] = size(x)
    // Polinomio n
    np = 1
    p = diferenciasDivididas(x(1), y(1))
    for k = 2:n
        np = poly(x(1:k-1), 'x', 'r')
        p = p + np * diferenciasDivididas(x(1:k), y(1:k))
    end
endfunction





// Ejercicio 1
// a)

// Cúbico:
// x = [0 .2 .4 .6];
// y = [1 1.2214 1.4918 1.8221];
// disp("Cúbicos")
// // Lagrange:
// p = interpolacionLagrange(x, y)
// disp("Lagrange", horner(p, 1/3))
// 1.3955494


// // Newton:
// p = interpolacionNewton(x, y)
// disp("Newton", horner(p, 1/3))
// 1.3955494

// perror = poly(x, 'x', 'roots') / factorial(4)
// // sabiendo que la derivada de e^x es e^x
// // y que e^x es creciente
// cotaErrorCubico = horner(perror, 1/3) * exp(.6)
// disp('Cota error cúbico: ' + string(cotaErrorCubico))

// // Lineal:
x = [.2 .4]
// y = [1.2214 1.4918]
// disp("Lineal")
// // Lagrange:
// p = interpolacionLagrange(x, y)
// disp("Lagrange", horner(p, 1/3))
// 1.4016667

// // Newton:
// p = interpolacionNewton(x, y)
// disp("Newton", horner(p, 1/3))
// 1.4016667

// perror = poly(x, 'x', 'roots') / 2
// // sabiendo que la derivada de e^x es e^x
// // y que e^x es creciente
// cotaErrorLineal = horner(perror, 1/3) * exp(.4)
// disp('Cota error lineal: ' + string(cotaErrorLineal))


// // b)
x = [0 .2 .4 .6];
y = [1 1.2214 1.4918 1.8221];
p = interpolacionLagrange(x, y)
errorExacto = horner(p, 1/3)-1.395612425
disp("Error cúbico: ", errorExacto)

// x = [.2 .4];
// y = [1.2214 1.4918];
// p = interpolacionNewton(x, y)
// errorExacto = horner(p, 1/3)-1.395612425
// disp("Error lineal: ", errorExacto)
// // TODO COTAS DE ERROR


// Ejercicio 2
// Por el Teorema 3 (Error de la Interpolación Polinómica)
// sabemos que para x0, x1,..., xn distintos en [a, b],
// para todo x en [a, b] existe ξ en (a, b) tal que

// f(x) - p(x) = (x - x0)(x - x1)...(x-xn) / (n+1)! * f^(n+1) (ξ) 

// donde p(x) es el polinomio interpolante de f para los puntos xi

// Sabiendo que f es un polinomio de orden menor o igual que n, 
// la derivada f^(n+1) (x) = 0 para todo x en [a, b]

// Por lo tanto f(x) = p(x) para todo x en [a, b]

// Dado que este razonamiento vale para cualquier conjunto [a, b]
// que contenga a los xi, vale para todo x en R

// Por lo tanto f(x) = p(x) para todo x

// Ejercicio 3
// Por el Teorema 3 (Error de la Interpolación Polinómica)
// sabemos que para x0, x1 distintos en [0, 1],
// existe ξ en (0, 1) tal que
// f(x) - p(x) = (x - x0)(x - x1) / 2! * f^(2) (ξ) 
// = (x - x0)(x - x1) / 2 * f^(2) (ξ)
//     como la derivada segunda de J0 es creciente en [0, 1]
//     J0'(x)= 1 / %pi * 
// <= (x - x0)(x - x1) / 2 * f^(2) (1)
// = (x - x0)(x - x1) / 2 * -13/40
// = (x - x0)(x - x1) * -13 / 80
//     y queremos
// (x - x0)(x - x1) * -13 / 80 < .5e-6
// (x - x0)(x - x1) * > 3.08e-6

// WTF TODO


// Ejercicio 4
// TODO PARECIDO ANTERIOR

// Ejercicio 6
// p(x) = 
//     f[x0] +                      
//     f[x0, x1] * (x - x0) + 
//     f[x0, x1, x2] * (x-x0) * (x-x1) +...+
//     f[x0,...,xn-1] * (x-x0) *...* (x-xn)

// x0 = -1
// x1 = 1
// x2 = 2
// x3 = 4

// p3(x) = 
//     f[x0] + 
//     f[x0, x1] * (x - x0) + 
//     f[x0, x1, x2] * (x-x0) * (x-x1) +
//     f[x0, x1, x2, x3] * (x-x0) * (x-x1) * (x-x2)
// =
//     2 + 
//     1 * (x + 1) +
//     -2 * (x+1) * (x-1) +
//     2 * (x+1) * (x-1) * (x-2)
// =
//     3 + x + -2 * (x+1) * (x-1) + 2 * (x+1) * (x-1) * (x-2)

// p3(0) = 3 + -2 * 1 * -1 + 2 * 1 * -1 * -2
//       = 3 + 2 + 2 * 2
//       = 11

// TODO C


function pol = minimosCuadrados(x, y, grado)
    [m, n] = size(x)
    // Definimos A tomando ɸ_k como x^k  
    A = zeros(n, grado + 1)
    for i = 1:n
        for j = 0:grado
            A(i, j+1) = x(i) ^ j
        end
    end
    a = inv(A' * A) * A' * y'
    pol = poly(a, 'x', 'c')
endfunction

// Ejercicio 5
// P_0,1 (x) = 2x + 1
// P_0,2 (x) = x + 1
// P_1,2,3 (2.5) = 3

// P_0,1 (x) = L_0(x)y_0 + L_1(x)y_1
// P_0,1 (x) = (x-1)/(0-1) * y_0 + L_1(x)y_1
// P_0,1 (x) = (x-1)/(-1) * y_0 + L_1(x)y_1
// P_0,1 (x) = -(x-1) * y_0 + (x-0)/(1-0) * y_1
// P_0,1 (x) = (1-x) * y_0 + x * y_1
// P_0,1 (x) =  y0 - x * y0 + x * y_1
// P_0,1 (x) =  y0 + x * (-y0 + y_1)    *1*

// P_0,1 (x) = 2x + 1                  *2*

// de *1* y *2* concluimos
y0 = 1
y1 = 3

// P_0,2 (x) = L_0(x)y_0 + L_2(x)y_2
// P_0,2 (x) = (x-2) / (0-2) * y_0 + x/2 * y_2
// P_0,2 (x) = (x-2) / (-2) * y_0 + x/2 * y_2
// P_0,2 (x) = (1 - x/2) * y_0 + x/2 * y_2
// P_0,2 (x) = y_0 - x/2 * y0 + x/2 * y_2
// P_0,2 (x) = y_0 + x * (-y0 + y2) / 2

// P_0,2 (x) = x + 1

// Por lo tanto
y2 = 3

// P_1,2,3 (2.5) = L1(2.5) * y1 + L2(2.5) * y2 + L3(2.5) * y3
// P_1,2,3 (2.5) = L1(2.5) * 3 + L2(2.5) * 3 + L3(2.5) * y3
// P_1,2,3 (2.5) = (2.5 - 2) * (2.5 - 3) / (1 - 2) / (1 - 3) * 3 + L2(2.5) * 3 + L3(2.5) * y3
// P_1,2,3 (2.5) = 
//     (2.5 - 2) * (2.5 - 3) / (1 - 2) / (1 - 3) * 3 + 
//     (2.5 - 1) * (2.5 - 3) / (2 - 1) / (2 - 3) * 3 + 
//     (2.5 - 1) * (2.5 - 2) / (3 - 1) / (3 - 2) * y3
// P_1,2,3 (2.5) = -3/8 + 9/4 + 3/8 * y3
// 3 = 1.875 + 3/8 * y3
// 1.125= 3/8 * y3
y3 = 3

horner(interpolacionLagrange([0, 1, 2, 3], [y0, y1, y2, y3]), 2.5)
2.8750000


// Ejercicio 7
x = [0  .15     .31     .5      .6      .75];
y = [1  1.004   1.31    1.117   1.223   1.422];

p1 = minimosCuadrados(x, y, 1)
p2 = minimosCuadrados(x, y, 2)
p3 = minimosCuadrados(x, y, 3)

// disp(norm(horner(p1, x) - y))
//    0.2315637

// disp(norm(horner(p2, x) - y))
//    0.2305248

// disp(norm(horner(p3, x) - y))
//    0.2012893


// Ejercicio 8
x = [4      4.2     4.5     4.7     5.1     5.5     5.9     6.3     6.8     7.1]
y = [102.56 113.18  130.11  142.05  167.53  195.14  224.87  256.73  299.5   326.72]



p1 = minimosCuadrados(x, y, 1)
p2 = minimosCuadrados(x, y, 2)
p3 = minimosCuadrados(x, y, 3)

// TODO REVISAR DAN MUY MAL
// disp(norm(horner(p1, x) - y))
// plot(x, y)
// plot(x, horner(p2, x))

// disp(norm(horner(p2, x) - y))

// disp(norm(horner(p3, x) - y))



// Ejercicio 9

deff('y = f(x)', 'y = 1 / 1 + x^2') 


// TODO ES CON POLINOMIO DE INTERPOLACION
//xgrid()
// for n=2:2:14
//     x = linspace(-5, 5, n)'
//     y = f(x)
//     p = interpolacionLagrange(x, y)

//     t = -5:.01:5
//     plot(t, abs (f(t) - horner(p, t)))
// end

// Los errores tienden a aumentar a medida que 
// nos alejamos del 0

// Vemos que a medida que aumenta el n
// obtenemos errores mayores

// Parece ser el fenómeno de Runge


function [p, r] = chebyshev(n)
    // Casos base
    c(1) = poly(1, 'x', 'c')
    c(2) = poly([0 1], 'x', 'c')
    // Caso inductivo
    for i=3:n+1
        c(i) = poly([0 2], 'x', 'c') * c(i-1) - c(i-2)
    end
    p = c(n+1)
    r = roots(p)'
endfunction

// Ejercicio 10
// a)
[p, x] = chebyshev(4)
y = exp(x)

p = interpolacionLagrange(x, y)

// b)
// t=[-1:.1:1]
// plot(t, horner(p, t)-exp(t))

// Ejercicio 11
function nc = nodosChebyshev(n, a, b)
    [p, r] = chebyshev(n)
    nc = ((b + a) + r * (b - a)) / 2
endfunction

a = 0
b = %pi/4

x = nodosChebyshev(4, a, b)
y = cos(x)

p = interpolacionLagrange(x, y)

// t=[a:.1:b]
// plot(t, horner(p, t) - cos(t))
