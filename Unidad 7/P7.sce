// TODAS LAS FUNCIONES TRABAJAN CON VECTORES FILA

// Calcula el polinomio L_k para la 
// interpolación de Lagrange
// a partir del vector x de preimágenes
function y = Lk(x,k)
    [m, n] = size(x)
    r = x
    r(k) = []
    p = poly(r,"x","roots")
    pk = prod(x(k) * ones(r) - r)
    y = p / pk
endfunction

// Calcula el polinomio de interpolación de Lagrange
// para los puntos x con imagen y
function p = interpolacionLagrange(x, y)
    [m,n] = size(x)
    p = 0
    for k = 1:n
        p = p + Lk(x, k) * y(k)
    end
endfunction

// Calcula las diferencias divididas de Newton
// para los puntos x con imagen y
function dd = diferenciasDivididas(x, y)
    [m, n] = size(x)
    if n == 1
        dd = y(1)
    else
        dd = (diferenciasDivididas(x(2:n), y(2:n)) - diferenciasDivididas(x(1:n-1), y(1:n-1))) / (x(n) - x(1))
    end
endfunction

// Versión sin encajar
// Calcula el polinomio de interpolación de Newton por diferencias divididas
// para los puntos x con imagen y
// function p = interpolacionNewton(x, y)
//     [m,n] = size(x)
//     // Polinomio resultante
//     p = diferenciasDivididas(x(1), y(1))
//     for k = 2:n
//         np = poly(x(1:k-1), 'x', 'r')
//         p = p + np * diferenciasDivididas(x(1:k), y(1:k))
//     end
// endfunction


// Versión con multiplicaciones encajadas
// Calcula el polinomio de interpolación de Newton por diferencias divididas
// para los puntos x con imagen y
function pn = interpolacionNewton(x, y)
    [m,n] = size(x)
    // Polinomio resultante
    pn = 0
    for k = n:-1:1
        dd = diferenciasDivididas(x(1:k), y(1:k))   // D_k
        p = poly(x(k), 'x', 'r')                    // (x-x_k)
        pn = dd + p * pn                            // Iteración
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
// x = [.2 .4]
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
// x = [0 .2 .4 .6];
// y = [1 1.2214 1.4918 1.8221];
// p = interpolacionLagrange(x, y)
// errorExacto = horner(p, 1/3)-1.395612425
// disp("Error cúbico: ", errorExacto)

// x = [.2 .4];
// y = [1.2214 1.4918];
// p = interpolacionNewton(x, y)
// errorExacto = horner(p, 1/3)-1.395612425
// disp("Error lineal: ", errorExacto)


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
// TODO en papel, son varias interpolaciones lineales


// Ejercicio 4
x = [2 2.1 2.2 2.3 2.4 2.5]
y = [.2239 .1666 .1104 .0555 .0025 -.0484]

function y = J0(x)
    deff('y=f(t)', 'y=cos(x*sin(t))')
    i = integrate('f(t)', 't', 0, %pi)
    y = i / %pi
endfunction

p = interpolacionNewton(x, y)
v1 = horner(p, 2.15)
v2 = horner(p, 2.35)

v1real = J0(2.15)
v2real = J0(2.35)

// Calculemos el error 
// (sabiendo por el ejercicio 3 que la derivada enésima está acotada por 0)
perror = poly(x, 'x', 'roots') / factorial(6)
v1cotaError = horner(perror, 2.15)
v2cotaError = horner(perror, 2.35)

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

// disp(horner(interpolacionLagrange([0 1 2 3], [y0 y1 y2 y3]), 2.5))
2.8750000



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
//       = 9

// c)
// f(x) - p3(x) = (x - x0)(x - x1)...(x-xn) / (n+1)! * f^(n+1) (ξ) 
// f(0) - p(0) = (0 - x0)(0 - x1)(0 - x2)(0 - x3) / 4! * f^(4) (ξ) 
// f(0) - p(0) = (-x0)(-x1)(-x2)(-x3) / 24 * f^(4) (ξ) 
// f(0) - p(0) = x0 * x1 * x2 * x3 / 24 * f^(4) (ξ) 
// f(0) - p(0) = -1 * 1 * 2 * 4 / 24 * f^(4) (ξ) 
// f(0) - p(0) = -1 / 3 * f^(4) (ξ) 
// |f(0) - p(0)| = 1 / 3 * |f^(4) (ξ)| <= 1 / 3 * 33.6 
// |f(0) - p(0)| <= 11.2 


function pol = minimosCuadrados(x, y, grado)
    [m, n] = size(x)
    // Definimos A tomando ɸ_k como x^k  
    A = zeros(n, grado + 1)
    for i = 1:n
        for j = 0:grado
            A(i, j+1) = x(i) ^ j
        end
    end
    // Calculamos los coeficientes utilizando el Teorema 2
    a = (A' * A) \ (A' * y')
    pol = poly(a, 'x', 'c')
endfunction


// Ejercicio 7
x = [0  .15     .31     .5      .6      .75];
y = [1  1.004   1.31    1.117   1.223   1.422];

p1 = minimosCuadrados(x, y, 1)
p2 = minimosCuadrados(x, y, 2)
p3 = minimosCuadrados(x, y, 3)

// disp(sum(abs((horner(p1, x) - y))))
// 0.4784433
// disp(sum(abs((horner(p2, x) - y))))
// 0.4665116
// disp(sum(abs((horner(p3, x) - y))))
// 0.3821412

// disp(norm(horner(p1, x) - y))
// 0.2315637
// disp(norm(horner(p2, x) - y))
// 0.2305248
// disp(norm(horner(p3, x) - y))
// 0.2012893

// La mejor aproximación de mínimos cuadrados es la de grado 3

// Ejercicio 8
x = [4      4.2     4.5     4.7     5.1     5.5     5.9     6.3     6.8     7.1]
y = [102.56 113.18  130.11  142.05  167.53  195.14  224.87  256.73  299.5   326.72]

// a)
p1 = minimosCuadrados(x, y, 1)
p2 = minimosCuadrados(x, y, 2)
p3 = minimosCuadrados(x, y, 3)

// Gráficos de funciones
// plot2d(x, y, style=color("red"))
// plot2d(x, horner(p1, x), style=color("green"))
// plot2d(x, horner(p2, x), style=color("blue"))
// plot2d(x, horner(p3, x), style=color("pink"))

// Errores
// disp(norm(horner(p1, x) - y))
// disp(norm(horner(p2, x) - y))
// disp(norm(horner(p3, x) - y))

// Podemos ver que el polinomio de grado 1 tiene un error grande
// mientras que los polinomios de grado 2 y 3 son casi idénticos a la función

// Gráficos de errores
// plot2d(x, abs(horner(p1, x) - y), style=color("green"))
// plot2d(x, abs(horner(p2, x) - y), style=color("blue"))
// plot2d(x, abs(horner(p3, x) - y), style=color("pink"))

// Ejercicio 9

deff('y = f(x)', 'y = 1 ./ (1 + x^2)') 

// xgrid()
function Ejercicio9(n, c)
    x = linspace(-5, 5, n)
    y = f(x)
    p = interpolacionLagrange(x, y)

    t = -5:.01:5
    plot2d(t, f(t) - horner(p, t), style = color(c))
endfunction

// Ejercicio9(2, 'red')
// Ejercicio9(4, 'green')
// Ejercicio9(6, 'blue')
// Ejercicio9(10, 'pink')
// Ejercicio9(14, 'magenta')


// Los errores tienden a aumentar a medida que 
// nos alejamos del 0

// Vemos que a medida que aumenta el n
// obtenemos errores menores cerca del 0
// pero mayores hacia los extremos -5 y 5

// Parece ser el fenómeno de Runge






// Calcula n nodos de Chebyshev
// para el intervalo [-1, 1]
// Retora el polinomio de Chebyshev p de grado n y sus raices r
function [p, r] = chebyshevConRaices(n)
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
[p, x] = chebyshevConRaices(4)
y = exp(x)

p = interpolacionLagrange(x, y)

// b)
// t=[-1:.1:1]
// plot(t, exp(t) - horner(p, t))

// Ejercicio 11

// Calcula n nodos de Chebyshev
// para el intervalo [a, b]
// a partir de un cambio lineal de variable 
// Retora el polinomio de Chebyshev p de grado n y sus raices r
function nc = nodosChebyshev(n, a, b)
    [p, r] = chebyshevConRaices(n)
    nc = ((b + a) + r * (b - a)) / 2
endfunction

a = 0
b = %pi/4

n = 4

x = nodosChebyshev(n, a, b)
y = cos(x)

p = interpolacionLagrange(x, y)

t=[a:.1:b]
// plot2d(t, horner(p, t) - cos(t), style=color("red"))
