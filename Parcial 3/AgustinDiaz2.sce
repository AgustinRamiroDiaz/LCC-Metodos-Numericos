// Practica 7
// Ejercicio 6
// Hice la resolucion con Scilab y a mano

// a)
function p = resolver_ejercicio6()
    x = [-1, 1, 2, 4]
    dd = [2, 1, -2, 2]
  
    p = dd(1)
    p_n = 1
    for i = 2 : 4
      p_n = poly(x(1:i-1), 'x', 'r')
      p = p + p_n * dd(i)
    end
endfunction

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


// b)
// p3(0) = 3 + -2 * 1 * -1 + 2 * 1 * -1 * -2
//       = 3 + 2 + 2 * 2
//       = 9

p3 = resolver_ejercicio6()
disp('p3(x):', p3)
disp('b:  p3(0) = ' + string(horner(p3, 0)))

// Resultado de Scilab
// "p3(x):"

// 9 -x -6x² +2x³

// "b:  p3(0) = 9"




// c)
// f(x) - p3(x) = (x - x0)(x - x1)...(x-xn) / (n+1)! * f^(n+1) (ξ) 
// f(0) - p(0) = (0 - x0)(0 - x1)(0 - x2)(0 - x3) / 4! * f^(4) (ξ) 
// f(0) - p(0) = (-x0)(-x1)(-x2)(-x3) / 24 * f^(4) (ξ) 
// f(0) - p(0) = x0 * x1 * x2 * x3 / 24 * f^(4) (ξ) 
// f(0) - p(0) = -1 * 1 * 2 * 4 / 24 * f^(4) (ξ) 
// f(0) - p(0) = -1 / 3 * f^(4) (ξ) 
// |f(0) - p(0)| = 1 / 3 * |f^(4) (ξ)| <= 1 / 3 * 33.6 
// |f(0) - p(0)| <= 11.2 