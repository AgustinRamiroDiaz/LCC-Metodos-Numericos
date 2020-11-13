// Dado un sitema de ecuaciones F y un punto inicial x,
// obtiene una solución aproximada al sistema F(v) = 0 con el método de newton multivariable
// con una tolerancia de epsilon y un máximo iter de iteraciones
function v = metodo_newton_multivariable(F, x, eps, iter)
    for i = 1: iter
        Jinv = inv(numderivative(F, x))
        x = x - (Jinv * F(x))
        if norm(F(x)) < eps then
            v = x; return
        end
    end
    v = %nan
endfunction


function y = F10(x)
  y(1) = x(1)^2 + (x(1) * x(2)^3) - 9
  y(2) = (3 * x(1)^2 * x(2)) - 4 - x(2)^3
endfunction

// Item a)
x = metodo_newton_multivariable(F10, [1.2, 2.5]', 1e-8, 100)
disp("Item a", x)
disp("Corroboramos que de aproximadamente 0", F10(x))


// Item b)
x = metodo_newton_multivariable(F10, [-2, 2.5]', 1e-8, 100)
disp("Item b", x)
disp("Corroboramos que de aproximadamente 0", F10(x))

// Item c)
x = metodo_newton_multivariable(F10, [-1.2, -2.5]', 1e-8, 100)
disp("Item c", x)
disp("Corroboramos que de aproximadamente 0", F10(x))


// Item d)
x = metodo_newton_multivariable(F10, [2, -2.5]', 1e-8, 100)
disp("Item d", x)
disp("Corroboramos que de aproximadamente 0", F10(x))
  

// Resultado de Scilab:
// "Item a"

// 1.3363554
// 1.7542352

// "Corroboramos que de aproximadamente 0"

// -1.776D-15
// -1.776D-15

// "Item b"

// -0.9012662
// -2.0865876

// "Corroboramos que de aproximadamente 0"

// 7.726D-10
// 1.821D-10

// "Item c"

// -0.9012662
// -2.0865876

// "Corroboramos que de aproximadamente 0"

// 1.776D-15
// 3.553D-15

// "Item d"

// -3.0016249
// 0.148108 

// "Corroboramos que de aproximadamente 0"

// -4.086D-14
// 5.381D-13
