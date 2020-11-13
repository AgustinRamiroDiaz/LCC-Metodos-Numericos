// Aplica el metodo de la bisección con la función f
// para reducir el intervalo [a, b]
// a la mitad p veces
function [a, b] = metodo_biseccion(f, a, b, p)
    for i=1:p 
        c = (a + b) / 2
        if (f(b) * f(c)) <= 0 then
            a = c
        else
            b = c
        end
    end
endfunction


// Aplica el método de newton con la función f
// partiendo desde x con una tolerancia de epsilon 
// y un cantidad máxima de iteraciones iter
function x = metodo_newton(f, x, eps, iter)
    for i = 1: iter
        d = numderivative(f, x)
        xNuevo = x - f(x) / d
        if abs(xNuevo - x) < eps then
            x = xNuevo; return
        end
        x = xNuevo
    end
    x = %nan    
endfunction



// Item a)

// Paso 1:
// Partiendo del intervalo [a, b], lo reduce p veces a la mitad
// con el método de la bisección
// Paso 2:
// Aplica el método de Newton dentro del intervalo obtenido anteriormente
// hasta que se cumpla una tolerancia de eps
// Si el método de newton se llegase a salir del intervalo, se vuleve al paso 1
function x = metodo_mixto(f, a, b, eps, p)
    while 1
        [a, b] = metodo_biseccion(f, a, b, p)
        x = (a + b) / 2
        while a < x & x < b
            d = numderivative(f, x)
            xNuevo = x - f(x) / d
            if abs(xNuevo - x) < eps then
                x = xNuevo; return
            end
            x = xNuevo
        end
    end
endfunction


// Item b)

resultado = metodo_mixto(tanh, -9, 14, 1e-30, 1)
disp(resultado)
// 4.139D-49