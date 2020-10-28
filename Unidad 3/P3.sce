


// 3

//Recibe una funcion y dos puntos de un intervalo
//Convergencia mas rapida que lineal
//Convergencia no asegurada
//Tratar de que F'(a) != 0
function output = secante(f, anterior, actual, epsilon)
    while abs (f (actual) - f (anterior)) >= epsilon
        siguiente = actual - f (actual) * (actual - anterior) / (f (actual) - f (anterior))
        anterior = actual        
        actual = siguiente 
    end
    output = actual
endfunction

// 1.9337538
function output = tres(x)
// Description of tres(x)
    output = x * x / 4 - sin(x)
endfunction


// 4
// Usamos el Teorema 2 con g = cos, y obtenemos que existe una única solución de g(x) = x.
// Por lo tanto, al aplicar reiteradas veces la función coseno obtenemos su punto fijo


// 5
// Usamos el Teorema 2
// g (x) = 2 ^ (x - 1)
// Buscamos que g' (x) < 1 => x < 1.53
// Definimos a = -inf, b = 1.53
// Corroboramos que si -inf <= x <= 1.53 => -inf <= g (x) <= 1.53
// Por lo tanto, converge si partimos de un x < 1.53
// Como la solución es única y sabemos que 1 es solución, entonces 1 es la única solución
// Por lo tanto converge a 1



// 6


// 7

function output = puntoFijo(f, x, epsilon)
    // Description of puntoFijo(f, x, epsilon)
    while abs (f (x) - x) > epsilon
        x = f (x)
    end
    output = x
endfunction


function output = longitudDeOnda(d)
    output = 4 * %pi ^ 2 / (5 ^ 2 * 9.8 * tanh (4 * d)) 
endfunction

// a    0.2835513

// b    

// function output = metodoNewton(f, x, epsilon)
//     output = puntoFijo(deff('y = newton(x)','y = x - f (x) / numderivative(f, x)'), x, epsilon)
// endfunction

function output = metodoNewton(f, x, epsilon)
    while abs (f (x) - x) > epsilon
        x = x - f (x) / numderivative(f, x)
    end
    output = x
endfunction

