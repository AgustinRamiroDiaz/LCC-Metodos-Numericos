// Pregunta 4

// Aplica el método de la potencia a la matriz A
// partiendo del vector inicial z
// con un máximo de iteraciones maxIter
// para obtener su autovalor de módulo máximo
// junto a su autovector asociado
//      lambda es el autovalor
//      z el autovector asociado a lambda
//      i la cantidad de iteraciones del método
// El método toma como regla de corte
// norm(z^(i+1) - z^i) < eps
function [lambda, z, i] = metodoDeLaPotencia(A, z, eps, maxIter)
    for i=1:maxIter
        w = A * z
        znuevo = w / norm(w, 'inf')
        if norm(znuevo - z) < eps   // Vemos si se cumple la regla de corte
            z = znuevo  // z representa nuestro z^n
            break
        end
        z = znuevo
    end
    [_, k] = max(abs(w)) // Encontramos posición de una componente no nula (en este caso la de mayor módulo)
    w = A*z        // conseguimos w^(n+1)
    lambda = w(k) / z(k)
endfunction


// a)
// Aplicamos el teorema de Gerschgorin para acotar los autovalores
// Sea r cualquier autovalor de A, sabemos que se cumple

// |r - A(i, i)| <= Σ|A(i, j)| - |A(i, i)|

// Acotemos por filas
// Sabemos que se cumple alguno de los siguientes casos

// |r - 0.4| <= 0.6
// |r - 0.7| <= 0.3
// |r - 0.6| <= 0.4

// Por lo tanto podemos acotar el modulo de r de la siguiente manera:
// r pertenece pertenece a un circulo centrado en c de radio z, 
// entonces estará acotado por |r| <= |c| + z

// Por lo tanto

// |r| <= 0.4 + 0.6 = 1
// o bien
// |r| <= 0.7 + 0.3 = 1
// o bien
// |r| <= 0.6 + 0.4 = 1

// De esto concluimos que todo r autovalor está acotado por 1 =>
// El autovalor maximo esta acotado por 1


// b)
// A = [0.4 0.2 0.4;
//      0.3 0.7 0  ;
//      0.3 0.1 0.6];

// [lambda, z, i] = metodoDeLaPotencia(A, [1 2 3]', 1e-12, 100)

// disp("Obtuve por el metodo de la potencia: ")
// disp("Autovalor de modulo maximo: " + string(lambda))
// disp("Autovector asociado:", z)
// disp("Cantidad de iteraciones: " + string(i))

// av = spec(A)
// [_, k] = max(abs(av))
// avmaxmod = av(k)

// disp("El autovalor de modulo maximo calculado por scilab es: " + string(avmaxmod))

// disp("Verificamos que nuestros resultados obtenidos por el metodo de la potencia: ")
// disp("A*v-lambda*v = ")
// disp(A*z - lambda * z)
// disp("Vemos que nos da practicamente 0 excepto por unos pequeños errores de calculo numerico")



// Resultado de Scilab
// "Obtuve por el metodo de la potencia: "

// "Autovalor de modulo maximo: 1"

// "Autovector asociado:"

//  1.
//  1.
//  1.

// "Cantidad de iteraciones: 52"

// "El autovalor de modulo maximo calculado por scilab es: 1"

// "Verificamos que nuestros resultados obtenidos por el metodo de la potencia: "

// "A*v-lambda*v = "

//  1.166D-13
//  5.821D-13
//  0.       

// "Vemos que nos da practicamente 0 excepto por unos pequeños errores de calculo numerico"