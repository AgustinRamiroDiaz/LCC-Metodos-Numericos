// Ejercicio 1
// SI LO PIDEN HACERLO A MANO
// SI NO ES SIMÉTRICA, DIBUJAR CIRCULOS
// SI ES SIMÉTRICA DAR INTERVALOS REALES

// a)
// |λ - A(i, i)| <= ΣA(i, j) - A(i, i)

// A = [1 0 0; 
//     -1 0 1; 
//     -1 -1 2];


// Filas
// λ = 1     
// |λ| <= 2     
// |λ - 2| <= 2

// Columnas
// |λ - 1| <= 2     
// |λ| <= 1    
// |λ - 2| <= 1


// // b)
// A = [1 0 0; 
//     -.1 0 .1; 
//     -.1 -.1 2]



// Filas
// λ = 1
// |λ| <= .2
// |λ - 2| <= .2

// Columnas
// |λ - 1| <= .2
// |λ| <= .1
// |λ - 2| <= .1


// // c)
// A = [1 0 0; 
//     -.25 0 .25; 
//     -.25 -.25 2]

// Filas
// λ = 1
// |λ| <= .5
// |λ - 2| <= .5

// Columnas
// |λ - 1| <= .5
// |λ| <= .25
// |λ - 2| <= .25

// // d)
// A = [4 -1 0; 
//     -1 4 -1; 
//     -1 -1 4]

// Filas
// |λ - 4| <= 1
// |λ - 4| <= 2
// |λ - 4| <= 2

// Columnas
// |λ - 4| <= 2
// |λ - 4| <= 2
// |λ - 4| <= 1

// Es decir
// |λ - 4| <= 2

// // e)
// A = [3 2 1; 
//     2 3 0; 
//     1 0 3]

// Filas
// |λ - 3| <= 3
// |λ - 3| <= 2
// |λ - 3| <= 1

// Columnas ídem a filas por ser simétrica

// Es decir
// |λ - 3| <= 3

// // f)
// A = [4.75 2.25 -.25; 2.25 4.75 1.25; −.25 1.25 4.75]

// Filas
// |λ - 4.75| <= 2.5
// |λ - 4.75| <= 3.5
// |λ - 4.75| <= 1.5

// Columnas ídem a filas por ser simétrica

// Es decir
// |λ - 4.75| <= 3.5


// Ejercicio 2
// TODO


// Ejercicio 3
// TODO REVISAR ACOTAR EL ERROR CON LOS RESULTADOS DEL EJERCICIO 2
// i)
// ii)
k = 0:10
eps = .1*k

raices = zeros(11, 3)
av = zeros(11, 3)

// for i=k+1
//     disp('Epsilon: ', eps(i))
//     A = [1 -1 0; -2 4 -2; 0 -1 1] + [0 0 0; 0 0 0; 0 0 eps(i)];
//     // disp('Matriz A: ', A)
//     p = poly(A, 'x')
//     disp('Polinomio Característico: ', p)
//     raices(i, :) = gsort(roots(p))
//     disp('Raices: ', raices(i, :))
//     av(i, :) = gsort(spec(A))
//     disp('Autovalores de A: ', av(i, :))
// end


// Podría ver de aplicar a partir del polinomio caracteístico
// p(λ) = λ^n + a_n-1 λ^(n-1) + ... + a_1 λ + a_0
// para calcular las raícs r
// |r| <= 1   o   |r + a_(n-1)| <= |a_0| + ... + |a_(n-2)|
// y
// |r| <= a_0   o   |r| <= 1 + |a_i| para i=1...n-2   o   |r + a_(n-1)| <= 1

// Ejercicio 4

// a)

// Dado un radio r y dos coordenadas x e y
// dibuja un circulo en la pantalla 
// centrado en la coordenada (x, y) de radio r
function circ(r, x, y)
    xarc(x-r, y+r, 2*r, 2*r, 0, 360*64)
endfunction

// b)
// Dada una matriz A
// dibuja todos sus círculos de Gerschgorin
function Gers(A)
    [m, n] = size(A)

    t = linspace(0, 2*%pi, 100);

    centros = diag(A)
    radios = sum(abs(A), 'c') - abs(centros)

    minx = round(min(centros-radios)-1)
    miny = round(min(-radios)-1)

    maxx = round(max(centros+radios)+1)
    maxy = round(max(radios)+1)

    rectangulo = [minx, miny, maxx, maxy]

    plot2d([], [], rect = rectangulo)
    xgrid(2);

    for i=1:m
        circ(radios(i), centros(i), 0)
    end
endfunction

// c)
// Dada una matriz A
// dibuja todos sus círculos de Gerschgorin
// y marca sus autovalores
function CircGersValor(A)
    Gers(A)
    autov = spec(A)
    plot2d(real(autov), imag(autov), -1)
endfunction

A = [1 0 0; -1 0 1; -1 -1 2];

// CircGersValor(A)

// Ejercicio 5
// a)

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

A=[-7 13 -16; 13 -10 13; -16 13 -7];

A1 = [6 4 4 1;
      4 6 1 4;
      4 1 6 4;
      1 4 4 6]

[lambda1, z1, i1] = metodoDeLaPotencia(A1, [1 2 3 4]', 1e-10, 100)

A2 = [12 1 3 4;
      1 -3 1 5;
      3 1 6 -2;
      4 5 -2 -1]

[lambda2, z2, i2] = metodoDeLaPotencia(A2, [1 0 0 0]', 1e-10, 100)

// Compara la diferencia entre el autovalor aproximado por el 
// método de la potencia y el mayor autovalor, 
// considerando el número de iteraciones realizadas
function [diff, i] = comparar(A, z, eps, maxIter)
    autoValores = spec(A)
    [_, k] = max(abs(autoValores))
    autoValor_max = autoValores(k)
    
    // lambda es el vector de todos los valores de las iteraciones
    // del autovalor que queremos calcular
    for i=1:maxIter
        [lambda(i), z, _] = metodoDeLaPotencia(A, z, 0, 1)
        if abs(lambda(i) - autoValor_max) <= eps
            break
        end
    end

    diff = lambda - autoValor_max
endfunction

[diff1, it1] = comparar(A1, rand(4, 1), 1e-10, 100)

[diff2, it2] = comparar(A2, rand(4, 1), 1e-10, 100)

// plot(diff1)
// plot(diff2)



function [lambda, zs, autoValor_max, i] = f(A, z, eps, maxIter)
    autoValores = spec(A)
    [_, k] = max(abs(autoValores))
    autoValor_max = autoValores(k)
    
    //zs = zeros(2, 1)
    // lambda es el vector de todos los valores de las iteraciones
    // del autovalor que queremos calcular
    for i=1:maxIter
        [lambda(i), zs(i, :), _] = metodoDeLaPotencia(A, z, 0, 1)
        z = zs(i, :)'
        if abs(lambda(i) - autoValor_max) <= eps
            break
        end
    end

endfunction