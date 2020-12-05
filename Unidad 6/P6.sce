// Ejercicio 1

// a)
// TODO POR COLUMNAS
// |λ - A(i, i)| <= ΣA(i, j) - A(i, i)

// A = 
//     [1 0 0; 
//     -1 0 1; 
//     -1 -1 2];

// λ = 1     
// o
// |λ| <= 2     
// o
// |λ - 2| <= 2

// y

// |λ - 1| <= 2     
// o
// |λ| <= 1    
// o
// |λ - 2| <= 1


// // b)
// A = [1 0 0; -.1 0 .1; -.1 -.1 2]
// |λ1 - 1| <= 0
// |λ2| <= .2
// |λ3 - 2| <= .2

// // c)
// A = [1 0 0; -.25 0 .25; -.25 -.25 2]
// |λ1 - 1| <= 0
// |λ2| <= .5
// |λ3 - 2| <= .5

// // d)
// A = [4 -1 0; -1 4 -1; -1 -1 4]
// |λ1 - 4| <= 1
// |λ2 - 4| <= 2
// |λ3 - 4| <= 2

// // e)
// A = [3 2 1; 2 3 0; 1 0 3]
// |λ1 - 3| <= 3
// |λ2 - 3| <= 2
// |λ3 - 3| <= 1

// // f)
// A = [4.75 2.25 -.25; 2.25 4.75 1.25; −.25 1.25 4.75]
// |λ1 - 4.75| <= 2.5
// |λ2 - 4.75| <= 3.5
// |λ3 - 4.75| <= 5

// Ejercicio 2
// TODO

// Ejercicio 3
// i)
// ii)
k = 0:10
eps = .1*k

raices = zeros(11, 3)
av = zeros(11, 3)

for i=k+1
    disp('Epsilon: ', eps(i))
    A = [1 -1 0; -2 4 -2; 0 -1 1] + [0 0 0; 0 0 0; 0 0 eps(i)];
    // disp('Matriz A: ', A)
    p = poly(A, 'x')
    disp('Polinomio Característico: ', p)
    raices(i, :) = gsort(roots(p))
    disp('Raices: ', raices(i, :))
    av(i, :) = gsort(spec(A))
    disp('Autovalores de A: ', av(i, :))
end

// Ejercicio 4

// a)
function circ(r, x, y)
    xarc(x-r, y+r, 2*r, 2*r, 0, 360*64)
endfunction

// b)
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
function CircGersValor(A)
    Gers(A)
    autov = spec(A)
    plot2d(real(autov), imag(autov), -1)
endfunction

// Ejercicio 5
// a)

// TODO CORREGIR Y TERMINAR
function [lambda, z, i] = metodoDeLaPotencia(A, z, eps, maxIter)
    for i=1:maxIter
        w = A * z
        znuevo = w / max(abs(w))
        if norm(znuevo - z) < eps
            k = find(w, 1)
            w = A*z
            lambda = w(k) / z(k)
            return
        end
        z = znuevo
    end
    k = find(w, 1)
    w = A*z
    lambda = w(k) / z(k)
endfunction

A1 = [6 4 4 1;
      4 6 1 4;
      4 1 6 4;
      1 4 4 6]

metodoDeLaPotencia(A1, rand(4, 1), 1e-12, 100)

A2 = [12 1 3 4;
      1 -2 1 5;
      3 1 6 -2;
      4 5 -2 -1]

metodoDeLaPotencia(A2, rand(4, 1), 1e-12, 100)
