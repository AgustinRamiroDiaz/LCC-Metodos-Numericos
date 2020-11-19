// Ejercicio 1

// a)

A = [1 0 0; -1 0 1; -1 -1 2];

|λ1 - 1| <= 0
|λ2 + 1| <= 0
|λ3 - 2| <= 2

// b)
A = [1 0 0; -.1 0 .1; -.1 -.1 2]
|λ1 - 1| <= 0
|λ2| <= .2
|λ3 - 2| <= .2

// c)
A = [1 0 0; -.25 0 .25; -.25 -.25 2]
|λ1 - 1| <= 0
|λ2| <= .5
|λ3 - 2| <= .5





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