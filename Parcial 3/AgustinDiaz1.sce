// b)
b = [1 0 0; 
    -.1 0 .1; 
    -.1 -.1 2]

disp('Autovalores de b', spec(b))

// d)
d = [4 -1 0; 
    -1 4 -1; 
    -1 -1 4]

disp('Autovalores de d', spec(d))

// f)
f = [4.75 2.25 -.25; 2.25 4.75 1.25; -.25 1.25 4.75]

disp('Autovalores de f', spec(f))


// Resultado de Scilab

// "Autovalores de b"

// 1.9949874 + 0.i
// 0.0050126 + 0.i
// 1.        + 0.i

// "Autovalores de d"

// 4.618034 + 0.i
// 2.381966 + 0.i
// 5.       + 0.i

// "Autovalores de f"

// 2.0646374
// 4.9616991
// 7.2236635





// Si desea graficarlo en Scilab, puede descomentar las ultimas lineas 
// (por partes para mantener coherencia)

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

// b
// CircGersValor(b)
// CircGersValor(b')
// d
// CircGersValor(d)
// CircGersValor(d')
// f
// CircGersValor(f)
// CircGersValor(f')