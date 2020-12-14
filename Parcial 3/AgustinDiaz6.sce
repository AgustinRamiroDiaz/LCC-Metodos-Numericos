// Pregunta 6


// Aplica f a todos los elementos
// de la matriz A
function B = matrixfun(f,A)
    B = A
    for i = 1 : size(A,'r')
        for j = 1 : size(A,'c')
            B(i,j) = f(A(i,j))
        end
    end
endfunction


// Aplica el metodo compuesto de Simpson
// para aproximar la integral definida
// de f entre a y b
// partida en n subintervalos
// siendo n par
function integral = metodoCompuestoSimpson(f, a, b, n)
    h = (b - a) / n
    ximpares = a+h:2*h:b-h
    xpares = a+2*h:2*h:b-2*h

    integral = f(a) + 4 * sum(matrixfun(f, ximpares)) + 2 * sum(matrixfun(f, xpares)) + f(b)
    integral = h / 3 * integral
endfunction


deff('y=f(x)', 'y=cos(2*acos(x))')

a=-1
b=1

i2 = metodoCompuestoSimpson(f, a, b, 2)
i4 = metodoCompuestoSimpson(f, a, b, 4)
i6 = metodoCompuestoSimpson(f, a, b, 6)

disp('Los valores con el metodo compuesto de simpson son:')
disp('Con 2 intervalos: ' + string(i2))
disp('Con 4 intervalos: ' + string(i4))
disp('Con 6 intervalos: ' + string(i6))
disp('La integral calculada por Scilab es: ' + string(integrate('f(x)', 'x', -1, 1))) 

t = -1:.1:1

deff('y=g(x)', 'y=2*x^2-1')
plot2d(t, g(t), style=color('blue'))

plot2d(t, f(t), style=color('magenta'))



// Podemos ver que la función coincide con un polinomio de segundo grado en el intervalo [-1, 1]

// Debido a que el metodo compuesto de Simpson aproxima mediante
// interpolacion de polinomios de segundo grado, 
// (y dado que el polinomio de interpolacion de segundo grado coincide con el polinomio de segundo grado que se interpola)
// la aproximacion es siempre exacta (obviando errores numericos)
// sin importar cuantos intervalos utilizemos (mientras sean par y mayor que 1)
// Por lo cual no cambian los valores de la integral aproximada a medida que aumentamos 
// la cantidad de subintervalos


// Resultado de Scilab

// "Los valores con el metodo compuesto de simpson son:"

// "Con 2 intervalos: -0.6666667"

// "Con 4 intervalos: -0.6666667"

// "Con 6 intervalos: -0.6666667"

// "La integral calculada por Scilab es: -0.6666667"

// Y un plot de la funcion donde coinciden las funciones f y g
