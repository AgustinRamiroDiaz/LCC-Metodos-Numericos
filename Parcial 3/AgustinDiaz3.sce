// Practica 8
// Ejercicio 2

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

// Aplica el metodo compuesto del trapecio
// para aproximar la integral definida
// de f entre a y b
// partida en n subintervalos
function integral = metodoCompuestoTrapecio(f, a, b, n)
    h = (b - a) / n
    x = a+h:h:b-h
    integral = f(a) / 2 + sum(matrixfun(f, x)) + f(b) / 2
    integral = h * integral
endfunction

function Ejercicio2(f, a, b, n)
    disp('a=' + string(a))
    disp('b=' + string(b))
    disp('n=' + string(n))
    
    t = metodoCompuestoTrapecio(f, a, b, n)
    disp('Metodo Compuesto del Trapecio: ' + string(t))
    
    i = integrate('f(x)', 'x', a, b, n)
    disp('Integral Real: ' + string(i))
    
    disp('Error absoluto: ' + string(abs(i-t)))
    disp('Error relativo: ' + string(abs(1-t/i)))
endfunction

// a)
deff('y=a(x)', 'y=1/x')
disp('f(x)=1/x')
Ejercicio2(a, 1, 3, 4)

// d)
deff('y=d(x)', 'y=sin(%pi*x)')
disp('f(x)=sin(%pi*x)')
Ejercicio2(d, 0, 1, 8)

// e)
deff('y=e(x)', 'y = x * sin(x)')
disp('f(x) =  x * sin(x)')
Ejercicio2(e, 0, 2*%pi, 8)

// Resultado Scilab

// "f(x)=1/x"

// "a=1"

// "b=3"

// "n=4"

// "Metodo Compuesto del Trapecio: 1.1166667"

// "Integral Real: 1.0986123"

// "Error absoluto: 0.0180544"

// "Error relativo: 0.0164338"

// "f(x)=sin(%pi*x)"

// "a=0"

// "b=1"

// "n=8"

// "Metodo Compuesto del Trapecio: 0.6284174"

// "Integral Real: 0.6366198"

// "Error absoluto: 0.0082023"

// "Error relativo: 0.0128842"

// "f(x) =  x * sin(x)"

// "a=0"

// "b=6.2831853"

// "n=8"

// "Metodo Compuesto del Trapecio: -5.9568332"

// "Integral Real: -6.2831853"

// "Error absoluto: 0.3263521"

// "Error relativo: 0.0519406"
