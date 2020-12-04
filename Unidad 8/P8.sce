
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


// Aplica la regla del trapecio
// para aproximar la integral definida
// de f entre a y b
function integral = reglaTrapecio(f, a, b)
    h = b - a
    integral = h / 2 * (f(a) + f(b))
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



// Aplica la regla de Simpson
// para aproximar la integral definida
// de f entre a y b
function integral = reglaSimpson(f, a, b)
    h = (b - a) / 2
    integral = h / 3 * (f(a) + 4 * f(a + h) + f(b))
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


// Ejercicio 1
// TODO ENCONTRAR ERROR
// // a)
// disp('f(x)=log(x)')
// disp('a=1, b=2')

// deff('y=f(x)', 'y=log(x)')
// a = 1
// b = 2

// disp('Regla del Trapecio', reglaTrapecio(f, a, b))
// 0.3465736

// disp('Regla de Simpson', reglaSimpson(f, a, b))
// 0.3465736

// disp('Integral Real', integrate('f(x)', 'x', a, b))
// 0.3862944

// // b)
// disp('f(x)=log(x)')
// disp('a=0, b=0.1')

// deff('y=f(x)', 'y=x^(1/3)')
// a = 0
// b = .1

// disp('Regla del Trapecio', reglaTrapecio(f, a, b))
// 0.0232079

// disp('Regla de Simpson', reglaSimpson(f, a, b))
// 0.0322962

// disp('Integral Real', integrate('f(x)', 'x', a, b))
// 0.0348119

// // c)
// disp('f(x)=sin^2(x)')
// disp('a=0, b=%pi/3')

// deff('y=f(x)', "y=sin(x)^2")
// a = 0
// b = %pi/3

// disp('Regla del Trapecio', reglaTrapecio(f, a, b))
// 0.3926991

// disp('Regla de Simpson', reglaSimpson(f, a, b))
// 0.3054326

// disp('Integral Real', integrate('f(x)', 'x', a, b))
// 0.3070924


// Ejercicio 2 
// Ejercicio 3
function Ejercicio2y3(f, a, b, n)
    disp('a=' + string(a))
    disp('b=' + string(b))
    disp('n=' + string(n))
    disp('Metodo Compuesto del Trapecio: ' + string(metodoCompuestoTrapecio(f, a, b, n)))
    disp('Metodo Compuesto de Simpson: ' + string(metodoCompuestoSimpson(f, a, b, n)))
    disp('Integral Real: ' + string(integrate('f(x)', 'x', a, b, n)))
endfunction

// // a)
// deff('y=f(x)', 'y=1/x')
// disp('f(x)=1/x')
// Ejercicio2y3(f, 1, 3, 4)

// // b)
// deff('y=f(x)', 'y=x^3')
// disp('f(x)=x^3')
// Ejercicio2y3(f, 0, 2, 4)

// // c)
// deff('y=f(x)', 'y = x * (1 + x^2)^.5')
// disp('f(x) = x * (1 + x^2)^.5')
// Ejercicio2y3(f, 0, 3, 6)

// // d)
// deff('y=f(x)', 'y=sin(%pi*x)')
// disp('f(x)=sin(%pi*x)')
// Ejercicio2y3(f, 0, 1, 8)

// // e)
// deff('y=f(x)', 'y = x * sin(x)')
// disp('f(x) =  x * sin(x)')
// Ejercicio2y3(f, 0, 2*%pi, 8)

// // f)
// deff('y=f(x)', 'y = x^2 * exp(x)')
// disp('f(x) =  x^2 * exp(x)')
// Ejercicio2y3(f, 0, 1, 8)


// Ejercicio 4
deff('y=f(x)', 'y=1/(x+1)')
a = 0
b = 1.5
n = 10

// a)
t = metodoCompuestoTrapecio(f, a, b, n)
disp("Trapecio: " + string(t))

// b)
s = metodoCompuestoSimpson(f, a, b, n)
disp("Simpson: " + string(s))

// c)
I = .9262907

disp("Error absoluto Trapecio: " + string(abs(I - t)))
disp("Error relativo Trapecio: " + string(abs((I - t) / I)))

disp("Error absoluto Simpson: " + string(abs(I - s)))
disp("Error relativo Simpson: " + string(abs((I - s) / I)))



// Ejercicio 5













