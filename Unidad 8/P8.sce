
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

function i = GTrapecio(f, c, d, x, n)
    a = c(x)
    b = d(x)

    h = (b - a) / n
    Y = a+h:h:b-h
    i = (f(x, a) + f(x, b)) / 2 
    for y=Y
        i = i + f(x, y)
    end
    i = i * h
endfunction

function I = integralBidimensionalTrapecio(f, a, b, c, d, n)
    h = (b - a) / n
    X = a+h:h:b-h
    I = (GTrapecio(f, c, d, a, n) + GTrapecio(f, c, d, b, n)) / 2 
    for x=X
        I = I + GTrapecio(f, c, d, x, n)
    end
    I = I * h
endfunction


deff('y=f(x, y)', 'y=sin(x+y)')
deff('y=c(x)', 'y=0')
deff('y=d(x)', 'y=1')

r = integralBidimensionalTrapecio(f, 0, 2, c, d, 2)



// Ejercicio 6
function i = GSimpson(f, c, d, x, n)
    a = c(x)
    b = d(x)

    h = (b - a) / n
    yimpares = a+h:2*h:b-h
    ypares = a+2*h:2*h:b-2*h
    i = f(x, a) + f(x, b)
    for y=yimpares
        i = i + 4 * f(x, y)
    end
    for y=ypares
        i = i + 2 * f(x, y)
    end
    i = i * h / 3
endfunction

function I = integralBidimensionalSimpson(f, a, b, c, d, n)
    h = (b - a) / n
    ximpares = a+h:2*h:b-h
    xpares = a+2*h:2*h:b-2*h
    I = GSimpson(f, c, d, a, n) + GSimpson(f, c, d, b, n) 
    for x=ximpares
        I = I + 4 * GSimpson(f, c, d, x, n)
    end
    for x=xpares
        I = I + 2 * GSimpson(f, c, d, x, n)
    end
    I = I * h / 3
endfunction


deff('y=f(x, y)', 'y=1')
deff('y=c(x)', 'y=-(x * (2 - x))^.5')
deff('y=d(x)', 'y=(x * (2 - x))^.5')

deff('y=c(x)', 'y=-sin(acos(x-1))')
deff('y=d(x)', 'y=sin(acos(x-1))')


t = integralBidimensionalTrapecio(f, 0, 2, c, d, 100)
s = integralBidimensionalSimpson(f, 0, 2, c, d, 100)

// TODO OTROS MÉTODOS
