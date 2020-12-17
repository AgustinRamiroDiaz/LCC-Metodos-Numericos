// Ejercicio 1

// Calcula las raíces de forma robusta para un polinomio de grado 2
// Retorna un vector columna con las soluciones de x- y x+
function x = raices(p)
    a = coeff(p, 2)
    b = coeff(p, 1)
    c = coeff(p, 0)

    determinante = b^2 - 4*a*c

    // Veo que el discriminante no sea negativo
    if determinante < 0 then x(1) = %nan; x(2) = %nan; return end

    sqrtDet = sqrt(determinante)

    if b < 0 
        then
        x(1) = 2*c / (-b + sqrtDet);
        x(2) = (-b + sqrtDet) / 2 / a;
    else
        if b == 0
            then 
                x(1) = sqrt(-c / a)
                x(2) = -x(1)
            else
                x(1) = (-b - sqrtDet) / 2 / a;
                x(2) = 2*c / (-b - sqrtDet);
        end
    end
endfunction

// Ejercicio 3
// a y c en pdf

// b
function v = hornerNew(p, x)
    v = hornerAux (p, x, 0)
    
endfunction

function v = hornerAux(p, x, n)
    if n == degree(p) 
        then v = coeff(p, n)
        else v = coeff(p, n) + x * hornerNewAux(p, x, n+1)
    end
endfunction

// Evalúa de forma eficiente el polinomio p
// aplicado a x con el algoritmo de Horner
function b = Horner(p, x)
    b = coeff(p, degree(p)) 
    for (i = degree(p) - 1: -1: 0)
        b = coeff(p, i) + x * b
    end
endfunction

//d

// Implementa una generalización del algoritmo Horner
// de tal forma que se pueda calcular p(x) y p'(x)
function x = HornerConDerivada(p, x)
    b = coeff(p, degree(p)) 
    c = b * x ^ (degree(p) - 1)

    for (i = degree(p) - 1: -1: 1)
        b = coeff(p, i) + x * b
        c = c + b * x ^ (i - 1)
    end

    // Hay que hacer una iteración por fuera ya que la derivada tiene un término menos
    b = coeff(p, 0) + x * b
    x(1) = b
    x(2) = c
endfunction


//Ejercicio 4

// Toma una función f, un valor v, un orden n y un paso h 
// y retorna el valor de evaluar la derivada de f de orden n en el punto v
// utilizando el cociente incremental
function r = derivar(f, v, n, h)
    if n == 0 
        then r = f (v)
        else r = (derivar (f, v + h, n-1) - derivar (f, v, n-1)) / h
    end
endfunction


function x = derivarIterativo(f,v,n,h)
    deff("y=D0F(x)", "y="+f)
    for i = 1:1:n-1
        deff("y=D"+string(i)+"F(x)", "y=(D"+string(i-1)+"F(x+h)-D"+string(i-1)+"F(x))/h")
    end
    deff ("y=DnF(x)", "y=(D"+string(n-1)+"F(x+h)-D"+string(n-1)+"F(x))/h")
    x = DnF(v)
endfunction

// a
// Estamos aproximando la función derviada a partir de un paso h
// También tenemos errores de redondeo

//b
// dividimos n veces por h, lo cual hace que dividimos por un número muy cercano a 0

//5
// Calcula el valor de un polinomio de Taylor de la función f aplicado en v 
// alrededor del punto a de grado n, aproximando sus derivadas con el cociente incremental con paso h
function r = Taylor(f, v, n, h, a)
    r = 0
    for (i = 0:n) 
        r = r + derivar (f, a, i, h) * (v - a) ^ i / factorial (i)
    end
endfunction

//6

// b

coeficientes = 1 ./ factorial(0:10)
t = poly(coeficientes, "x", "c")
    // a
a = Horner(t, -2)
    // b 
b = 1 / Horner(t, 2)

valorCorrecto = 0.135

aerror = abs(a - valorCorrecto) / valorCorrecto
berror = abs(b - valorCorrecto) / valorCorrecto

// 1/e^2 aproxima mejor ya que hace solo suma de términos positivos
// mientras que e^-2 suma términos alternantes, lo cual aumenta el error