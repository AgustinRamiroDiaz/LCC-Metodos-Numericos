
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
    integral = f(a) / 2 + sum(f(x)) + f(b) / 2
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

    integral = f(a) + 4 * sum(f(ximpares)) + 2 * sum(f(xpares)) + f(b)
    integral = h / 3 * integral
endfunction
