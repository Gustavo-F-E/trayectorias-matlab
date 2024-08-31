export function obtenerEnteroPositivo(mensaje) {
    while (true) {
        let valor = parseInt(prompt(mensaje));
        if (!isNaN(valor) && valor > 0) {
            return valor;
        } else {
            console.log("Error: Por favor, ingrese un número entero positivo.");
        }
    }
}

export function obtenerFlotantePositivo(mensaje) {
    while (true) {
        let valor = parseFloat(prompt(mensaje));
        if (!isNaN(valor) && valor > 0) {
            return valor;
        } else {
            console.log("Error: Por favor, ingrese un número positivo.");
        }
    }
}

export function obtenerAnguloBobinado(mensaje) {
    while (true) {
        let valor = parseInt(prompt(mensaje));
        if (isNaN(valor) || valor <= 0) {
            console.log("Error: Por favor, ingrese un número entero positivo.");
        } else if (valor < 5) {
            console.log("No es posible realizar bobinados de menos de 5º");
        } else if (valor > 89) {
            console.log("No es posible realizar bobinados mayores a 89º");
        } else {
            return valor;
        }
    }
}
