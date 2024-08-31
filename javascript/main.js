//Módulos externos
//import * as THREE from "./three.js";
//import { OrbitControls } from "./orbitControls.js";

//Módulos propios
import { primeros_calculos } from "./primeros_calculos.js";
import { graficoMandril } from "./grafico_mandril.js";
import { graficoTrayectoria } from "./grafico_trayectoria.js";

// Función principal que se ejecutará cuando se haga clic en el botón
function calcularYGraficar() {
    const ancho_de_la_fibra = parseInt(
        document.getElementById("anchoFibra").value
    );
    const angulo_de_bobinado = parseInt(
        document.getElementById("anguloBobinado").value
    );
    const diametro_del_mandril = parseFloat(
        document.getElementById("diametroMandril").value
    );
    const longitud_a_bobinar = parseFloat(
        document.getElementById("longitud_a_bobinar").value
    );

    console.log("Ancho de la fibra:", ancho_de_la_fibra);
    console.log("Ángulo de bobinado:", angulo_de_bobinado);
    console.log("Diámetro del mandril:", diametro_del_mandril);
    console.log("Longitud a bobinar:", longitud_a_bobinar);

    // Primeros cálculos
    let { alfa, b_eff, numero_de_divisiones, X_lim, Y_lim, D, L, b, vector_y } =
        primeros_calculos(
            ancho_de_la_fibra,
            diametro_del_mandril,
            longitud_a_bobinar,
            angulo_de_bobinado
        );

    // Limpiar gráficos anteriores
    document.getElementById("mandrilGrafico").innerHTML = "";
    document.getElementById("trayectoriaGrafico").innerHTML = "";

    // Graficar el mandril
    graficoMandril(D, L, "mandrilGrafico");

    // Gráfico en los ejes X e Y
    for (let i = 0; i < numero_de_divisiones; i++) {
        graficoTrayectoria(
            i,
            vector_y,
            b,
            b_eff,
            alfa,
            D,
            L,
            X_lim,
            Y_lim,
            numero_de_divisiones,
            vector_y,
            "trayectoriaGrafico"
        );
    }
}

// Agregar el event listener al botón de cálculo
document.addEventListener("DOMContentLoaded", function () {
    document
        .getElementById("calcularBtn")
        .addEventListener("click", calcularYGraficar);
});