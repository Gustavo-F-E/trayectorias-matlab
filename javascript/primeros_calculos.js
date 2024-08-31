export function primeros_calculos(
    ancho_de_la_fibra,
    diametro_del_mandril,
    longitud_a_bobinar,
    angulo_de_bobinado
) {
    let b = ancho_de_la_fibra;
    let D = diametro_del_mandril;
    let L = longitud_a_bobinar;
    let X_lim = L;

    let calculo_1 =
        (2 * Math.PI * diametro_del_mandril) /
        (ancho_de_la_fibra / Math.cos((angulo_de_bobinado * Math.PI) / 180));
    let calculo_2 = Math.floor(
        (2 * Math.PI * diametro_del_mandril) /
            (ancho_de_la_fibra / Math.cos((angulo_de_bobinado * Math.PI) / 180))
    );
    let calculo_3 = calculo_1 - calculo_2;

    let alfa, calculo_4;

    if (0 < calculo_3 && calculo_3 < 0.5) {
        calculo_4 =
            (calculo_2 * ancho_de_la_fibra) /
            (2 * Math.PI * diametro_del_mandril);
        alfa = Math.acos(calculo_4);
        console.log("Se ha aplicado el caso donde 0 < calculo_3 < 0.5");
    } else if (calculo_3 >= 0.5) {
        calculo_2 = Math.ceil(calculo_1);
        calculo_4 =
            (calculo_2 * ancho_de_la_fibra) /
            (2 * Math.PI * diametro_del_mandril);
        alfa = Math.acos(calculo_4);
        console.log("Se ha aplicado el caso donde calculo_3 >= 0.5");
    } else if (calculo_3 === 0) {
        alfa = (angulo_de_bobinado * Math.PI) / 180;
        console.log("Se ha aplicado el caso donde calculo_3 == 0");
    } else {
        console.log("No se ha aplicado ningún caso");
    }

    let b_eff = b / Math.cos(alfa);

    // ... (resto del código con console.log para imprimir resultados)

    let Y_lim = 2 * Math.PI * D;
    let numero_de_divisiones = Math.floor(Y_lim / b_eff);
    let paso_en_y = Y_lim / numero_de_divisiones;
    let vector_y = Array.from(
        { length: numero_de_divisiones },
        (_, i) => i * paso_en_y
    );

    return {
        alfa,
        b_eff,
        numero_de_divisiones,
        X_lim,
        Y_lim,
        D,
        L,
        b,
        vector_y,
    };
}