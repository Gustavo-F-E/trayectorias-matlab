//import * as THREE from "https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.module.js";
//import { OrbitControls } from "https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/examples/js/controls/OrbitControls.js";

export function evaluarXY(yOrigen, XLim, YLim, b, alfa) {
    const bEffY = b / Math.cos(alfa);
    const bEffX = b / Math.sin(alfa);
    let n = 0;
    let xIni = 0;
    const paso = XLim / 1000;
    const xValues = [];
    const yValues = [];
    let contador = 0;

    while (xIni < XLim) {
        contador++;
        let x = xIni;
        let y = x * Math.tan(alfa) - n * YLim + yOrigen;
        if (y <= YLim && x < XLim) {
            xIni = xIni + paso;
            xValues.push(x);
            yValues.push(y);
        } else if (y <= YLim && x >= XLim) {
            xIni = XLim;
            xValues.push(xIni);
            yValues.push(x * Math.tan(alfa) - n * YLim + yOrigen);
        } else if (y > YLim && x < XLim) {
            n = n + 1;
            xIni = (n * YLim) / Math.tan(alfa);
            xValues.push(xIni);
            yValues.push(xIni * Math.tan(alfa) - n * YLim + yOrigen);
            xIni = xIni + paso;
        } else {
            console.log("Caso 4: problemas");
        }
    }
    xValues.push(XLim);
    yValues.push(XLim * Math.tan(alfa) - n * YLim + yOrigen);
    return [xValues, yValues];
}
/*
export function graficoTrayectoria(
    i,
    yOrigen,
    b,
    bEff,
    alfa,
    D,
    L,
    XLim,
    YLim,
    numeroDeDivisiones,
    vectorY,
    containerId
) {
    const [xValues, yValues] = evaluarXY(yOrigen[i], XLim, YLim, b, alfa);
    const yParalela = yValues.map((y) => y + bEff);

    // Configuración de la escena
    const scene = new THREE.Scene();
    const camera = new THREE.PerspectiveCamera(
        75,
        window.innerWidth / window.innerHeight,
        0.1,
        1000
    );
    const renderer = new THREE.WebGLRenderer();
    renderer.setSize(window.innerWidth, window.innerHeight);
    document.getElementById(containerId).appendChild(renderer.domElement);

    // Crear las líneas de trayectoria
    const material = new THREE.LineBasicMaterial({ color: 0x0000ff });
    const points = xValues.map(
        (x, index) => new THREE.Vector3(x, yValues[index], 0)
    );
    const geometry = new THREE.BufferGeometry().setFromPoints(points);
    const line = new THREE.Line(geometry, material);
    scene.add(line);

    // Crear la línea paralela
    const parallelPoints = xValues.map(
        (x, index) => new THREE.Vector3(x, yParalela[index], 0)
    );
    const parallelGeometry = new THREE.BufferGeometry().setFromPoints(
        parallelPoints
    );
    const dashedMaterial = new THREE.LineDashedMaterial({
        color: 0x0000ff,
        dashSize: 3,
        gapSize: 1,
    });
    const parallelLine = new THREE.Line(parallelGeometry, dashedMaterial);
    parallelLine.computeLineDistances();
    scene.add(parallelLine);

    // Configurar la cámara
    camera.position.set(XLim / 2, YLim / 2, Math.max(XLim, YLim));
    camera.lookAt(XLim / 2, YLim / 2, 0);

    // Añadir controles de órbita
    const controls = new OrbitControls(camera, renderer.domElement);

    // Añadir ejes
    const axesHelper = new THREE.AxesHelper(Math.max(XLim, YLim));
    scene.add(axesHelper);

    // Añadir texto con información
    const canvas = document.createElement("canvas");
    const context = canvas.getContext("2d");
    context.font = "12px Arial";
    context.fillStyle = "white";
    context.fillText(
        `Ángulo (alfa) = ${((alfa * 180) / Math.PI).toFixed(2)}°`,
        10,
        20
    );
    context.fillText(`Ancho efectivo (b_eff) = ${bEff.toFixed(2)}`, 10, 40);
    context.fillText(`Diámetro = ${D}`, 10, 60);
    context.fillText(`Longitud = ${L}`, 10, 80);
    context.fillText(`Longitud del cilindro extendida = ${YLim}`, 10, 100);
    context.fillText(`Número_de_divisiones = ${numeroDeDivisiones}`, 10, 120);

    const texture = new THREE.CanvasTexture(canvas);
    const spriteMaterial = new THREE.SpriteMaterial({ map: texture });
    const sprite = new THREE.Sprite(spriteMaterial);
    sprite.scale.set(100, 50, 1);
    sprite.position.set(XLim / 2, YLim, 0);
    scene.add(sprite);

    // Función de animación
    function animate() {
        requestAnimationFrame(animate);
        controls.update();
        renderer.render(scene, camera);
    }

    // Iniciar la animación
    animate();

    // Manejar el redimensionamiento de la ventana
    window.addEventListener("resize", onWindowResize, false);
    function onWindowResize() {
        camera.aspect = window.innerWidth / window.innerHeight;
        camera.updateProjectionMatrix();
        renderer.setSize(window.innerWidth, window.innerHeight);
    }
}
    */
export function graficoTrayectoria(
    i,
    yOrigen,
    b,
    bEff,
    alfa,
    D,
    L,
    XLim,
    YLim,
    numeroDeDivisiones,
    vectorY,
    containerId
) {
    const [xValues, yValues] = evaluarXY(yOrigen[i], XLim, YLim, b, alfa);
    const yParalela = yValues.map((y) => y + bEff);

    const data = [
        {
            type: "scatter3d",
            mode: "lines",
            x: xValues,
            y: yValues,
            z: Array(xValues.length).fill(0),
            line: {
                width: 6,
                color: "blue",
                reversescale: false,
            },
            name: "Línea principal",
        },
        {
            type: "scatter3d",
            mode: "lines",
            x: xValues,
            y: yParalela,
            z: Array(xValues.length).fill(0),
            line: {
                width: 6,
                color: "red",
                dash: "dash",
                reversescale: false,
            },
            name: "Línea paralela",
        },
    ];

    const layout = {
        title: `Trayectoria ${i}`,
        scene: {
            xaxis: { title: "X", range: [0, XLim] },
            yaxis: { title: "Y", range: [0, YLim] },
            zaxis: { title: "Z", range: [-1, 1] },
            aspectratio: { x: 1, y: 1, z: 0.1 },
        },
        annotations: [
            {
                x: XLim / 2,
                y: YLim,
                z: 0,
                text:
                    `Ángulo (alfa) = ${((alfa * 180) / Math.PI).toFixed(
                        2
                    )}°<br>` +
                    `Ancho efectivo (b_eff) = ${bEff.toFixed(2)}<br>` +
                    `Diámetro = ${D}<br>` +
                    `Longitud = ${L}<br>` +
                    `Longitud del cilindro extendida = ${YLim}<br>` +
                    `Número_de_divisiones = ${numeroDeDivisiones}`,
                showarrow: false,
            },
        ],
    };

    Plotly.newPlot(containerId, data, layout);
}