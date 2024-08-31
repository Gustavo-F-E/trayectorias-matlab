/*

export function graficoMandril(D, L, containerId) {
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

    // Crear el cilindro
    const geometry = new THREE.CylinderGeometry(D / 2, D / 2, L, 32, 1, true);
    const material = new THREE.MeshPhongMaterial({
        color: 0x00ffff,
        transparent: true,
        opacity: 0.8,
        side: THREE.DoubleSide,
    });
    const cylinder = new THREE.Mesh(geometry, material);
    cylinder.rotation.z = Math.PI / 2; // Rotar para que coincida con la orientación del gráfico original
    scene.add(cylinder);

    // Configurar la cámara
    camera.position.set(L, D, L);
    camera.lookAt(L / 2, 0, 0);

    // Añadir luz
    const light = new THREE.PointLight(0xffffff, 1, 100);
    light.position.set(L, D, L);
    scene.add(light);

    // Añadir controles de órbita
    const controls = new OrbitControls(camera, renderer.domElement);

    // Añadir ejes
    const axesHelper = new THREE.AxesHelper(Math.max(L, D));
    scene.add(axesHelper);

    // Añadir texto
    const canvas = document.createElement("canvas");
    const context = canvas.getContext("2d");
    context.font = "24px Arial";
    context.fillStyle = "white";
    context.fillText(`Diámetro (D) = ${D}`, 10, 30);
    context.fillText(`Longitud (L) = ${L}`, 10, 60);
    const texture = new THREE.CanvasTexture(canvas);
    const spriteMaterial = new THREE.SpriteMaterial({ map: texture });
    const sprite = new THREE.Sprite(spriteMaterial);
    sprite.scale.set(0.5, 0.25, 1);
    sprite.position.set(0, D, 0);
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
}*/
export function graficoMandril(D, L, containerId) {
    const theta = Array.from({ length: 100 }, (_, i) => (i * 2 * Math.PI) / 99);
    const x = theta.map((t) => (D / 2) * Math.cos(t));
    const y = theta.map((t) => (D / 2) * Math.sin(t));
    const z = Array.from({ length: 100 }, () => L / 2);

    const data = [
        {
            type: "scatter3d",
            mode: "lines",
            x: x,
            y: y,
            z: z,
            opacity: 0.8,
            line: {
                width: 6,
                color: "cyan",
                reversescale: false,
            },
        },
        {
            type: "scatter3d",
            mode: "lines",
            x: x,
            y: y,
            z: z.map((v) => -v),
            opacity: 0.8,
            line: {
                width: 6,
                color: "cyan",
                reversescale: false,
            },
        },
    ];

    const layout = {
        title: "Mandril",
        scene: {
            xaxis: { title: "X" },
            yaxis: { title: "Y" },
            zaxis: { title: "Z" },
            aspectratio: { x: 1, y: 1, z: 1 },
        },
        annotations: [
            {
                x: 0,
                y: 0,
                z: L / 2,
                text: `Diámetro (D) = ${D}<br>Longitud (L) = ${L}`,
                showarrow: false,
            },
        ],
    };

    Plotly.newPlot(containerId, data, layout);
}