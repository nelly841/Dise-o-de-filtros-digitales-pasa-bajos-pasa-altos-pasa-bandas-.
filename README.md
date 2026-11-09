# Filtros digitales para una señal de vibración

Proyecto reproducible para comparar filtros pasa bajos, pasa altos y pasa bandas sobre una lectura sintética de un sensor de vibración.

## Contenido

- `analisis_filtros.py`: genera la señal, diseña y aplica filtros, calcula amplitudes y crea las gráficas.
- `figuras/`: visualizaciones exportadas en PNG.
- `resultados/amplitudes_componentes.csv`: amplitud medida de cada tono antes y después del filtrado.
- `resultados/resumen.json`: parámetros, valores RMS y respuesta teórica.

## Requisitos y ejecución

Se recomienda Python 3.10 o posterior.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python analisis_filtros.py
```

En Windows, activar el entorno con `.venv\Scripts\activate`.

La ejecución usa una semilla fija para que el ruido sea reproducible. Las figuras y los resultados se guardan en las carpetas del proyecto. Los filtros se aplican en ambas direcciones para compensar el desfase en este análisis sin conexión; para procesamiento en tiempo real se requiere una implementación causal.

## Publicar el proyecto en GitHub

Después de crear un repositorio vacío en GitHub, abre una terminal en esta carpeta y ejecuta:

```bash
git init
git add .
git commit -m "Implementa análisis de filtros digitales"
git branch -M main
git remote add origin https://github.com/USUARIO/NOMBRE-DEL-REPOSITORIO.git
git push -u origin main
```

Sustituye `USUARIO/NOMBRE-DEL-REPOSITORIO` por la dirección del repositorio que hayas creado.

## Diseños incluidos

- Butterworth IIR pasa bajos, orden 6, frecuencia de corte 35 Hz.
- Chebyshev tipo I IIR pasa altos, orden 5, rizado de paso de 0.5 dB y corte de 80 Hz.
- FIR pasa bandas de 241 coeficientes con ventana Blackman, bandas de 55 a 90 Hz.

Los tonos de prueba están en 18, 72 y 145 Hz, con ruido blanco gaussiano de desviación estándar 0.42.
