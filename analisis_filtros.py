"""Diseño y evaluación de filtros para una señal sintética de vibración."""

from pathlib import Path
import csv
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, cheby1, sosfiltfilt, firwin, filtfilt, sosfreqz, freqz


ROOT = Path(__file__).resolve().parent
FIG = ROOT / "figuras"
OUT = ROOT / "resultados"
FIG.mkdir(exist_ok=True)
OUT.mkdir(exist_ok=True)


def amplitud_en_frecuencia(x, t, frecuencia):
    """Estima amplitud pico ajustando seno, coseno y componente continua."""
    matriz = np.column_stack((np.sin(2 * np.pi * frecuencia * t),
                              np.cos(2 * np.pi * frecuencia * t),
                              np.ones_like(t)))
    coef, *_ = np.linalg.lstsq(matriz, x, rcond=None)
    return float(np.hypot(coef[0], coef[1]))


def db_seguro(valor):
    return float(20 * np.log10(max(valor, 1e-12)))


def main():
    # Señal sintética de un sensor de vibración: tres componentes y ruido blanco.
    fs = 1000.0
    duracion = 6.0
    t = np.arange(int(fs * duracion)) / fs
    rng = np.random.default_rng(2026)
    amplitudes = {18.0: 1.0, 72.0: 0.70, 145.0: 0.45}
    limpia = sum(a * np.sin(2 * np.pi * f * t) for f, a in amplitudes.items())
    ruido = rng.normal(0.0, 0.42, size=t.size)
    entrada = limpia + ruido

    # IIR pasa bajos: conserva la oscilación lenta de 18 Hz.
    sos_lp = butter(6, 35.0, btype="lowpass", fs=fs, output="sos")
    salida_lp = sosfiltfilt(sos_lp, entrada)

    # IIR pasa altos: conserva la componente rápida de 145 Hz.
    sos_hp = cheby1(5, 0.5, 80.0, btype="highpass", fs=fs, output="sos")
    salida_hp = sosfiltfilt(sos_hp, entrada)

    # FIR pasa bandas: aísla la vibración intermedia alrededor de 72 Hz.
    taps_bp = firwin(241, [55.0, 90.0], pass_zero=False,
                     fs=fs, window="blackman")
    salida_bp = filtfilt(taps_bp, [1.0], entrada)

    resultados = {
        "Pasa bajos Butterworth": salida_lp,
        "Pasa altos Chebyshev I": salida_hp,
        "Pasa bandas FIR Blackman": salida_bp,
    }

    # Se descartan extremos para reducir la influencia del relleno de filtfilt.
    margen = 0.5
    valido = (t >= margen) & (t <= duracion - margen)
    tv = t[valido]
    xv = entrada[valido]
    filas = []
    for frecuencia in amplitudes:
        amp_in = amplitud_en_frecuencia(xv, tv, frecuencia)
        for nombre, y in resultados.items():
            amp_out = amplitud_en_frecuencia(y[valido], tv, frecuencia)
            filas.append({
                "frecuencia_hz": frecuencia,
                "filtro": nombre,
                "amplitud_entrada": amp_in,
                "amplitud_salida": amp_out,
                "ganancia_medida_db": db_seguro(amp_out / amp_in),
            })

    with (OUT / "amplitudes_componentes.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=filas[0].keys())
        writer.writeheader()
        writer.writerows(filas)

    resumen = {
        "frecuencia_muestreo_hz": fs,
        "duracion_s": duracion,
        "semilla_ruido": 2026,
        "rms_entrada": float(np.sqrt(np.mean(entrada[valido] ** 2))),
        "rms_salida": {name: float(np.sqrt(np.mean(y[valido] ** 2)))
                       for name, y in resultados.items()},
        "ganancia_teorica_db": {},
    }
    f_lp, h_lp = sosfreqz(sos_lp, worN=8192, fs=fs)
    f_hp, h_hp = sosfreqz(sos_hp, worN=8192, fs=fs)
    f_bp, h_bp = freqz(taps_bp, worN=8192, fs=fs)
    for name, frequencies, response in [
        ("Pasa bajos Butterworth", f_lp, h_lp),
        ("Pasa altos Chebyshev I", f_hp, h_hp),
        ("Pasa bandas FIR Blackman", f_bp, h_bp),
    ]:
        resumen["ganancia_teorica_db"][name] = {
            str(int(freq)): db_seguro(float(np.interp(freq, frequencies, np.abs(response))))
            for freq in [18, 72, 145]
        }
    (OUT / "resumen.json").write_text(json.dumps(resumen, indent=2, ensure_ascii=False), encoding="utf-8")

    # Figura 1: tramo temporal y espectro de la medición original.
    fig, ax = plt.subplots(2, 1, figsize=(10, 6.8), constrained_layout=True)
    tramo = t < 0.35
    ax[0].plot(t[tramo], entrada[tramo], color="#1e5c65", lw=0.8)
    ax[0].set(title="Registro de vibración con interferencia aleatoria", xlabel="Tiempo (s)", ylabel="Amplitud")
    ax[0].grid(alpha=0.25)
    freqs = np.fft.rfftfreq(t.size, 1 / fs)
    espectro = 2 * np.abs(np.fft.rfft(entrada)) / t.size
    ax[1].plot(freqs, espectro, color="#c56b3f", lw=0.8)
    ax[1].set(xlim=(0, 220), title="Espectro de amplitud antes del filtrado", xlabel="Frecuencia (Hz)", ylabel="Amplitud estimada")
    ax[1].grid(alpha=0.25)
    for f in amplitudes:
        ax[1].axvline(f, color="#59636e", ls="--", lw=0.8)
    fig.savefig(FIG / "entrada_tiempo_espectro.png", dpi=180)
    plt.close(fig)

    # Figura 2: respuestas de magnitud calculadas a partir de los coeficientes.
    fig, ax = plt.subplots(figsize=(10, 5.3), constrained_layout=True)
    for name, frequencies, response, color in [
        ("Pasa bajos Butterworth", f_lp, h_lp, "#147d74"),
        ("Pasa altos Chebyshev I", f_hp, h_hp, "#bd5b36"),
        ("Pasa bandas FIR Blackman", f_bp, h_bp, "#415f9b"),
    ]:
        ax.plot(frequencies, 20 * np.log10(np.maximum(np.abs(response), 1e-8)), label=name, color=color, lw=1.8)
    for f in amplitudes:
        ax.axvline(f, color="#888888", ls=":", lw=0.8)
    ax.set(xlim=(0, 220), ylim=(-100, 5), title="Respuesta en frecuencia de los tres diseños", xlabel="Frecuencia (Hz)", ylabel="Magnitud (dB)")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, loc="lower right")
    fig.savefig(FIG / "respuesta_frecuencia.png", dpi=180)
    plt.close(fig)

    # Figura 3: comparación temporal en el intervalo de observación.
    fig, axes = plt.subplots(4, 1, figsize=(10, 8), sharex=True, constrained_layout=True)
    signals = [("Entrada", entrada, "#333f4f"),
               ("Pasa bajos", salida_lp, "#147d74"),
               ("Pasa altos", salida_hp, "#bd5b36"),
               ("Pasa bandas", salida_bp, "#415f9b")]
    for axis, (label, signal, color) in zip(axes, signals):
        axis.plot(t[tramo], signal[tramo], color=color, lw=0.8)
        axis.set_ylabel(label)
        axis.grid(alpha=0.22)
    axes[-1].set_xlabel("Tiempo (s)")
    fig.suptitle("Efecto temporal de cada filtro sobre la misma medición")
    fig.savefig(FIG / "comparacion_temporal.png", dpi=180)
    plt.close(fig)

    # Figura 4: amplitud de cada componente al pasar por cada filtro.
    fig, ax = plt.subplots(figsize=(9.5, 5.2), constrained_layout=True)
    labels = ["18 Hz", "72 Hz", "145 Hz"]
    x = np.arange(len(labels))
    width = 0.2
    input_amps = [amplitud_en_frecuencia(xv, tv, f) for f in amplitudes]
    ax.bar(x - 1.5 * width, input_amps, width, label="Entrada", color="#697583")
    palette = ["#147d74", "#bd5b36", "#415f9b"]
    for j, (name, y) in enumerate(resultados.items()):
        values = [amplitud_en_frecuencia(y[valido], tv, f) for f in amplitudes]
        ax.bar(x + (j - 0.5) * width, values, width, label=name, color=palette[j])
    ax.set_xticks(x, labels)
    ax.set(title="Amplitud estimada de los tonos antes y después", xlabel="Componente de la señal", ylabel="Amplitud pico")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(FIG / "amplitudes_componentes.png", dpi=180)
    plt.close(fig)

    print(json.dumps(resumen, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
