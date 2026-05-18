import numpy as np
import matplotlib.pyplot as plt
import os

# Rutas
csv_in = r"c:\astro\Codigos practica + docs + papers\codigos\Oka et al 2011 fig 7\oka et al 2011 dataset fig7 eje x invertido.csv"
csv_out = r"c:\astro\Codigos practica + docs + papers\codigos\Oka et al 2011 fig 7\oka_et_al_2011_data_corrected.csv"
plot_out = r"c:\astro\Codigos practica + docs + papers\codigos\oka_fig7_reproduced.png"

# Leer el CSV manualmente
with open(csv_in, 'r') as f:
    lines = f.readlines()

x_vals = []
y_vals = []

# No hay encabezado, empezamos desde 0
for line in lines:
    line = line.strip()
    if not line or 'x;y' in line:
        continue
    # Reemplazar coma por punto
    line = line.replace(',', '.')
    parts = line.split(';')
    if len(parts) == 2:
        x = float(parts[0])
        y = float(parts[1])
        x_vals.append(x)
        y_vals.append(y)

x_vals = np.array(x_vals)
y_vals = np.array(y_vals)

# El usuario indicó: "invertir la lista del eje x, es decir el primero es el ulimto"
# Además sacamos el valor absoluto porque eran negativos.
x_positive = np.abs(x_vals)
x_corrected = x_positive[::-1]

# Guardar el nuevo CSV manualmente sin pandas
with open(csv_out, 'w') as f:
    f.write("Mdot_Msun_yr;Rsnow_AU\n")
    for i in range(len(x_corrected)):
        f.write(f"{x_corrected[i]};{y_vals[i]}\n")
print(f"Datos corregidos guardados en: {csv_out}")

# Plotear igual a la figura original
plt.figure(figsize=(8, 6))

# Dibujar la línea de la curva extraída
plt.plot(x_corrected, y_vals, 'k-', linewidth=2, label='our calculation')

# Configurar escalas
plt.xscale('log')
plt.yscale('log')

# Invertir el eje X para que vaya de 10^-7 a 10^-12 igual que en el paper
plt.xlim(1e-7, 1e-12)
plt.ylim(0.1, 10)

# Personalizar los ticks para que coincidan con la figura
# Eje X: 10^-7, 10^-8, 10^-9, 10^-10, 10^-11, 10^-12
plt.xticks([1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12])
# Eje Y: 0.1, 1, 10
plt.yticks([0.1, 1, 10], ['0.1', '1', '10'])

# Etiquetas
plt.xlabel(r'$\dot{M} \ [M_\odot \ \mathrm{yr}^{-1}]$', fontsize=16)
plt.ylabel(r'$R_{\mathrm{snow}} \ [\mathrm{AU}]$', fontsize=16)



# Estilo de los ticks
plt.tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize=14)

plt.tight_layout()
plt.savefig(plot_out, dpi=300)
print(f"Gráfico guardado en: {plot_out}")
