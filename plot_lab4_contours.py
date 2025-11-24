#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Skrypt do generowania wykresów poziomic dla Lab 4
6 wykresów z trajektoriami optymalizacji naniesionych na wykres poziomic funkcji testowej
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import cm

# Funkcja testowa Lab 4
def f(x1, x2):
    """f(x1, x2) = (1/6)*x1^6 - 1.05*x1^4 + 2*x1^2 + x2^2 + x1*x2"""
    return (1/6) * x1**6 - 1.05 * x1**4 + 2 * x1**2 + x2**2 + x1 * x2

# Siatka dla wykresu poziomic
x1_range = np.linspace(-2, 2, 200)
x2_range = np.linspace(-2, 2, 200)
X1, X2 = np.meshgrid(x1_range, x2_range)
Z = f(X1, X2)

# Wczytanie danych z CSV
try:
    df = pd.read_csv('data/lab4_wykresy.csv', header=None)
    
    # Kolumny CSV (18 kolumn):
    # 0-1: SD s=0.05 (x1, x2)
    # 2-3: SD s=0.25 (x1, x2)
    # 4-5: SD zmiennokrokowy (x1, x2)
    # 6-7: CG s=0.05 (x1, x2)
    # 8-9: CG s=0.25 (x1, x2)
    # 10-11: CG zmiennokrokowy (x1, x2)
    # 12-13: Newton s=0.05 (x1, x2)
    # 14-15: Newton s=0.25 (x1, x2)
    # 16-17: Newton zmiennokrokowy (x1, x2)
    
    col_names = ['SD_x1_005', 'SD_x2_005', 'SD_x1_025', 'SD_x2_025', 'SD_x1_var', 'SD_x2_var',
                 'CG_x1_005', 'CG_x2_005', 'CG_x1_025', 'CG_x2_025', 'CG_x1_var', 'CG_x2_var',
                 'Newton_x1_005', 'Newton_x2_005', 'Newton_x1_025', 'Newton_x2_025', 'Newton_x1_var', 'Newton_x2_var']
    
    df.columns = col_names
    
except FileNotFoundError:
    print("BŁĄD: Nie znaleziono pliku data/lab4_wykresy.csv")
    print("Najpierw uruchom program C++ aby wygenerować dane.")
    exit(1)

# Funkcja pomocnicza do rysowania wykresu
def plot_contour_with_path(ax, paths, labels, colors, title):
    """
    Rysuje wykres poziomic z trajektoriami optymalizacji
    
    ax: matplotlib axis
    paths: lista ścieżek [(x1_array, x2_array), ...]
    labels: lista etykiet dla każdej ścieżki
    colors: lista kolorów dla każdej ścieżki
    title: tytuł wykresu
    """
    # Wykres poziomic z wypełnieniem kolorowym
    contourf = ax.contourf(X1, X2, Z, levels=50, cmap='jet', alpha=0.8)
    contour = ax.contour(X1, X2, Z, levels=20, colors='black', alpha=0.4, linewidths=0.5)
    ax.clabel(contour, inline=True, fontsize=7, fmt='%.2f')
    
    # Naniesienie trajektorii
    for (x1, x2), label, color in zip(paths, labels, colors):
        # Usuń wartości NaN
        valid = ~(np.isnan(x1) | np.isnan(x2))
        x1_clean = x1[valid]
        x2_clean = x2[valid]
        
        if len(x1_clean) > 0:
            # Ścieżka optymalizacji
            ax.plot(x1_clean, x2_clean, '-o', label=label, color=color, 
                   markersize=4, linewidth=1.5, alpha=0.8)
            # Punkt startowy (większy marker)
            ax.plot(x1_clean[0], x2_clean[0], 'o', color=color, 
                   markersize=10, markeredgecolor='black', markeredgewidth=1.5)
            # Punkt końcowy (gwiazdka)
            ax.plot(x1_clean[-1], x2_clean[-1], '*', color=color, 
                   markersize=15, markeredgecolor='black', markeredgewidth=1.5)
    
    ax.set_xlabel('x1')
    ax.set_ylabel('x2')
    ax.set_title(title)
    ax.legend(loc='best', fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)
    ax.set_aspect('equal')

# ===== WYKRES 1: s=0.05, wszystkie metody =====
fig1, ax1 = plt.subplots(figsize=(10, 8))
paths1 = [
    (df['SD_x1_005'].values, df['SD_x2_005'].values),
    (df['CG_x1_005'].values, df['CG_x2_005'].values),
    (df['Newton_x1_005'].values, df['Newton_x2_005'].values)
]
labels1 = ['Najszybszy spadek', 'Gradienty sprzężone', 'Newton']
colors1 = ['red', 'magenta', 'green']
plot_contour_with_path(ax1, paths1, labels1, colors1, 'Wykres 1: s = 0.05 - wszystkie metody')
plt.tight_layout()
plt.savefig('data/lab4_wykres1.png', dpi=300, bbox_inches='tight')
print("Zapisano: data/lab4_wykres1.png")

# ===== WYKRES 2: s=0.25, wszystkie metody =====
fig2, ax2 = plt.subplots(figsize=(10, 8))
paths2 = [
    (df['SD_x1_025'].values, df['SD_x2_025'].values),
    (df['CG_x1_025'].values, df['CG_x2_025'].values),
    (df['Newton_x1_025'].values, df['Newton_x2_025'].values)
]
labels2 = ['Najszybszy spadek', 'Gradienty sprzężone', 'Newton']
colors2 = ['red', 'magenta', 'green']
plot_contour_with_path(ax2, paths2, labels2, colors2, 'Wykres 2: s = 0.25 - wszystkie metody')
plt.tight_layout()
plt.savefig('data/lab4_wykres2.png', dpi=300, bbox_inches='tight')
print("Zapisano: data/lab4_wykres2.png")

# ===== WYKRES 3: zmiennokrokowy, wszystkie metody =====
fig3, ax3 = plt.subplots(figsize=(10, 8))
paths3 = [
    (df['SD_x1_var'].values, df['SD_x2_var'].values),
    (df['CG_x1_var'].values, df['CG_x2_var'].values),
    (df['Newton_x1_var'].values, df['Newton_x2_var'].values)
]
labels3 = ['Najszybszy spadek', 'Gradienty sprzężone', 'Newton']
colors3 = ['red', 'magenta', 'green']
plot_contour_with_path(ax3, paths3, labels3, colors3, 'Wykres 3: zmiennokrokowy (złoty podział) - wszystkie metody')
plt.tight_layout()
plt.savefig('data/lab4_wykres3.png', dpi=300, bbox_inches='tight')
print("Zapisano: data/lab4_wykres3.png")

# ===== WYKRES 4: Najszybszy spadek, wszystkie długości kroku =====
fig4, ax4 = plt.subplots(figsize=(10, 8))
paths4 = [
    (df['SD_x1_005'].values, df['SD_x2_005'].values),
    (df['SD_x1_025'].values, df['SD_x2_025'].values),
    (df['SD_x1_var'].values, df['SD_x2_var'].values)
]
labels4 = ['s = 0.05', 's = 0.25', 'zmiennokrokowy']
colors4 = ['magenta', 'red', 'green']
plot_contour_with_path(ax4, paths4, labels4, colors4, 'Wykres 4: Metoda najszybszego spadku - różne długości kroku')
plt.tight_layout()
plt.savefig('data/lab4_wykres4.png', dpi=300, bbox_inches='tight')
print("Zapisano: data/lab4_wykres4.png")

# ===== WYKRES 5: Gradienty sprzężone, wszystkie długości kroku =====
fig5, ax5 = plt.subplots(figsize=(10, 8))
paths5 = [
    (df['CG_x1_005'].values, df['CG_x2_005'].values),
    (df['CG_x1_025'].values, df['CG_x2_025'].values),
    (df['CG_x1_var'].values, df['CG_x2_var'].values)
]
labels5 = ['s = 0.05', 's = 0.25', 'zmiennokrokowy']
colors5 = ['magenta', 'red', 'green']
plot_contour_with_path(ax5, paths5, labels5, colors5, 'Wykres 5: Metoda gradientów sprzężonych - różne długości kroku')
plt.tight_layout()
plt.savefig('data/lab4_wykres5.png', dpi=300, bbox_inches='tight')
print("Zapisano: data/lab4_wykres5.png")

# ===== WYKRES 6: Newton, wszystkie długości kroku =====
fig6, ax6 = plt.subplots(figsize=(10, 8))
paths6 = [
    (df['Newton_x1_005'].values, df['Newton_x2_005'].values),
    (df['Newton_x1_025'].values, df['Newton_x2_025'].values),
    (df['Newton_x1_var'].values, df['Newton_x2_var'].values)
]
labels6 = ['s = 0.05', 's = 0.25', 'zmiennokrokowy']
colors6 = ['magenta', 'red', 'green']
plot_contour_with_path(ax6, paths6, labels6, colors6, 'Wykres 6: Metoda Newtona - różne długości kroku')
plt.tight_layout()
plt.savefig('data/lab4_wykres6.png', dpi=300, bbox_inches='tight')
print("Zapisano: data/lab4_wykres6.png")

print("\nWszystkie wykresy zostały wygenerowane!")
print("Pliki PNG zapisane w katalogu data/")
