#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Wizualizacja wyników klasyfikacji logistycznej - Lab 4 Problem Rzeczywisty
Wykres przedstawia dane uczące (przyjęci/nieprzyjęci) oraz granicę klasyfikacji
"""

import numpy as np
import matplotlib.pyplot as plt
import csv

# Wczytanie danych
data = []
with open('data/lab4_classification_data.csv', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        data.append(row)

# Odczytanie parametrów najlepszego klasyfikatora (zapisane tylko w pierwszym wierszu)
theta0 = float(data[0]['best_theta0'])
theta1 = float(data[0]['best_theta1'])
theta2 = float(data[0]['best_theta2'])

print(f"Parametry najlepszego klasyfikatora:")
print(f"  θ₀ = {theta0:.6f}")
print(f"  θ₁ = {theta1:.6f}")
print(f"  θ₂ = {theta2:.6f}")

# Podział danych na przyjętych (y=1) i nieprzyjętych (y=0)
x1_przyjeci = []
x2_przyjeci = []
x1_nieprzyjeci = []
x2_nieprzyjeci = []

for row in data:
    x1 = float(row['x1'])
    x2 = float(row['x2'])
    y = float(row['y'])
    
    if y == 1:
        x1_przyjeci.append(x1)
        x2_przyjeci.append(x2)
    else:
        x1_nieprzyjeci.append(x1)
        x2_nieprzyjeci.append(x2)

# Tworzenie wykresu
plt.figure(figsize=(10, 8))

# Punkty danych
plt.scatter(x1_przyjeci, x2_przyjeci, 
           c='green', marker='o', s=100, alpha=0.7, 
           edgecolors='darkgreen', linewidth=1.5,
           label='Przyjęci (y=1)')

plt.scatter(x1_nieprzyjeci, x2_nieprzyjeci, 
           c='red', marker='x', s=100, alpha=0.7, 
           linewidth=2,
           label='Nieprzyjęci (y=0)')

# Granica klasyfikacji: h_theta(x) = 0.5
# Dla hipotezy sigmoidalnej: h_theta(x) = 0.5 ⟺ θ₀ + θ₁·x₁ + θ₂·x₂ = 0
# Stąd: x₂ = -(θ₀ + θ₁·x₁) / θ₂

if abs(theta2) > 1e-10:  # Sprawdzenie czy theta2 != 0
    # Zakres wartości x1
    all_x1 = [float(row['x1']) for row in data]
    all_x2 = [float(row['x2']) for row in data]
    
    x1_min = min(all_x1) - 5
    x1_max = max(all_x1) + 5
    x1_line = np.linspace(x1_min, x1_max, 200)
    
    # Obliczenie x2 z równania granicy
    x2_line = -(theta0 + theta1 * x1_line) / theta2
    
    # Rysowanie granicy klasyfikacji
    plt.plot(x1_line, x2_line, 'b-', linewidth=2.5, 
            label=f'Granica klasyfikacji\n(θ₀={theta0:.3f}, θ₁={theta1:.3f}, θ₂={theta2:.3f})')
    
    # Zacieniowanie obszarów
    y_min = min(all_x2) - 5
    y_max = max(all_x2) + 5
    
    # Obszar przyjęcia (h >= 0.5)
    plt.fill_between(x1_line, x2_line, y_max, alpha=0.1, color='green', 
                     label='Obszar przyjęcia (h≥0.5)')
    
    # Obszar odrzucenia (h < 0.5)
    plt.fill_between(x1_line, x2_line, y_min, alpha=0.1, color='red',
                     label='Obszar odrzucenia (h<0.5)')
    
    # Dostosowanie zakresów osi
    margin = 5
    plt.xlim(min(all_x1) - margin, max(all_x1) + margin)
    plt.ylim(min(all_x2) - margin, max(all_x2) + margin)
else:
    print("UWAGA: θ₂ ≈ 0, nie można narysować granicy w postaci x₂(x₁)")

# Opis wykresu
plt.xlabel('Ocena z przedmiotu 1 (x₁)', fontsize=12)
plt.ylabel('Ocena z przedmiotu 2 (x₂)', fontsize=12)
plt.title('Klasyfikator logistyczny - Przyjęcie na uczelnię\nGranica decyzyjna dla najlepszego modelu', 
         fontsize=14, fontweight='bold')
plt.legend(loc='best', fontsize=10)
plt.grid(True, alpha=0.3, linestyle='--')

plt.tight_layout()

# Zapis wykresu
plt.savefig('data/lab4_classification_plot.png', dpi=300, bbox_inches='tight')
print("\nWykres zapisany do: data/lab4_classification_plot.png")

plt.show()

# Wyświetlenie statystyk
print(f"\nStatystyki danych:")
print(f"  Liczba przyjętych: {len(x1_przyjeci)}")
print(f"  Liczba nieprzyjętych: {len(x1_nieprzyjeci)}")
print(f"  Łącznie przykładów: {len(data)}")
