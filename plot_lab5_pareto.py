import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Wczytaj dane z CSV
# Kolumny: l0[mm], d0[mm], l*[mm], d*[mm], masa*[kg], ugięcie*[m], naprężenie*[MPa], f_calls
data = pd.read_csv('data/lab5_real.csv', header=None, 
                   names=['l0_mm', 'd0_mm', 'l_star_mm', 'd_star_mm', 
                          'masa_kg', 'ugiecie_m', 'naprezenie_MPa', 'f_calls'])

print(f"Wczytano {len(data)} wierszy danych\n")

# Filtruj nieprawidłowe wyniki - BARDZO restrykcyjne filtry
# l* ∈ [200, 1000] mm, d* ∈ [10, 50] mm
# u_max = 2.5 mm, sigma_max = 300 MPa
# Wszystkie wartości muszą być dodatnie i w sensownych zakresach
valid_mask = (
    (data['l_star_mm'] >= 200) & (data['l_star_mm'] <= 1000) &
    (data['d_star_mm'] >= 10) & (data['d_star_mm'] <= 50) &
    (data['ugiecie_m'] > 0) & (data['ugiecie_m'] <= 0.1) &  # ugięcie w [0, 0.1] m = [0, 100] mm
    (data['masa_kg'] > 0) & (data['masa_kg'] <= 100) &  # masa w sensownym zakresie
    (data['naprezenie_MPa'] > 0) & (data['naprezenie_MPa'] <= 1000)  # naprężenie dodatnie i rozsądne
)

data_valid = data[valid_mask].copy()
print(f"Prawidłowe rozwiązania: {len(data_valid)} / {len(data)}")

if len(data_valid) == 0:
    print("\nBRAK PRAWIDŁOWYCH ROZWIĄZAŃ!")
    print("Przykładowe nieprawidłowe dane:")
    print(data.head(10))
    exit(1)

# Usuń wartości odstające używając metody IQR (Interquartile Range)
def remove_outliers_iqr(df, columns, multiplier=3.0):
    """
    Usuwa wartości odstające używając metody IQR
    multiplier: im większy, tym mniej restrykcyjne (1.5 = standardowe, 3.0 = tylko ekstremalne)
    """
    df_filtered = df.copy()
    for col in columns:
        Q1 = df_filtered[col].quantile(0.25)
        Q3 = df_filtered[col].quantile(0.75)
        IQR = Q3 - Q1
        lower_bound = Q1 - multiplier * IQR
        upper_bound = Q3 + multiplier * IQR
        
        before_count = len(df_filtered)
        df_filtered = df_filtered[(df_filtered[col] >= lower_bound) & (df_filtered[col] <= upper_bound)]
        removed = before_count - len(df_filtered)
        if removed > 0:
            print(f"  {col}: usunięto {removed} wartości odstających (zakres: [{lower_bound:.3f}, {upper_bound:.3f}])")
    
    return df_filtered

print("\nUsuwanie wartości odstających (metoda IQR, multiplier=3.0):")
data_valid = remove_outliers_iqr(data_valid, ['masa_kg', 'ugiecie_m', 'naprezenie_MPa'], multiplier=3.0)
print(f"Po usunięciu wartości odstających: {len(data_valid)} rozwiązań")

if len(data_valid) == 0:
    print("\nBRAK ROZWIĄZAŃ PO FILTROWANIU!")
    exit(1)

# Funkcje celu to masa i ugięcie
f1 = data_valid['masa_kg'].values  # masa [kg]
f2 = data_valid['ugiecie_m'].values * 1000  # ugięcie [mm] - konwersja dla lepszej czytelności

# Znajdź rozwiązania Pareto-optymalne
def is_pareto_efficient(costs):
    """
    Znajduje punkty Pareto-optymalne
    costs: array (n_points, n_costs) - niższe wartości są lepsze
    """
    is_efficient = np.ones(costs.shape[0], dtype=bool)
    for i, c in enumerate(costs):
        if is_efficient[i]:
            # Punkt i jest zdominowany jeśli istnieje punkt lepszy w obu kryteriach
            is_efficient[is_efficient] = np.any(costs[is_efficient] < c, axis=1)
            is_efficient[i] = True
    return is_efficient

# Stwórz tablicę z funkcjami celu
costs = np.column_stack([f1, f2])

# Znajdź punkty Pareto-optymalne
pareto_mask = is_pareto_efficient(costs)

# Wyodrębnij punkty Pareto i nie-Pareto
f1_pareto = f1[pareto_mask]
f2_pareto = f2[pareto_mask]
f1_non_pareto = f1[~pareto_mask]
f2_non_pareto = f2[~pareto_mask]

# Posortuj punkty Pareto dla lepszej wizualizacji frontu
pareto_points = np.column_stack([f1_pareto, f2_pareto])
pareto_sorted = pareto_points[np.argsort(pareto_points[:, 0])]

# Tworzenie wykresu
plt.figure(figsize=(12, 8))

# Wszystkie rozwiązania (szare punkty)
if len(f1_non_pareto) > 0:
    plt.scatter(f1_non_pareto, f2_non_pareto, c='lightgray', alpha=0.6, 
                s=50, label='Rozwiązania nieoptymaline', zorder=1)

# Rozwiązania Pareto-optymalne (czerwone punkty)
plt.scatter(f1_pareto, f2_pareto, c='red', alpha=0.9, 
            s=100, label='Rozwiązania Pareto-optymalne', 
            edgecolors='darkred', linewidth=2, zorder=3)

# Linia łącząca front Pareto
plt.plot(pareto_sorted[:, 0], pareto_sorted[:, 1], 'r--', 
         alpha=0.7, linewidth=2.5, label='Front Pareto', zorder=2)

plt.xlabel('$f_1$ - Masa belki [kg]', fontsize=14, fontweight='bold')
plt.ylabel('$f_2$ - Ugięcie belki [mm]', fontsize=14, fontweight='bold')
plt.title('Wykres rozwiązań Pareto-optymalnych\n(Optymalizacja wielokryterialna belki - Lab 5)', 
          fontsize=16, fontweight='bold', pad=20)
plt.legend(fontsize=11, loc='best')
plt.grid(True, alpha=0.4, linestyle='--')
plt.tight_layout()

# Zapisz wykres
plt.savefig('data/lab5_pareto_front.png', dpi=300, bbox_inches='tight')
print(f"Wykres zapisany do: data/lab5_pareto_front.png")

# Wyświetl wykres
plt.show()

# Wypisz statystyki
print(f"\nStatystyki:")
print(f"Liczba prawidłowych rozwiązań: {len(f1)}")
print(f"Liczba rozwiązań Pareto-optymalnych: {np.sum(pareto_mask)}")
print(f"\nZakres f1 (masa): [{f1.min():.3f}, {f1.max():.3f}] kg")
print(f"Zakres f2 (ugięcie): [{f2.min():.3f}, {f2.max():.3f}] mm")

# Sprawdź ograniczenia
u_violations = np.sum(data_valid['ugiecie_m'] > 0.0025)
sigma_violations = np.sum(data_valid['naprezenie_MPa'] > 300)
print(f"\nNaruszenia ograniczeń:")
print(f"Ugięcie > 2.5 mm: {u_violations} rozwiązań")
print(f"Naprężenie > 300 MPa: {sigma_violations} rozwiązań")

# Wypisz kilka przykładowych rozwiązań Pareto
print(f"\nPrzykładowe rozwiązania Pareto-optymalne:")
print(f"{'Masa [kg]':<12} {'Ugięcie [mm]':<15} {'Naprężenie [MPa]':<20} {'l* [mm]':<12} {'d* [mm]':<12}")
print("-" * 80)
pareto_indices = np.where(pareto_mask)[0]
for idx in pareto_indices[:5]:  # Pierwsze 5 rozwiązań
    valid_idx = data_valid.index[idx]
    print(f"{data_valid['masa_kg'].iloc[idx]:<12.3f} {data_valid['ugiecie_m'].iloc[idx]*1000:<15.3f} {data_valid['naprezenie_MPa'].iloc[idx]:<20.1f} {data_valid['l_star_mm'].iloc[idx]:<12.1f} {data_valid['d_star_mm'].iloc[idx]:<12.1f}")

if len(pareto_indices) > 5:
    print("...")
    for idx in pareto_indices[-2:]:  # Ostatnie 2 rozwiązania
        valid_idx = data_valid.index[idx]
        print(f"{data_valid['masa_kg'].iloc[idx]:<12.3f} {data_valid['ugiecie_m'].iloc[idx]*1000:<15.3f} {data_valid['naprezenie_MPa'].iloc[idx]:<20.1f} {data_valid['l_star_mm'].iloc[idx]:<12.1f} {data_valid['d_star_mm'].iloc[idx]:<12.1f}")
