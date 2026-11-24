# Guia d'usuari

## 1. Requisits

Cal tenir instal·lat Python i les biblioteques següents, que es poden instal·lar amb:

```bash
pip install numpy matplotlib pillow
```

## 2. Ordre d'execució

Els fitxers auxiliars (`_euler_explicit_.py`, `_gauss_seidel_.py` i `solucio_analitica.py`) no s'han d'executar directament, ja que són importats pels altres programes.

Per executar els diferents apartats de la pràctica, cal utilitzar:

1. `Euler_explicit.py`: mètode d'Euler explícit.
2. `Euler_implicit.py`: mètode d'Euler implícit.
3. `Crank_Nicolson.py`: mètode de Crank-Nicolson.
4. `comparacio_errors.py`: comparació dels errors dels tres mètodes.
5. `solucio_problema.py`: resolució del problema d'ablació.
6. `animacions.py`: generació de les animacions d'Euler explícit.

Per exemple, des del terminal:

```bash
python comparacio_errors.py
```

No cal executar prèviament els altres mètodes per fer la comparació, ja que el programa importa automàticament els fitxers necessaris.

## 3. Observacions

- Les animacions es generen en format GIF i s'han guarden a la carpeta Animacions.
- És important mantenir tots els fitxers Python al mateix directori perquè les importacions funcionin correctament.