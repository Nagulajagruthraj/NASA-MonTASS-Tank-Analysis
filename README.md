# OpenMonTASS — Open-source Monocoque Tank Analysis & Sizing System

A lightweight Python tool for thin-shell stress and buckling analysis of rocket propellant tank structures (domes, cylinders, cones).

## Quick Start

```bash
pip install -r requirements.txt
python open_montass.py
```

## What It Does

| Capability | Method |
|---|---|
| **Membrane stress** | Thin-wall pressure-vessel theory (meridional & hoop) |
| **Yield margin** | Von-Mises equivalent stress vs. material yield strength |
| **Buckling load** | NASA SP-8007 knockdown factors for axial compression |
| **Excel export** | One-click `.xlsx` report matching the original MonTASS layout |

## Usage

```python
from open_montass import OpenMonTASS, al_2219

tass = OpenMonTASS()
tass.add_dome(al_2219, R=60, t=0.12)                              # forward dome
tass.add_cone(al_2219, R1=60, R2=45, length=120, t=0.10, alpha=15) # transition cone
results = tass.analyze(P_int=25)                                    # 25 psi internal
print(results)
```

## Pre-defined Materials

| Name | E (psi) | ν | σ_y (psi) |
|---|---|---|---|
| `al_2219` — 2219-T87 Aluminum | 10.5 × 10⁶ | 0.33 | 52 000 |
| `ti_6al4v` — Ti-6Al-4V Titanium | 16.0 × 10⁶ | 0.31 | 130 000 |

Custom materials can be created with the `Material` class:

```python
from open_montass import Material
steel = Material("Steel 304L", 28e6, 0.29, 25000, 9.6e-6, 0.290)
```

## File Overview

| File | Purpose |
|---|---|
| `open_montass.py` | Main analysis module — all classes and demo |
| `requirements.txt` | Python dependencies |
| `README.md` | This file |
