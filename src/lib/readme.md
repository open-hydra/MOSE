# Gestione della Configurazione – Proposta Architetturale

## Logica generale

Il flusso che sto considerando per la gestione della configurazione è il seguente:
             ┌──────────────┐
             │ input_reader │
             └──────┬───────┘
                    │
             ┌──────▼───────┐
             │ config_types │
             └──────┬───────┘
                    │
             ┌──────▼───────┐
             │    setup     │
             └──────┬───────┘
                    │
             ┌──────▼───────┐
             │    solver    │
             └──────────────┘

## Ruolo dei moduli

- Input reader (per ora INI, domani chissà...): legge il file di input e popola le variabili di configurazione (config_types).
- Config types: definisce i tipi di dati per la configurazione (es. input_chemistry_t, input_thermo_t, etc.). Validazione dei dati di configurazione (warnings + errors).
- Setup: utilizza i dati di configurazione per inizializzare il solver. Fa da ponte tra i dati di configurazione e i moduli di fisica/numerica.
- Solver: esegue la simulazione utilizzando i dati di configurazione.

## Dipendenze

Le dipendenze sono:
config_types  → usato da input e setup
input_reader  → usa config_types
setup         → usa config_types e tutti i moduli di fisica e numerica

## Vantaggi rispetto alla struttura precedente

- Le routine di lettura sono completamente separate da setup.
- I moduli di fisica/numerica non si occupano della configurazione.
- La validazione è centralizzata in `config_types`.
- Il codice è più modulare.
- È più semplice passare da INI → YAML (o altro).
- Migliore manutenibilità ed estendibilità futura.

# To-Do

- Definire nomi coerenti e leggibili per:
  - file
  - tipi
  - moduli
- Strutturare le cartelle:
  - abbastanza per chiarire la struttura
  - non troppe per evitare frammentazione
- Definire i tipi di configurazione (`config_types`)
- Implementare validazione (warnings + errors)
- Implementare `input_reader`
- Implementare `setup`
- Disegnare un DAG delle dipendenze



# Bonus - da INI a ...

## 🏆 Le 3 Alternative

### 1️⃣ YAML  → 🔥 migliore per CFD

Esempio:

```yaml
simulation:
  final_time: 1.0
  dt: 1e-3

physics:
  type: navier_stokes
  gamma: 1.4
  viscous: true
  turbulence:
    model: komega
    beta: 0.075

scheme:
  flux: roe
  reconstruction_order: 2
```

✅ Pro

* nesting naturale
* molto leggibile
* ottimo per configurazioni scientifiche
* facile da estendere
* supporta liste e oggetti complessi

❌ Contro

* parsing un po’ più complesso
* richiede libreria (ma ne esistono per Fortran)

👉 Per CFD serio, YAML è spesso la scelta migliore.

---

### 2️⃣ JSON  → più rigido ma robusto

```json
{
  "simulation": {
    "final_time": 1.0,
    "dt": 1e-3
  },
  "physics": {
    "type": "navier_stokes",
    "gamma": 1.4,
    "viscous": true
  }
}
```

✅ Pro

* formalmente ben definito
* facile da validare
* molte librerie
* molto usato in HPC

❌ Contro

* meno leggibile di YAML
* niente commenti (ufficialmente)

👉 Ottimo se vuoi rigore + validazione.

---

### 3️⃣ TOML  → compromesso elegante

```toml
[simulation]
final_time = 1.0
dt = 1e-3

[physics]
type = "navier_stokes"
gamma = 1.4

[physics.turbulence]
model = "komega"
beta = 0.075
```

✅ Pro

* supporta nesting
* più semplice di YAML
* più strutturato di INI
* molto leggibile

❌ Contro

* meno diffuso in HPC

👉 Se ti piace INI ma vuoi nesting → TOML è naturale evoluzione.
