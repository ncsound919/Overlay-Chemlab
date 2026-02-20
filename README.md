# Overlay-ChemLab

**Professional-grade digital chemistry laboratory** — a full-stack scientific platform
built on methodologies from leading research institutions including Pfizer, GSK,
the ACS Green Chemistry Institute, Astex Pharmaceuticals, and NIST.

---

## Features

### Core Science Engine

| Module | Description |
|--------|-------------|
| **Arrhenius Kinetics** | Fischer-Speier esterification ODE with RK4 integration, temperature ramp, catalyst modelling |
| **SMILES Parser** | Full organic-subset tokenizer; molecular formula (Hill system), atom/bond counts, substructure matching |
| **Morgan Fingerprints** | Circular FNV fingerprints with configurable radius and bit-length; Tanimoto and cosine similarity; k-NN search |
| **Bayesian Risk Scorer** | Logistic regression toxicity scorer with bootstrap confidence intervals and batch scoring |
| **Reaction Parser** | Arrow-notation parser; condition extraction (catalyst, temperature, time, solvent); balance checker; reaction-type classifier |

---

### Professional Scientific Modules

#### Drug-Likeness Assessment *(Pfizer · GSK · Astex)*

Industry-standard screening rules used in pharmaceutical hit-identification and lead optimisation:

- **Lipinski Rule of Five** (Pfizer, 1997) — MW ≤ 500, logP ≤ 5, HBD ≤ 5, HBA ≤ 10 with one-violation tolerance
- **Veber Oral Bioavailability Rules** (GSK, 2002) — rotatable bonds ≤ 10, TPSA ≤ 140 Å²
- **Rule of Three for FBDD** (Astex, 2003) — fragment-based drug discovery screening
- **TPSA** — atom-contribution method (Ertl *et al.* 2000, validated: aspirin = 63.6 Å²)
- **Estimated logP** — simplified Wildman–Crippen atom contributions (Wildman & Crippen 1999)
- **Per-atom HBD/HBA** — derived from SMILES connectivity and implicit hydrogen rules
- **Rotatable bonds** — bridge-detection algorithm (Tarjan DFS), matches PubChem/RDKit definition

#### Green Chemistry Metrics *(ACS GCI · Pfizer · CHEM21)*

ACS Green Chemistry Institute Pharmaceutical Roundtable toolkit:

- **Atom Economy** (Trost, 1991) — stoichiometric efficiency: MW(product) / Σ MW(reactants)
- **E-Factor** (Sheldon, 1992) — kg waste per kg product; benchmark: APIs = 25–100
- **Process Mass Intensity (PMI)** (ACS GCI, 2011) — total mass per kg product; PMI = E + 1
- **Reaction Mass Efficiency (RME)** (Curzons *et al.*, 2001) — combined yield × AE × stoichiometric factor
- **Carbon Efficiency** — moles C in product / moles C in feedstock
- **CHEM21 Solvent Guide** (Henderson *et al.*, 2011) — classifies 35+ solvents as *recommended / usable / problematic / hazardous*
- **Composite green score** — weighted A–F grade (AE 30 %, E-factor 25 %, yield 30 %, solvent 15 %)

#### PubChem Integration *(NIH/NLM)*

- Live CAS number, IUPAC name, XLogP3, TPSA, HBD/HBA from the PubChem REST API
- Compound autocomplete for editor suggestions
- LRU cache with configurable TTL; no API key required

---

### Monaco-Style Web IDE

- Syntax-highlighted `.chem` editor with reaction-arrow rendering
- Real-time 2D structure preview via SmilesDrawer
- Kinetics simulation panel with time-series chart
- AI chemistry assistant with topic chips (mechanism, NMR, IR, yield optimisation)
- Molecule property panel

---

## API Reference

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/api/health` | Server status and version |
| GET | `/api/molecule/properties?smiles=` | MW, formula, atom/bond counts |
| GET | `/api/molecule/fingerprint?smiles=` | Morgan fingerprint vector |
| GET | `/api/molecule/similarity?smiles1=&smiles2=` | Tanimoto similarity |
| POST | `/api/molecule/risk` | Bayesian toxicity risk score |
| **GET** | **`/api/molecule/druglikeness?smiles=`** | **Lipinski Ro5, Veber rules, TPSA, logP, HBD/HBA, RotBonds** |
| POST | `/api/simulate/esterification` | Single-point Arrhenius kinetics |
| POST | `/api/simulate/kinetics` | ODE time-series (RK4) |
| POST | `/api/reaction/parse` | Parse reaction string to AST |
| GET | `/api/pubchem/compound?smiles=` | Live PubChem compound data |
| GET | `/api/pubchem/autocomplete?query=` | Compound name autocomplete |
| **POST** | **`/api/green/metrics`** | **Atom Economy, E-factor, PMI, RME, green score** |
| **GET** | **`/api/green/solvent?name=`** | **CHEM21 solvent classification** |
| POST | `/api/chem-assistant` | Chemistry Q&A knowledge base |

---

## Quick Start

```bash
npm install
npm start          # http://localhost:3000
npm test           # 82 tests across 7 suites
```

---

## Scientific References

The modules in this project implement peer-reviewed methodologies:

| Reference | Module |
|-----------|--------|
| Lipinski *et al.* *Adv. Drug Deliv. Rev.* **23**, 3–25 (1997) | Lipinski Rule of Five |
| Veber *et al.* *J. Med. Chem.* **45**, 2615–2623 (2002) | Oral bioavailability rules |
| Congreve *et al.* *Drug Discov. Today* **8**, 876–877 (2003) | Rule of Three (FBDD) |
| Ertl *et al.* *J. Med. Chem.* **43**, 3714–3717 (2000) | TPSA atom contributions |
| Wildman & Crippen *J. Chem. Inf. Comput. Sci.* **39**, 868–873 (1999) | logP estimation |
| Trost *Science* **254**, 1471–1477 (1991) | Atom Economy |
| Sheldon *Chem. Ind.* (London) 903–906 (1992) | E-Factor |
| Curzons *et al.* *Green Chem.* **3**, 7–9 (2001) | Reaction Mass Efficiency |
| Jiménez-González *et al.* *Org. Process Res. Dev.* **15**, 912–917 (2011) | PMI benchmarks |
| Henderson *et al.* *Green Chem.* **13**, 854–862 (2011) | CHEM21 Solvent Guide |

---

## License

MIT
