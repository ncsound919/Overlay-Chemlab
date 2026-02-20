'use strict';

/**
 * Drug-likeness assessment module.
 *
 * Implements industry-standard screening rules used in pharmaceutical
 * hit identification and lead optimisation:
 *
 *   - Lipinski Rule of Five  (Pfizer, 1997)
 *     Lipinski et al. Adv. Drug Deliv. Rev. 23, 3–25 (1997)
 *
 *   - Veber Oral Bioavailability Rules  (GSK, 2002)
 *     Veber et al. J. Med. Chem. 45, 2615–2623 (2002)
 *
 *   - Rule of Three for FBDD  (Astex Pharmaceuticals, 2003)
 *     Congreve et al. Drug Discov. Today 8, 876–877 (2003)
 *
 *   - TPSA atom contributions
 *     Ertl et al. J. Med. Chem. 43, 3714–3717 (2000)
 *
 *   - logP atom contributions (simplified Wildman–Crippen)
 *     Wildman & Crippen J. Chem. Inf. Comput. Sci. 39, 868–873 (1999)
 *     NOTE: accuracy ±1–2 log units; use PubChem XLogP3 for precise values.
 */

const { parseSMILES, molecularWeight } = require('./molecular-embeddings.js');

// ─────────────────────────────────────────────────────────────
// Internal helpers
// ─────────────────────────────────────────────────────────────

const VALENCES = { C: 4, N: 3, O: 2, S: 2, P: 3, F: 1, Cl: 1, Br: 1, I: 1 };

/**
 * Compute implicit hydrogen count for a heavy atom.
 */
function _implicitH(symbol, aromatic, bondOrderSum) {
  const v = VALENCES[symbol];
  if (v === undefined) return 0;
  return Math.max(0, (aromatic ? v - 1 : v) - bondOrderSum);
}

/**
 * Build per-atom descriptor array: {symbol, aromatic, bondOrderSum, degree, implicitH}.
 * @param {string} smiles
 * @returns {{ atoms: object[], bonds: object[] }}
 */
function _atomDetails(smiles) {
  const parsed = parseSMILES(smiles);
  const bondSum = new Array(parsed.atoms.length).fill(0);
  const degree  = new Array(parsed.atoms.length).fill(0);

  for (const b of parsed.bonds) {
    bondSum[b.from] += b.order;
    bondSum[b.to]   += b.order;
    degree[b.from]++;
    degree[b.to]++;
  }

  return {
    atoms: parsed.atoms.map((a, i) => ({
      symbol:       a.symbol,
      aromatic:     a.aromatic,
      bondOrderSum: bondSum[i],
      degree:       degree[i],
      implicitH:    _implicitH(a.symbol, a.aromatic, bondSum[i]),
    })),
    bonds: parsed.bonds,
  };
}

// ─────────────────────────────────────────────────────────────
// H-bond Donors and Acceptors
// ─────────────────────────────────────────────────────────────

/**
 * Count H-bond donors: O–H and N–H groups.
 * PubChem/Lipinski definition: number of NH and OH groups (not individual H atoms).
 *
 * @param {string} smiles
 * @returns {number}
 */
function countHBD(smiles) {
  const { atoms } = _atomDetails(smiles);
  return atoms.filter(a =>
    (a.symbol === 'O' || a.symbol === 'N') && a.implicitH > 0
  ).length;
}

/**
 * Count H-bond acceptors: all N and O heavy atoms.
 * Uses the simplified Lipinski/PubChem definition.
 *
 * @param {string} smiles
 * @returns {number}
 */
function countHBA(smiles) {
  const parsed = parseSMILES(smiles);
  return parsed.atoms.filter(a => a.symbol === 'O' || a.symbol === 'N').length;
}

// ─────────────────────────────────────────────────────────────
// Rotatable Bonds
// ─────────────────────────────────────────────────────────────

/**
 * Find the set of bond indices that are part of at least one ring.
 * Uses Tarjan bridge-finding (DFS) — O(V+E).
 * A bond is a ring bond iff it is NOT a bridge (i.e., its removal
 * leaves the graph connected).
 *
 * @param {object[]} bonds  - Parsed bond array {from, to, order}
 * @param {number}   n      - Number of atoms
 * @returns {Set<number>} Indices of ring bonds
 */
function _findRingBonds(bonds, n) {
  // Adjacency list: each entry is {to, idx}
  const adj = Array.from({ length: n }, () => []);
  for (let i = 0; i < bonds.length; i++) {
    adj[bonds[i].from].push({ to: bonds[i].to,   idx: i });
    adj[bonds[i].to  ].push({ to: bonds[i].from, idx: i });
  }

  const bridges = new Set();
  const visited = new Array(n).fill(false);
  const disc    = new Array(n).fill(-1);
  const low     = new Array(n).fill(-1);
  let timer = 0;

  function dfs(u, parentEdgeIdx) {
    visited[u] = true;
    disc[u] = low[u] = timer++;
    for (const { to: v, idx } of adj[u]) {
      if (idx === parentEdgeIdx) continue;
      if (!visited[v]) {
        dfs(v, idx);
        low[u] = Math.min(low[u], low[v]);
        if (low[v] > disc[u]) bridges.add(idx);  // bridge condition
      } else {
        low[u] = Math.min(low[u], disc[v]);
      }
    }
  }

  for (let i = 0; i < n; i++) {
    if (!visited[i]) dfs(i, -1);
  }

  // Ring bonds = bonds that are NOT bridges
  const ringBonds = new Set();
  for (let i = 0; i < bonds.length; i++) {
    if (!bridges.has(i)) ringBonds.add(i);
  }
  return ringBonds;
}

/**
 * Count rotatable single bonds.
 *
 * A bond is counted as rotatable if it:
 *   - is a single bond between two heavy atoms
 *   - is not in a ring
 *   - is not terminal (both endpoints have heavy-atom degree > 1)
 *
 * This matches the PubChem/RDKit rotatable bond definition
 * (amide and ester bonds are included).
 *
 * @param {string} smiles
 * @returns {number}
 */
function countRotatableBonds(smiles) {
  const parsed = parseSMILES(smiles);
  const { bonds, atoms } = parsed;

  const degree = new Array(atoms.length).fill(0);
  for (const b of bonds) {
    degree[b.from]++;
    degree[b.to]++;
  }

  const ringBonds = _findRingBonds(bonds, atoms.length);

  let count = 0;
  for (let i = 0; i < bonds.length; i++) {
    const b = bonds[i];
    if (b.order !== 1)           continue;  // single bonds only
    if (ringBonds.has(i))        continue;  // skip ring bonds
    if (degree[b.from] <= 1 || degree[b.to] <= 1) continue; // skip terminal
    count++;
  }
  return count;
}

// ─────────────────────────────────────────────────────────────
// logP Estimation (simplified Wildman–Crippen)
// ─────────────────────────────────────────────────────────────

/**
 * Estimate logP using simplified atom-type contributions.
 *
 * Atom types are determined from SMILES connectivity:
 *   - Aromatic C       +0.355
 *   - sp3 C (all single bonds)  +0.175
 *   - sp2 C (double bond or high bond-order sum)  +0.050
 *   - O–H (hydroxyl, carboxyl)  –0.590
 *   - C=O oxygen (carbonyl)     –0.520
 *   - Ether / ester O           –0.290
 *   - N (NH₂)   –0.960  |  N (NH)  –0.450  |  N–  –0.080
 *   - Aromatic N  –0.550
 *   - S, halogens, P  — published Crippen values
 *
 * Accuracy: ±1–2 log units for typical drug-like molecules.
 * For more precise XLogP3 values use GET /api/pubchem/compound?smiles=
 *
 * @param {string} smiles
 * @returns {number} Estimated logP
 */
function estimateLogP(smiles) {
  const { atoms, bonds } = parseSMILES(smiles);
  const bondSum = new Array(atoms.length).fill(0);
  const degree  = new Array(atoms.length).fill(0);

  for (const b of bonds) {
    bondSum[b.from] += b.order;
    bondSum[b.to]   += b.order;
    degree[b.from]++;
    degree[b.to]++;
  }

  let logP = 0;

  for (let i = 0; i < atoms.length; i++) {
    const a  = atoms[i];
    const bs = bondSum[i];
    const deg = degree[i];
    const h = _implicitH(a.symbol, a.aromatic, bs);

    switch (a.symbol) {
      case 'C':
        if (a.aromatic)                logP += 0.355;  // aromatic
        else if (bs >= 4 || deg < bs)  logP += 0.050;  // sp2 (double bond)
        else                           logP += 0.175;  // sp3
        break;

      case 'O':
        if (h > 0)       logP += -0.590;  // –OH (alcohol / carboxyl)
        else if (deg <= 1) logP += -0.520; // C=O (carbonyl oxygen)
        else               logP += -0.290; // ether / ester oxygen
        break;

      case 'N':
        if (a.aromatic)    logP += -0.550;
        else if (h >= 2)   logP += -0.960;  // –NH₂
        else if (h === 1)  logP += -0.450;  // –NH–
        else               logP += -0.080;  // –N< (tertiary)
        break;

      case 'S':
        logP += (h > 0) ? 0.148 : 0.347;
        break;

      case 'P':  logP += -0.228; break;
      case 'F':  logP +=  0.375; break;
      case 'Cl': logP +=  0.530; break;
      case 'Br': logP +=  0.876; break;
      case 'I':  logP +=  1.296; break;
      default:   break;
    }
  }

  return Math.round(logP * 100) / 100;
}

// ─────────────────────────────────────────────────────────────
// TPSA Estimation  (Ertl 2000)
// ─────────────────────────────────────────────────────────────

/**
 * Estimate Topological Polar Surface Area (TPSA) in Å².
 *
 * Uses atom-surface contributions from:
 *   Ertl et al. J. Med. Chem. 43, 3714–3717 (2000).
 *
 * Atom type assignment is based on element, aromaticity, degree,
 * and implicit hydrogen count derived from SMILES connectivity.
 *
 * Accuracy: typically within ±5 Å² for drug-like molecules.
 * Validated: aspirin (CC(=O)Oc1ccccc1C(=O)O) → 63.6 Å² (PubChem: 63.6 Å²).
 *
 * @param {string} smiles
 * @returns {number} TPSA in Å²
 */
function estimateTPSA(smiles) {
  const { atoms, bonds } = parseSMILES(smiles);
  const bondSum = new Array(atoms.length).fill(0);
  const degree  = new Array(atoms.length).fill(0);

  for (const b of bonds) {
    bondSum[b.from] += b.order;
    bondSum[b.to]   += b.order;
    degree[b.from]++;
    degree[b.to]++;
  }

  let tpsa = 0;

  for (let i = 0; i < atoms.length; i++) {
    const a   = atoms[i];
    const bs  = bondSum[i];
    const deg = degree[i];
    const h   = _implicitH(a.symbol, a.aromatic, bs);

    if (a.symbol === 'O') {
      if (h > 0)            tpsa += 20.23;  // –OH (alcohol, phenol, carboxyl OH)
      else if (deg <= 1)    tpsa += 17.07;  // C=O (carbonyl oxygen)
      else if (a.aromatic)  tpsa += 13.14;  // aromatic O (furan-type)
      else                  tpsa +=  9.23;  // ether / ester O

    } else if (a.symbol === 'N') {
      if (a.aromatic && h > 0)  tpsa += 15.79;  // pyrrole-like NH
      else if (a.aromatic)      tpsa += 12.89;  // pyridine-like N
      else if (h >= 2)          tpsa += 26.02;  // –NH₂
      else if (h === 1)         tpsa += 17.07;  // –NH–
      else                      tpsa +=  3.24;  // tertiary –N<

    } else if (a.symbol === 'S') {
      tpsa += (h > 0) ? 38.80 : 25.30;

    } else if (a.symbol === 'P') {
      tpsa += 40.74;
    }
  }

  return Math.round(tpsa * 10) / 10;
}

// ─────────────────────────────────────────────────────────────
// Drug-likeness rules
// ─────────────────────────────────────────────────────────────

/**
 * Evaluate Lipinski Rule of Five (Pfizer, 1997).
 *
 * Predicts oral drug absorption / permeability for drug-like small molecules:
 *   1. MW  ≤ 500 Da
 *   2. logP ≤ 5
 *   3. HBD ≤ 5
 *   4. HBA ≤ 10
 *
 * One violation is acceptable (the original paper tolerates one failure,
 * e.g. for substrates of active transporters or efflux pumps).
 *
 * @param {string} smiles
 * @returns {{ pass: boolean, violations: number, criteria: object[], mw, logP, hbd, hba }}
 */
function lipinskiRuleOfFive(smiles) {
  const mw   = molecularWeight(smiles);
  const logP = estimateLogP(smiles);
  const hbd  = countHBD(smiles);
  const hba  = countHBA(smiles);

  const criteria = [
    { name: 'MW ≤ 500',  value: mw,   threshold: 500, unit: 'Da', pass: mw   <= 500 },
    { name: 'logP ≤ 5',  value: logP, threshold:   5, unit: '',   pass: logP <=   5 },
    { name: 'HBD ≤ 5',   value: hbd,  threshold:   5, unit: '',   pass: hbd  <=   5 },
    { name: 'HBA ≤ 10',  value: hba,  threshold:  10, unit: '',   pass: hba  <=  10 },
  ];

  const violations = criteria.filter(c => !c.pass).length;

  return {
    pass: violations <= 1,
    violations,
    criteria,
    mw,
    logP,
    hbd,
    hba,
  };
}

/**
 * Evaluate Veber rules for oral bioavailability (GSK, 2002).
 *
 * Derived from a retrospective study of 1100 drug candidates in rats.
 * Compounds satisfying both criteria show good oral bioavailability:
 *   1. Rotatable bonds ≤ 10
 *   2. TPSA ≤ 140 Å²  OR  (HBD + HBA) ≤ 12
 *
 * @param {string} smiles
 * @returns {{ pass: boolean, rotBonds, rotPass, tpsa, tpsaPass, hbSum, hbSumPass }}
 */
function veberRules(smiles) {
  const rotBonds = countRotatableBonds(smiles);
  const tpsa     = estimateTPSA(smiles);
  const hbd      = countHBD(smiles);
  const hba      = countHBA(smiles);

  const rotPass   = rotBonds <= 10;
  const tpsaPass  = tpsa     <= 140;
  const hbSumPass = (hbd + hba) <= 12;

  return {
    pass: rotPass && (tpsaPass || hbSumPass),
    rotBonds, rotPass,
    tpsa,     tpsaPass,
    hbSum: hbd + hba, hbSumPass,
  };
}

/**
 * Evaluate Rule of Three for fragment-based drug discovery (Astex, 2003).
 *
 * Fragments are lower-MW starting points for structure-based drug design.
 * One violation is tolerated:
 *   1. MW  ≤ 300 Da
 *   2. logP ≤ 3
 *   3. HBD ≤ 3
 *   4. HBA ≤ 3
 *
 * @param {string} smiles
 * @returns {{ pass: boolean, violations: number, criteria: object[] }}
 */
function ruleOfThree(smiles) {
  const mw   = molecularWeight(smiles);
  const logP = estimateLogP(smiles);
  const hbd  = countHBD(smiles);
  const hba  = countHBA(smiles);

  const criteria = [
    { name: 'MW ≤ 300',  value: mw,   threshold: 300, pass: mw   <= 300 },
    { name: 'logP ≤ 3',  value: logP, threshold:   3, pass: logP <=   3 },
    { name: 'HBD ≤ 3',   value: hbd,  threshold:   3, pass: hbd  <=   3 },
    { name: 'HBA ≤ 3',   value: hba,  threshold:   3, pass: hba  <=   3 },
  ];

  const violations = criteria.filter(c => !c.pass).length;
  return { pass: violations <= 1, violations, criteria };
}

/**
 * Comprehensive drug-likeness assessment report.
 * Combines all screening rules and descriptors into a single object.
 *
 * @param {string} smiles
 * @returns {object} Full drug-likeness report
 */
function assessDruglikeness(smiles) {
  const mw       = molecularWeight(smiles);
  const logP     = estimateLogP(smiles);
  const tpsa     = estimateTPSA(smiles);
  const hbd      = countHBD(smiles);
  const hba      = countHBA(smiles);
  const rotBonds = countRotatableBonds(smiles);

  const lipinski = lipinskiRuleOfFive(smiles);
  const veber    = veberRules(smiles);
  const ro3      = ruleOfThree(smiles);

  return {
    smiles,
    descriptors: { mw, logP, tpsa, hbd, hba, rotBonds },
    lipinski,
    veber,
    ruleOfThree: ro3,
    overallDruglike: lipinski.pass && veber.pass,
  };
}

module.exports = {
  countHBD,
  countHBA,
  countRotatableBonds,
  estimateLogP,
  estimateTPSA,
  lipinskiRuleOfFive,
  veberRules,
  ruleOfThree,
  assessDruglikeness,
};
