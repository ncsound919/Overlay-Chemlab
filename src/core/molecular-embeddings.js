'use strict';

const ATOMIC_WEIGHTS = {
  H: 1.008, C: 12.011, N: 14.007, O: 15.999,
  S: 32.065, P: 30.974, F: 18.998, Cl: 35.453,
  Br: 79.904, I: 126.904,
};

const AROMATIC_SET = new Set(['c', 'n', 'o', 's', 'p']);
const TWO_LETTER = new Set(['Cl', 'Br']);

/**
 * Simple deterministic hash of a string to a 32-bit integer.
 */
function hashString(str) {
  let h = 0x811c9dc5; // FNV offset basis
  for (let i = 0; i < str.length; i++) {
    h ^= str.charCodeAt(i);
    h = Math.imul(h, 0x01000193); // FNV prime
  }
  return h >>> 0; // unsigned
}

/**
 * Parse a SMILES string into an atom/bond graph.
 *
 * @param {string} smiles
 * @returns {{ atoms: Array<{symbol: string, index: number, aromatic: boolean}>,
 *             bonds: Array<{from: number, to: number, order: number}>,
 *             rings: number }}
 */
function parseSMILES(smiles) {
  const atoms = [];
  const bonds = [];
  const branchStack = [];
  const ringOpenings = new Map();
  let prevAtomIdx = -1;
  let pendingBondOrder = null;
  let ringCount = 0;

  let i = 0;
  while (i < smiles.length) {
    const ch = smiles[i];

    // Branch
    if (ch === '(') {
      branchStack.push(prevAtomIdx);
      i++;
      continue;
    }
    if (ch === ')') {
      prevAtomIdx = branchStack.pop();
      i++;
      continue;
    }

    // Bonds
    if (ch === '-') { pendingBondOrder = 1; i++; continue; }
    if (ch === '=') { pendingBondOrder = 2; i++; continue; }
    if (ch === '#') { pendingBondOrder = 3; i++; continue; }
    if (ch === ':') { pendingBondOrder = 1; i++; continue; } // aromatic

    // Stereochemistry markers — skip
    if (ch === '/' || ch === '\\') { i++; continue; }

    // Ring closure
    if (ch >= '0' && ch <= '9') {
      const digit = parseInt(ch, 10);
      handleRing(digit);
      i++;
      continue;
    }
    if (ch === '%' && i + 2 < smiles.length) {
      const num = parseInt(smiles.substring(i + 1, i + 3), 10);
      handleRing(num);
      i += 3;
      continue;
    }

    // Bracket atom
    if (ch === '[') {
      const close = smiles.indexOf(']', i);
      const inner = smiles.substring(i + 1, close === -1 ? smiles.length : close);
      const element = parseBracketAtom(inner);
      addAtom(element.symbol, element.aromatic);
      i = (close === -1 ? smiles.length : close + 1);
      continue;
    }

    // Two-letter organic atoms
    if (i + 1 < smiles.length && TWO_LETTER.has(smiles.substring(i, i + 2))) {
      addAtom(smiles.substring(i, i + 2), false);
      i += 2;
      continue;
    }

    // Single-letter atoms
    if (/[A-Za-z]/.test(ch)) {
      const aromatic = AROMATIC_SET.has(ch);
      const symbol = aromatic ? ch.toUpperCase() : ch;
      addAtom(symbol, aromatic);
      i++;
      continue;
    }

    // Dot (fragment separator) and other characters
    if (ch === '.') { prevAtomIdx = -1; i++; continue; }
    i++;
  }

  return { atoms, bonds, rings: ringCount };

  function addAtom(symbol, aromatic) {
    const idx = atoms.length;
    atoms.push({ symbol, index: idx, aromatic });

    if (prevAtomIdx >= 0) {
      const order = pendingBondOrder !== null ? pendingBondOrder
        : (aromatic && atoms[prevAtomIdx].aromatic ? 1 : 1);
      bonds.push({ from: prevAtomIdx, to: idx, order });
    }
    pendingBondOrder = null;
    prevAtomIdx = idx;
  }

  function handleRing(digit) {
    if (ringOpenings.has(digit)) {
      const openIdx = ringOpenings.get(digit);
      ringOpenings.delete(digit);
      const order = pendingBondOrder !== null ? pendingBondOrder : 1;
      bonds.push({ from: openIdx, to: prevAtomIdx, order });
      pendingBondOrder = null;
      ringCount++;
    } else {
      ringOpenings.set(digit, prevAtomIdx);
    }
  }

  function parseBracketAtom(inner) {
    let s = inner.replace(/^[0-9]+/, ''); // isotope
    s = s.replace(/@+/g, '');             // chirality
    s = s.replace(/[+\-][0-9]*/g, '');    // charge
    s = s.replace(/H[0-9]*/g, '');        // H count

    const aromatic = s.length >= 1 && s[0] === s[0].toLowerCase() && /[a-z]/.test(s[0]);

    // Extract element
    let symbol;
    if (s.length >= 2 && ATOMIC_WEIGHTS[s[0].toUpperCase() + s[1].toLowerCase()]) {
      symbol = s[0].toUpperCase() + s[1].toLowerCase();
    } else if (s.length >= 1) {
      symbol = s[0].toUpperCase();
    } else {
      symbol = 'H';
    }
    return { symbol, aromatic };
  }
}

/**
 * Build adjacency list from parsed SMILES graph.
 */
function buildAdjacency(parsed) {
  const adj = new Array(parsed.atoms.length).fill(null).map(() => []);
  for (const bond of parsed.bonds) {
    adj[bond.from].push({ neighbor: bond.to, order: bond.order });
    adj[bond.to].push({ neighbor: bond.from, order: bond.order });
  }
  return adj;
}

/**
 * Generate Morgan-style circular fingerprint.
 *
 * Algorithm: For each atom, at each radius level (0 to radius), collect the
 * sorted list of neighbour environment hashes, concatenate with the atom's own
 * hash, re-hash, and set the corresponding bit in the fingerprint.
 *
 * NOTE: Uses FNV-1a (32-bit) hash mapped via modulo to fingerprint bits.
 * The default 128-bit length may produce significant hash collisions for
 * large molecule sets. For production similarity searches, use nBits=2048
 * or higher to reduce collision rates.
 *
 * @param {string} smiles
 * @param {object} [options]
 * @param {number} [options.radius=2]
 * @param {number} [options.nBits=128]
 * @returns {Uint8Array} Bit vector of length nBits
 */
function morganFingerprint(smiles, { radius = 2, nBits = 128 } = {}) {
  const parsed = parseSMILES(smiles);
  const adj = buildAdjacency(parsed);
  const fp = new Uint8Array(nBits);

  // Initialise atom invariants (radius 0): element + degree + aromatic
  let identifiers = parsed.atoms.map((atom, i) => {
    const deg = adj[i].length;
    return hashString(`${atom.symbol}:${deg}:${atom.aromatic ? 1 : 0}`);
  });

  // Set bits for radius 0
  for (const id of identifiers) {
    fp[id % nBits] = 1;
  }

  // Iterate radii
  for (let r = 1; r <= radius; r++) {
    const nextIds = new Array(identifiers.length);
    for (let a = 0; a < parsed.atoms.length; a++) {
      const neighborHashes = adj[a]
        .map(e => identifiers[e.neighbor] + e.order * 7)
        .sort((x, y) => x - y);
      const combined = `${identifiers[a]}|${neighborHashes.join(',')}|r${r}`;
      nextIds[a] = hashString(combined);
      fp[nextIds[a] % nBits] = 1;
    }
    identifiers = nextIds;
  }

  return fp;
}

/**
 * Cosine similarity between two fingerprint vectors.
 */
function cosineSimilarity(fp1, fp2) {
  const len = Math.min(fp1.length, fp2.length);
  let dot = 0, mag1 = 0, mag2 = 0;
  for (let i = 0; i < len; i++) {
    dot += fp1[i] * fp2[i];
    mag1 += fp1[i] * fp1[i];
    mag2 += fp2[i] * fp2[i];
  }
  const denom = Math.sqrt(mag1) * Math.sqrt(mag2);
  return denom === 0 ? 0 : dot / denom;
}

/**
 * Tanimoto coefficient for binary fingerprints: |A∩B| / |A∪B|
 */
function tanimotoSimilarity(fp1, fp2) {
  const len = Math.min(fp1.length, fp2.length);
  let andCount = 0, orCount = 0;
  for (let i = 0; i < len; i++) {
    const a = fp1[i] ? 1 : 0;
    const b = fp2[i] ? 1 : 0;
    andCount += a & b;
    orCount += a | b;
  }
  return orCount === 0 ? 0 : andCount / orCount;
}

/**
 * k-Nearest Neighbours search using Tanimoto similarity.
 *
 * @param {Uint8Array} queryFp - Query fingerprint
 * @param {Array<{id: string, fingerprint: Uint8Array}>} database
 * @param {number} [k=5]
 * @returns {Array<{id: string, similarity: number}>}
 */
function knnSearch(queryFp, database, k = 5) {
  const scored = database.map(entry => ({
    id: entry.id,
    similarity: tanimotoSimilarity(queryFp, entry.fingerprint),
  }));

  scored.sort((a, b) => b.similarity - a.similarity);
  return scored.slice(0, k);
}

/**
 * Calculate molecular weight from SMILES.
 * Accounts for implicit hydrogens using simple valence rules.
 */
function molecularWeight(smiles) {
  const parsed = parseSMILES(smiles);
  const adj = buildAdjacency(parsed);

  const valences = { C: 4, N: 3, O: 2, S: 2, P: 3, F: 1, Cl: 1, Br: 1, I: 1 };

  let mw = 0;
  for (let i = 0; i < parsed.atoms.length; i++) {
    const atom = parsed.atoms[i];
    const weight = ATOMIC_WEIGHTS[atom.symbol] || 0;
    mw += weight;

    // Implicit hydrogens
    const v = valences[atom.symbol];
    if (v !== undefined) {
      let bondOrderSum = 0;
      for (const edge of adj[i]) bondOrderSum += edge.order;
      const effective = atom.aromatic ? v - 1 : v;
      const implicitH = Math.max(0, effective - bondOrderSum);
      mw += implicitH * ATOMIC_WEIGHTS.H;
    }
  }

  return Math.round(mw * 1000) / 1000;
}

module.exports = {
  parseSMILES,
  morganFingerprint,
  cosineSimilarity,
  tanimotoSimilarity,
  knnSearch,
  molecularWeight,
};
