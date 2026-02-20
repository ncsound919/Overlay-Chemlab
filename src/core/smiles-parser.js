'use strict';

const ATOMIC_WEIGHTS = {
  H: 1.008, C: 12.011, N: 14.007, O: 15.999,
  S: 32.065, P: 30.974, F: 18.998, Cl: 35.453,
  Br: 79.904, I: 126.904,
};

// Two-letter element symbols that can appear in SMILES (organic subset + brackets)
const TWO_LETTER = new Set(['Cl', 'Br']);

/**
 * Tokenize a SMILES string into an array of token objects.
 * Token types: 'atom', 'bond', 'branch_open', 'branch_close', 'ring'
 */
function tokenize(smiles) {
  const tokens = [];
  let i = 0;
  while (i < smiles.length) {
    const ch = smiles[i];

    // Branch delimiters
    if (ch === '(') {
      tokens.push({ type: 'branch_open', value: '(' });
      i++;
      continue;
    }
    if (ch === ')') {
      tokens.push({ type: 'branch_close', value: ')' });
      i++;
      continue;
    }

    // Explicit bonds
    if (ch === '-' || ch === '=' || ch === '#' || ch === ':') {
      tokens.push({ type: 'bond', value: ch });
      i++;
      continue;
    }

    // Ring closure digits (single digit or %nn)
    if (ch >= '0' && ch <= '9') {
      tokens.push({ type: 'ring', value: parseInt(ch, 10) });
      i++;
      continue;
    }
    if (ch === '%' && i + 2 < smiles.length) {
      const num = parseInt(smiles.substring(i + 1, i + 3), 10);
      tokens.push({ type: 'ring', value: num });
      i += 3;
      continue;
    }

    // Bracket atoms [...]
    if (ch === '[') {
      const close = smiles.indexOf(']', i);
      if (close === -1) {
        tokens.push({ type: 'atom', value: smiles.substring(i), bracket: true });
        i = smiles.length;
      } else {
        tokens.push({ type: 'atom', value: smiles.substring(i + 1, close), bracket: true });
        i = close + 1;
      }
      continue;
    }

    // Two-letter organic atoms
    if (i + 1 < smiles.length && TWO_LETTER.has(smiles.substring(i, i + 2))) {
      tokens.push({ type: 'atom', value: smiles.substring(i, i + 2), bracket: false });
      i += 2;
      continue;
    }

    // Single-letter atoms (organic subset + aromatic)
    if (/[A-Za-z]/.test(ch)) {
      tokens.push({ type: 'atom', value: ch, bracket: false });
      i++;
      continue;
    }

    // Charge and other characters inside brackets are already consumed above.
    // Skip unrecognised characters (e.g., '.', '/', '\\' for stereochemistry).
    i++;
  }
  return tokens;
}

/**
 * Resolve an atom token to its element symbol (title-case).
 */
function resolveElement(token) {
  let raw = token.value;
  if (token.bracket) {
    // Strip isotope prefix, H-count, charge, etc.  Keep element letters.
    raw = raw.replace(/^[0-9]*/, '');          // isotope
    raw = raw.replace(/[H][0-9]*/g, '');       // explicit H (simple)
    raw = raw.replace(/[+\-][0-9]*/g, '');     // charge
    raw = raw.replace(/@+/g, '');               // chirality
    if (raw.length === 0) return 'H';           // e.g. [H]
  }
  // Aromatic atoms: lowercase → title-case
  const aromaticMap = { c: 'C', n: 'N', o: 'O', s: 'S', p: 'P' };
  if (aromaticMap[raw]) return aromaticMap[raw];
  if (raw.length >= 2) {
    const two = raw[0].toUpperCase() + raw[1].toLowerCase();
    if (ATOMIC_WEIGHTS[two]) return two;
  }
  const one = raw[0].toUpperCase();
  if (ATOMIC_WEIGHTS[one]) return one;
  return raw[0].toUpperCase() + raw.substring(1).toLowerCase();
}

/**
 * Count implicit hydrogens for an organic-subset atom using simple valence rules.
 */
function implicitHydrogens(element, bondOrderSum, aromatic) {
  const valences = { C: 4, N: 3, O: 2, S: 2, P: 3, F: 1, Cl: 1, Br: 1, I: 1 };
  const v = valences[element];
  if (v === undefined) return 0;
  const effective = aromatic ? v - 1 : v;
  return Math.max(0, effective - bondOrderSum);
}

/**
 * Validate SMILES syntax.
 * @returns {{ valid: boolean, error?: string }}
 */
function validate(smiles) {
  if (typeof smiles !== 'string' || smiles.length === 0) {
    return { valid: false, error: 'Empty SMILES string' };
  }

  // Check balanced parentheses
  let depth = 0;
  for (const ch of smiles) {
    if (ch === '(') depth++;
    if (ch === ')') depth--;
    if (depth < 0) return { valid: false, error: 'Unmatched closing parenthesis' };
  }
  if (depth !== 0) return { valid: false, error: 'Unmatched opening parenthesis' };

  // Check balanced brackets
  let bDepth = 0;
  for (const ch of smiles) {
    if (ch === '[') bDepth++;
    if (ch === ']') bDepth--;
    if (bDepth < 0) return { valid: false, error: 'Unmatched closing bracket' };
  }
  if (bDepth !== 0) return { valid: false, error: 'Unmatched opening bracket' };

  // Check ring closure pairing
  const tokens = tokenize(smiles);
  const ringDigits = new Map();
  for (const tok of tokens) {
    if (tok.type === 'ring') {
      if (ringDigits.has(tok.value)) {
        ringDigits.delete(tok.value);
      } else {
        ringDigits.set(tok.value, true);
      }
    }
  }
  if (ringDigits.size > 0) {
    return { valid: false, error: 'Unclosed ring digit(s): ' + [...ringDigits.keys()].join(', ') };
  }

  // Must contain at least one atom
  if (!tokens.some(t => t.type === 'atom')) {
    return { valid: false, error: 'No atoms found' };
  }

  return { valid: true };
}

/**
 * Basic canonicalization: normalise aromatic notation to Kekulé-style uppercase.
 * This is a simplified canonicalizer — it normalises atom case and removes
 * stereochemistry markers (/, \\, @).
 */
function canonicalize(smiles) {
  let s = smiles;
  // Remove stereochemistry
  s = s.replace(/[/\\]/g, '');
  s = s.replace(/@+/g, '');

  // Convert lowercase aromatic atoms to uppercase (simple approach)
  // Only convert standalone aromatic atoms, not inside brackets
  let result = '';
  let inBracket = false;
  for (let i = 0; i < s.length; i++) {
    const ch = s[i];
    if (ch === '[') { inBracket = true; result += ch; continue; }
    if (ch === ']') { inBracket = false; result += ch; continue; }
    if (!inBracket && 'cnops'.includes(ch)) {
      result += ch.toUpperCase();
    } else {
      result += ch;
    }
  }
  return result;
}

/**
 * Count heavy atoms by element.
 * @returns {Object} e.g. { C: 6, H: 6, O: 1 }
 */
function atomCount(smiles) {
  const tokens = tokenize(smiles);

  // Build connectivity to calculate implicit H
  const atoms = [];   // [{element, aromatic, bondOrderSum}]
  const atomStack = [];
  let prevAtomIdx = -1;
  let pendingBond = null;

  for (const tok of tokens) {
    if (tok.type === 'atom') {
      const element = resolveElement(tok);
      const aromatic = tok.value.length === 1 && tok.value === tok.value.toLowerCase() && /[a-z]/.test(tok.value);
      const idx = atoms.length;
      atoms.push({ element, aromatic, bondOrderSum: 0, bracket: tok.bracket || false, raw: tok.value });

      // Bond to previous atom
      if (prevAtomIdx >= 0) {
        let order = 1;
        if (pendingBond === '=') order = 2;
        else if (pendingBond === '#') order = 3;
        else if (pendingBond === ':') order = 1; // aromatic
        else if (aromatic && atoms[prevAtomIdx].aromatic && !pendingBond) order = 1; // aromatic bond
        atoms[prevAtomIdx].bondOrderSum += order;
        atoms[idx].bondOrderSum += order;
      }
      pendingBond = null;
      prevAtomIdx = idx;
    } else if (tok.type === 'bond') {
      pendingBond = tok.value;
    } else if (tok.type === 'branch_open') {
      atomStack.push(prevAtomIdx);
    } else if (tok.type === 'branch_close') {
      prevAtomIdx = atomStack.pop();
    } else if (tok.type === 'ring') {
      // Ring closure adds a bond to a previous atom — handled below
    }
  }

  // Handle ring closures: second pass
  const ringOpens = new Map();
  let atomIdx = -1;
  for (const tok of tokens) {
    if (tok.type === 'atom') atomIdx++;
    if (tok.type === 'ring') {
      if (ringOpens.has(tok.value)) {
        const openIdx = ringOpens.get(tok.value);
        ringOpens.delete(tok.value);
        atoms[openIdx].bondOrderSum += 1;
        atoms[atomIdx].bondOrderSum += 1;
      } else {
        ringOpens.set(tok.value, atomIdx);
      }
    }
  }

  // Tally elements
  const counts = {};
  let totalH = 0;

  for (const atom of atoms) {
    const el = atom.element;
    counts[el] = (counts[el] || 0) + 1;

    // Implicit hydrogens for non-bracket atoms
    if (!atom.bracket) {
      const hCount = implicitHydrogens(el, atom.bondOrderSum, atom.aromatic);
      totalH += hCount;
    } else {
      // Bracket atoms: check for explicit H count in raw value
      const hMatch = atom.raw.match(/H([0-9]*)/);
      if (hMatch) {
        totalH += hMatch[1] ? parseInt(hMatch[1], 10) : 1;
      }
    }
  }

  if (totalH > 0) counts.H = (counts.H || 0) + totalH;

  return counts;
}

/**
 * Calculate molecular formula string in Hill system order.
 */
function molecularFormula(smiles) {
  const counts = atomCount(smiles);

  // Hill system: C first, then H, then alphabetical
  const parts = [];
  if (counts.C) {
    parts.push('C' + (counts.C > 1 ? counts.C : ''));
    if (counts.H) parts.push('H' + (counts.H > 1 ? counts.H : ''));
  }

  const remaining = Object.keys(counts)
    .filter(e => e !== 'C' && e !== 'H')
    .sort();

  for (const el of remaining) {
    parts.push(el + (counts[el] > 1 ? counts[el] : ''));
  }

  // If no carbon, just alphabetical
  if (!counts.C && counts.H) {
    // H wasn't added above
    const all = Object.keys(counts).sort();
    parts.length = 0;
    for (const el of all) {
      parts.push(el + (counts[el] > 1 ? counts[el] : ''));
    }
  }

  return parts.join('');
}

/**
 * Count bonds by type.
 * @returns {{ single: number, double: number, triple: number, aromatic: number }}
 */
function bondCount(smiles) {
  const tokens = tokenize(smiles);
  const result = { single: 0, double: 0, triple: 0, aromatic: 0 };

  let prevAtom = null;
  let pendingBond = null;
  const atomStack = [];
  const ringOpens = new Map();
  let atomIdx = -1;

  for (const tok of tokens) {
    if (tok.type === 'atom') {
      atomIdx++;
      const aromatic = tok.value.length === 1 && tok.value === tok.value.toLowerCase() && /[a-z]/.test(tok.value);

      if (prevAtom !== null) {
        if (pendingBond === '=') result.double++;
        else if (pendingBond === '#') result.triple++;
        else if (pendingBond === ':') result.aromatic++;
        else if (aromatic && prevAtom.aromatic) result.aromatic++;
        else result.single++;
      }
      pendingBond = null;
      prevAtom = { index: atomIdx, aromatic };
    } else if (tok.type === 'bond') {
      pendingBond = tok.value;
    } else if (tok.type === 'branch_open') {
      atomStack.push(prevAtom);
    } else if (tok.type === 'branch_close') {
      prevAtom = atomStack.pop();
    } else if (tok.type === 'ring') {
      if (ringOpens.has(tok.value)) {
        const openAtom = ringOpens.get(tok.value);
        ringOpens.delete(tok.value);
        // Ring closure bond — default single unless there's a pending bond
        if (pendingBond === '=') result.double++;
        else if (pendingBond === '#') result.triple++;
        else if (pendingBond === ':') result.aromatic++;
        else if (prevAtom && prevAtom.aromatic && openAtom.aromatic) result.aromatic++;
        else result.single++;
        pendingBond = null;
      } else {
        ringOpens.set(tok.value, prevAtom);
      }
    }
  }

  return result;
}

/**
 * Simple substructure check via canonical substring matching.
 *
 * WARNING: This uses naive substring matching on canonicalized SMILES and can
 * produce false positives (e.g., "C(=O)O" matches esters as well as free
 * carboxylic acids). For accurate results, use a graph-based substructure
 * matching algorithm such as the one provided by RDKit.
 */
function hasSubstructure(smiles, pattern) {
  const canonical = canonicalize(smiles);
  const canonicalPattern = canonicalize(pattern);
  return canonical.includes(canonicalPattern);
}

module.exports = {
  ATOMIC_WEIGHTS,
  tokenize,
  validate,
  canonicalize,
  molecularFormula,
  atomCount,
  bondCount,
  hasSubstructure,
};
