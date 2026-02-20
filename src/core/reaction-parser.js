'use strict';

const { atomCount } = require('./smiles-parser.js');

const ARROW_PATTERNS = [
  { re: /⟶/, symbol: '⟶' },
  { re: /→/, symbol: '→' },
  { re: /=>/, symbol: '=>' },
  { re: /->/, symbol: '->' },
];

/**
 * Parse a reaction string with arrow notation.
 * Supports optional conditions in square brackets after the arrow.
 * Example: "A + B ->[H3PO4, 85°C] C + D"
 */
function parseReaction(reactionStr) {
  if (typeof reactionStr !== 'string' || !reactionStr.trim()) {
    throw new Error('Reaction string must be a non-empty string');
  }

  const str = reactionStr.trim();
  let arrowSymbol = null;
  let arrowIndex = -1;
  let arrowLength = 0;

  for (const { re, symbol } of ARROW_PATTERNS) {
    const match = re.exec(str);
    if (match) {
      arrowSymbol = symbol;
      arrowIndex = match.index;
      arrowLength = match[0].length;
      break;
    }
  }

  if (arrowIndex === -1) {
    throw new Error('No reaction arrow found. Use ->, =>, →, or ⟶');
  }

  const lhs = str.slice(0, arrowIndex).trim();
  let rhs = str.slice(arrowIndex + arrowLength).trim();

  // Extract optional conditions in square brackets immediately after the arrow
  let conditions = null;
  if (rhs.startsWith('[')) {
    const closeBracket = rhs.indexOf(']');
    if (closeBracket !== -1) {
      conditions = rhs.slice(1, closeBracket).trim();
      rhs = rhs.slice(closeBracket + 1).trim();
    }
  }

  const reactants = splitComponents(lhs);
  const products = splitComponents(rhs);

  if (reactants.length === 0) throw new Error('No reactants found');
  if (products.length === 0) throw new Error('No products found');

  return { reactants, products, conditions, arrow: arrowSymbol };
}

/** Split "A + B + C" into trimmed, non-empty SMILES strings. */
function splitComponents(str) {
  return str
    .split('+')
    .map(s => s.trim())
    .filter(s => s.length > 0);
}

/**
 * Basic atom-balance check between reactants and products.
 * Uses atomCount from smiles-parser.
 */
function balanceCheck(reaction) {
  const reactantAtoms = sumAtoms(reaction.reactants);
  const productAtoms = sumAtoms(reaction.products);

  const allElements = new Set([
    ...Object.keys(reactantAtoms),
    ...Object.keys(productAtoms),
  ]);

  let balanced = true;
  for (const el of allElements) {
    if ((reactantAtoms[el] || 0) !== (productAtoms[el] || 0)) {
      balanced = false;
      break;
    }
  }

  return { balanced, reactantAtoms, productAtoms };
}

/** Sum atom counts across an array of SMILES strings. */
function sumAtoms(smilesList) {
  const total = {};
  for (const smi of smilesList) {
    const counts = atomCount(smi);
    for (const [el, n] of Object.entries(counts)) {
      total[el] = (total[el] || 0) + n;
    }
  }
  return total;
}

/**
 * Classify a reaction type by heuristic analysis of reactant/product SMILES.
 */
function reactionType(reaction) {
  const rSmiles = reaction.reactants.join('.');
  const pSmiles = reaction.products.join('.');

  // Esterification: carboxylic acid + alcohol → ester + water
  if (hasEsterFormation(rSmiles, pSmiles)) return 'esterification';

  // Hydrolysis: ester + water → acid + alcohol
  if (hasHydrolysis(rSmiles, pSmiles)) return 'hydrolysis';

  // Condensation: two molecules combine, small molecule (H2O) lost
  if (isCondensation(reaction)) return 'condensation';

  // Oxidation/reduction heuristics
  const rAtoms = sumAtoms(reaction.reactants);
  const pAtoms = sumAtoms(reaction.products);
  const oGain = (pAtoms.O || 0) - (rAtoms.O || 0);
  const hLoss = (rAtoms.H || 0) - (pAtoms.H || 0);

  if (oGain > 0 && hLoss >= 0) return 'oxidation';
  if (oGain < 0 || (hLoss < 0 && oGain <= 0)) return 'reduction';

  // Elimination: fewer products atoms in total, loss of small molecule
  if (reaction.products.length > reaction.reactants.length && hasSmallMolecule(pSmiles)) {
    return 'elimination';
  }

  // Addition: two reactants combine into one product, no small molecule lost
  if (reaction.reactants.length >= 2 && reaction.products.length === 1) return 'addition';

  // Substitution: same number of reactants and products, atom swap
  if (reaction.reactants.length === reaction.products.length) return 'substitution';

  return 'unknown';
}

function hasEsterFormation(rSmiles, pSmiles) {
  // Reactants contain C(=O)O (acid) or OC(=O) and alcohol (O), product contains ester linkage
  const hasAcid = /C\(=O\)O/.test(rSmiles) || /OC\(=O\)/.test(rSmiles);
  const hasEster = /C\(=O\)O[^H)\s]/.test(pSmiles) || /OC\(=O\)[^O\s)]/.test(pSmiles);
  const hasWater = pSmiles.includes('O') && hasAcid;
  return hasAcid && hasEster && hasWater;
}

function hasHydrolysis(rSmiles, pSmiles) {
  // Reactant has ester, products have acid + alcohol
  const hasEster = /C\(=O\)O[^H)\s]/.test(rSmiles) || /OC\(=O\)[^O\s)]/.test(rSmiles);
  const hasAcid = /C\(=O\)O/.test(pSmiles) || /OC\(=O\)/.test(pSmiles);
  return hasEster && hasAcid;
}

function isCondensation(reaction) {
  if (reaction.reactants.length < 2) return false;
  const pSmiles = reaction.products.join('.');
  // Water or other small molecule in products
  return reaction.products.length >= 2 && hasSmallMolecule(pSmiles);
}

function hasSmallMolecule(smiles) {
  const parts = smiles.split('.');
  for (const part of parts) {
    const trimmed = part.trim();
    if (trimmed === 'O' || trimmed === '[H]O[H]' || trimmed === 'N' ||
        trimmed === 'OC(O)=O' || trimmed === 'Cl' || trimmed === '[HCl]') {
      return true;
    }
  }
  return false;
}

/**
 * Format a parsed reaction back to a string.
 * @param {object} parsed - Parsed reaction object from parseReaction
 * @param {object} [options]
 * @param {boolean} [options.unicode=true] - Use unicode arrow (→) vs ASCII (->)
 */
function formatReaction(parsed, { unicode = true } = {}) {
  const arrow = unicode ? '→' : '->';
  const lhs = parsed.reactants.join(' + ');
  const rhs = parsed.products.join(' + ');

  let arrowPart = arrow;
  if (parsed.conditions) {
    arrowPart = `${arrow}[${parsed.conditions}]`;
  }

  return `${lhs} ${arrowPart} ${rhs}`;
}

/**
 * Parse a condition string into structured data.
 * Example: "H3PO4, 85°C, 15 min" →
 *   { catalyst: 'H3PO4', temperature: 85, time: 15, solvent: null }
 */
function extractConditions(conditionStr) {
  if (typeof conditionStr !== 'string' || !conditionStr.trim()) {
    return { catalyst: null, temperature: null, time: null, solvent: null };
  }

  const parts = conditionStr.split(',').map(s => s.trim()).filter(Boolean);

  let catalyst = null;
  let temperature = null;
  let time = null;
  let solvent = null;

  const KNOWN_SOLVENTS = [
    'H2O', 'water', 'EtOH', 'ethanol', 'MeOH', 'methanol',
    'THF', 'DCM', 'DMF', 'DMSO', 'acetone', 'toluene',
    'hexane', 'diethyl ether', 'Et2O', 'CHCl3', 'CCl4',
    'acetonitrile', 'MeCN', 'AcOH', 'benzene', 'dioxane',
  ];

  for (const part of parts) {
    // Temperature: match patterns like 85°C, 85 °C, 85C, 300K, 85 degC
    const tempMatch = part.match(/(-?\d+(?:\.\d+)?)\s*(?:°\s*C|degC|℃)/i);
    if (tempMatch) {
      temperature = parseFloat(tempMatch[1]);
      if (temperature < -273.15) temperature = -273.15;
      continue;
    }
    const tempK = part.match(/(\d+(?:\.\d+)?)\s*K\b/);
    if (tempK) {
      temperature = parseFloat(tempK[1]) - 273.15;
      temperature = Math.round(temperature * 100) / 100;
      continue;
    }

    // Time: match patterns like 15 min, 2 h, 30 s, 1.5 hr
    const timeMatch = part.match(/(\d+(?:\.\d+)?)\s*(min|mins|minute|minutes|h|hr|hrs|hour|hours|s|sec|secs|second|seconds|d|day|days)\b/i);
    if (timeMatch) {
      let val = parseFloat(timeMatch[1]);
      const unit = timeMatch[2].toLowerCase();
      // Normalize to minutes
      if (unit.startsWith('h')) val *= 60;
      else if (unit.startsWith('s')) val /= 60;
      else if (unit.startsWith('d')) val *= 1440;
      time = Math.round(val * 1000) / 1000;
      continue;
    }

    // Solvent check
    const lowerPart = part.toLowerCase();
    const isSolvent = KNOWN_SOLVENTS.some(
      s => lowerPart === s.toLowerCase() || lowerPart === `in ${s.toLowerCase()}`
    );
    if (isSolvent) {
      solvent = part.replace(/^in\s+/i, '');
      continue;
    }

    // Default: treat as catalyst if we don't have one yet
    if (!catalyst) {
      catalyst = part;
    } else if (!solvent) {
      // Second unknown part could be a solvent
      solvent = part;
    }
  }

  return { catalyst, temperature, time, solvent };
}

module.exports = {
  parseReaction,
  balanceCheck,
  reactionType,
  formatReaction,
  extractConditions,
};
