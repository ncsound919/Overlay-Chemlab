'use strict';

/**
 * Green Chemistry Metrics module.
 *
 * Implements the ACS Green Chemistry Institute (ACS GCI) Pharmaceutical
 * Roundtable toolkit and CHEM21 metrics for evaluating reaction sustainability.
 *
 * References:
 *   - Trost, B.M. Science 254, 1471–1477 (1991)             — Atom Economy
 *   - Sheldon, R.A. Chem. Ind. (London) 903–906 (1992)       — E-Factor
 *   - Curzons, A.D. et al. Green Chem. 3, 7–9 (2001)         — RME
 *   - Jiménez-González, C. et al. Org. Process Res. Dev. 15, 912–917 (2011) — PMI
 *   - Henderson, R.K. et al. Green Chem. 13, 854–862 (2011)  — CHEM21 Solvent Guide
 */

// ─────────────────────────────────────────────────────────────
// Atom Economy  (Trost, 1991)
// ─────────────────────────────────────────────────────────────

/**
 * Calculate Atom Economy.
 *
 * AE = MW(desired product) / Σ MW(all reactants) × 100
 *
 * Measures the fraction of reactant atoms that end up in the desired product.
 * A 100% AE means no atoms are wasted (e.g. addition reactions).
 *
 * @param {number} productMW   - Molecular weight of the desired product (g/mol)
 * @param {number} reactantsMW - Sum of MW of all stoichiometric reactants (g/mol)
 *                               (excludes catalysts, solvents)
 * @returns {number} Atom economy as a percentage (0–100)
 */
function atomEconomy(productMW, reactantsMW) {
  if (typeof productMW  !== 'number' || !isFinite(productMW))  throw new Error('productMW must be a finite number');
  if (typeof reactantsMW !== 'number' || !isFinite(reactantsMW)) throw new Error('reactantsMW must be a finite number');
  if (reactantsMW <= 0) throw new Error('reactantsMW must be > 0');
  if (productMW   <  0) throw new Error('productMW must be >= 0');
  return Math.min(100, Math.round((productMW / reactantsMW) * 10000) / 100);
}

// ─────────────────────────────────────────────────────────────
// E-Factor  (Sheldon, 1992)
// ─────────────────────────────────────────────────────────────

/**
 * Calculate the E-Factor (Environmental Factor).
 *
 * E = mass_waste_kg / mass_product_kg
 *
 * Sheldon benchmark values:
 *   Bulk chemicals       E <  1–5
 *   Fine chemicals       E =  5–50
 *   Active pharma APIs   E = 25–100  (can exceed 100)
 *
 * @param {number} wasteKg   - Total mass of all waste streams (kg)
 * @param {number} productKg - Mass of isolated product (kg)
 * @returns {number} E-factor (dimensionless)
 */
function eFactor(wasteKg, productKg) {
  if (productKg <= 0) throw new Error('productKg must be > 0');
  if (wasteKg   <  0) throw new Error('wasteKg must be >= 0');
  return Math.round((wasteKg / productKg) * 1e6) / 1e6;
}

/**
 * Calculate E-factor from total inputs and product mass.
 *
 * E = (Σ mass_inputs − mass_product) / mass_product  =  PMI − 1
 *
 * @param {number} totalInputsKg - Sum of all inputs: reactants + solvents + catalysts (kg)
 * @param {number} productKg     - Mass of isolated product (kg)
 * @returns {number} E-factor
 */
function eFactorFromInputs(totalInputsKg, productKg) {
  if (productKg     <= 0) throw new Error('productKg must be > 0');
  if (totalInputsKg <= 0) throw new Error('totalInputsKg must be > 0');
  if (totalInputsKg < productKg) throw new Error('totalInputsKg cannot be less than productKg');
  return Math.round(((totalInputsKg - productKg) / productKg) * 1e6) / 1e6;
}

// ─────────────────────────────────────────────────────────────
// Process Mass Intensity  (ACS GCI Pharmaceutical Roundtable, 2011)
// ─────────────────────────────────────────────────────────────

/**
 * Calculate Process Mass Intensity (PMI).
 *
 * PMI = total_mass_of_all_materials / mass_product
 *     = E-factor + 1    (minimum possible PMI = 1, i.e. perfect reaction)
 *
 * ACS GCI Pharmaceutical Roundtable 2011 benchmark:
 *   New molecular entities (full synthesis): PMI median ≈ 200 kg/kg
 *   Target for API synthesis step: PMI < 100 kg/kg
 *
 * @param {number} totalMassKg - Total mass of all inputs incl. solvents (kg)
 * @param {number} productKg   - Mass of isolated product (kg)
 * @returns {number} PMI (dimensionless; minimum = 1)
 */
function pmi(totalMassKg, productKg) {
  if (productKg  <= 0) throw new Error('productKg must be > 0');
  if (totalMassKg <= 0) throw new Error('totalMassKg must be > 0');
  if (totalMassKg < productKg) throw new Error('totalMassKg cannot be less than productKg');
  return Math.round((totalMassKg / productKg) * 1e6) / 1e6;
}

// ─────────────────────────────────────────────────────────────
// Reaction Mass Efficiency  (Curzons et al., 2001)
// ─────────────────────────────────────────────────────────────

/**
 * Calculate Reaction Mass Efficiency (RME).
 *
 * Combined metric that incorporates yield, atom economy and reagent excess:
 *   RME = yield_fraction × (AE / 100) / stoichiometricFactor
 *
 * Ideal RME = 1 (100 % yield, 100 % AE, no excess reagents).
 * Typical pharmaceutical reactions: RME = 0.2–0.5.
 *
 * @param {number} yieldFraction         - Fractional isolated yield (0–1)
 * @param {number} atomEconomyPct        - Atom economy in % (0–100)
 * @param {number} [stoichiometricFactor=1] - Ratio of actual reagent mass used
 *                                           to the theoretical stoichiometric mass
 * @returns {number} RME as a fraction (0–1)
 */
function reactionMassEfficiency(yieldFraction, atomEconomyPct, stoichiometricFactor = 1) {
  if (yieldFraction < 0 || yieldFraction > 1)   throw new Error('yieldFraction must be in [0, 1]');
  if (atomEconomyPct < 0 || atomEconomyPct > 100) throw new Error('atomEconomyPct must be in [0, 100]');
  if (stoichiometricFactor <= 0) throw new Error('stoichiometricFactor must be > 0');
  return Math.round((yieldFraction * atomEconomyPct / 100 / stoichiometricFactor) * 1e6) / 1e6;
}

// ─────────────────────────────────────────────────────────────
// Carbon Efficiency
// ─────────────────────────────────────────────────────────────

/**
 * Calculate Carbon Efficiency.
 *
 * CE = (moles C in product / moles C in all C-containing reactants) × 100
 *
 * Measures how effectively carbon feedstock atoms end up in the product.
 *
 * @param {number} carbonInProduct   - C atoms (or moles) in the product
 * @param {number} carbonInFeedstock - Total C atoms in all C-containing reactants
 * @returns {number} Carbon efficiency as a percentage (0–100)
 */
function carbonEfficiency(carbonInProduct, carbonInFeedstock) {
  if (carbonInFeedstock <= 0) throw new Error('carbonInFeedstock must be > 0');
  if (carbonInProduct   <  0) throw new Error('carbonInProduct must be >= 0');
  return Math.min(100, Math.round((carbonInProduct / carbonInFeedstock) * 10000) / 100);
}

// ─────────────────────────────────────────────────────────────
// CHEM21 Solvent Classification
// ─────────────────────────────────────────────────────────────

/**
 * CHEM21 / ACS GCI Pharmaceutical Roundtable solvent selection guide.
 *
 * Based on: Henderson et al. Green Chem. 13, 854–862 (2011) and the
 * updated CHEM21 Green Metrics Toolkit solvent selection guide.
 *
 * Lookup key is the solvent common name, lower-cased.
 */
const SOLVENT_CLASSIFICATION = {
  // ── Recommended ─────────────────────────────────────────────
  'water':                    { score: 'recommended', concern: null },
  'ethanol':                  { score: 'recommended', concern: null },
  'isopropanol':              { score: 'recommended', concern: null },
  'n-butanol':                { score: 'recommended', concern: null },
  'ethyl acetate':            { score: 'recommended', concern: null },
  'acetone':                  { score: 'recommended', concern: null },
  'methyl ethyl ketone':      { score: 'recommended', concern: null },
  '2-methylthf':              { score: 'recommended', concern: null },
  'heptane':                  { score: 'recommended', concern: null },
  'anisole':                  { score: 'recommended', concern: null },
  'dimethyl isosorbide':      { score: 'recommended', concern: null },
  'cyrene':                   { score: 'recommended', concern: null },

  // ── Usable ───────────────────────────────────────────────────
  'methanol':                 { score: 'usable', concern: 'toxic' },
  'tert-butanol':             { score: 'usable', concern: null },
  'isobutanol':               { score: 'usable', concern: null },
  'tert-amyl alcohol':        { score: 'usable', concern: null },
  'methyl isobutyl ketone':   { score: 'usable', concern: null },
  'thf':                      { score: 'usable', concern: 'peroxide_formation' },
  'diisopropyl ether':        { score: 'usable', concern: 'peroxide_formation' },
  'diethyl ether':            { score: 'usable', concern: 'peroxide_formation' },
  'toluene':                  { score: 'usable', concern: 'reproductive_hazard' },
  'xylene':                   { score: 'usable', concern: 'reproductive_hazard' },
  'cyclohexane':              { score: 'usable', concern: 'voc' },
  'hexane':                   { score: 'usable', concern: 'neurotoxic' },
  'acetic acid':              { score: 'usable', concern: 'corrosive' },
  'propylene glycol':         { score: 'usable', concern: null },

  // ── Problematic ──────────────────────────────────────────────
  'acetonitrile':             { score: 'problematic', concern: 'toxic_manufacturing' },
  'dcm':                      { score: 'problematic', concern: 'carcinogen_suspect' },
  'dichloromethane':          { score: 'problematic', concern: 'carcinogen_suspect' },
  'dmf':                      { score: 'problematic', concern: 'reproductive_toxin' },
  'n,n-dimethylformamide':    { score: 'problematic', concern: 'reproductive_toxin' },
  'dmso':                     { score: 'problematic', concern: null },
  'dimethyl sulfoxide':       { score: 'problematic', concern: null },
  'dioxane':                  { score: 'problematic', concern: 'carcinogen' },
  '1,4-dioxane':              { score: 'problematic', concern: 'carcinogen' },
  'pyridine':                 { score: 'problematic', concern: 'toxic' },
  'nitrobenzene':             { score: 'problematic', concern: 'toxic' },

  // ── Hazardous ────────────────────────────────────────────────
  'benzene':                  { score: 'hazardous', concern: 'carcinogen' },
  'chloroform':               { score: 'hazardous', concern: 'carcinogen_suspect' },
  'chcl3':                    { score: 'hazardous', concern: 'carcinogen_suspect' },
  'ccl4':                     { score: 'hazardous', concern: 'carcinogen' },
  'carbon tetrachloride':     { score: 'hazardous', concern: 'carcinogen' },
  '1,2-dichloroethane':       { score: 'hazardous', concern: 'carcinogen' },
  'nitromethane':             { score: 'hazardous', concern: 'explosive' },
  'hexamethylphosphoramide':  { score: 'hazardous', concern: 'carcinogen' },
  'hmpa':                     { score: 'hazardous', concern: 'carcinogen' },
  'dmac':                     { score: 'hazardous', concern: 'reproductive_toxin' },
  'n,n-dimethylacetamide':    { score: 'hazardous', concern: 'reproductive_toxin' },
};

/**
 * Classify a solvent according to the CHEM21 / ACS GCI solvent selection guide.
 *
 * @param {string} solventName - Common name of the solvent (case-insensitive)
 * @returns {{ name: string, score: 'recommended'|'usable'|'problematic'|'hazardous'|'unknown', concern: string|null }}
 */
function classifySolvent(solventName) {
  if (typeof solventName !== 'string' || !solventName.trim()) {
    throw new Error('solventName must be a non-empty string');
  }
  const key   = solventName.trim().toLowerCase();
  const entry = SOLVENT_CLASSIFICATION[key] || { score: 'unknown', concern: null };
  return { name: solventName, score: entry.score, concern: entry.concern };
}

// ─────────────────────────────────────────────────────────────
// Composite Green Score
// ─────────────────────────────────────────────────────────────

/**
 * Calculate a composite green chemistry score for a reaction step.
 *
 * Each sub-score is normalised to 0–100 and combined with weights:
 *   Atom Economy  30 %
 *   E-Factor      25 %  (inverse: lower E = higher score)
 *   Yield         30 %
 *   Solvent       15 %  (CHEM21 classification)
 *
 * Grade scale: A ≥ 80, B ≥ 60, C ≥ 40, D ≥ 20, F < 20.
 *
 * @param {object} params
 * @param {number} params.atomEconomyPct  - Atom economy (0–100)
 * @param {number} params.eFact          - E-factor (lower = greener)
 * @param {number} params.yieldPct       - Isolated yield percentage (0–100)
 * @param {string} [params.solvent]      - Solvent name for CHEM21 lookup
 * @returns {{ score: number, grade: string, details: object }}
 */
function greenScore({ atomEconomyPct = 0, eFact = 0, yieldPct = 0, solvent } = {}) {
  // Atom economy: direct %
  const aeScore = Math.min(100, Math.max(0, atomEconomyPct));

  // E-factor: 0 → 100 pts; ≥ 100 → 0 pts (linear interpolation)
  const efScore = Math.max(0, 100 - eFact);

  // Yield: direct %
  const yScore = Math.min(100, Math.max(0, yieldPct));

  // Solvent score
  let solventScore = 50; // neutral when not specified
  if (solvent) {
    const sc = classifySolvent(solvent).score;
    solventScore = sc === 'recommended' ? 100
      : sc === 'usable'       ? 70
      : sc === 'problematic'  ? 40
      : sc === 'hazardous'    ? 10
      : 50;
  }

  const score   = 0.30 * aeScore + 0.25 * efScore + 0.30 * yScore + 0.15 * solventScore;
  const rounded = Math.round(score * 10) / 10;
  const grade   = rounded >= 80 ? 'A' : rounded >= 60 ? 'B' : rounded >= 40 ? 'C' : rounded >= 20 ? 'D' : 'F';

  return {
    score: rounded,
    grade,
    details: { aeScore, efScore, yScore, solventScore },
  };
}

module.exports = {
  atomEconomy,
  eFactor,
  eFactorFromInputs,
  pmi,
  reactionMassEfficiency,
  carbonEfficiency,
  classifySolvent,
  greenScore,
  SOLVENT_CLASSIFICATION,
};
