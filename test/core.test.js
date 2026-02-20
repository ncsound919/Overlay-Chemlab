'use strict';

const { describe, it } = require('node:test');
const assert = require('node:assert/strict');

const arrhenius = require('../src/core/arrhenius');
const smilesParser = require('../src/core/smiles-parser');
const molEmbed = require('../src/core/molecular-embeddings');
const bayesian = require('../src/core/bayesian-risk');
const reactionParser = require('../src/core/reaction-parser');
const drugLikeness = require('../src/core/drug-likeness');
const greenChem = require('../src/core/green-chemistry');

// ────────────────────── Arrhenius ──────────────────────

describe('arrhenius', () => {
  it('arrheniusRate returns correct k at 25°C', () => {
    const k = arrhenius.arrheniusRate(65300, 4.2e9, 25);
    assert.ok(k > 0, 'rate constant must be positive');
    // k = 4.2e9 * exp(-65300 / (8.314 * 298.15))
    const expected = 4.2e9 * Math.exp(-65300 / (8.314 * 298.15));
    assert.ok(Math.abs(k - expected) < 1e-15);
  });

  it('arrheniusRate increases with temperature', () => {
    const k25 = arrhenius.arrheniusRate(65300, 4.2e9, 25);
    const k80 = arrhenius.arrheniusRate(65300, 4.2e9, 80);
    assert.ok(k80 > k25);
  });

  it('simulateEsterification returns valid result', () => {
    const res = arrhenius.simulateEsterification({
      temperature: 80, time: 3600, equiv: 1.5, catalystDrops: 3,
    });
    assert.ok(res.conversion >= 0 && res.conversion <= 1);
    assert.ok(res.yield >= 0 && res.yield <= 1);
    assert.ok(res.pH > 0 && res.pH < 14);
    assert.ok(Array.isArray(res.sideProducts));
    assert.ok(res.rateConstant > 0);
  });

  it('simulateEsterification produces side products at high T', () => {
    const res = arrhenius.simulateEsterification({
      temperature: 200, time: 3600, equiv: 1.5, catalystDrops: 3,
    });
    assert.ok(res.sideProducts.length > 0);
  });

  it('runKineticsODE returns time-series data', () => {
    const data = arrhenius.runKineticsODE({
      temperature: 120, totalTime: 7200, steps: 50, equiv: 1.5, catalystDrops: 5,
    });
    assert.equal(data.length, 51); // steps + 1
    assert.equal(data[0].time, 0);
    assert.equal(data[0].conversion, 0);
    assert.ok(data[data.length - 1].conversion > 0);
    // Conversion should be monotonically increasing
    for (let i = 1; i < data.length; i++) {
      assert.ok(data[i].conversion >= data[i - 1].conversion - 1e-12);
    }
  });

  it('calculateYield accounts for losses', () => {
    const res = arrhenius.calculateYield({
      conversion: 0.8, equiv: 1.5, crystallizationLoss: 0.05, filtrationLoss: 0.02,
    });
    assert.ok(res.isolatedYield < res.theoreticalYield);
    assert.ok(res.isolatedYield > 0);
    assert.equal(res.loss, 0.07);
  });
});

// ────────────────────── SMILES Parser ──────────────────────

describe('smiles-parser', () => {
  it('tokenize parses ethanol', () => {
    const tokens = smilesParser.tokenize('CCO');
    const atoms = tokens.filter(t => t.type === 'atom');
    assert.equal(atoms.length, 3);
  });

  it('validate accepts valid SMILES', () => {
    assert.ok(smilesParser.validate('CCO').valid);
    assert.ok(smilesParser.validate('c1ccccc1').valid);
    assert.ok(smilesParser.validate('CC(=O)O').valid);
  });

  it('validate rejects invalid SMILES', () => {
    assert.ok(!smilesParser.validate('').valid);
    assert.ok(!smilesParser.validate('CC(O').valid);
    assert.ok(!smilesParser.validate('C1CC').valid);
  });

  it('molecularFormula for aspirin', () => {
    // Aspirin: CC(=O)Oc1ccccc1C(=O)O → C9H8O4
    const formula = smilesParser.molecularFormula('CC(=O)Oc1ccccc1C(=O)O');
    assert.equal(formula, 'C9H8O4');
  });

  it('atomCount for ethanol', () => {
    const counts = smilesParser.atomCount('CCO');
    assert.equal(counts.C, 2);
    assert.equal(counts.O, 1);
    assert.equal(counts.H, 6);
  });

  it('bondCount for ethene', () => {
    const bonds = smilesParser.bondCount('C=C');
    assert.equal(bonds.double, 1);
    assert.equal(bonds.single, 0);
  });

  it('canonicalize normalises aromatic atoms', () => {
    const canon = smilesParser.canonicalize('c1ccccc1');
    assert.equal(canon, 'C1CCCCC1');
  });

  it('hasSubstructure finds pattern', () => {
    assert.ok(smilesParser.hasSubstructure('CC(=O)O', 'C(=O)O'));
    assert.ok(!smilesParser.hasSubstructure('CCO', 'N'));
  });
});

// ────────────────────── Molecular Embeddings ──────────────────────

describe('molecular-embeddings', () => {
  it('parseSMILES parses benzene', () => {
    const parsed = molEmbed.parseSMILES('c1ccccc1');
    assert.equal(parsed.atoms.length, 6);
    assert.equal(parsed.rings, 1);
    assert.ok(parsed.atoms[0].aromatic);
  });

  it('parseSMILES parses ethanol', () => {
    const parsed = molEmbed.parseSMILES('CCO');
    assert.equal(parsed.atoms.length, 3);
    assert.equal(parsed.bonds.length, 2);
  });

  it('morganFingerprint returns Uint8Array', () => {
    const fp = molEmbed.morganFingerprint('CCO');
    assert.ok(fp instanceof Uint8Array);
    assert.equal(fp.length, 128);
    // Should have some bits set
    const setBits = fp.reduce((s, v) => s + v, 0);
    assert.ok(setBits > 0);
  });

  it('identical molecules have Tanimoto = 1', () => {
    const fp1 = molEmbed.morganFingerprint('CCO');
    const fp2 = molEmbed.morganFingerprint('CCO');
    assert.equal(molEmbed.tanimotoSimilarity(fp1, fp2), 1);
  });

  it('different molecules have Tanimoto < 1', () => {
    const fp1 = molEmbed.morganFingerprint('CCO');
    const fp2 = molEmbed.morganFingerprint('c1ccccc1');
    const sim = molEmbed.tanimotoSimilarity(fp1, fp2);
    assert.ok(sim >= 0 && sim < 1);
  });

  it('cosineSimilarity for identical fingerprints is 1', () => {
    const fp = molEmbed.morganFingerprint('CCO');
    assert.ok(Math.abs(molEmbed.cosineSimilarity(fp, fp) - 1) < 1e-10);
  });

  it('knnSearch returns sorted results', () => {
    const db = [
      { id: 'ethanol', fingerprint: molEmbed.morganFingerprint('CCO') },
      { id: 'methanol', fingerprint: molEmbed.morganFingerprint('CO') },
      { id: 'benzene', fingerprint: molEmbed.morganFingerprint('c1ccccc1') },
    ];
    const query = molEmbed.morganFingerprint('CCO');
    const results = molEmbed.knnSearch(query, db, 2);
    assert.equal(results.length, 2);
    assert.equal(results[0].id, 'ethanol');
    assert.equal(results[0].similarity, 1);
  });

  it('molecularWeight for water', () => {
    // H2O: O with 2 implicit H → 15.999 + 2*1.008 = 18.015
    const mw = molEmbed.molecularWeight('O');
    assert.ok(Math.abs(mw - 18.015) < 0.01);
  });

  it('molecularWeight for methane', () => {
    // CH4: C with 4 implicit H → 12.011 + 4*1.008 = 16.043
    const mw = molEmbed.molecularWeight('C');
    assert.ok(Math.abs(mw - 16.043) < 0.01);
  });
});

// ────────────────────── Bayesian Risk ──────────────────────

describe('bayesian-risk', () => {
  it('logistic(0) = 0.5', () => {
    assert.ok(Math.abs(bayesian.logistic(0) - 0.5) < 1e-10);
  });

  it('logistic is bounded [0, 1]', () => {
    assert.ok(bayesian.logistic(-100) > 0);
    assert.ok(bayesian.logistic(100) <= 1);
    assert.ok(bayesian.logistic(-100) < 0.01);
    assert.ok(bayesian.logistic(100) > 0.99);
  });

  it('riskScore returns probability in [0,1] and matches expected value', () => {
    const emb = new Array(128).fill(0).map((_, i) => (i % 2 ? 1 : 0));
    const linear = emb.reduce(
      (sum, val, idx) => sum + val * bayesian.DEFAULT_WEIGHTS[idx],
      bayesian.DEFAULT_BIAS,
    );
    const expected = bayesian.logistic(linear);
    const score = bayesian.riskScore(emb, bayesian.DEFAULT_WEIGHTS, bayesian.DEFAULT_BIAS);
    assert.ok(score >= 0 && score <= 1);
    assert.ok(Math.abs(score - expected) < 1e-12);
  });

  it('classifyRisk correctly classifies', () => {
    assert.equal(bayesian.classifyRisk(0.1).level, 'low');
    assert.equal(bayesian.classifyRisk(0.5).level, 'medium');
    assert.equal(bayesian.classifyRisk(0.8).level, 'high');
    assert.ok(bayesian.classifyRisk(0.8).exceeds);
    assert.ok(!bayesian.classifyRisk(0.3).exceeds);
  });

  it('batchRiskScore returns array of scores', () => {
    const embs = [
      new Array(128).fill(1),
      new Array(128).fill(0),
    ];
    const scores = bayesian.batchRiskScore(embs, bayesian.DEFAULT_WEIGHTS, bayesian.DEFAULT_BIAS);
    assert.equal(scores.length, 2);
    scores.forEach(s => assert.ok(s >= 0 && s <= 1));
  });

  it('confidenceInterval returns valid CI', () => {
    const scores = [0.3, 0.4, 0.35, 0.38, 0.42, 0.31, 0.39, 0.36];
    const ci = bayesian.confidenceInterval(scores);
    assert.ok(ci.lower <= ci.mean);
    assert.ok(ci.mean <= ci.upper);
    assert.ok(ci.lower >= 0);
    assert.ok(ci.upper <= 1);
  });

  it('DEFAULT_WEIGHTS has 128 dimensions', () => {
    assert.equal(bayesian.DEFAULT_WEIGHTS.length, 128);
  });

  it('DEFAULT_BIAS is a number', () => {
    assert.equal(typeof bayesian.DEFAULT_BIAS, 'number');
  });
});

// ────────────────────── Reaction Parser ──────────────────────

describe('reaction-parser', () => {
  it('parseReaction splits reactants and products', () => {
    const r = reactionParser.parseReaction('A + B -> C + D');
    assert.deepEqual(r.reactants, ['A', 'B']);
    assert.deepEqual(r.products, ['C', 'D']);
    assert.equal(r.arrow, '->');
    assert.equal(r.conditions, null);
  });

  it('parseReaction extracts conditions in brackets', () => {
    const r = reactionParser.parseReaction('A + B ->[H3PO4, 85°C] C');
    assert.equal(r.conditions, 'H3PO4, 85°C');
    assert.deepEqual(r.reactants, ['A', 'B']);
    assert.deepEqual(r.products, ['C']);
  });

  it('parseReaction supports unicode arrows', () => {
    const r = reactionParser.parseReaction('A → B');
    assert.equal(r.arrow, '→');
    assert.deepEqual(r.reactants, ['A']);
    assert.deepEqual(r.products, ['B']);
  });

  it('parseReaction throws on missing arrow', () => {
    assert.throws(() => reactionParser.parseReaction('A B C'));
  });

  it('balanceCheck detects balanced reactions', () => {
    const parsed = reactionParser.parseReaction('CCO + O -> CC=O + O');
    const check = reactionParser.balanceCheck(parsed);
    assert.equal(typeof check.balanced, 'boolean');
    assert.ok(check.reactantAtoms);
    assert.ok(check.productAtoms);
  });

  it('formatReaction produces correct string', () => {
    const parsed = { reactants: ['A', 'B'], products: ['C'], conditions: null };
    assert.equal(reactionParser.formatReaction(parsed), 'A + B → C');
    assert.equal(reactionParser.formatReaction(parsed, { unicode: false }), 'A + B -> C');
  });

  it('formatReaction includes conditions', () => {
    const parsed = { reactants: ['A'], products: ['B'], conditions: 'H3PO4' };
    assert.equal(reactionParser.formatReaction(parsed), 'A →[H3PO4] B');
  });

  it('extractConditions parses temperature and time', () => {
    const conds = reactionParser.extractConditions('H3PO4, 85°C, 15 min');
    assert.equal(conds.catalyst, 'H3PO4');
    assert.equal(conds.temperature, 85);
    assert.equal(conds.time, 15);
  });

  it('extractConditions returns nulls for empty string', () => {
    const conds = reactionParser.extractConditions('');
    assert.equal(conds.catalyst, null);
    assert.equal(conds.temperature, null);
  });

  it('reactionType returns a string', () => {
    const parsed = reactionParser.parseReaction('CC(=O)OC(=O)C + OC(=O)c1ccccc1O -> CC(=O)Oc1ccccc1C(=O)O + CC(O)=O');
    const rtype = reactionParser.reactionType(parsed);
    assert.equal(typeof rtype, 'string');
    assert.ok(rtype.length > 0);
  });
});

// ────────────────────── Drug-Likeness ──────────────────────

describe('drug-likeness', () => {
  // Aspirin: CC(=O)Oc1ccccc1C(=O)O — a well-characterised drug molecule
  const aspirinSMILES = 'CC(=O)Oc1ccccc1C(=O)O';

  it('countHBD for aspirin returns 1', () => {
    assert.equal(drugLikeness.countHBD(aspirinSMILES), 1);
  });

  it('countHBA for aspirin returns 4', () => {
    assert.equal(drugLikeness.countHBA(aspirinSMILES), 4);
  });

  it('countHBD for ethanol returns 1', () => {
    assert.equal(drugLikeness.countHBD('CCO'), 1);
  });

  it('countHBA for ethanol returns 1', () => {
    assert.equal(drugLikeness.countHBA('CCO'), 1);
  });

  it('countHBD for benzene returns 0', () => {
    assert.equal(drugLikeness.countHBD('c1ccccc1'), 0);
  });

  it('countRotatableBonds for aspirin returns 3', () => {
    assert.equal(drugLikeness.countRotatableBonds(aspirinSMILES), 3);
  });

  it('countRotatableBonds for benzene returns 0', () => {
    assert.equal(drugLikeness.countRotatableBonds('c1ccccc1'), 0);
  });

  it('countRotatableBonds for n-butane returns 1', () => {
    // CCCC: bond C1-C2 only (C0-C1 and C2-C3 are terminal)
    assert.equal(drugLikeness.countRotatableBonds('CCCC'), 1);
  });

  it('estimateTPSA for aspirin matches PubChem (63.6 Å²)', () => {
    const tpsa = drugLikeness.estimateTPSA(aspirinSMILES);
    assert.ok(Math.abs(tpsa - 63.6) < 1, `Expected ~63.6, got ${tpsa}`);
  });

  it('estimateTPSA for benzene returns 0', () => {
    assert.equal(drugLikeness.estimateTPSA('c1ccccc1'), 0);
  });

  it('estimateTPSA for ethanol is > 0', () => {
    assert.ok(drugLikeness.estimateTPSA('CCO') > 0);
  });

  it('estimateLogP for benzene is positive', () => {
    // Benzene logP = 2.13; our simplified model should give a positive value
    assert.ok(drugLikeness.estimateLogP('c1ccccc1') > 0);
  });

  it('estimateLogP for ethanol is negative', () => {
    // Ethanol logP = −0.31; our model should give a negative value
    assert.ok(drugLikeness.estimateLogP('CCO') < 0);
  });

  it('lipinskiRuleOfFive: aspirin passes all criteria', () => {
    const result = drugLikeness.lipinskiRuleOfFive(aspirinSMILES);
    assert.ok(result.pass);
    assert.equal(result.violations, 0);
    assert.equal(result.criteria.length, 4);
    result.criteria.forEach(c => assert.ok(c.pass));
  });

  it('lipinskiRuleOfFive: large lipophilic molecule fails', () => {
    // Cyclosporin A-like large MW lipophilic peptide — very high MW + logP
    // Use a simple heavy large molecule: C60 fullerene skeleton (very high MW)
    // Instead use a synthetic high-MW molecule
    const largeSMILES = 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'; // MW ~450, very high logP
    const result = drugLikeness.lipinskiRuleOfFive(largeSMILES);
    assert.ok(!result.pass);
    assert.ok(result.violations > 0);
  });

  it('veberRules: aspirin passes', () => {
    const result = drugLikeness.veberRules(aspirinSMILES);
    assert.ok(result.pass);
    assert.ok(result.rotBonds <= 10);
    assert.ok(result.tpsa <= 140);
  });

  it('ruleOfThree: aspirin passes with one tolerated HBA violation', () => {
    // Aspirin MW=180.16, logP≈0.45, HBD=1, HBA=4 → HBA violates (>3) but one violation is tolerated
    const result = drugLikeness.ruleOfThree(aspirinSMILES);
    assert.ok(result.pass);
    assert.equal(result.violations, 1);
  });

  it('assessDruglikeness returns complete report', () => {
    const report = drugLikeness.assessDruglikeness(aspirinSMILES);
    assert.equal(report.smiles, aspirinSMILES);
    assert.ok(typeof report.descriptors.mw === 'number');
    assert.ok(typeof report.descriptors.logP === 'number');
    assert.ok(typeof report.descriptors.tpsa === 'number');
    assert.ok(typeof report.descriptors.hbd === 'number');
    assert.ok(typeof report.descriptors.hba === 'number');
    assert.ok(typeof report.descriptors.rotBonds === 'number');
    assert.ok(typeof report.overallDruglike === 'boolean');
    assert.ok(report.lipinski);
    assert.ok(report.veber);
    assert.ok(report.ruleOfThree);
  });
});

// ────────────────────── Green Chemistry ──────────────────────

describe('green-chemistry', () => {
  // Aspirin synthesis: salicylic acid (MW 138.12) + acetic anhydride (MW 102.09)
  // → aspirin (MW 180.16) + acetic acid (MW 60.05)
  const reactantsMW = 138.12 + 102.09; // 240.21
  const productMW   = 180.16;

  it('atomEconomy for aspirin synthesis', () => {
    const ae = greenChem.atomEconomy(productMW, reactantsMW);
    assert.ok(ae > 0 && ae <= 100);
    // 180.16/240.21 = 75.0%
    assert.ok(Math.abs(ae - 75.0) < 1);
  });

  it('atomEconomy for addition reaction is 100%', () => {
    // ethene + H2 → ethane: all atoms conserved
    const ae = greenChem.atomEconomy(30, 30);
    assert.equal(ae, 100);
  });

  it('atomEconomy throws on zero reactantsMW', () => {
    assert.throws(() => greenChem.atomEconomy(100, 0));
  });

  it('eFactor calculation', () => {
    const ef = greenChem.eFactor(4.5, 1.0);
    assert.equal(ef, 4.5);
  });

  it('eFactor throws on zero productKg', () => {
    assert.throws(() => greenChem.eFactor(1, 0));
  });

  it('eFactorFromInputs matches eFactor definition', () => {
    // E = (inputs - product) / product
    const ef = greenChem.eFactorFromInputs(10, 1);
    assert.equal(ef, 9);
  });

  it('eFactorFromInputs = PMI - 1', () => {
    const ef  = greenChem.eFactorFromInputs(15, 3);
    const p   = greenChem.pmi(15, 3);
    assert.ok(Math.abs((p - 1) - ef) < 1e-9);
  });

  it('pmi for perfect reaction is 1', () => {
    assert.equal(greenChem.pmi(1, 1), 1);
  });

  it('pmi throws when inputs < product', () => {
    assert.throws(() => greenChem.pmi(0.5, 1));
  });

  it('reactionMassEfficiency for ideal reaction is 1', () => {
    const rme = greenChem.reactionMassEfficiency(1.0, 100, 1.0);
    assert.equal(rme, 1);
  });

  it('reactionMassEfficiency combines yield and AE', () => {
    // 80% yield, 75% AE, SF=1 → RME = 0.80 × 0.75 = 0.60
    const rme = greenChem.reactionMassEfficiency(0.80, 75, 1.0);
    assert.ok(Math.abs(rme - 0.60) < 1e-6);
  });

  it('reactionMassEfficiency throws on invalid yield', () => {
    assert.throws(() => greenChem.reactionMassEfficiency(1.5, 80));
  });

  it('carbonEfficiency for complete conversion is 100', () => {
    assert.equal(greenChem.carbonEfficiency(9, 9), 100);
  });

  it('carbonEfficiency calculates correctly', () => {
    // Aspirin (9 C) from salicylic acid (7 C) + acetic anhydride (4 C of 2×Ac)
    // CE = 9 / (7 + 4) × 100 = 81.8%
    const ce = greenChem.carbonEfficiency(9, 11);
    assert.ok(Math.abs(ce - 81.82) < 0.1);
  });

  it('classifySolvent: water is recommended', () => {
    const result = greenChem.classifySolvent('water');
    assert.equal(result.score, 'recommended');
    assert.equal(result.concern, null);
  });

  it('classifySolvent: benzene is hazardous', () => {
    const result = greenChem.classifySolvent('benzene');
    assert.equal(result.score, 'hazardous');
    assert.equal(result.concern, 'carcinogen');
  });

  it('classifySolvent: DCM is problematic', () => {
    const result = greenChem.classifySolvent('dcm');
    assert.equal(result.score, 'problematic');
  });

  it('classifySolvent: case-insensitive', () => {
    assert.equal(greenChem.classifySolvent('Ethanol').score, 'recommended');
    assert.equal(greenChem.classifySolvent('BENZENE').score, 'hazardous');
  });

  it('classifySolvent: unknown solvent returns unknown', () => {
    const result = greenChem.classifySolvent('unobtanium');
    assert.equal(result.score, 'unknown');
    assert.equal(result.concern, null);
  });

  it('classifySolvent throws on empty string', () => {
    assert.throws(() => greenChem.classifySolvent(''));
  });

  it('greenScore returns score and grade', () => {
    const result = greenChem.greenScore({
      atomEconomyPct: 75, eFact: 10, yieldPct: 80, solvent: 'ethyl acetate',
    });
    assert.ok(result.score >= 0 && result.score <= 100);
    assert.ok(['A', 'B', 'C', 'D', 'F'].includes(result.grade));
    assert.ok(typeof result.details.aeScore === 'number');
  });

  it('greenScore: ideal reaction scores A', () => {
    const result = greenChem.greenScore({
      atomEconomyPct: 100, eFact: 0, yieldPct: 100, solvent: 'water',
    });
    assert.equal(result.grade, 'A');
    assert.ok(result.score >= 80);
  });

  it('greenScore: poor reaction scores low', () => {
    const result = greenChem.greenScore({
      atomEconomyPct: 10, eFact: 150, yieldPct: 10, solvent: 'benzene',
    });
    assert.ok(result.score < 40);
  });
});
