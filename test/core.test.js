'use strict';

const { describe, it } = require('node:test');
const assert = require('node:assert/strict');

const arrhenius = require('../src/core/arrhenius');
const smilesParser = require('../src/core/smiles-parser');
const molEmbed = require('../src/core/molecular-embeddings');
const bayesian = require('../src/core/bayesian-risk');

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
      assert.ok(data[i].conversion >= data[i - 1].conversion - 1e-9);
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

  it('riskScore returns probability in [0,1]', () => {
    const emb = new Array(128).fill(0).map((_, i) => i % 2 ? 1 : 0);
    const score = bayesian.riskScore(emb, bayesian.DEFAULT_WEIGHTS, bayesian.DEFAULT_BIAS);
    assert.ok(score >= 0 && score <= 1);
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
