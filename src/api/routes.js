'use strict';

const { Router } = require('express');
const { molecularWeight, morganFingerprint, tanimotoSimilarity } = require('../core/molecular-embeddings.js');
const { molecularFormula, atomCount, bondCount } = require('../core/smiles-parser.js');
const { riskScore, classifyRisk, DEFAULT_WEIGHTS, DEFAULT_BIAS } = require('../core/bayesian-risk.js');
const { simulateEsterification, runKineticsODE } = require('../core/arrhenius.js');
const { parseReaction } = require('../core/reaction-parser.js');
const pubchem = require('./pubchem.js');
const { version } = require('../../package.json');

const router = Router();
const startTime = Date.now();

// GET /api/health
router.get('/health', (_req, res) => {
  res.json({
    status: 'ok',
    version,
    uptime: Math.floor((Date.now() - startTime) / 1000),
  });
});

// GET /api/molecule/properties?smiles=...
router.get('/molecule/properties', (req, res) => {
  try {
    const { smiles } = req.query;
    if (!smiles) return res.status(400).json({ error: 'smiles query parameter is required' });

    const mw = molecularWeight(smiles);
    const formula = molecularFormula(smiles);
    const atoms = atomCount(smiles);
    const bonds = bondCount(smiles);

    res.json({ smiles, molecularWeight: mw, formula, atoms, bonds });
  } catch (err) {
    res.status(400).json({ error: err.message });
  }
});

// GET /api/molecule/fingerprint?smiles=...&radius=2&nBits=128
router.get('/molecule/fingerprint', (req, res) => {
  try {
    const { smiles } = req.query;
    if (!smiles) return res.status(400).json({ error: 'smiles query parameter is required' });

    const radius = parseInt(req.query.radius, 10) || 2;
    const nBits = parseInt(req.query.nBits, 10) || 128;

    const fp = morganFingerprint(smiles, { radius, nBits });
    res.json({ smiles, radius, nBits, fingerprint: Array.from(fp) });
  } catch (err) {
    res.status(400).json({ error: err.message });
  }
});

// GET /api/molecule/similarity?smiles1=...&smiles2=...
router.get('/molecule/similarity', (req, res) => {
  try {
    const { smiles1, smiles2 } = req.query;
    if (!smiles1 || !smiles2) {
      return res.status(400).json({ error: 'Both smiles1 and smiles2 query parameters are required' });
    }

    const fp1 = morganFingerprint(smiles1);
    const fp2 = morganFingerprint(smiles2);
    const similarity = tanimotoSimilarity(fp1, fp2);

    res.json({ smiles1, smiles2, tanimoto: similarity });
  } catch (err) {
    res.status(400).json({ error: err.message });
  }
});

// POST /api/molecule/risk
router.post('/molecule/risk', (req, res) => {
  try {
    const { smiles, embedding } = req.body;
    if (!smiles && !embedding) {
      return res.status(400).json({ error: 'Either smiles or embedding is required in request body' });
    }

    let emb = embedding;
    if (!emb) {
      const fp = morganFingerprint(smiles);
      emb = Array.from(fp);
    }

    const score = riskScore(emb, DEFAULT_WEIGHTS, DEFAULT_BIAS);
    const classification = classifyRisk(score);

    res.json({ score, ...classification });
  } catch (err) {
    res.status(400).json({ error: err.message });
  }
});

// POST /api/simulate/esterification
router.post('/simulate/esterification', (req, res) => {
  try {
    const { temperature, time, equiv, catalystDrops } = req.body;
    if (temperature == null || time == null) {
      return res.status(400).json({ error: 'temperature and time are required' });
    }

    const result = simulateEsterification({
      temperature: Number(temperature),
      time: Number(time),
      equiv: equiv != null ? Number(equiv) : 1.0,
      catalystDrops: catalystDrops != null ? Number(catalystDrops) : 3,
    });

    res.json(result);
  } catch (err) {
    res.status(400).json({ error: err.message });
  }
});

// POST /api/simulate/kinetics
router.post('/simulate/kinetics', (req, res) => {
  try {
    const { temperature, totalTime, steps, equiv, catalystDrops } = req.body;
    if (temperature == null || totalTime == null) {
      return res.status(400).json({ error: 'temperature and totalTime are required' });
    }

    const result = runKineticsODE({
      temperature: Number(temperature),
      totalTime: Number(totalTime),
      steps: steps != null ? Number(steps) : 200,
      equiv: equiv != null ? Number(equiv) : 1.0,
      catalystDrops: catalystDrops != null ? Number(catalystDrops) : 3,
    });

    res.json(result);
  } catch (err) {
    res.status(400).json({ error: err.message });
  }
});

// POST /api/reaction/parse
router.post('/reaction/parse', (req, res) => {
  try {
    const { reaction } = req.body;
    if (!reaction) {
      return res.status(400).json({ error: 'reaction string is required in request body' });
    }

    const result = parseReaction(reaction);
    res.json(result);
  } catch (err) {
    res.status(400).json({ error: err.message });
  }
});

// GET /api/pubchem/compound?smiles=...
router.get('/pubchem/compound', async (req, res) => {
  try {
    const { smiles, name } = req.query;
    if (!smiles && !name) {
      return res.status(400).json({ error: 'smiles or name query parameter is required' });
    }

    const result = smiles
      ? await pubchem.fetchBySMILES(smiles)
      : await pubchem.fetchByName(name);

    res.json(result);
  } catch (err) {
    const status = err.statusCode || 400;
    res.status(status).json({ error: err.message });
  }
});

// GET /api/pubchem/autocomplete?query=...
router.get('/pubchem/autocomplete', async (req, res) => {
  try {
    const { query, limit } = req.query;
    if (!query) {
      return res.status(400).json({ error: 'query parameter is required' });
    }

    const suggestions = await pubchem.autocomplete(query, limit ? parseInt(limit, 10) : 10);
    res.json({ suggestions });
  } catch (err) {
    const status = err.statusCode || 400;
    res.status(status).json({ error: err.message });
  }
});

module.exports = router;
