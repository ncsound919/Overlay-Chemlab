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

const MAX_SMILES_LENGTH = 1000;
const MAX_SIM_STEPS = 10000;

function validateSmiles(smiles) {
  if (!smiles) return 'smiles query parameter is required';
  if (smiles.length > MAX_SMILES_LENGTH) return `smiles exceeds maximum length of ${MAX_SMILES_LENGTH} characters`;
  return null;
}

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
    const err = validateSmiles(smiles);
    if (err) return res.status(400).json({ error: err });

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
    const err = validateSmiles(smiles);
    if (err) return res.status(400).json({ error: err });

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
    const err1 = validateSmiles(smiles1);
    if (err1) return res.status(400).json({ error: err1 });
    const err2 = validateSmiles(smiles2);
    if (err2) return res.status(400).json({ error: err2 });

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

    const stepsVal = steps != null ? Number(steps) : 200;
    if (stepsVal > MAX_SIM_STEPS) {
      return res.status(400).json({ error: `steps exceeds maximum of ${MAX_SIM_STEPS}` });
    }

    const result = runKineticsODE({
      temperature: Number(temperature),
      totalTime: Number(totalTime),
      steps: stepsVal,
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

// POST /api/chem-assistant — server-side chemistry knowledge base
const CHEM_KNOWLEDGE = {
  'h3po4': '**H₃PO₄ vs H₂SO₄**: Phosphoric acid is preferred because it\'s a weaker acid (pKₐ₁ = 2.15 vs −3 for H₂SO₄), reducing charring and side-product formation. H₂SO₄ can sulphonate the aromatic ring at >100°C, yielding unwanted aryl sulfonic acids. H₃PO₄ provides sufficient catalytic protons without degrading the salicylic acid substrate.\n[chips: "What is the mechanism?" | "Compare yields with each catalyst" | "Temperature effect on selectivity"]',
  'green': '**Green alternatives for aspirin synthesis**: Replace acetic anhydride with acetic acid (lower atom economy but less corrosive). Use **microwave-assisted** synthesis (2 min vs 15 min, 90% yield). Consider **solvent-free mechanochemical** grinding. Water is already the quench/recrystallization solvent — E-factor ≈ 4.3. Biocatalytic routes using lipases (CAL-B) show promise at 40°C.\n[chips: "Microwave conditions?" | "Calculate E-factor" | "Enzyme catalysis details"]',
  'nmr': '**¹H NMR of Aspirin** (CDCl₃, 400 MHz):\n• **δ 2.36** (s, 3H, OCOCH₃) — acetyl methyl\n• **δ 7.12** (dd, 1H, H-3, J=8.1, 1.1 Hz)\n• **δ 7.33** (td, 1H, H-5, J=7.6, 1.2 Hz)\n• **δ 7.61** (td, 1H, H-4, J=7.8, 1.8 Hz)\n• **δ 8.11** (dd, 1H, H-6, J=7.8, 1.8 Hz)\n• **δ 11.0** (br s, 1H, COOH — exchangeable)\nThe methyl singlet at δ 2.36 confirms acetylation vs salicylic acid (no methyl).\n[chips: "¹³C NMR shifts?" | "How to confirm purity by NMR?" | "Compare with salicylic acid NMR"]',
  'side products': '**Side products above 95°C**:\n1. **Acetylsalicylic anhydride** — over-acetylation at high [Ac₂O]\n2. **Salicylic acid dimer** — dehydration condensation\n3. **Polymeric tar** — charring from acid-catalyzed degradation >110°C\n4. **Acetic acid** (expected co-product, but excess at high T)\nAt 85°C: ~3% side products. At 100°C: ~8%. At 120°C: >15% with visible darkening.\n[chips: "How to minimize side products?" | "TLC monitoring protocol" | "Purification methods"]',
  'mechanism': '**Fischer esterification mechanism** (acid-catalyzed):\n1. Protonation of Ac₂O carbonyl oxygen by H₃PO₄\n2. Nucleophilic attack by phenolic −OH of salicylic acid\n3. Tetrahedral intermediate forms\n4. Proton transfer\n5. Loss of acetic acid (AcOH) as leaving group\n6. Deprotonation yields aspirin\nRate-determining step: nucleophilic addition (step 2). Ea ≈ 65.3 kJ/mol.\n[chips: "Draw the transition state" | "Why is step 2 rate-limiting?" | "Compare with Schotten-Baumann"]',
  'yield': '**Optimizing yield**: Current conditions give ~82% yield. To improve:\n• **Increase Ac₂O to 2.0 equiv** (drives equilibrium, ~87% yield)\n• **Extend time to 20 min** at 85°C (+3% conversion)\n• **Use dry glassware** — moisture hydrolyzes Ac₂O\n• **Recrystallize from ethanol/water** (3:1) instead of pure water\n• **Theoretical max**: ~92% (limited by crystallization losses)\n[chips: "Calculate theoretical yield" | "Effect of excess reagent" | "Recrystallization tips"]',
  'melting point': '**Melting point analysis**: Pure aspirin mp = **135–136°C**. If your crystals melt at 128–132°C, impurities (salicylic acid, mp 159°C) are present — use mixed mp test. Broad melting range (>2°C) indicates need for recrystallization. DSC shows sharp endotherm at 141°C (decomposition begins).\n[chips: "How to do mixed melting point?" | "DSC vs mp apparatus" | "Common impurities"]',
  'ir': '**IR spectrum of aspirin**:\n• **1754 cm⁻¹** — C=O ester stretch (confirms acetylation)\n• **1689 cm⁻¹** — C=O carboxylic acid\n• **2500–3300 cm⁻¹** — broad O−H stretch (COOH)\n• **1185 cm⁻¹** — C−O ester stretch\n• **No broad 3200–3550 cm⁻¹** peak = no free phenolic OH (salicylic acid gone)\n[chips: "Compare with salicylic acid IR" | "How to identify ester vs acid C=O" | "Sample preparation for IR"]',
  'solubility': '**Aspirin solubility**: 3.3 g/L in water at 20°C, 10 g/L at 37°C. Freely soluble in ethanol (200 g/L), acetone, chloroform. The low water solubility enables precipitation upon quenching with ice water. pKₐ = 3.49 means at stomach pH (~2), aspirin is mostly unionized → crosses gastric membrane → GI absorption.\n[chips: "Why does it dissolve in blood?" | "Buffered aspirin formulation" | "Henderson-Hasselbalch calculation"]',
  'default': 'I can help with aspirin synthesis analysis. Try asking about reaction mechanisms, NMR/IR spectra, optimizing yield, green chemistry alternatives, or side product formation.\n[chips: "Explain the mechanism" | "Predict NMR shifts" | "Green alternatives"]'
};

function getChemAnswer(msg) {
  const q = msg.toLowerCase();
  if (q.includes('h3po4') || q.includes('h₃po₄') || q.includes('phosphoric') || q.includes('h2so4') || q.includes('catalyst')) return CHEM_KNOWLEDGE['h3po4'];
  if (q.includes('green') || q.includes('solvent') || q.includes('sustainable') || q.includes('eco')) return CHEM_KNOWLEDGE['green'];
  if (q.includes('nmr') || q.includes('shift') || q.includes('spectrum') || q.includes('spectra')) return CHEM_KNOWLEDGE['nmr'];
  if (q.includes('side product') || q.includes('temperature') || q.includes('above 95') || q.includes('high-temp') || q.includes('byproduct')) return CHEM_KNOWLEDGE['side products'];
  if (q.includes('mechanism') || q.includes('how does') || q.includes('step')) return CHEM_KNOWLEDGE['mechanism'];
  if (q.includes('yield') || q.includes('optim') || q.includes('improve') || q.includes('condition')) return CHEM_KNOWLEDGE['yield'];
  if (q.includes('melting') || q.includes('purity') || q.includes('mp')) return CHEM_KNOWLEDGE['melting point'];
  if (q.includes('ir') || q.includes('infrared') || q.includes('1754') || q.includes('stretch')) return CHEM_KNOWLEDGE['ir'];
  if (q.includes('solub') || q.includes('dissolv') || q.includes('pka') || q.includes('water')) return CHEM_KNOWLEDGE['solubility'];
  return CHEM_KNOWLEDGE['default'];
}

router.post('/chem-assistant', (req, res) => {
  try {
    const { message } = req.body;
    if (!message || typeof message !== 'string') {
      return res.status(400).json({ error: 'message is required in request body' });
    }
    const answer = getChemAnswer(message);
    res.json({ answer });
  } catch (err) {
    res.status(500).json({ error: err.message });
  }
});

module.exports = router;
