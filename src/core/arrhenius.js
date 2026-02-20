'use strict';

const R = 8.314;           // J/(mol·K)
const DEFAULT_Ea = 65300;  // J/mol — Fischer esterification activation energy
const DEFAULT_A = 4.2e9;   // Pre-exponential factor (s⁻¹)

/**
 * Calculate rate constant k from the Arrhenius equation: k = A * exp(-Ea / (R * T))
 * @param {number} Ea  - Activation energy in J/mol
 * @param {number} A   - Pre-exponential factor (s⁻¹)
 * @param {number} T_celsius - Temperature in °C
 * @returns {number} Rate constant k (s⁻¹)
 */
function arrheniusRate(Ea, A, T_celsius) {
  const T = T_celsius + 273.15;
  if (T <= 0) return 0;
  return A * Math.exp(-Ea / (R * T));
}

/**
 * Simulate Fischer-Speier esterification.
 * Alcohol + carboxylic acid  ⇌  ester + water  (acid-catalysed)
 *
 * @param {object} params
 * @param {number} params.temperature - °C
 * @param {number} params.time        - reaction time in seconds
 * @param {number} params.equiv       - molar equivalents of alcohol to acid (≥1)
 * @param {number} params.catalystDrops - drops of H₂SO₄ catalyst (0–10)
 * @param {number} [params.Ea]        - activation energy (J/mol)
 * @param {number} [params.A]         - pre-exponential factor
 * @returns {{ conversion: number, yield: number, pH: number, sideProducts: string[], rateConstant: number }}
 */
function simulateEsterification({
  temperature,
  time,
  equiv = 1.0,
  catalystDrops = 3,
  Ea = DEFAULT_Ea,
  A = DEFAULT_A,
} = {}) {
  const k = arrheniusRate(Ea, A, temperature);

  // Catalyst effectiveness factor (each drop roughly doubles rate, diminishing)
  const catalystFactor = 1 + 0.8 * catalystDrops / (1 + 0.15 * catalystDrops);

  const kEff = k * catalystFactor;

  // Second-order reversible kinetics simplified to pseudo-first-order
  // with excess alcohol acting as solvent. Equilibrium constant Keq ≈ 4 for
  // Fischer esterification. conversion = Keq*equiv / (Keq*equiv + 1) at equilibrium.
  const Keq = 4.0;
  const eqConversion = (Keq * equiv) / (Keq * equiv + 1);

  // Approach to equilibrium: x(t) = x_eq * (1 - exp(-kEff * t * (1 + 1/Keq)))
  const expTerm = Math.exp(-kEff * time * (1 + 1 / Keq));
  const conversion = eqConversion * (1 - expTerm);

  // Side products increase at high temperature
  const sideProducts = [];
  if (temperature > 140) sideProducts.push('diethyl_ether');
  if (temperature > 160) sideProducts.push('alkene_elimination');
  if (temperature > 180) sideProducts.push('charring');

  // Side-product loss
  const sideLoss = sideProducts.length * 0.03;

  const rawYield = Math.max(0, conversion * (1 - sideLoss));

  // pH estimate: strong acid catalyst diluted into reaction mixture
  // More catalyst → lower pH. Base line ~2–4.
  const pH = catalystDrops > 0
    ? Math.max(0.5, 4.0 - Math.log10(catalystDrops + 1) * 2)
    : 7.0;

  return {
    conversion: Math.min(1, Math.max(0, conversion)),
    yield: Math.min(1, Math.max(0, rawYield)),
    pH: Math.round(pH * 100) / 100,
    sideProducts,
    rateConstant: kEff,
  };
}

/**
 * Run an ODE-style kinetics simulation returning time-series data.
 * Uses Euler integration of  dx/dt = kEff*(1-x)*(equiv - x) - (kEff/Keq)*x^2
 *
 * @param {object} params
 * @param {number} params.temperature   - °C
 * @param {number} params.totalTime     - total time in seconds
 * @param {number} [params.steps=200]   - number of integration steps
 * @param {number} [params.Ea]          - activation energy (J/mol)
 * @param {number} [params.A]           - pre-exponential factor
 * @param {number} [params.equiv=1]     - molar equivalents
 * @param {number} [params.catalystDrops=3]
 * @returns {Array<{time: number, conversion: number, concentration: number}>}
 */
function runKineticsODE({
  temperature,
  totalTime,
  steps = 200,
  Ea = DEFAULT_Ea,
  A = DEFAULT_A,
  equiv = 1.0,
  catalystDrops = 3,
} = {}) {
  const k = arrheniusRate(Ea, A, temperature);
  const catalystFactor = 1 + 0.8 * catalystDrops / (1 + 0.15 * catalystDrops);
  const kEff = k * catalystFactor;
  const Keq = 4.0;

  const dt = totalTime / steps;
  let x = 0; // conversion fraction

  const data = [];

  // RK4 derivative
  function dxdt(xVal) {
    const forward = kEff * (1 - xVal) * (equiv - xVal);
    const reverse = (kEff / Keq) * xVal * xVal;
    return forward - reverse;
  }

  for (let i = 0; i <= steps; i++) {
    const t = i * dt;
    const conc = 1 - x; // normalised reactant concentration
    data.push({
      time: Math.round(t * 1e6) / 1e6,
      conversion: Math.round(x * 1e6) / 1e6,
      concentration: Math.round(conc * 1e6) / 1e6,
    });

    // Adaptive sub-stepping with RK4 to handle stiff kinetics
    const maxStep = 0.5 / (kEff * (equiv + 1) + 1e-30);
    const nSub = Math.max(1, Math.ceil(dt / maxStep));
    const h = dt / nSub;

    for (let s = 0; s < nSub; s++) {
      const k1 = dxdt(x);
      const k2 = dxdt(x + 0.5 * h * k1);
      const k3 = dxdt(x + 0.5 * h * k2);
      const k4 = dxdt(x + h * k3);
      x += (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
      x = Math.min(1, Math.max(0, x));
    }
  }

  return data;
}

/**
 * Calculate final isolated yield accounting for work-up losses.
 *
 * @param {object} params
 * @param {number} params.conversion          - fractional conversion (0–1)
 * @param {number} [params.equiv=1]           - molar equivalents
 * @param {number} [params.crystallizationLoss=0.05] - fractional loss during crystallization
 * @param {number} [params.filtrationLoss=0.02]      - fractional loss during filtration
 * @returns {{ theoreticalYield: number, isolatedYield: number, loss: number }}
 */
function calculateYield({
  conversion,
  equiv = 1.0,
  crystallizationLoss = 0.05,
  filtrationLoss = 0.02,
} = {}) {
  // Theoretical yield is limited by the limiting reagent
  const limitingFactor = Math.min(1, equiv);
  const theoreticalYield = conversion * limitingFactor;

  const totalLoss = crystallizationLoss + filtrationLoss;
  const isolatedYield = Math.max(0, theoreticalYield * (1 - totalLoss));

  return {
    theoreticalYield: Math.round(theoreticalYield * 1e6) / 1e6,
    isolatedYield: Math.round(isolatedYield * 1e6) / 1e6,
    loss: Math.round(totalLoss * 1e6) / 1e6,
  };
}

module.exports = {
  R,
  DEFAULT_Ea,
  DEFAULT_A,
  arrheniusRate,
  simulateEsterification,
  runKineticsODE,
  calculateYield,
};
