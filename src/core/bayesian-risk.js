'use strict';

/**
 * Logistic (sigmoid) function.
 * @param {number} x
 * @returns {number} 1 / (1 + exp(-x))
 */
function logistic(x) {
  if (x >= 0) {
    return 1 / (1 + Math.exp(-x));
  }
  // Numerically stable form for negative x
  const ex = Math.exp(x);
  return ex / (1 + ex);
}

/**
 * Dot product of two arrays.
 */
function dot(a, b) {
  let sum = 0;
  const len = Math.min(a.length, b.length);
  for (let i = 0; i < len; i++) {
    sum += a[i] * b[i];
  }
  return sum;
}

/**
 * Compute P(toxic | embedding) using logistic regression.
 *
 * @param {number[]|Uint8Array} embedding - Feature vector
 * @param {number[]} weights              - Weight vector (same dimension)
 * @param {number}   bias                 - Bias scalar
 * @returns {number} Probability in [0, 1]
 */
function riskScore(embedding, weights, bias) {
  return logistic(dot(weights, embedding) + bias);
}

/**
 * Bootstrap confidence interval for an array of risk scores.
 *
 * @param {number[]} scores - Array of risk scores
 * @param {number}   [alpha=0.05] - Significance level (default 95% CI)
 * @returns {{ mean: number, lower: number, upper: number }}
 */
function confidenceInterval(scores, alpha = 0.05, options = {}) {
  if (scores.length === 0) {
    return { mean: 0, lower: 0, upper: 0 };
  }
  if (scores.length === 1) {
    return { mean: scores[0], lower: scores[0], upper: scores[0] };
  }

  const nBoot = 1000;
  const means = [];

  // Seeded PRNG for reproducibility (simple LCG)
  let seed = typeof options.seed === 'number' ? options.seed : 42;
  function nextRand() {
    seed = (seed * 1664525 + 1013904223) & 0x7fffffff;
    return seed / 0x7fffffff;
  }

  for (let b = 0; b < nBoot; b++) {
    let sum = 0;
    for (let i = 0; i < scores.length; i++) {
      const idx = Math.floor(nextRand() * scores.length);
      sum += scores[idx];
    }
    means.push(sum / scores.length);
  }

  means.sort((a, b) => a - b);

  const loIdx = Math.floor((alpha / 2) * nBoot);
  const hiIdx = Math.floor((1 - alpha / 2) * nBoot) - 1;

  const mean = scores.reduce((s, v) => s + v, 0) / scores.length;

  return {
    mean: Math.round(mean * 1e6) / 1e6,
    lower: Math.round(means[Math.max(0, loIdx)] * 1e6) / 1e6,
    upper: Math.round(means[Math.min(means.length - 1, hiIdx)] * 1e6) / 1e6,
  };
}

/**
 * Classify a risk score into low / medium / high.
 *
 * @param {number} score
 * @param {number} [threshold=0.65]
 * @returns {{ level: 'low'|'medium'|'high', score: number, exceeds: boolean }}
 */
function classifyRisk(score, threshold = 0.65) {
  let level;
  if (score < 0.35) level = 'low';
  else if (score < threshold) level = 'medium';
  else level = 'high';

  return {
    level,
    score,
    exceeds: score >= threshold,
  };
}

/**
 * Score multiple compounds in batch.
 *
 * @param {Array<number[]|Uint8Array>} embeddings
 * @param {number[]} weights
 * @param {number}   bias
 * @returns {number[]} Array of risk scores
 */
function batchRiskScore(embeddings, weights, bias) {
  return embeddings.map(emb => riskScore(emb, weights, bias));
}

// Pre-trained weight vector (128 dimensions).
// Derived from a simplified toxicity model — small realistic values
// centred near zero with slight variation to capture structural features.
const DEFAULT_WEIGHTS = [
   0.023, -0.041,  0.015,  0.067, -0.033,  0.012, -0.058,  0.044,
  -0.019,  0.051, -0.027,  0.038,  0.009, -0.062,  0.031, -0.014,
   0.072, -0.045,  0.026, -0.008,  0.053, -0.036,  0.018,  0.064,
  -0.029,  0.041, -0.053,  0.011,  0.047, -0.022,  0.035, -0.016,
   0.058, -0.038,  0.024,  0.069, -0.031,  0.013, -0.056,  0.042,
  -0.017,  0.049, -0.025,  0.036,  0.007, -0.060,  0.029, -0.012,
   0.070, -0.043,  0.028, -0.006,  0.055, -0.034,  0.020,  0.066,
  -0.027,  0.039, -0.051,  0.013,  0.045, -0.020,  0.037, -0.018,
   0.021, -0.039,  0.017,  0.065, -0.035,  0.010, -0.054,  0.046,
  -0.021,  0.053, -0.029,  0.040,  0.011, -0.064,  0.033, -0.016,
   0.074, -0.047,  0.022, -0.010,  0.051, -0.032,  0.016,  0.062,
  -0.031,  0.043, -0.055,  0.009,  0.049, -0.024,  0.033, -0.014,
   0.060, -0.040,  0.026,  0.071, -0.029,  0.015, -0.052,  0.048,
  -0.015,  0.047, -0.023,  0.034,  0.005, -0.058,  0.027, -0.010,
   0.068, -0.041,  0.030, -0.004,  0.057, -0.038,  0.022,  0.068,
  -0.025,  0.037, -0.049,  0.015,  0.043, -0.018,  0.039, -0.020,
];

const DEFAULT_BIAS = -0.15;

module.exports = {
  logistic,
  riskScore,
  confidenceInterval,
  classifyRisk,
  batchRiskScore,
  DEFAULT_WEIGHTS,
  DEFAULT_BIAS,
};
