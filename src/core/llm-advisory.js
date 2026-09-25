'use strict';

/**
 * Fail-soft Jev + MiniCPM advisories for Chemlab risk calls.
 *
 * classifyRisk stays authoritative (deterministic). These helpers add a
 * calibrated second signal (LocalJev System One) and a plain-language
 * explanation (local MiniCPM, OpenAI-compatible). Both return offline on any
 * failure — never a fabricated recommendation. Zero dependencies (global fetch).
 *
 * Env: JEV_LOCAL_BASE_URL (default http://127.0.0.1:8080), JEV_LOCAL_MODEL
 *      (default localjev-latest), MINICPM_BASE_URL (default
 *      http://127.0.0.1:11434), MINICPM_MODEL (auto-discovered), timeouts via
 *      JEV_TIMEOUT_MS / MINICPM_TIMEOUT_MS.
 */

const JEV_BASE_URL_DEFAULT = 'http://127.0.0.1:8080';
const JEV_MODEL_DEFAULT = 'localjev-latest';
const MINICPM_BASE_URL_DEFAULT = 'http://127.0.0.1:11434';

let cachedMinicpmModel = null;
let minicpmDownUntil = 0;

function baseUrl(value, fallback) {
  return String(value || fallback || '').replace(/\/+$/, '');
}

/**
 * One Jev System One round trip.
 * @param {object} state
 * @param {object} questions
 * @param {{baseUrl?: string, model?: string, timeoutMs?: number, fetchImpl?: Function}} [opts]
 * @returns {Promise<{ok: boolean, source: 'localjev'|'offline', answers?: object, error?: string}>}
 */
async function systemOne(state, questions, opts = {}) {
  const base = baseUrl(opts.baseUrl || process.env.JEV_LOCAL_BASE_URL, JEV_BASE_URL_DEFAULT);
  const model = opts.model || process.env.JEV_LOCAL_MODEL || JEV_MODEL_DEFAULT;
  const timeoutMs = opts.timeoutMs || Number(process.env.JEV_TIMEOUT_MS) || 20000;
  const runFetch = opts.fetchImpl || fetch;
  try {
    const res = await runFetch(`${base}/v1/systemone`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ model, state, questions }),
      signal: AbortSignal.timeout(timeoutMs),
    });
    if (!res.ok) return { ok: false, source: 'offline', error: `JEV /v1/systemone -> HTTP ${res.status}` };
    const data = await res.json();
    if (!data || typeof data.answers !== 'object') {
      return { ok: false, source: 'offline', error: 'JEV response missing answers' };
    }
    return { ok: true, source: 'localjev', answers: data.answers };
  } catch (err) {
    return { ok: false, source: 'offline', error: err instanceof Error ? err.message : String(err) };
  }
}

async function discoverMinicpmModel(base, runFetch) {
  const res = await runFetch(`${base}/v1/models`, { signal: AbortSignal.timeout(10000) });
  if (!res.ok) return null;
  const data = await res.json();
  const first = Array.isArray(data.models) ? data.models[0] : null;
  if (first && typeof first.model === 'string' && first.model.length > 0) return first.model;
  if (first && typeof first.name === 'string' && first.name.length > 0) return first.name;
  return null;
}

/**
 * Local MiniCPM prose completion.
 * @param {{system: string, user: string}} input
 * @param {{baseUrl?: string, model?: string, maxTokens?: number, timeoutMs?: number, fetchImpl?: Function}} [opts]
 * @returns {Promise<{ok: boolean, source: 'local'|'offline', content?: string, error?: string}>}
 */
async function minicpmChat(input, opts = {}) {
  const base = baseUrl(opts.baseUrl || process.env.MINICPM_BASE_URL, MINICPM_BASE_URL_DEFAULT);
  const timeoutMs = opts.timeoutMs || Number(process.env.MINICPM_TIMEOUT_MS) || 60000;
  const runFetch = opts.fetchImpl || fetch;
  let model = opts.model || process.env.MINICPM_MODEL;
  if (!model && cachedMinicpmModel === null && Date.now() < minicpmDownUntil) {
    return { ok: false, source: 'offline', error: 'MiniCPM model discovery failed' };
  }
  if (!model && cachedMinicpmModel === null) {
    try {
      cachedMinicpmModel = await discoverMinicpmModel(base, runFetch);
    } catch {
      cachedMinicpmModel = null;
    }
    minicpmDownUntil = cachedMinicpmModel ? 0 : Date.now() + 30000;
  }
  model = model || cachedMinicpmModel || '';
  if (!model) return { ok: false, source: 'offline', error: 'MiniCPM model discovery failed' };
  try {
    const res = await runFetch(`${base}/v1/chat/completions`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        model,
        max_tokens: opts.maxTokens || 256,
        temperature: 0.2,
        messages: [
          { role: 'system', content: input.system },
          { role: 'user', content: input.user },
        ],
        chat_template_kwargs: { enable_thinking: false },
      }),
      signal: AbortSignal.timeout(timeoutMs),
    });
    if (!res.ok) return { ok: false, source: 'offline', error: `MiniCPM /v1/chat/completions -> HTTP ${res.status}` };
    const data = await res.json();
    const raw = data && data.choices && data.choices[0] && data.choices[0].message ? data.choices[0].message.content : '';
    const content = typeof raw === 'string' ? raw : raw == null ? '' : JSON.stringify(raw);
    if (!content) return { ok: false, source: 'offline', error: 'MiniCPM returned empty content' };
    return { ok: true, source: 'local', content };
  } catch (err) {
    return { ok: false, source: 'offline', error: err instanceof Error ? err.message : String(err) };
  }
}

/**
 * Jev screening advisory for a classifyRisk verdict.
 * @param {{score: number, threshold?: number, level: string}} verdict
 * @param {object} [opts]
 * @returns {Promise<{ok: boolean, source: string, holdForReview?: boolean, noul?: number, riskScore?: number, error?: string}>}
 */
async function riskAdvisory(verdict, opts = {}) {
  const state = {
    action: 'molecule_risk_screening',
    riskScore: Math.round(verdict.score * 10000) / 10000,
    threshold: verdict.threshold != null ? verdict.threshold : 0.65,
    deterministicLevel: verdict.level,
  };
  const questions = {
    hold_for_review: {
      type: 'noul',
      instructions:
        'Should this molecule be flagged for human review before progressing, independent of the deterministic low/medium/high bucket?',
      criteria: { true: 'Flag for review', false: 'No additional review needed' },
    },
    severity: {
      type: 'score',
      instructions: 'Rate the toxicological concern of this molecule.',
      criteria: ['Minimal', 'Mild', 'Moderate', 'Severe'],
    },
  };
  const result = await systemOne(state, questions, opts);
  if (!result.ok) return { ok: false, source: result.source, error: result.error };
  const hold = result.answers.hold_for_review;
  const severity = result.answers.severity;
  const noul = hold && hold.type === 'noul' ? hold.noul : undefined;
  return {
    ok: true,
    source: result.source,
    holdForReview: noul === undefined ? undefined : noul >= 0.5,
    noul,
    riskScore: severity && severity.type === 'score' ? severity.score : undefined,
  };
}

/**
 * MiniCPM plain-language explanation of a risk verdict.
 * @param {{score: number, level: string, smiles?: string}} verdict
 * @param {object} [opts]
 * @returns {Promise<{ok: boolean, source: string, prose?: string, error?: string}>}
 */
async function explainRisk(verdict, opts = {}) {
  const system =
    'You explain a deterministic Bayesian toxicity risk score to a bench chemist. ' +
    'Do not invent specific toxicity endpoints not stated. Do not change the score or its bucket. ' +
    'Under 120 words, plain language, no markdown headers.';
  const user = [
    `Risk score: ${Math.round(verdict.score * 10000) / 10000} (0-1).`,
    `Deterministic bucket: ${verdict.level}.`,
    verdict.smiles ? `SMILES: ${verdict.smiles}` : 'SMILES: not provided.',
    'Explain what this score means for triage and what would raise confidence.',
  ].join('\n');
  const result = await minicpmChat({ system, user }, opts);
  if (!result.ok) return { ok: false, source: result.source, error: result.error };
  return { ok: true, source: result.source, prose: result.content };
}

module.exports = {
  systemOne,
  minicpmChat,
  riskAdvisory,
  explainRisk,
  _resetMinicpmModelCache() {
    cachedMinicpmModel = null;
    minicpmDownUntil = 0;
  },
};