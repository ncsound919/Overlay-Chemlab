'use strict';

const { test } = require('node:test');
const assert = require('node:assert/strict');
const advisory = require('../src/core/llm-advisory');

function jsonFetch(routes) {
  return async (url, init) => {
    const u = String(url);
    const route = routes[u];
    if (!route) return new Response('not found', { status: 404 });
    return new Response(JSON.stringify(route.body), { status: route.status || 200 });
  };
}

test('riskAdvisory returns Jev screening signal', async () => {
  const fetchImpl = jsonFetch({
    'http://127.0.0.1:8080/v1/systemone': {
      body: {
        answers: {
          hold_for_review: { type: 'noul', noul: 0.8 },
          severity: { type: 'score', score: 3 },
        },
      },
    },
  });
  const result = await advisory.riskAdvisory({ score: 0.7, level: 'high' }, { fetchImpl });
  assert.equal(result.ok, true);
  assert.equal(result.source, 'localjev');
  assert.equal(result.holdForReview, true);
  assert.equal(result.noul, 0.8);
  assert.equal(result.riskScore, 3);
});

test('riskAdvisory fails soft offline', async () => {
  const result = await advisory.riskAdvisory({ score: 0.5, level: 'medium' }, { fetchImpl: jsonFetch({}) });
  assert.equal(result.ok, false);
  assert.equal(result.source, 'offline');
  assert.ok(result.error);
});

test('explainRisk returns prose from MiniCPM', async () => {
  advisory._resetMinicpmModelCache();
  const fetchImpl = jsonFetch({
    'http://127.0.0.1:11434/v1/models': { body: { models: [{ model: 'minicpm5.gguf' }] } },
    'http://127.0.0.1:11434/v1/chat/completions': {
      body: { choices: [{ message: { content: 'This score is moderate; review before advancing.' } }] },
    },
  });
  const result = await advisory.explainRisk({ score: 0.5, level: 'medium' }, { fetchImpl });
  assert.equal(result.ok, true);
  assert.equal(result.source, 'local');
  assert.match(result.prose, /moderate/);
});

test('explainRisk fails soft when MiniCPM is offline', async () => {
  advisory._resetMinicpmModelCache();
  const result = await advisory.explainRisk({ score: 0.5, level: 'medium' }, { fetchImpl: jsonFetch({}) });
  assert.equal(result.ok, false);
  assert.equal(result.source, 'offline');
});

test('deterministic classifyRisk is unchanged by advisory module', () => {
  const bayesian = require('../src/core/bayesian-risk');
  assert.equal(bayesian.classifyRisk(0.1).level, 'low');
  assert.equal(bayesian.classifyRisk(0.5).level, 'medium');
  assert.equal(bayesian.classifyRisk(0.8).level, 'high');
});