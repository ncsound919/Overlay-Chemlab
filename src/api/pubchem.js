'use strict';

const https = require('https');

const DEFAULT_TTL = 5 * 60 * 1000; // 5 minutes
const MAX_CACHE_SIZE = 1000;
const REQUEST_TIMEOUT_MS = 15000; // 15 seconds

const cache = new Map();

function cacheGet(key) {
  const entry = cache.get(key);
  if (!entry) return undefined;
  if (Date.now() - entry.ts > entry.ttl) {
    cache.delete(key);
    return undefined;
  }
  // LRU: move to end by reinserting
  cache.delete(key);
  cache.set(key, entry);
  return entry.value;
}

function cacheSet(key, value, ttl = DEFAULT_TTL) {
  if (!cache.has(key) && cache.size >= MAX_CACHE_SIZE) {
    const oldestKey = cache.keys().next().value;
    if (oldestKey !== undefined) cache.delete(oldestKey);
  }
  cache.set(key, { value, ts: Date.now(), ttl });
}

/**
 * Make an HTTPS GET request and return parsed JSON.
 */
function httpsGetJSON(url) {
  return new Promise((resolve, reject) => {
    const req = https.get(url, (res) => {
      const chunks = [];
      res.on('data', (chunk) => chunks.push(chunk));
      res.on('end', () => {
        const body = Buffer.concat(chunks).toString();
        if (res.statusCode < 200 || res.statusCode >= 300) {
          const err = new Error(`PubChem API returned ${res.statusCode}`);
          err.statusCode = res.statusCode;
          err.body = body;
          return reject(err);
        }
        try {
          resolve(JSON.parse(body));
        } catch (e) {
          reject(new Error('Invalid JSON from PubChem'));
        }
      });
      res.on('error', reject);
    });
    req.on('error', reject);
    req.setTimeout(REQUEST_TIMEOUT_MS, () => {
      req.destroy(new Error('PubChem request timed out'));
    });
  });
}

const PROPERTY_LIST = [
  'MolecularFormula', 'MolecularWeight', 'IUPACName', 'InChIKey',
  'CanonicalSMILES', 'XLogP', 'TPSA', 'HBondDonorCount',
  'HBondAcceptorCount', 'RotatableBondCount',
].join(',');

/**
 * Fetch compound data by SMILES from PubChem.
 */
async function fetchBySMILES(smiles) {
  if (!smiles || typeof smiles !== 'string') {
    throw new Error('SMILES string is required');
  }

  const cacheKey = `smiles:${smiles}`;
  const cached = cacheGet(cacheKey);
  if (cached) return cached;

  const encoded = encodeURIComponent(smiles);
  const cidData = await httpsGetJSON(
    `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${encoded}/cids/JSON`
  );

  const cid = cidData.IdentifierList && cidData.IdentifierList.CIDs && cidData.IdentifierList.CIDs[0];
  if (!cid) throw new Error('No CID found for the given SMILES');

  const propData = await httpsGetJSON(
    `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/property/${PROPERTY_LIST}/JSON`
  );

  const props = propData.PropertyTable && propData.PropertyTable.Properties && propData.PropertyTable.Properties[0];
  const result = { cid, ...(props || {}) };

  cacheSet(cacheKey, result);
  return result;
}

/**
 * Fetch compound data by common name from PubChem.
 */
async function fetchByName(name) {
  if (!name || typeof name !== 'string') {
    throw new Error('Compound name is required');
  }

  const cacheKey = `name:${name}`;
  const cached = cacheGet(cacheKey);
  if (cached) return cached;

  const encoded = encodeURIComponent(name);
  const cidData = await httpsGetJSON(
    `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${encoded}/cids/JSON`
  );

  const cid = cidData.IdentifierList && cidData.IdentifierList.CIDs && cidData.IdentifierList.CIDs[0];
  if (!cid) throw new Error('No CID found for the given name');

  const propData = await httpsGetJSON(
    `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/property/${PROPERTY_LIST}/JSON`
  );

  const props = propData.PropertyTable && propData.PropertyTable.Properties && propData.PropertyTable.Properties[0];
  const result = { cid, ...(props || {}) };

  cacheSet(cacheKey, result);
  return result;
}

/**
 * Get autocomplete suggestions from PubChem.
 */
async function autocomplete(query, limit = 10) {
  if (!query || typeof query !== 'string') {
    throw new Error('Query string is required');
  }

  const cacheKey = `ac:${query}:${limit}`;
  const cached = cacheGet(cacheKey);
  if (cached) return cached;

  const encoded = encodeURIComponent(query);
  const data = await httpsGetJSON(
    `https://pubchem.ncbi.nlm.nih.gov/rest/autocomplete/compound/${encoded}/JSON?limit=${limit}`
  );

  const suggestions = (data.dictionary_terms && data.dictionary_terms.compound) || [];
  cacheSet(cacheKey, suggestions);
  return suggestions;
}

module.exports = {
  fetchBySMILES,
  fetchByName,
  autocomplete,
};
