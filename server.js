'use strict';

const express = require('express');
const path = require('path');
const apiRoutes = require('./src/api/routes.js');

const app = express();
const PORT = process.env.PORT || 3000;

// Configure CORS origins via CORS_ALLOWED_ORIGINS env var in production.
// WARNING: wildcard origin is used in development only.
const ALLOWED_ORIGINS = (process.env.CORS_ALLOWED_ORIGINS || '')
  .split(',').map(o => o.trim()).filter(Boolean);

app.use((req, res, next) => {
  const origin = req.headers.origin;
  if (process.env.NODE_ENV === 'production' && ALLOWED_ORIGINS.length > 0) {
    if (origin && ALLOWED_ORIGINS.includes(origin)) {
      res.setHeader('Access-Control-Allow-Origin', origin);
    }
  } else {
    res.setHeader('Access-Control-Allow-Origin', '*');
  }
  res.setHeader('Access-Control-Allow-Methods', 'GET, POST, OPTIONS');
  res.setHeader('Access-Control-Allow-Headers', 'Content-Type');
  if (req.method === 'OPTIONS') return res.sendStatus(204);
  next();
});

// JSON body parsing
app.use(express.json());

// Serve static files from public/
app.use(express.static(path.join(__dirname, 'public')));

// Mount API routes
app.use('/api', apiRoutes);

// Start server when run directly
if (require.main === module) {
  app.listen(PORT, () => {
    console.log(`Overlay-Chemlab server running on http://localhost:${PORT}`);
  });
}

module.exports = app;
