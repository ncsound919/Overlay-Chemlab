'use strict';

const express = require('express');
const path = require('path');
const apiRoutes = require('./src/api/routes.js');

const app = express();
const PORT = process.env.PORT || 3000;

// CORS headers for development
app.use((_req, res, next) => {
  res.setHeader('Access-Control-Allow-Origin', '*');
  res.setHeader('Access-Control-Allow-Methods', 'GET, POST, OPTIONS');
  res.setHeader('Access-Control-Allow-Headers', 'Content-Type');
  if (_req.method === 'OPTIONS') return res.sendStatus(204);
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
