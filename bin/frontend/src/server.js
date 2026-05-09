// server.js — SnakeLauncher backend
// Run with: node server.js
// Dependencies: npm install express cors js-yaml

import express from 'express'
import cors from 'cors'
import fs from 'fs'
import path from 'path'
import { spawn } from 'child_process'
import { parse as parseYaml, dump as dumpYaml } from 'js-yaml'

const app = express()
const PORT = 3001

// ── CONFIG ────────────────────────────────────────────────────────────────────
// Change this to the folder where your .yaml config files live
let CONFIGS_DIR = process.env.CONFIGS_DIR || path.join(process.env.HOME, 'snakemake-configs')

// The snakemake command template.
// {configFile} is replaced with the full path to the chosen config.
const SNAKEMAKE_CMD = process.env.SNAKEMAKE_CMD || 'snakemake'
const SNAKEMAKE_ARGS = ['--configfile', '{configFile}', '--cores', 'all']
// ─────────────────────────────────────────────────────────────────────────────

app.use(cors())
app.use(express.json())

// POST /set-dir — update the configs directory at runtime
app.post('/set-dir', (req, res) => {
  const { dir } = req.body
  if (!dir) return res.status(400).json({ error: 'dir required' })
  const expanded = dir.replace(/^~/, process.env.HOME)
  if (!fs.existsSync(expanded)) {
    try { fs.mkdirSync(expanded, { recursive: true }) } catch (e) {
      return res.status(400).json({ error: `Cannot create directory: ${e.message}` })
    }
  }
  CONFIGS_DIR = expanded
  res.json({ dir: CONFIGS_DIR })
})

// GET /configs — list all .yaml/.yml files in CONFIGS_DIR
app.get('/configs', (req, res) => {
  try {
    if (!fs.existsSync(CONFIGS_DIR)) {
      fs.mkdirSync(CONFIGS_DIR, { recursive: true })
    }
    const files = fs.readdirSync(CONFIGS_DIR)
      .filter(f => f.endsWith('.yaml') || f.endsWith('.yml'))
      .map(f => {
        const fullPath = path.join(CONFIGS_DIR, f)
        const stat = fs.statSync(fullPath)
        return { name: f, size: stat.size, mtime: stat.mtime }
      })
      .sort((a, b) => new Date(b.mtime) - new Date(a.mtime))
    res.json({ configs: files, dir: CONFIGS_DIR })
  } catch (err) {
    res.status(500).json({ error: err.message })
  }
})

// GET /configs/:filename — read and parse a config file
app.get('/configs/:filename', (req, res) => {
  try {
    const safeName = path.basename(req.params.filename)
    const fullPath = path.join(CONFIGS_DIR, safeName)
    if (!fs.existsSync(fullPath)) return res.status(404).json({ error: 'Not found' })
    const raw = fs.readFileSync(fullPath, 'utf8')
    const parsed = parseYaml(raw)
    res.json({ name: safeName, raw, parsed })
  } catch (err) {
    res.status(500).json({ error: err.message })
  }
})

// POST /configs — save a new config file
// Body: { filename: 'my_config.yaml', params: { key: value, ... } }
app.post('/configs', (req, res) => {
  try {
    const { filename, params } = req.body
    if (!filename || !params) return res.status(400).json({ error: 'filename and params required' })
    const safeName = path.basename(filename)
    if (!safeName.endsWith('.yaml') && !safeName.endsWith('.yml')) {
      return res.status(400).json({ error: 'Filename must end in .yaml or .yml' })
    }
    const fullPath = path.join(CONFIGS_DIR, safeName)
    const yamlStr = dumpYaml(params, { lineWidth: -1 })
    fs.writeFileSync(fullPath, yamlStr, 'utf8')
    res.json({ saved: safeName })
  } catch (err) {
    res.status(500).json({ error: err.message })
  }
})

// GET /run?config=filename — stream snakemake stdout/stderr via SSE
app.get('/run', (req, res) => {
  const configName = req.query.config
  if (!configName) return res.status(400).end('config query param required')

  const safeName = path.basename(configName)
  const configPath = path.join(CONFIGS_DIR, safeName)

  if (!fs.existsSync(configPath)) {
    return res.status(404).end('Config file not found')
  }

  res.setHeader('Content-Type', 'text/event-stream')
  res.setHeader('Cache-Control', 'no-cache')
  res.setHeader('Connection', 'keep-alive')
  res.flushHeaders()

  const sendEvent = (type, data) => {
    res.write(`data: ${JSON.stringify({ type, data })}\n\n`)
  }

  const args = SNAKEMAKE_ARGS.map(a => a.replace('{configFile}', configPath))
  sendEvent('info', `$ ${SNAKEMAKE_CMD} ${args.join(' ')}`)

  const proc = spawn(SNAKEMAKE_CMD, args, {
    env: { ...process.env },
    shell: false,
  })

  proc.stdout.on('data', chunk => sendEvent('stdout', chunk.toString()))
  proc.stderr.on('data', chunk => sendEvent('stderr', chunk.toString()))

  proc.on('close', code => {
    sendEvent('exit', code === 0 ? '✓ Pipeline finished successfully.' : `✗ Exited with code ${code}`)
    res.end()
  })

  proc.on('error', err => {
    sendEvent('error', `Failed to start process: ${err.message}`)
    res.end()
  })

  req.on('close', () => {
    proc.kill()
  })
})

app.listen(PORT, () => {
  console.log(`SnakeLauncher backend running on http://localhost:${PORT}`)
  console.log(`Config directory: ${CONFIGS_DIR}`)
})