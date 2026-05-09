// App.jsx — SnakeLauncher frontend
// Drop this into src/App.jsx in your Vite + React project

import { useState, useEffect, useRef, useCallback } from 'react'
import './App.css'

import FileExplorer from './components/FileExplorer'
import './App.css'

const API = 'http://localhost:8000'


function App() {

  return (
    <FileExplorer />
  )
}

export default App



// ── tiny helpers ──────────────────────────────────────────────────────────────
function timeAgo(dateStr) {
  const diff = (Date.now() - new Date(dateStr)) / 1000
  if (diff < 60) return 'just now'
  if (diff < 3600) return `${Math.floor(diff / 60)}m ago`
  if (diff < 86400) return `${Math.floor(diff / 3600)}h ago`
  return `${Math.floor(diff / 86400)}d ago`
}

function fmt(bytes) {
  return bytes < 1024 ? `${bytes} B` : `${(bytes / 1024).toFixed(1)} KB`
}

// ── Step indicator ────────────────────────────────────────────────────────────
function StepBadge({ n, active, done }) {
  return (
    <div
      className={`step-badge ${active ? 'active' : ''} ${done ? 'done' : ''}`}
      aria-label={`Step ${n}`}
    >
      {done ? '✓' : n}
    </div>
  )
}

// ── Directory picker ──────────────────────────────────────────────────────────
function DirPicker({ currentDir, onDirSet }) {
  const [input, setInput] = useState(currentDir || '')
  const [status, setStatus] = useState('')
  const [loading, setLoading] = useState(false)

  useEffect(() => { setInput(currentDir || '') }, [currentDir])

  const handleSet = async () => {
    if (!input.trim()) return
    setLoading(true)
    setStatus('')
    try {
      const res = await fetch(`${API}/set-dir`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ dir: input.trim() }),
      })
      const data = await res.json()
      if (!res.ok) throw new Error(data.error)
      setStatus('ok')
      onDirSet(data.dir)
    } catch (e) {
      setStatus(e.message)
    } finally {
      setLoading(false)
    }
  }

  const handleKey = (e) => { if (e.key === 'Enter') handleSet() }

  return (
    <div className="card dir-picker-card">
      <label className="step-label" htmlFor="dir-input">Configs directory</label>
      <div className="dir-row">
        <span className="dir-icon mono">📁</span>
        <input
          id="dir-input"
          className="mono param-input dir-input"
          value={input}
          onChange={e => { setInput(e.target.value); setStatus('') }}
          onKeyDown={handleKey}
          placeholder="/home/you/snakemake-configs"
          spellCheck={false}
        />
        <button className="btn primary" onClick={handleSet} disabled={loading}>
          {loading ? '…' : 'Set'}
        </button>
      </div>
      {status === 'ok' && <p className="save-msg">✓ Directory set</p>}
      {status && status !== 'ok' && <p className="save-msg" style={{color:'var(--red,#E24B4A)'}}>{status}</p>}
    </div>
  )
}


function ConfigList({ configs, selected, onSelect, configsDir, onRefresh, loading }) {
  return (
    <div className="card">
      <div className="card-header">
        <StepBadge n={1} active={!selected} done={!!selected} />
        <span className="step-label">Choose a config file</span>
        <button className="icon-btn" onClick={onRefresh} title="Refresh" aria-label="Refresh config list">
          ↺
        </button>
      </div>
      {configsDir && (
        <div className="config-dir">
          <span className="mono dim">📁 {configsDir}</span>
        </div>
      )}
      {loading && <p className="dim">Loading configs…</p>}
      {!loading && configs.length === 0 && (
        <p className="dim">No .yaml / .yml files found in that directory.</p>
      )}
      <div className="config-list">
        {configs.map(c => (
          <button
            key={c.name}
            className={`config-item ${selected?.name === c.name ? 'selected' : ''}`}
            onClick={() => onSelect(c)}
          >
            <span className="config-icon">⚙</span>
            <span className="config-info">
              <span
                className="mono config-name"
                title={c.name}
              >
                {c.name.replace(".yaml", "")}
              </span>
              <span className="config-meta dim">{timeAgo(c.mtime)} · {fmt(c.size)}</span>
            </span>
            {selected?.name === c.name && <span className="selected-dot" aria-hidden="true" />}
          </button>
        ))}
      </div>
    </div>
  )
}

// ── Param editor ──────────────────────────────────────────────────────────────
function ParamEditor({ configName, params, onChange, onSaveNew, onReset, isDirty }) {
  const [newName, setNewName] = useState('')
  const [saving, setSaving] = useState(false)
  const [saveMsg, setSaveMsg] = useState('')

  const handleSave = async () => {
    const target = newName.trim() || configName
    if (!target.endsWith('.yaml') && !target.endsWith('.yml')) {
      setSaveMsg('Filename must end in .yaml or .yml')
      return
    }
    setSaving(true)
    setSaveMsg('')
    try {
      await onSaveNew(target, params)
      setSaveMsg(`Saved as ${target}`)
    } catch (e) {
      setSaveMsg(`Error: ${e.message}`)
    } finally {
      setSaving(false)
    }
  }

  return (
    <div className="card">
      <div className="card-header">
        <StepBadge n={2} active done={false} />
        <span className="step-label">Review & edit parameters</span>
        {isDirty && <span className="dirty-badge">edited</span>}
      </div>
      <div className="param-grid">
        {Object.entries(params).map(([k, v]) => (
          <div className="param-row" key={k}>
            <label className="mono param-key" htmlFor={`param-${k}`}>{k}</label>
            <input
              id={`param-${k}`}
              className="mono param-input"
              value={String(v)}
              onChange={e => onChange(k, e.target.value)}
            />
          </div>
        ))}
      </div>
      {isDirty && (
        <div className="save-row">
          <input
            className="mono new-name-input"
            placeholder={`New name (default: ${configName})`}
            value={newName}
            onChange={e => setNewName(e.target.value)}
          />
          <button className="btn primary" onClick={handleSave} disabled={saving}>
            {saving ? 'Saving…' : 'Save as new config'}
          </button>
          <button className="btn" onClick={onReset}>Reset</button>
        </div>
      )}
      {saveMsg && <p className="save-msg">{saveMsg}</p>}
    </div>
  )
}

// ── Log viewer ────────────────────────────────────────────────────────────────
function LogViewer({ lines, running }) {
  const bottomRef = useRef(null)

  useEffect(() => {
    bottomRef.current?.scrollIntoView({ behavior: 'smooth' })
  }, [lines])

  return (
    <div className="log-panel" aria-label="Pipeline output" aria-live="polite">
      {lines.map((l, i) => (
        <div key={i} className={`log-line ${l.type}`}>
          {l.text}
        </div>
      ))}
      {running && <div className="log-line info blink">▌</div>}
      <div ref={bottomRef} />
    </div>
  )
}

// ── Main App ──────────────────────────────────────────────────────────────────
export default function App() {
  const [configs, setConfigs] = useState([])
  const [configsDir, setConfigsDir] = useState('')
  const [loadingConfigs, setLoadingConfigs] = useState(false)

  const [selectedConfig, setSelectedConfig] = useState(null)
  const [params, setParams] = useState(null)
  const [originalParams, setOriginalParams] = useState(null)
  const [loadingParams, setLoadingParams] = useState(false)

  const [running, setRunning] = useState(false)
  const [logLines, setLogLines] = useState([])
  const [activeConfig, setActiveConfig] = useState(null) // the config actually used to run

  const fetchConfigs = useCallback(async () => {
    setLoadingConfigs(true)
    try {
      const res = await fetch(`${API}/configs`)
      const data = await res.json()
      setConfigs(data.configs || [])
      setConfigsDir(data.dir || '')
    } catch {
      setConfigs([])
    } finally {
      setLoadingConfigs(false)
    }
  }, [])

  useEffect(() => { fetchConfigs() }, [fetchConfigs])

  const handleSelectConfig = async (c) => {
    setSelectedConfig(c)
    setParams(null)
    setOriginalParams(null)
    setLogLines([])
    setLoadingParams(true)
    try {
      const res = await fetch(`${API}/configs/${encodeURIComponent(c.name)}`)
      const data = await res.json()
      setParams(data.parsed || {})
      setOriginalParams(data.parsed || {})
    } catch (e) {
      setParams({ error: e.message })
    } finally {
      setLoadingParams(false)
    }
  }

  const handleParamChange = (key, value) => {
    setParams(p => ({ ...p, [key]: value }))
  }

  const handleReset = () => {
    setParams({ ...originalParams })
  }

  const isDirty = params && originalParams &&
    JSON.stringify(params) !== JSON.stringify(originalParams)

  const handleSaveNew = async (filename, newParams) => {
    const res = await fetch(`${API}/configs`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ filename, params: newParams }),
    })
    if (!res.ok) {
      const err = await res.json()
      throw new Error(err.error)
    }
    await fetchConfigs()
    // switch selection to the newly saved file
    const saved = (await (await fetch(`${API}/configs`)).json()).configs
      .find(c => c.name === filename)
    if (saved) {
      setSelectedConfig(saved)
      setOriginalParams({ ...newParams })
    }
  }

  const handleRun = () => {
    if (!selectedConfig || running) return
    setRunning(true)
    setLogLines([])
    setActiveConfig(selectedConfig.name)

    const es = new EventSource(
      `${API}/run?config=${encodeURIComponent(selectedConfig.name)}`
    )

    const push = (type, text) => {
      setLogLines(l => [...l, { type, text }])
    }

    es.onmessage = (e) => {
      const { type, data } = JSON.parse(e.data)
      if (type === 'exit') {
        push('exit', data)
        setRunning(false)
        es.close()
      } else {
        // split multiline chunks into individual lines
        String(data).split('\n').filter(Boolean).forEach(line => push(type, line))
      }
    }

    es.onerror = () => {
      push('error', 'Connection to backend lost.')
      setRunning(false)
      es.close()
    }
  }

  return (
    <div className="sk-app">
      <header className="sk-header">
        <span className="header-snake" aria-hidden="true">🐍</span>
        <div>
          <h1 className="header-title">SnakeLauncher</h1>
          <p className="mono header-sub dim">
            $ snakemake --configfile [config] --cores all
          </p>
        </div>
      </header>

      <main>
        {/* Step 1: Config list */}
        <ConfigList
          configs={configs}
          selected={selectedConfig}
          onSelect={handleSelectConfig}
          configsDir={configsDir}
          onRefresh={fetchConfigs}
          loading={loadingConfigs}
        />

        {/* Step 2: Param editor */}
        {loadingParams && <p className="dim loading-params">Loading parameters…</p>}
        {params && selectedConfig && (
          <ParamEditor
            configName={selectedConfig.name}
            params={params}
            onChange={handleParamChange}
            onSaveNew={handleSaveNew}
            onReset={handleReset}
            isDirty={isDirty}
          />
        )}

        {/* Step 3: Run */}
        {params && selectedConfig && (
          <div className="card">
            <div className="card-header">
              <StepBadge n={3} active done={false} />
              <span className="step-label">Run the pipeline</span>
            </div>
            <div className="run-row">
              <button
                className={`btn run-btn ${running ? 'running' : ''}`}
                onClick={handleRun}
                disabled={running}
              >
                {running ? '⏳ Running…' : '▶ Run pipeline'}
              </button>
              <span className="mono dim run-config">
                using <strong>{selectedConfig.name}</strong>
              </span>
            </div>

            {logLines.length > 0 && (
              <LogViewer lines={logLines} running={running} />
            )}
          </div>
        )}
      </main>
    </div>
  )
}