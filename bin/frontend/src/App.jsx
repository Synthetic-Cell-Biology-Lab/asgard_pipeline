// App.jsx — Asgard Pipeline
import { useState, useEffect, useRef, useCallback, useMemo } from 'react'
import FileExplorer from './components/FileExplorer'
import DagView from './components/DagView'
import './App.css'
import logoImg from './assets/logo.png'

const API = 'http://localhost:8000'

// ── helpers ───────────────────────────────────────────────────────────────────
function timeAgo(ts) {
  const diff = (Date.now() - ts * 1000) / 1000
  if (diff < 60) return 'just now'
  if (diff < 3600) return `${Math.floor(diff / 60)}m ago`
  if (diff < 86400) return `${Math.floor(diff / 3600)}h ago`
  return `${Math.floor(diff / 86400)}d ago`
}

function fmt(bytes) {
  return bytes < 1024 ? `${bytes} B` : `${(bytes / 1024).toFixed(1)} KB`
}

function coerceValue(value, previous) {
  if (typeof previous === 'number') return Number(value)
  if (typeof previous === 'boolean') return value === 'true'
  return value
}

function setNestedValue(target, path, nextValue) {
  if (path.length === 0) return nextValue
  const [head, ...rest] = path
  if (Array.isArray(target)) {
    const copy = [...target]
    copy[head] = setNestedValue(copy[head], rest, nextValue)
    return copy
  }
  return {
    ...(target || {}),
    [head]: setNestedValue(target?.[head], rest, nextValue),
  }
}

// ── StepBadge ─────────────────────────────────────────────────────────────────
function StepBadge({ n, active, done }) {
  return (
    <div className={`step-badge ${active ? 'active' : ''} ${done ? 'done' : ''}`}>
      {done ? '✓' : n}
    </div>
  )
}

// ── ConfigList ────────────────────────────────────────────────────────────────
function ConfigList({ configs, selected, onSelect, configsDir, onRefresh, loading }) {
  return (
    <div className="card">
      <div className="card-header">
        <StepBadge n={1} active={!selected} done={!!selected} />
        <span className="step-label">Choose a config file</span>
        <button className="icon-btn" onClick={onRefresh} title="Refresh">↺</button>
      </div>
      {configsDir && <div className="config-dir mono dim">📁 {configsDir}</div>}
      {loading && <p className="dim">Loading…</p>}
      {!loading && configs.length === 0 && (
        <p className="dim">No .yaml / .yml files found.</p>
      )}
      <div className="config-list">
        {configs.map(c => (
          <button
            key={c.name}
            className={`config-item ${selected === c.name ? 'selected' : ''}`}
            onClick={() => onSelect(c.name)}
          >
            <span className="config-icon">⚙</span>
            <span className="config-info">
              <span className="mono config-name">{c.name}</span>
              <span className="config-meta dim">{timeAgo(c.mtime)} · {fmt(c.size)}</span>
            </span>
            {selected === c.name && <span className="selected-dot" />}
          </button>
        ))}
      </div>
    </div>
  )
}

// ── ParamEditor ───────────────────────────────────────────────────────────────
function ParamField({ name, value, path, onChange, depth = 0 }) {
  const inputId = `p-${path.join('-')}`
  const keyLabel = depth === 0 ? name : `↳ ${name}`

  if (value && typeof value === 'object' && !Array.isArray(value)) {
    return (
      <div className="param-group">
        <div className="mono param-group-label">{keyLabel}</div>
        <div className="param-group-body">
          {Object.entries(value).map(([k, v]) => (
            <ParamField
              key={k}
              name={k}
              value={v}
              path={[...path, k]}
              onChange={onChange}
              depth={depth + 1}
            />
          ))}
        </div>
      </div>
    )
  }

  if (Array.isArray(value)) {
    return (
      <div className="param-group">
        <div className="mono param-group-label">{keyLabel}</div>
        <div className="param-group-body">
          {value.map((item, idx) => (
            <ParamField
              key={`${name}-${idx}`}
              name={`[${idx}]`}
              value={item}
              path={[...path, idx]}
              onChange={onChange}
              depth={depth + 1}
            />
          ))}
        </div>
      </div>
    )
  }

  return (
    <div className="param-row">
      <label className="mono param-key" htmlFor={inputId}>{keyLabel}</label>
      <input
        id={inputId}
        className="mono param-input"
        value={String(value ?? '')}
        onChange={e => onChange(path, coerceValue(e.target.value, value))}
      />
    </div>
  )
}

function ParamEditor({ configName, params, onChange, onSaveNew, onReset, isDirty }) {
  const [saveAs, setSaveAs] = useState('')
  const [saving, setSaving] = useState(false)
  const [saveMsg, setSaveMsg] = useState('')

  const handleSave = async () => {
    const target = saveAs.trim() || configName
    if (!target.endsWith('.yaml') && !target.endsWith('.yml')) {
      setSaveMsg('Filename must end in .yaml or .yml')
      return
    }
    setSaving(true)
    setSaveMsg('')
    try {
      await onSaveNew(target, params)
      setSaveMsg(`✓ Saved as ${target}`)
      setSaveAs('')
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
          <ParamField key={k} name={k} value={v} path={[k]} onChange={onChange} />
        ))}
      </div>
      {isDirty && (
        <div className="save-row">
          <input
            className="mono new-name-input"
            placeholder={`New name (default: ${configName})`}
            value={saveAs}
            onChange={e => setSaveAs(e.target.value)}
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

// ── LogViewer ─────────────────────────────────────────────────────────────────
function LogViewer({ lines, running }) {
  const bottomRef = useRef(null)
  useEffect(() => { bottomRef.current?.scrollIntoView({ behavior: 'smooth' }) }, [lines])

  return (
    <div className="log-panel" aria-label="Pipeline output" aria-live="polite">
      {lines.map((l, i) => (
        <div key={i} className={`log-line ${l.type}`}>{l.text}</div>
      ))}
      {running && <div className="log-line info blink">▌</div>}
      <div ref={bottomRef} />
    </div>
  )
}

// ── PipelinePanel ─────────────────────────────────────────────────────────────
function PipelinePanel() {
  const [configs, setConfigs] = useState([])
  const [configsDir, setConfigsDir] = useState('')
  const [loadingConfigs, setLoadingConfigs] = useState(false)

  const [selected, setSelected] = useState(null)
  const [mode, setMode] = useState(null)
  const [params, setParams] = useState(null)
  const [original, setOriginal] = useState(null)
  const [loadingParams, setLoadingParams] = useState(false)

  const [running, setRunning] = useState(false)
  const [logs, setLogs] = useState([])
  const [runId, setRunId] = useState(null)
  const [runStatus, setRunStatus] = useState('idle')

  const isDirty = useMemo(
    () => params && original && JSON.stringify(params) !== JSON.stringify(original),
    [params, original]
  )
  const templateCandidates = useMemo(
    () => configs.filter(c => c.name.endsWith('.template.yaml')).map(c => c.name),
    [configs]
  )

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

  const handleSelectConfig = async (name) => {
    setMode('existing')
    setSelected(name)
    setParams(null)
    setOriginal(null)
    setLogs([])
    setRunStatus('idle')
    setRunId(null)
    setLoadingParams(true)
    try {
      const res = await fetch(`${API}/configs/${encodeURIComponent(name)}`)
      const data = await res.json()
      setParams(data.parsed || {})
      setOriginal(data.parsed || {})
    } catch (e) {
      setParams({ error: e.message })
    } finally {
      setLoadingParams(false)
    }
  }

  const handleParamChange = (path, value) => setParams(p => setNestedValue(p, path, value))
  const handleReset = () => setParams({ ...original })

  const handleStartFromTemplate = async (name) => {
    const templateConfigName = `${name}.template.yaml`
    await handleSelectConfig(templateConfigName)
    setMode('template')
  }

  const handleSaveNew = async (filename, newParams) => {
    const res = await fetch(`${API}/configs`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ filename, params: newParams }),
    })
    if (!res.ok) {
      const err = await res.json().catch(() => ({}))
      throw new Error(err.detail || err.error || res.statusText)
    }
    await fetchConfigs()
    await handleSelectConfig(filename)
  }

  const handleRun = async () => {
    if (!selected || running) return
    setRunning(true)
    setLogs([])
    setRunStatus('running')

    const push = (type, text) => setLogs(l => [...l, { type, text }])

    // Step 1: create the run
    let runData
    try {
      const res = await fetch(`${API}/runs`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ config: selected }),
      })
      if (!res.ok) {
        const err = await res.json().catch(() => ({ detail: res.statusText }))
        throw new Error(err.detail || JSON.stringify(err))
      }
      runData = await res.json()
    } catch (e) {
      push('error', `Failed to start run: ${e.message}`)
      setRunning(false)
      setRunStatus('error')
      return
    }

    setRunId(runData.id)

    // Step 2: stream logs via SSE
    const es = new EventSource(`${API}/runs/${runData.id}/stream`)
    es.onmessage = (msg) => {
      const payload = JSON.parse(msg.data)
      const text = payload.line ?? payload.data ?? JSON.stringify(payload)
      const type = payload.type || 'stdout'
      if (type === 'exit') {
        push('exit', text)
        setRunStatus(text)
        setRunning(false)
        es.close()
      } else {
        String(text).split('\n').filter(Boolean).forEach(t => push(type, t))
      }
    }
    es.onerror = () => {
      push('error', 'Lost connection to backend.')
      setRunning(false)
      setRunStatus('error')
      es.close()
    }
  }

  const showConfigPicker = mode === 'existing' || (!!selected && mode === 'template')

  return (
    <div className="pipeline-panel">
      <header className="sk-header">
        <img src={logoImg} alt="Asgard Pipeline" className="header-logo" />
        <div>
          <h1 className="header-title">Asgard Pipeline</h1>
          <p className="mono header-sub dim">Discover Asgard</p>
        </div>
      </header>

      {!mode && (
        <div className="card">
          <div className="card-header">
            <StepBadge n={1} active done={false} />
            <span className="step-label">Choose how to start</span>
          </div>
          <div className="save-row">
            <button className="btn primary" onClick={() => setMode('existing')}>
              Use existing config
            </button>
            <button
              className="btn"
              onClick={() => {
                if (templateCandidates.length > 0) {
                  handleSelectConfig(templateCandidates[0])
                  setMode('template')
                }
              }}
              disabled={templateCandidates.length === 0}
            >
              Create from template
            </button>
          </div>
          <p className="dim" style={{ marginTop: '8px', fontSize: '12px' }}>
            Tip: "Create from template" opens templates/configs/*.template.yaml for editing and saving as a new config.
          </p>
        </div>
      )}

      {showConfigPicker && (
        <ConfigList
          configs={configs}
          selected={selected}
          onSelect={handleSelectConfig}
          configsDir={configsDir}
          onRefresh={fetchConfigs}
          loading={loadingConfigs}
        />
      )}

      {!mode && (
        <div className="card">
          <div className="card-header">
            <StepBadge n={2} active done={false} />
            <span className="step-label">Choose how to start</span>
          </div>
          <div className="save-row">
            <button className="btn primary" onClick={() => setMode('existing')}>
              Use existing config
            </button>
            <button
              className="btn"
              onClick={() => {
                if (templateCandidates.length > 0) {
                  handleSelectConfig(templateCandidates[0])
                  setMode('template')
                }
              }}
              disabled={templateCandidates.length === 0}
            >
              Create from template
            </button>
          </div>
          <p className="dim" style={{ marginTop: '8px', fontSize: '12px' }}>
            Tip: "Create from template" opens templates/configs/*.template.yaml for editing and saving as a new config.
          </p>
        </div>
      )}

      {loadingParams && <p className="dim loading-params">Loading parameters…</p>}

      {params && selected && mode && (
        <ParamEditor
          configName={selected}
          params={params}
          onChange={handleParamChange}
          onSaveNew={handleSaveNew}
          onReset={handleReset}
          isDirty={isDirty}
        />
      )}

      {params && selected && (
        <div className="card">
          <div className="card-header">
            <StepBadge n={3} active done={false} />
            <span className="step-label">Run the pipeline</span>
            {runId && (
              <span className="mono dim" style={{ marginLeft: 'auto', fontSize: '11px' }}>
                run #{runId.slice(0, 8)}
              </span>
            )}
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
              using <strong>{selected}</strong>
            </span>
            {runStatus !== 'idle' && runStatus !== 'running' && (
              <span className="mono dim" style={{ fontSize: '11px' }}>{runStatus}</span>
            )}
          </div>
          {logs.length > 0 && <LogViewer lines={logs} running={running} />}
          {(runId || running) && <DagView runId={runId} running={running} />}
        </div>
      )}
    </div>
  )
}

// ── App root ──────────────────────────────────────────────────────────────────
export default function App() {
  return (
    <div className="workspace-grid">
      <FileExplorer />
      <PipelinePanel />
    </div>
  )
}
