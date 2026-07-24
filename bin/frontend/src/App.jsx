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

// Convert snake_case / UPPER_CASE key to a readable label
function toLabel(key) {
  return key
    .replace(/_/g, ' ')
    .replace(/([a-z])([A-Z])/g, '$1 $2')
    .toLowerCase()
    .replace(/^\w/, c => c.toUpperCase())
}


// ── LandingPage ──────────────────────────────────────────────────────────────
function LandingPage({ onChoose }) {
  return (
    <main className="landing-page">
      <section className="landing-hero">
        <img src={logoImg} alt="Asgard Pipeline" className="landing-logo" />
        <p className="mono eyebrow">Asgard discovery workspace</p>
        <h1>Where would you like to start?</h1>
        <p className="landing-copy">
          Keep the first screen focused on intent: run a reproducible workflow or search
          sequences quickly, then progressively reveal the detailed tools each path needs.
        </p>
        <div className="landing-actions" aria-label="Choose a workflow">
          <button className="choice-card primary-choice" onClick={() => onChoose('pipeline')}>
            <span className="choice-icon" aria-hidden="true">⚙</span>
            <span>
              <strong>Pipeline</strong>
              <small>Configure templates, edit YAML, run Snakemake, and inspect outputs.</small>
            </span>
          </button>
          <button className="choice-card" onClick={() => onChoose('blast')}>
            <span className="choice-icon" aria-hidden="true">🔎</span>
            <span>
              <strong>BLAST search</strong>
              <small>Paste a sequence, choose a target database, and review hits.</small>
            </span>
          </button>
        </div>
      </section>
    </main>
  )
}

function BlastPanel({ onBack }) {
  const [databases, setDatabases] = useState([])
  const [blastDir, setBlastDir] = useState('')
  const [program, setProgram] = useState('blastp')
  const [database, setDatabase] = useState('')
  const [query, setQuery] = useState('')
  const [evalue, setEvalue] = useState('1e-5')
  const [maxTargets, setMaxTargets] = useState(50)
  const [search, setSearch] = useState(null)
  const [logs, setLogs] = useState([])
  const [hits, setHits] = useState([])
  const [loadingDatabases, setLoadingDatabases] = useState(false)
  const [running, setRunning] = useState(false)
  const [error, setError] = useState('')

  const fetchBlastDatabases = useCallback(async () => {
    setLoadingDatabases(true)
    setError('')
    try {
      const res = await fetch(`${API}/blast/databases`)
      const data = await res.json()
      setDatabases(data.databases || [])
      setBlastDir(data.dir || '')
      setDatabase(current => current || data.databases?.[0]?.name || '')
    } catch (e) {
      setError(`Failed to load BLAST databases: ${e.message}`)
      setDatabases([])
    } finally {
      setLoadingDatabases(false)
    }
  }, [])

  useEffect(() => {
    const loadBlastDatabases = async () => {
      await fetchBlastDatabases()
    }
    loadBlastDatabases()
  }, [fetchBlastDatabases])

  const loadResults = async (searchId) => {
    const res = await fetch(`${API}/blast/searches/${searchId}/results`)
    if (!res.ok) return
    const data = await res.json()
    setHits(data.hits || [])
  }

  const handleSubmit = async (event) => {
    event.preventDefault()
    setError('')
    setLogs([])
    setHits([])
    setSearch(null)

    if (!query.trim()) {
      setError('Paste a FASTA sequence before starting a BLAST search.')
      return
    }
    if (!database) {
      setError('Select a BLAST database before starting a search.')
      return
    }

    setRunning(true)
    try {
      const res = await fetch(`${API}/blast/searches`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          query,
          database,
          program,
          evalue: Number(evalue),
          max_targets: Number(maxTargets),
        }),
      })
      const data = await res.json()
      if (!res.ok) throw new Error(data.detail || res.statusText)
      setSearch(data)

      const es = new EventSource(`${API}/blast/searches/${data.id}/stream`)
      es.onmessage = (msg) => {
        const payload = JSON.parse(msg.data)
        const type = payload.type || 'stdout'
        const text = payload.line ?? payload.data ?? JSON.stringify(payload)
        if (type === 'exit') {
          setLogs(current => [...current, { type: 'exit', text }])
          setRunning(false)
          setSearch(current => current ? { ...current, status: text } : current)
          es.close()
          if (text === 'succeeded') loadResults(data.id)
        } else {
          setLogs(current => [...current, { type, text }])
        }
      }
      es.onerror = () => {
        setLogs(current => [...current, { type: 'error', text: 'Lost connection to backend.' }])
        setRunning(false)
        es.close()
      }
    } catch (e) {
      setError(e.message)
      setRunning(false)
    }
  }

  return (
    <main className="pipeline-panel focused-panel">
      <header className="sk-header">
        <img src={logoImg} alt="Asgard Pipeline" className="header-logo" />
        <div>
          <h1 className="header-title">BLAST search</h1>
          <p className="mono header-sub dim">Sequence search workspace</p>
        </div>
      </header>

      <form className="card blast-form" onSubmit={handleSubmit}>
        <div className="card-header">
          <StepBadge n={1} active done={false} />
          <span className="step-label">Query setup</span>
          <button type="button" className="icon-btn" onClick={onBack} title="Back to landing">←</button>
          <button type="button" className="icon-btn" onClick={fetchBlastDatabases} title="Refresh databases">↺</button>
        </div>

        {blastDir && <div className="config-dir mono dim">📁 {blastDir}</div>}
        {loadingDatabases && <p className="dim">Loading BLAST databases…</p>}
        {!loadingDatabases && databases.length === 0 && (
          <p className="dim blast-empty">
            No BLAST databases found. Add formatted BLAST database files under the directory above.
          </p>
        )}

        <div className="blast-grid">
          <label className="blast-field">
            <span className="mono param-key">Program</span>
            <select className="param-input" value={program} onChange={e => setProgram(e.target.value)}>
              <option value="blastp">blastp</option>
              <option value="blastn">blastn</option>
              <option value="blastx">blastx</option>
              <option value="tblastn">tblastn</option>
              <option value="tblastx">tblastx</option>
            </select>
          </label>
          <label className="blast-field">
            <span className="mono param-key">Database</span>
            <select className="param-input" value={database} onChange={e => setDatabase(e.target.value)}>
              <option value="">Select a database</option>
              {databases.map(db => (
                <option key={db.name} value={db.name}>{db.name} ({db.type})</option>
              ))}
            </select>
          </label>
          <label className="blast-field">
            <span className="mono param-key">E-value</span>
            <input className="param-input" value={evalue} onChange={e => setEvalue(e.target.value)} />
          </label>
          <label className="blast-field">
            <span className="mono param-key">Max hits</span>
            <input className="param-input" type="number" min="1" value={maxTargets} onChange={e => setMaxTargets(e.target.value)} />
          </label>
        </div>

        <label className="blast-field blast-query">
          <span className="mono param-key">FASTA query</span>
          <textarea
            className="mono blast-textarea"
            value={query}
            onChange={e => setQuery(e.target.value)}
            placeholder={">query\nMSEQUENCE..."}
            rows={9}
          />
        </label>

        {error && <p className="save-msg err-msg">{error}</p>}
        <div className="save-row">
          <button className="btn primary" type="submit" disabled={running || databases.length === 0}>
            {running ? 'Searching…' : 'Run BLAST'}
          </button>
          <button className="btn" type="button" onClick={onBack}>Back to landing</button>
          {search && <span className="mono dim run-config">search #{search.id.slice(0, 8)} · {search.status}</span>}
        </div>
      </form>

      {logs.length > 0 && (
        <div className="card">
          <div className="card-header">
            <StepBadge n={2} active={running} done={!running && search?.status === 'succeeded'} />
            <span className="step-label">Search output</span>
          </div>
          <LogViewer lines={logs} running={running} />
        </div>
      )}

      {hits.length > 0 && (
        <div className="card">
          <div className="card-header">
            <StepBadge n={3} active done />
            <span className="step-label">Results</span>
            {search && (
              <a className="btn" href={`${API}/blast/searches/${search.id}/download`}>
                Download JSON
              </a>
            )}
          </div>
          <div className="blast-results-wrap">
            <table className="blast-results">
              <thead>
                <tr>
                  <th>Accession</th>
                  <th>Title</th>
                  <th>E-value</th>
                  <th>Bitscore</th>
                  <th>Identity</th>
                  <th>Length</th>
                </tr>
              </thead>
              <tbody>
                {hits.map((hit, idx) => (
                  <tr key={`${hit.accession}-${idx}`}>
                    <td className="mono">{hit.accession}</td>
                    <td>{hit.title}</td>
                    <td className="mono">{hit.evalue}</td>
                    <td className="mono">{hit.bitscore}</td>
                    <td className="mono">{hit.identity_pct ?? '—'}%</td>
                    <td className="mono">{hit.alignment_length}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>
      )}
    </main>
  )
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
function ConfigList({ configs, selected, onSelect, configsDir, onRefresh, onBack, loading, isTemplate }) {
  return (
    <div className="card">
      <div className="card-header">
        <StepBadge n={1} active={!selected} done={!!selected} />
        <span className="step-label">
          {isTemplate ? 'Choose a template' : 'Choose a config file'}
        </span>
        <button className="icon-btn" onClick={onBack} title="Back">←</button>
        <button className="icon-btn" onClick={onRefresh} title="Refresh">↺</button>
      </div>
      {configsDir && <div className="config-dir mono dim">📁 {configsDir}</div>}
      {loading && <p className="dim">Loading…</p>}
      {!loading && configs.length === 0 && (
        <p className="dim">
          {isTemplate ? 'No .template.yaml / .template.yml files found.' : 'No .yaml / .yml files found.'}
        </p>
      )}
      <div className="config-list">
        {configs.map(c => (
          <button
            key={c.name}
            className={`config-item ${selected === c.name ? 'selected' : ''}`}
            onClick={() => onSelect(c.name)}
          >
            <span className="config-icon">{isTemplate ? '📋' : '⚙'}</span>
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

// ── Questionnaire ─────────────────────────────────────────────────────────────
function QuestionnaireField({ fieldKey, value, path, onChange }) {
  const label = toLabel(fieldKey)
  const inputId = `q-${path.join('-')}`

  if (value !== null && typeof value === 'object' && !Array.isArray(value)) {
    return (
      <div className="param-group">
        <div className="mono param-group-label">{label}</div>
        <div className="param-group-body">
          {Object.entries(value).map(([k, v]) => (
            <QuestionnaireField
              key={k}
              fieldKey={k}
              value={v}
              path={[...path, k]}
              onChange={onChange}
            />
          ))}
        </div>
      </div>
    )
  }

  if (Array.isArray(value)) {
    return (
      <div className="param-row">
        <label className="mono param-key" htmlFor={inputId}>{label}</label>
        <input
          id={inputId}
          className="mono param-input"
          value={value.join(', ')}
          onChange={e =>
            onChange(path, e.target.value.split(',').map(s => s.trim()).filter(Boolean))
          }
        />
      </div>
    )
  }

  if (typeof value === 'boolean') {
    return (
      <div className="param-row">
        <label className="mono param-key" htmlFor={inputId}>{label}</label>
        <select
          id={inputId}
          className="mono param-input"
          value={String(value)}
          onChange={e => onChange(path, e.target.value === 'true')}
        >
          <option value="false">false</option>
          <option value="true">true</option>
        </select>
      </div>
    )
  }

  if (typeof value === 'number') {
    return (
      <div className="param-row">
        <label className="mono param-key" htmlFor={inputId}>{label}</label>
        <input
          id={inputId}
          className="mono param-input"
          type="number"
          value={value}
          onChange={e => onChange(path, Number(e.target.value))}
        />
      </div>
    )
  }

  return (
    <div className="param-row">
      <label className="mono param-key" htmlFor={inputId}>{label}</label>
      <input
        id={inputId}
        className="mono param-input"
        value={value ?? ''}
        onChange={e => onChange(path, e.target.value)}
      />
    </div>
  )
}

function Questionnaire({ params, onSubmit, onBack }) {
  const [answers, setAnswers] = useState(() => JSON.parse(JSON.stringify(params)))

  const handleChange = (path, value) =>
    setAnswers(a => setNestedValue(a, path, value))

  return (
    <div className="card">
      <div className="card-header">
        <StepBadge n={2} active done={false} />
        <span className="step-label">Configure your pipeline</span>
      </div>
      <div className="param-grid">
        {Object.entries(answers).map(([k, v]) => (
          <QuestionnaireField
            key={k}
            fieldKey={k}
            value={v}
            path={[k]}
            onChange={handleChange}
          />
        ))}
      </div>
      <div className="save-row" style={{ marginTop: '1rem' }}>
        <button className="btn primary" onClick={() => onSubmit(answers)}>
          Continue to full editor →
        </button>
        <button className="btn" onClick={onBack}>Back</button>
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
        <StepBadge n={3} active done={false} />
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
  // modes: null | 'existing' | 'template' | 'questionnaire' | 'template-edit'
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

  // ── fetch configs (existing flow) ──────────────────────────────────────────
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

  // ── fetch templates (template flow) ────────────────────────────────────────
  const fetchTemplates = useCallback(async () => {
    setLoadingConfigs(true)
    try {
      const res = await fetch(`${API}/templates`)
      const data = await res.json()
      setConfigs(data.configs || [])
      setConfigsDir(data.dir || '')
    } catch {
      setConfigs([])
    } finally {
      setLoadingConfigs(false)
    }
  }, [])

  // Load the right list whenever mode switches to a picker state
  useEffect(() => {
    const loadConfigsForMode = async () => {
      if (mode === 'existing') await fetchConfigs()
      else if (mode === 'template') await fetchTemplates()
    }
    loadConfigsForMode()
  }, [mode, fetchConfigs, fetchTemplates])

  // ── select a config or template from the list ──────────────────────────────
  // FIX: capture mode before any async work; don't overwrite it at the top of the fn
  const handleSelectConfig = async (name, explicitMode) => {
    // explicitMode lets callers pass the intended mode when state hasn't updated yet
    const currentMode = explicitMode ?? mode

    setSelected(name)
    setParams(null)
    setOriginal(null)
    setLogs([])
    setRunStatus('idle')
    setRunId(null)
    setLoadingParams(true)

    // Choose the right endpoint
    const endpoint =
      currentMode === 'template'
        ? `${API}/templates/${encodeURIComponent(name)}`
        : `${API}/configs/${encodeURIComponent(name)}`

    try {
      const res = await fetch(endpoint)
      const data = await res.json()
      const parsed = data.parsed || {}
      setParams(parsed)
      setOriginal(parsed)
      // Template → questionnaire step; existing config → straight to editor
      setMode(currentMode === 'template' ? 'questionnaire' : 'existing')
    } catch (e) {
      setParams({ error: e.message })
      setMode('existing')
    } finally {
      setLoadingParams(false)
    }
  }

  // FIX: single declaration of handleParamChange
  const handleParamChange = (path, value) => setParams(p => setNestedValue(p, path, value))

  // FIX: single declaration of handleReset
  const handleReset = () => setParams(JSON.parse(JSON.stringify(original)))

  const handleBack = () => {
    setMode(null)
    setSelected(null)
    setParams(null)
    setOriginal(null)
    setLogs([])
    setRunId(null)
    setRunStatus('idle')
    setRunning(false)
  }

  const handleQuestionnaireSubmit = (answers) => {
    setParams(answers)
    setMode('template-edit')
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
    // After saving, switch to existing mode and load the new config
    setMode('existing')
    await handleSelectConfigAsExisting(filename)
  }

  // Separate helper to load a config in existing mode (used after save)
  const handleSelectConfigAsExisting = async (name) => {
    setSelected(name)
    setLoadingParams(true)
    try {
      const res = await fetch(`${API}/configs/${encodeURIComponent(name)}`)
      const data = await res.json()
      const parsed = data.parsed || {}
      setParams(parsed)
      setOriginal(parsed)
    } catch (e) {
      setParams({ error: e.message })
    } finally {
      setLoadingParams(false)
    }
  }

  const handleRun = async () => {
    if (!selected || running) return
    setRunning(true)
    setLogs([])
    setRunStatus('running')

    const push = (type, text) => setLogs(l => [...l, { type, text }])

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

  // ── derived display flags ──────────────────────────────────────────────────
  const showConfigPicker  = mode === 'existing' || mode === 'template'
  const showQuestionnaire = mode === 'questionnaire' && !!params
  const showParamEditor   = !!params && !!selected && (mode === 'existing' || mode === 'template-edit')
  const showRunPanel      = !!params && !!selected && (mode === 'existing' || mode === 'template-edit')

  return (
    <div className="pipeline-panel">
      <header className="sk-header">
        <img src={logoImg} alt="Asgard Pipeline" className="header-logo" />
        <div>
          <h1 className="header-title">Asgard Pipeline</h1>
          <p className="mono header-sub dim">Discover Asgard</p>
        </div>
      </header>

      {/* Step 1 — choose how to start */}
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
            {/* FIX: single onClick, cleanly sets mode then fetches templates */}
            <button
              className="btn"
              onClick={() => {
                setSelected(null)
                setParams(null)
                setOriginal(null)
                setMode('template')
              }}
            >
              Create from template
            </button>
          </div>
          <p className="dim" style={{ marginTop: '8px', fontSize: '12px' }}>
            Tip: "Create from template" opens templates/configs/*.template.yaml for editing and saving as a new config.
          </p>
        </div>
      )}

      {/* Step 1 (continued) — pick a config / template from the list */}
      {showConfigPicker && (
        <ConfigList
          configs={configs}
          selected={selected}
          onSelect={handleSelectConfig}
          configsDir={configsDir}
          onRefresh={mode === 'template' ? fetchTemplates : fetchConfigs}
          onBack={handleBack}
          loading={loadingConfigs}
          isTemplate={mode === 'template'}
        />
      )}

      {loadingParams && <p className="dim loading-params">Loading parameters…</p>}

      {/* Step 2 (template flow) — questionnaire generated from template fields */}
      {showQuestionnaire && (
        <Questionnaire
          params={params}
          onSubmit={handleQuestionnaireSubmit}
          onBack={() => {
            setMode('template')
            setSelected(null)
            setParams(null)
            setOriginal(null)
          }}
        />
      )}

      {/* Step 2 (existing) / Step 3 (template) — full param editor */}
      {showParamEditor && (
        <ParamEditor
          configName={selected}
          params={params}
          onChange={handleParamChange}
          onSaveNew={handleSaveNew}
          onReset={handleReset}
          isDirty={isDirty}
        />
      )}

      {/* Final step — run */}
      {showRunPanel && (
        <div className="card">
          <div className="card-header">
            <StepBadge n={mode === 'existing' ? 3 : 4} active done={false} />
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
  const [appMode, setAppMode] = useState(null)

  if (!appMode) return <LandingPage onChoose={setAppMode} />
  if (appMode === 'blast') return <BlastPanel onBack={() => setAppMode(null)} />

  return (
    <div className="workspace-grid">
      <FileExplorer />
      <PipelinePanel />
    </div>
  )
}