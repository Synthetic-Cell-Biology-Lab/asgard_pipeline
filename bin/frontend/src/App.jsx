import { useEffect, useMemo, useState } from 'react'
import FileExplorer from './components/FileExplorer'
import './App.css'

const API = 'http://localhost:8000'

function coerceValue(value, previous) {
  if (typeof previous === 'number') return Number(value)
  if (typeof previous === 'boolean') return value === 'true'
  return value
}

function PipelinePanel() {
  const [configs, setConfigs] = useState([])
  const [selected, setSelected] = useState('')
  const [params, setParams] = useState({})
  const [original, setOriginal] = useState({})
  const [saveAs, setSaveAs] = useState('')
  const [runId, setRunId] = useState(null)
  const [runStatus, setRunStatus] = useState('idle')
  const [logs, setLogs] = useState([])

  const dirty = useMemo(() => JSON.stringify(params) !== JSON.stringify(original), [params, original])

  async function loadConfigs() {
    const res = await fetch(`${API}/configs`)
    const data = await res.json()
    setConfigs(data.configs || [])
  }

  async function loadConfig(name) {
    const res = await fetch(`${API}/configs/${encodeURIComponent(name)}`)
    const data = await res.json()
    setSelected(name)
    setParams(data.parsed || {})
    setOriginal(data.parsed || {})
  }

  async function saveConfig() {
    const filename = (saveAs || selected).trim()
    if (!filename) return
    await fetch(`${API}/configs`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ filename, params })
    })
    setSaveAs('')
    await loadConfigs()
    await loadConfig(filename)
  }

  async function runPipeline() {
    const res = await fetch(`${API}/runs`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ config: selected })
    })
    const data = await res.json()
    setRunId(data.id)
    setRunStatus(data.status)
    setLogs([])

    const evt = new EventSource(`${API}/runs/${data.id}/stream`)
    evt.onmessage = (msg) => {
      const payload = JSON.parse(msg.data)
      setLogs((prev) => [...prev, payload])
      if (payload.type === 'exit') {
        setRunStatus(payload.line)
        evt.close()
      }
    }
    evt.onerror = () => evt.close()
  }

  useEffect(() => { loadConfigs() }, [])

  return (
    <div className="panel-shell">
      <div className="panel-header">
        <h2>Pipeline Launcher</h2>
        <button onClick={loadConfigs}>Refresh</button>
      </div>

      <div className="field">
        <label>Config file</label>
        <select value={selected} onChange={(e) => loadConfig(e.target.value)}>
          <option value="">Select config…</option>
          {configs.map((c) => <option key={c.name} value={c.name}>{c.name}</option>)}
        </select>
      </div>

      {selected && (
        <>
          <div className="params-card">
            {Object.entries(params).map(([k, v]) => (
              <div className="param-row" key={k}>
                <label>{k}</label>
                <input value={String(v)} onChange={(e) => setParams((prev) => ({ ...prev, [k]: coerceValue(e.target.value, prev[k]) }))} />
              </div>
            ))}
          </div>

          <div className="action-row">
            <input placeholder="save as (optional)" value={saveAs} onChange={(e) => setSaveAs(e.target.value)} />
            <button onClick={saveConfig} disabled={!dirty && !saveAs}>Save config</button>
            <button className="primary" onClick={runPipeline}>Run pipeline</button>
          </div>
        </>
      )}

      <div className="run-meta">Run ID: {runId || '—'} · Status: {runStatus}</div>
      <div className="log-panel">
        {logs.map((l, idx) => <div key={idx} className={`log-line ${l.type || ''}`}>{l.line || JSON.stringify(l)}</div>)}
      </div>
    </div>
  )
}

function App() {
  return (
    <div className="workspace-grid">
      <FileExplorer />
      <PipelinePanel />
    </div>
  )
}

export default App
