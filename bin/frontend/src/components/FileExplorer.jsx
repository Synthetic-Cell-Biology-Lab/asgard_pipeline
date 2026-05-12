import { useEffect, useState } from 'react'
import Breadcrumbs from './explorer/Breadcrumbs'
import FileRow from './explorer/FileRow'
import PreviewPane from './preview/PreviewPane'

const API_BASE = 'http://localhost:8000'

export default function FileExplorer() {
  const [currentPath, setCurrentPath] = useState('')
  const [parentPath, setParentPath] = useState(null)
  const [entries, setEntries] = useState([])
  const [selected, setSelected] = useState(null)
  const [fileData, setFileData] = useState(null)
  const [loading, setLoading] = useState(false)
  const [previewLoading, setPreviewLoading] = useState(false)

  async function loadDirectory(path = '') {
    setLoading(true)
    try {
      const res = await fetch(`${API_BASE}/browse?path=${encodeURIComponent(path)}`)
      const data = await res.json()
      setCurrentPath(data.current_path)
      setParentPath(data.parent_path)
      setEntries(data.entries)
      setSelected(null)
      setFileData(null)
    } finally {
      setLoading(false)
    }
  }

  async function openEntry(entry) {
    setSelected(entry.path)
    if (entry.type === 'directory') {
      loadDirectory(entry.path)
      return
    }
    setPreviewLoading(true)
    setFileData(null)
    try {
      const res = await fetch(`${API_BASE}/file?path=${encodeURIComponent(entry.path)}`)
      setFileData(await res.json())
    } finally {
      setPreviewLoading(false)
    }
  }

  useEffect(() => { loadDirectory('database') }, [])

  return (
    <div className="explorer-shell">
      {/* Left: file tree */}
      <div className="explorer">
        <div className="explorer-header">
          <button
            className="nav-btn"
            onClick={() => parentPath !== null && loadDirectory(parentPath)}
            disabled={parentPath === null}
            title="Go up"
          >
            ↑
          </button>
          <Breadcrumbs currentPath={currentPath} onNavigate={loadDirectory} />
        </div>
        <div className="file-list">
          {loading && <p className="explorer-empty dim">Loading…</p>}
          {!loading && entries.length === 0 && (
            <p className="explorer-empty dim">Empty directory</p>
          )}
          {!loading && entries.map((entry) => (
            <FileRow
              key={entry.path}
              entry={entry}
              selected={selected === entry.path}
              onClick={() => openEntry(entry)}
            />
          ))}
        </div>
      </div>

      {/* Right: preview */}
      <PreviewPane
        fileData={fileData}
        loading={previewLoading}
        apiBase={API_BASE}
      />
    </div>
  )
}