import { useEffect, useState } from 'react'

const API_BASE = "http://localhost:8000"

export default function FileExplorer() {

  const [currentPath, setCurrentPath] = useState("")
  const [entries, setEntries] = useState([])

  async function loadDirectory(path = "") {

    const res = await fetch(
      `${API_BASE}/browse?path=${encodeURIComponent(path)}`
    )

    const data = await res.json()

    setCurrentPath(data.current_path)
    setEntries(data.entries)
  }

  useEffect(() => {
    loadDirectory("protein_sets")
  }, [])

  function openEntry(entry) {

    if (entry.type !== "directory") {
      return
    }

    const nextPath =
      currentPath === ""
        ? entry.name
        : `${currentPath}/${entry.name}`

    loadDirectory(nextPath)
  }

  function goUp() {

    if (!currentPath) return

    const split = currentPath.split("/")

    split.pop()

    loadDirectory(split.join("/"))
  }

  return (

    <div className="explorer">

      <div className="explorer-header">

        <button
          className="nav-btn"
          onClick={goUp}
        >
          ⬅ Up
        </button>

        <div className="path-display">
          /{currentPath}
        </div>

      </div>

      <div className="file-list">

        {entries.map(entry => (

          <div
            key={entry.name}
            className="file-row"
            onClick={() => openEntry(entry)}
          >

            <div className="file-left">

              <span className="file-icon">
                {entry.type === "directory"
                  ? "📁"
                  : "📄"}
              </span>

              <span className="file-name">
                {entry.name}
              </span>

            </div>

            <div className="file-size">

              {entry.size
                ? `${(entry.size / 1024).toFixed(1)} KB`
                : ""}

            </div>

          </div>

        ))}

      </div>

    </div>
  )
}