export default function PreviewPane({ fileData }) {
  if (!fileData) return <div className="preview-pane">Select a file to preview.</div>
  if (fileData.kind === 'binary') {
    return (
      <div className="preview-pane">
        <p>Binary file ({fileData.metadata.mime})</p>
        <a href={`http://localhost:8000${fileData.download_url}`} target="_blank" rel="noreferrer">Download / Open</a>
      </div>
    )
  }

  const ext = fileData.metadata.extension
  return (
    <div className="preview-pane">
      <h3>{fileData.metadata.name}</h3>
      <p className="preview-meta">{ext} · {fileData.metadata.mime}</p>
      <pre>{fileData.content}</pre>
    </div>
  )
}
