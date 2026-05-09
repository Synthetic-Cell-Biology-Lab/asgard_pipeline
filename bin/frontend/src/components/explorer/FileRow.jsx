export default function FileRow({ entry, selected, onClick }) {
  return (
    <div className={`file-row ${selected ? 'selected' : ''}`} onClick={onClick}>
      <span>{entry.type === 'directory' ? '📁' : '📄'}</span>
      <span className="file-name">{entry.name}</span>
      <span className="file-meta">{entry.type === 'file' && entry.size ? `${(entry.size / 1024).toFixed(1)} KB` : ''}</span>
    </div>
  )
}
