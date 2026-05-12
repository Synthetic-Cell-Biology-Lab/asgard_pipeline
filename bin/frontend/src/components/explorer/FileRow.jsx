// Maps file extensions to icons and colour hints
const EXT_META = {
  // images
  png:  { icon: '🖼', color: '#7F77DD' },
  jpg:  { icon: '🖼', color: '#7F77DD' },
  jpeg: { icon: '🖼', color: '#7F77DD' },
  gif:  { icon: '🖼', color: '#7F77DD' },
  svg:  { icon: '🖼', color: '#7F77DD' },
  webp: { icon: '🖼', color: '#7F77DD' },
  // text / data
  txt:  { icon: '📝', color: '#888' },
  md:   { icon: '📝', color: '#888' },
  log:  { icon: '📋', color: '#EF9F27' },
  csv:  { icon: '📊', color: '#1D9E75' },
  tsv:  { icon: '📊', color: '#1D9E75' },
  // bioinformatics
  fasta:{ icon: '🧬', color: '#E24B4A' },
  fa:   { icon: '🧬', color: '#E24B4A' },
  aln:  { icon: '🧬', color: '#E24B4A' },
  afa:  { icon: '🧬', color: '#E24B4A' },
  nwk:  { icon: '🌿', color: '#1D9E75' },
  nex:  { icon: '🌿', color: '#1D9E75' },
  // config
  yaml: { icon: '⚙', color: '#534AB7' },
  yml:  { icon: '⚙', color: '#534AB7' },
  json: { icon: '{ }', color: '#534AB7' },
  // shell
  sh:   { icon: '>', color: '#EF9F27' },
  py:   { icon: '🐍', color: '#EF9F27' },
}

function fmtSize(bytes) {
  if (bytes == null) return ''
  if (bytes < 1024) return `${bytes} B`
  if (bytes < 1024 * 1024) return `${(bytes / 1024).toFixed(1)} KB`
  return `${(bytes / 1024 / 1024).toFixed(1)} MB`
}

export default function FileRow({ entry, selected, onClick }) {
  const ext = (entry.extension || '').replace('.', '').toLowerCase()
  const meta = EXT_META[ext]

  const icon = entry.type === 'directory' ? '📁' : (meta?.icon ?? '📄')
  const iconColor = entry.type === 'directory' ? '#EF9F27' : (meta?.color ?? '#888')

  return (
    <div
      className={`file-row ${selected ? 'selected' : ''}`}
      onClick={onClick}
      title={entry.name}
    >
      <span className="file-icon" style={{ color: iconColor }}>{icon}</span>
      <span className="file-name">{entry.name}</span>
      {entry.type === 'file' && (
        <span className="file-meta">{fmtSize(entry.size)}</span>
      )}
      {entry.type === 'directory' && (
        <span className="file-meta">›</span>
      )}
    </div>
  )
}