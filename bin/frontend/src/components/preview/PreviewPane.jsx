const IMAGE_EXTS = new Set(['.png', '.jpg', '.jpeg', '.gif', '.webp', '.svg'])
const WRAP_EXTS  = new Set(['.md', '.txt', '.log'])
const BIO_EXTS   = new Set(['.fasta', '.fa', '.aln', '.afa', '.nwk', '.nex'])

function fmtSize(bytes) {
  if (bytes == null) return ''
  if (bytes < 1024) return `${bytes} B`
  if (bytes < 1024 * 1024) return `${(bytes / 1024).toFixed(1)} KB`
  return `${(bytes / 1024 / 1024).toFixed(1)} MB`
}

function fmtDate(ts) {
  return new Date(ts * 1000).toLocaleString()
}

// Colourises FASTA sequences: header lines purple, bases green
function FastaView({ content }) {
  const lines = content.split('\n')
  return (
    <pre className="preview-code fasta-view">
      {lines.map((line, i) => {
        if (line.startsWith('>')) {
          return <span key={i} className="fasta-header">{line}{'\n'}</span>
        }
        return <span key={i} className="fasta-seq">{line}{'\n'}</span>
      })}
    </pre>
  )
}

// Colourises YAML keys
function YamlView({ content }) {
  const lines = content.split('\n')
  return (
    <pre className="preview-code yaml-view">
      {lines.map((line, i) => {
        const m = line.match(/^(\s*)([\w-]+)(\s*:)(.*)$/)
        if (m) {
          return (
            <span key={i}>
              {m[1]}
              <span className="yaml-key">{m[2]}</span>
              <span className="yaml-colon">{m[3]}</span>
              <span className="yaml-val">{m[4]}</span>
              {'\n'}
            </span>
          )
        }
        return <span key={i}>{line}{'\n'}</span>
      })}
    </pre>
  )
}

// Plain text, word-wrapped for .md/.txt/.log
function TextPreview({ content, wrap }) {
  return (
    <pre className={`preview-code ${wrap ? 'wrapped' : ''}`}>
      {content}
    </pre>
  )
}

export default function PreviewPane({ fileData, loading, apiBase }) {
  if (loading) {
    return (
      <div className="preview-pane preview-empty">
        <span className="preview-spinner" aria-label="Loading">⏳</span>
      </div>
    )
  }

  if (!fileData) {
    return (
      <div className="preview-pane preview-empty">
        <span className="preview-empty-icon">👁</span>
        <p>Select a file to preview</p>
      </div>
    )
  }

  const { metadata } = fileData
  const ext = (metadata.extension || '').toLowerCase()

  // ── Binary / image ────────────────────────────────────
  if (fileData.kind === 'binary') {
    const isImage = IMAGE_EXTS.has(ext)
    const src = `${apiBase}/download?path=${encodeURIComponent(metadata.path)}`

    return (
      <div className="preview-pane">
        <div className="preview-header">
          <span className="preview-filename">{metadata.name}</span>
          <span className="preview-filemeta dim">{fmtSize(metadata.size)} · {metadata.mime}</span>
        </div>
        {isImage ? (
          <div className="preview-image-wrap">
            <img
              src={src}
              alt={metadata.name}
              className="preview-image"
            />
          </div>
        ) : (
          <div className="preview-binary">
            <p className="dim">Binary file — {metadata.mime}</p>
            <a
              className="btn primary"
              href={src}
              target="_blank"
              rel="noreferrer"
              style={{ display: 'inline-block', marginTop: '12px' }}
            >
              ↓ Download
            </a>
          </div>
        )}
        <div className="preview-footer dim">
          Modified {fmtDate(metadata.mtime)}
        </div>
      </div>
    )
  }

  // ── Text content ──────────────────────────────────────
  const { content } = fileData
  const isFasta = BIO_EXTS.has(ext)
  const isYaml  = ext === '.yaml' || ext === '.yml'
  const isWrap  = WRAP_EXTS.has(ext)

  let body
  if (isFasta) {
    body = <FastaView content={content} />
  } else if (isYaml) {
    body = <YamlView content={content} />
  } else {
    body = <TextPreview content={content} wrap={isWrap} />
  }

  return (
    <div className="preview-pane">
      <div className="preview-header">
        <span className="preview-filename">{metadata.name}</span>
        <span className="preview-filemeta dim">
          {fmtSize(metadata.size)} · {ext || metadata.mime}
        </span>
      </div>
      <div className="preview-body">
        {body}
      </div>
      <div className="preview-footer dim">
        Modified {fmtDate(metadata.mtime)}
      </div>
    </div>
  )
}