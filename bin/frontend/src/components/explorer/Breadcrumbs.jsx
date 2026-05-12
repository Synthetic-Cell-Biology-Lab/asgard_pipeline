export default function Breadcrumbs({ currentPath, onNavigate }) {
  const parts = currentPath ? currentPath.split('/') : []
  return (
    <div className="breadcrumbs">
      <button onClick={() => onNavigate('')}>~</button>
      {parts.map((part, idx) => {
        const path = parts.slice(0, idx + 1).join('/')
        return (
          <span key={path} className="breadcrumb-segment">
            <span className="breadcrumb-sep">/</span>
            <button onClick={() => onNavigate(path)}>{part}</button>
          </span>
        )
      })}
    </div>
  )
}