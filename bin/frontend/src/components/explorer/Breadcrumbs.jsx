export default function Breadcrumbs({ currentPath, onNavigate }) {
  const parts = currentPath ? currentPath.split('/') : []
  return (
    <div className="breadcrumbs">
      <button onClick={() => onNavigate('')}>root</button>
      {parts.map((part, idx) => {
        const path = parts.slice(0, idx + 1).join('/')
        return (
          <span key={path}>
            {' / '}
            <button onClick={() => onNavigate(path)}>{part}</button>
          </span>
        )
      })}
    </div>
  )
}
