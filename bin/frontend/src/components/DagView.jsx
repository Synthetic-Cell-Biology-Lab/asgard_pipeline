import { useEffect, useRef, useState, useCallback } from 'react'
import * as d3 from 'd3'

const API_BASE = 'http://localhost:8000'

const STATUS_COLOR = {
  pending: '#444466',
  running: '#C8960C',
  done:    '#1D9E75',
  error:   '#E24B4A',
  unknown: '#8A8AA6',
}

const STATUS_LABEL = {
  pending: 'Pending',
  running: 'Running',
  done:    'Done',
  error:   'Error',
  unknown: 'Unknown',
}

export default function DagView({ runId, running }) {
  const svgRef = useRef(null)
  const simRef = useRef(null)
  const [dag, setDag] = useState(null)
  const [error, setError] = useState(null)

  const fetchDag = useCallback(async () => {
    if (!runId) return
    try {
      const res = await fetch(`${API_BASE}/runs/${runId}/dag`)
      if (!res.ok) throw new Error(`${res.status} ${await res.text()}`)
      setDag(await res.json())
      setError(null)
    } catch (e) {
      setError(e.message)
    }
  }, [runId])

  // Initial fetch + poll while running
  useEffect(() => {
    fetchDag()
    if (!running) return
    const id = setInterval(fetchDag, 3000)
    return () => clearInterval(id)
  }, [fetchDag, running])

  // D3 render
  useEffect(() => {
    if (!dag || !svgRef.current) return
    const { nodes, edges } = dag
    if (!nodes.length) return

    const el   = svgRef.current
    const W    = el.clientWidth  || 700
    const H    = el.clientHeight || 420

    d3.select(el).selectAll('*').remove()

    const svg = d3.select(el)
      .attr('viewBox', `0 0 ${W} ${H}`)

    // Arrow marker
    svg.append('defs').append('marker')
      .attr('id', 'arrow')
      .attr('viewBox', '0 -5 10 10')
      .attr('refX', 28)
      .attr('refY', 0)
      .attr('markerWidth', 6)
      .attr('markerHeight', 6)
      .attr('orient', 'auto')
      .append('path')
      .attr('d', 'M0,-5L10,0L0,5')
      .attr('fill', '#555')

    const g = svg.append('g')

    // Zoom
    svg.call(d3.zoom().scaleExtent([0.3, 3]).on('zoom', e => {
      g.attr('transform', e.transform)
    }))

    // Build id→node index map
    const nodeById = Object.fromEntries(nodes.map(n => [n.id, n]))

    const links = edges.map(e => ({
      source: nodeById[e.source],
      target: nodeById[e.target],
    })).filter(l => l.source && l.target)

    // Simulation
    const sim = d3.forceSimulation(nodes)
      .force('link', d3.forceLink(links).distance(110).strength(0.8))
      .force('charge', d3.forceManyBody().strength(-320))
      .force('center', d3.forceCenter(W / 2, H / 2))
      .force('collision', d3.forceCollide(52))

    simRef.current = sim

    // Edges
    const link = g.append('g').selectAll('line')
      .data(links).join('line')
      .attr('stroke', '#555')
      .attr('stroke-width', 1.5)
      .attr('marker-end', 'url(#arrow)')

    // Nodes
    const node = g.append('g').selectAll('g')
      .data(nodes).join('g')
      .attr('cursor', 'grab')
      .call(d3.drag()
        .on('start', (e, d) => { if (!e.active) sim.alphaTarget(0.3).restart(); d.fx = d.x; d.fy = d.y })
        .on('drag',  (e, d) => { d.fx = e.x; d.fy = e.y })
        .on('end',   (e, d) => { if (!e.active) sim.alphaTarget(0); d.fx = null; d.fy = null })
      )

    // Node pill
    node.append('rect')
      .attr('x', -52).attr('y', -18)
      .attr('width', 104).attr('height', 36)
      .attr('rx', 18)
      .attr('fill', d => STATUS_COLOR[d.status] || STATUS_COLOR.pending)
      .attr('stroke', '#fff')
      .attr('stroke-width', 1.5)
      .attr('opacity', d => d.status === 'pending' ? 0.55 : 1)

    // Pulse ring for running nodes
    node.filter(d => d.status === 'running')
      .append('rect')
      .attr('x', -56).attr('y', -22)
      .attr('width', 112).attr('height', 44)
      .attr('rx', 22)
      .attr('fill', 'none')
      .attr('stroke', STATUS_COLOR.running)
      .attr('stroke-width', 2)
      .attr('opacity', 0.6)
      .attr('class', 'dag-pulse')

    // Label
    node.append('text')
      .attr('text-anchor', 'middle')
      .attr('dy', '0.35em')
      .attr('fill', '#fff')
      .attr('font-family', "'JetBrains Mono', monospace")
      .attr('font-size', 11)
      .attr('font-weight', 600)
      .text(d => d.label.length > 14 ? d.label.slice(0, 13) + '…' : d.label)

    // Tick
    sim.on('tick', () => {
      link
        .attr('x1', d => d.source.x).attr('y1', d => d.source.y)
        .attr('x2', d => d.target.x).attr('y2', d => d.target.y)
      node.attr('transform', d => `translate(${d.x},${d.y})`)
    })

    return () => sim.stop()
  }, [dag])

  if (!runId) return null

  const counts = dag ? {
    pending: dag.nodes.filter(n => n.status === 'pending').length,
    running: dag.nodes.filter(n => n.status === 'running').length,
    done:    dag.nodes.filter(n => n.status === 'done').length,
    error:   dag.nodes.filter(n => n.status === 'error').length,
    unknown: dag.nodes.filter(n => !STATUS_COLOR[n.status]).length,
  } : null

  return (
    <div className="dag-card card">
      <div className="card-header">
        <span className="step-label">Pipeline DAG</span>
        <button className="icon-btn" onClick={fetchDag} title="Refresh DAG">↺</button>
        {counts && (
          <div className="dag-legend">
            {Object.entries(STATUS_LABEL).map(([s, label]) => (
              counts[s] > 0 && (
                <span key={s} className="dag-legend-item">
                  <span className="dag-dot" style={{ background: STATUS_COLOR[s] }} />
                  {label} ({counts[s]})
                </span>
              )
            ))}
          </div>
        )}
      </div>

      {error && <p className="save-msg err-msg">{error}</p>}
      {dag?.debug?.unmatched_nodes?.length > 0 && (
        <p className="save-msg err-msg">
          Warning: {dag.debug.unmatched_nodes.length} DAG node(s) did not map to run logs.
        </p>
      )}

      {!dag && !error && (
        <p className="dim" style={{ fontSize: 12, padding: '8px 0' }}>
          Generating DAG…
        </p>
      )}

      <svg ref={svgRef} className="dag-svg" />
    </div>
  )
}
