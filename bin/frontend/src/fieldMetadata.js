// fieldMeta.js — Human-friendly labels and hints for pipeline config params.
// Add an entry for any top-level YAML key you want to surface in the questionnaire.
// Unknown keys fall back to showing the raw key name with a plain text input.
//
// Field shape:
//   label   — friendly display name (required)
//   hint    — short description shown below the input (optional)
//   type    — 'text' (default) | 'select' | 'boolean'
//   options — array of strings, required when type === 'select'

const FIELD_META = {
  pipeline: {
    label: 'Pipeline name',
    hint:  'Unique identifier for this pipeline run',
  },
  environment: {
    label:   'Environment',
    hint:    'Target deployment environment',
    type:    'select',
    options: ['dev', 'staging', 'prod'],
  },
  input_path: {
    label: 'Input path',
    hint:  'Path or URI to the source data',
  },
  output_path: {
    label: 'Output path',
    hint:  'Where results will be written',
  },
  schedule: {
    label: 'Schedule',
    hint:  'Cron expression or interval (e.g. @daily)',
  },
  retries: {
    label: 'Retries',
    hint:  'Number of retry attempts on failure',
  },
  timeout: {
    label: 'Timeout (s)',
    hint:  'Maximum run duration in seconds',
  },
  enabled: {
    label: 'Enabled',
    hint:  'Whether the pipeline is active',
    type:  'boolean',
  },
  // Add more keys here as your templates grow
}

export default FIELD_META