## vexp - Context-Aware AI Coding <!-- vexp v2.5.3 -->
vexp runs entirely on this machine: local daemon, index in `.vexp/`.
`run_pipeline` transmits nothing to any external service.
- `run_pipeline({ "task": "..." })` - orientation in one call (ranked pivot
  files with line ranges + blast radius + session notes) when a task does NOT
  name the files/symbols to touch. If it does, SKIP vexp - use your normal tools.
- `get_skeleton` - file structure at 70-90% token savings for files you only
  need to understand, not edit.
- `verify_done` - call once BEFORE declaring a multi-file task complete:
  returns mechanically broken references (imports of removed names, parse
  errors) and untouched dependents of the files you changed, with file:line.
- vexp may append a one-line hint to a prompt when orientation would help;
  otherwise it stays silent.

### Query shape (do this)
Anchor the task on real identifiers (ClassName, functionName) or file paths:
`run_pipeline({ "task": "fix JWT expiry in AuthService.validateToken" })`
<!-- /vexp -->